// sdc-bench campaign (Algorithm 1): fault injection into scalar GEMM with every checker
// evaluated on the same corrupted output.
//
// Incremental evaluation with exact operation order: an injection recomputes only the
// outputs its fault reaches and only the affected row/column checksums, but each of these
// is recomputed with exactly the same sequence of rounded operations as a full
// recomputation would use. The independent verifier (verify.cpp) therefore must agree bit
// for bit; any disagreement is an implementation error, not a property of the method.
//
// Outputs (in --out DIR):
//   injections.bin  InjRec per injection      clean_checks.bin  CleanRec per fault-free check
//   problems.csv    one row per generated problem (scale, forward error, normaliser ratios)
//   campaign_meta.json
#include "common.hpp"
#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>

using namespace sdc;

struct Opts { std::string out = "results"; std::vector<int> precs{P64, P32, P16, PB16}; int K = 2000; int clean = 50; int threads = 0; int fwd_probs = 5; };

struct ProbRow { int cfg, prob; long double maxC, frobC; double fwd_err, ratio_med, ratio_max; int n_big; };
struct CfgOut { std::vector<InjRec> inj; std::vector<CleanRec> clean; std::vector<ProbRow> probs; };

template <class T> struct Prob {
  int n, prec; std::vector<T> A, B, C;
  std::vector<T> rowsum, colsum, rref, cref, d;     // d: clean discrepancies, rows then columns
  std::vector<long double> norm, absd, a, c;        // per check (rows then columns)
  std::vector<T> ref;                               // per check reference value
  long double maxC = 0, frobC = 0, u, g2n, kr;      // kr = rigorous/classical threshold ratio
  std::array<int, 9> top9{};                        // indices of the 9 largest |C|
  int n_big = 0;
};

template <class T> static void load(Prob<T>& P, int n, int prec, int dist, int seed, int prob) {
  std::vector<double> Ar, Br; gen_real(n, dist, seed, prob, Ar, Br);
  P.n = n; P.prec = prec; P.A.resize((size_t)n * n); P.B.resize((size_t)n * n);
  for (size_t t = 0; t < Ar.size(); ++t) {
    if constexpr (std::is_same<T, double>::value) { P.A[t] = Ar[t]; P.B[t] = Br[t]; }
    else { P.A[t] = round_in_f(Ar[t], prec); P.B[t] = round_in_f(Br[t], prec); }
  }
}

// Dot product of row i of A and column j of B, sequential order k = 0..n-1.
template <class T> static inline T dot(const Prob<T>& P, int i, int j) {
  const int n = P.n; T acc = 0;
  for (int k = 0; k < n; ++k) { T p = P.A[(size_t)i * n + k] * P.B[(size_t)k * n + j]; acc = acc + p; }
  return acc;
}
template <class T> static inline T rowsum_of(const std::vector<T>& C, int n, int i) { T s = 0; for (int j = 0; j < n; ++j) s = s + C[(size_t)i * n + j]; return s; }
template <class T> static inline T colsum_of(const std::vector<T>& C, int n, int j) { T s = 0; for (int i = 0; i < n; ++i) s = s + C[(size_t)i * n + j]; return s; }

template <class T> static void compute_clean(Prob<T>& P) {
  const int n = P.n; P.C.resize((size_t)n * n);
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) P.C[(size_t)i * n + j] = dot(P, i, j);
  P.rowsum.resize(n); P.colsum.resize(n); P.rref.resize(n); P.cref.resize(n);
  for (int i = 0; i < n; ++i) P.rowsum[i] = rowsum_of(P.C, n, i);
  for (int j = 0; j < n; ++j) P.colsum[j] = colsum_of(P.C, n, j);
  // references from clean operands, working precision: A(B1) and (1^T A)B
  std::vector<T> B1(n), A1(n);
  for (int k = 0; k < n; ++k) { T s = 0; for (int j = 0; j < n; ++j) s = s + P.B[(size_t)k * n + j]; B1[k] = s; }
  for (int k = 0; k < n; ++k) { T s = 0; for (int i = 0; i < n; ++i) s = s + P.A[(size_t)i * n + k]; A1[k] = s; }
  for (int i = 0; i < n; ++i) { T s = 0; for (int k = 0; k < n; ++k) { T p = P.A[(size_t)i * n + k] * B1[k]; s = s + p; } P.rref[i] = s; }
  for (int j = 0; j < n; ++j) { T s = 0; for (int k = 0; k < n; ++k) { T p = A1[k] * P.B[(size_t)k * n + j]; s = s + p; } P.cref[j] = s; }
  // normalisers |A|(|B|1) and (1^T|A|)|B| in extended precision
  std::vector<long double> nB1(n, 0), nA1(n, 0);
  for (int k = 0; k < n; ++k) for (int j = 0; j < n; ++j) nB1[k] += std::fabs((long double)P.B[(size_t)k * n + j]);
  for (int k = 0; k < n; ++k) for (int i = 0; i < n; ++i) nA1[k] += std::fabs((long double)P.A[(size_t)i * n + k]);
  P.norm.assign(2 * n, 0); P.ref.resize(2 * n); P.d.resize(2 * n);
  for (int i = 0; i < n; ++i) { long double s = 0; for (int k = 0; k < n; ++k) s += std::fabs((long double)P.A[(size_t)i * n + k]) * nB1[k]; P.norm[i] = s; }
  for (int j = 0; j < n; ++j) { long double s = 0; for (int k = 0; k < n; ++k) s += nA1[k] * std::fabs((long double)P.B[(size_t)k * n + j]); P.norm[n + j] = s; }
  for (int i = 0; i < n; ++i) { P.ref[i] = P.rref[i]; P.d[i] = P.rowsum[i] - P.rref[i]; }
  for (int j = 0; j < n; ++j) { P.ref[n + j] = P.cref[j]; P.d[n + j] = P.colsum[j] - P.cref[j]; }
  P.u = unit_roundoff(P.prec); P.g2n = gamma_k(2.0L * n, P.u);
  P.kr = 2.0L * gamma_k(2.0L * n - 1, P.u) * (1.0L + P.u) / P.g2n;
  // scale of the output
  long double mx = 0, fr = 0; P.n_big = 0;
  std::vector<std::pair<long double, int>> mags((size_t)n * n);
  for (size_t t = 0; t < P.C.size(); ++t) {
    long double v = std::fabs((long double)P.C[t]); mx = std::max(mx, v); fr += v * v; mags[t] = {v, (int)t};
    if (v >= RANGE_BOUND) P.n_big++;
  }
  P.maxC = mx; P.frobC = std::sqrt(fr);
  std::partial_sort(mags.begin(), mags.begin() + 9, mags.end(), [](auto& x, auto& y) { return x.first > y.first || (x.first == y.first && x.second < y.second); });
  for (int t = 0; t < 9; ++t) P.top9[t] = mags[t].second;
  // clean check statistics
  P.absd.resize(2 * n); P.a.resize(2 * n); P.c.resize(2 * n);
  for (int t = 0; t < 2 * n; ++t) {
    long double ad = std::fabs((long double)P.d[t]), Tr = P.g2n * P.norm[t];
    P.absd[t] = ad;
    P.a[t] = Tr > 0 ? ad / Tr : (ad > 0 ? INFINITY : 0);
    P.c[t] = ad / (FIXED_ATOL + FIXED_RTOL * std::fabs((long double)P.ref[t]));
  }
}

static inline long double safe_ratio(long double num, long double den) { return den > 0 ? num / den : (num > 0 ? INFINITY : 0); }

template <class T> static void run_config(int cfg, const Opts& o, CfgOut& out) {
  Config cf = cfg_of(cfg); const int n = cf.n;
  // ---- fault-free problems (problem 0 also receives the injections)
  Prob<T> P0;
  for (int prob = 0; prob < o.clean; ++prob) {
    Prob<T> Pp; load(Pp, n, cf.prec, cf.dist, cf.seed, prob); compute_clean(Pp);
    for (int t = 0; t < 2 * n; ++t) {
      CleanRec r{(uint32_t)cfg, (uint16_t)prob, (uint16_t)t, (double)Pp.a[t], (double)safe_ratio(Pp.absd[t], Pp.maxC), (double)Pp.c[t]};
      out.clean.push_back(r);
    }
    ProbRow pr{cfg, prob, Pp.maxC, Pp.frobC, NAN, 0, 0, Pp.n_big};
    std::vector<long double> rat(n); for (int i = 0; i < n; ++i) rat[i] = safe_ratio(Pp.norm[i], Pp.maxC);
    std::sort(rat.begin(), rat.end()); pr.ratio_med = (double)rat[n / 2]; pr.ratio_max = (double)rat[n - 1];
    if (prob < o.fwd_probs) {   // forward error of the fault-free kernel against an extended-precision product
      long double em = 0, cm = 0;
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        long double s = 0; for (int k = 0; k < n; ++k) s += (long double)Pp.A[(size_t)i * n + k] * (long double)Pp.B[(size_t)k * n + j];
        em = std::max(em, std::fabs((long double)Pp.C[(size_t)i * n + j] - s)); cm = std::max(cm, std::fabs(s));
      }
      pr.fwd_err = (double)safe_ratio(em, cm);
    }
    out.probs.push_back(pr);
    if (prob == 0) P0 = std::move(Pp);
  }
  Prob<T>& P = P0;
  const Fmt fa = acc_fmt(cf.prec), fi = in_fmt(cf.prec);
  std::vector<T> rs_new(n), cs_new(n);
  for (int fm = 0; fm < 3; ++fm) {
    SiteStream ss(n, cf.dist, cf.seed, fm, cf.prec);
    for (int t = 0; t < o.K; ++t) {
      Site s = ss.next();
      // ---- affected outputs, recomputed with the fault (same operation order as dot())
      int na = 0; std::array<int, REUSE> pos{}; std::array<T, REUSE> nv{}, ov{};
      if (fm == ACC) {
        T acc = 0;
        for (int k = 0; k < n; ++k) { T p = P.A[(size_t)s.i * n + k] * P.B[(size_t)k * n + s.j]; acc = acc + p; if (k == s.k) acc = flip_acc(acc, s.bit); }
        pos[0] = s.i * n + s.j; nv[0] = acc; na = 1;
      } else if (fm == REUSE_A) {
        T ap = flip_in(P.A[(size_t)s.i * n + s.k], s.bit, cf.prec);
        for (int q = 0; q < REUSE; ++q) {
          int j = s.j + q; T acc = 0;
          for (int k = 0; k < n; ++k) { T av = (k == s.k) ? ap : P.A[(size_t)s.i * n + k]; T p = av * P.B[(size_t)k * n + j]; acc = acc + p; }
          pos[q] = s.i * n + j; nv[q] = acc;
        }
        na = REUSE;
      } else {
        T bp = flip_in(P.B[(size_t)s.k * n + s.j], s.bit, cf.prec);
        for (int q = 0; q < REUSE; ++q) {
          int i = s.i + q; T acc = 0;
          for (int k = 0; k < n; ++k) { T bv = (k == s.k) ? bp : P.B[(size_t)k * n + s.j]; T p = P.A[(size_t)i * n + k] * bv; acc = acc + p; }
          pos[q] = i * n + s.j; nv[q] = acc;
        }
        na = REUSE;
      }
      for (int q = 0; q < na; ++q) ov[q] = P.C[pos[q]];
      // ---- harm, flags, output scale
      InjRec r{}; r.cfg = cfg; r.idx = (uint32_t)(fm * o.K + t); r.fm = (uint8_t)fm; r.bit = (uint8_t)s.bit;
      r.bitclass = (uint8_t)bit_class(fm == ACC ? fa : fi, s.bit);
      bool nonfin = false, any_nan = false, any_inf = false, changed = false, big_aff = false;
      long double dmax = 0, dss = 0, cw = 0, Mnew = 0; int big_clean_aff = 0;
      for (int q = 0; q < na; ++q) {
        T v = nv[q], w = ov[q];
        if (std::memcmp(&v, &w, sizeof(T)) != 0) changed = true;
        if (std::isnan(v)) { any_nan = true; nonfin = true; } else if (std::isinf(v)) { any_inf = true; nonfin = true; }
        if (std::fabs((long double)w) >= RANGE_BOUND) big_clean_aff++;
        if (!std::isfinite(v) || std::fabs((long double)v) >= RANGE_BOUND) big_aff = true;
        if (std::isfinite(v)) {
          Mnew = std::max(Mnew, std::fabs((long double)v));
          long double dl = std::fabs((long double)v - (long double)w);
          dmax = std::max(dmax, dl); dss += dl * dl;
          long double c1 = std::fabs((long double)w) > 0 ? dl / std::fabs((long double)w) : (dl > 0 ? INFINITY : 0);
          cw = std::max(cw, c1);
        }
      }
      for (int q2 = 0; q2 < 9; ++q2) {          // largest unaffected clean magnitude
        int id = P.top9[q2]; bool aff = false; for (int q = 0; q < na; ++q) aff |= (pos[q] == id);
        if (!aff) { Mnew = std::max(Mnew, std::fabs((long double)P.C[id])); break; }
      }
      if (nonfin) { r.r_max = r.r_fro = r.r_cw = INFINITY; }
      else { r.r_max = (double)safe_ratio(dmax, P.maxC); r.r_fro = (double)safe_ratio(std::sqrt(dss), P.frobC); r.r_cw = (double)cw; }
      uint8_t fl = 0;
      if (changed) fl |= F_BITWISE; if (any_nan) fl |= F_NAN; if (any_inf) fl |= F_INF;
      if (big_aff || (P.n_big - big_clean_aff) > 0) fl |= F_RANGE;
      // ---- affected checksums, same summation order as compute_clean()
      for (int q = 0; q < na; ++q) P.C[pos[q]] = nv[q];
      std::array<int, REUSE + 1> rows{}, cols{}; int nr = 0, nc = 0;
      if (fm == ACC) { rows[nr++] = s.i; cols[nc++] = s.j; }
      else if (fm == REUSE_A) { rows[nr++] = s.i; for (int q = 0; q < REUSE; ++q) cols[nc++] = s.j + q; }
      else { cols[nc++] = s.j; for (int q = 0; q < REUSE; ++q) rows[nr++] = s.i + q; }
      for (int q = 0; q < nr; ++q) rs_new[q] = rowsum_of(P.C, n, rows[q]);
      for (int q = 0; q < nc; ++q) cs_new[q] = colsum_of(P.C, n, cols[q]);
      for (int q = 0; q < na; ++q) P.C[pos[q]] = ov[q];   // restore
      // naive (unguarded) max|C'|: IEEE comparisons, NaN ignored, Inf taken
      long double Mnaive = 0;
      { for (int q2 = 0; q2 < 9; ++q2) { int id = P.top9[q2]; bool aff = false; for (int q = 0; q < na; ++q) aff |= (pos[q] == id); if (!aff) { Mnaive = std::fabs((long double)P.C[id]); break; } }
        for (int q = 0; q < na; ++q) { long double v = std::fabs((long double)nv[q]); if (v > Mnaive) Mnaive = v; } }
      // ---- all 2n checks
      long double amax = 0, absdmax = 0, cmax = 0, cbd = -1, cbd_rig = -1;
      bool ung_rb = false, ung_fx = false, ung_c6 = false, ung_c3 = false;
      auto visit = [&](long double ad, long double a, long double c, long double Tr, long double Tf, long double rawd_abs) {
        amax = std::max(amax, a); absdmax = std::max(absdmax, ad); cmax = std::max(cmax, c);
        if (a > 1) cbd = std::max(cbd, ad);
        if (a > P.kr) cbd_rig = std::max(cbd_rig, ad);
        // unguarded comparisons on the raw |d| (NaN compares false, Inf true)
        if (rawd_abs > Tr) ung_rb = true;
        if (rawd_abs > Tf) ung_fx = true;
        if (rawd_abs > Tr && rawd_abs > 1e-6L * Mnaive) ung_c6 = true;
        if (rawd_abs > Tr && rawd_abs > 1e-3L * Mnaive) ung_c3 = true;
      };
      for (int t2 = 0; t2 < 2 * n; ++t2) {
        bool isrow = t2 < n; int idx = isrow ? t2 : t2 - n; int q = -1;
        if (isrow) { for (int z = 0; z < nr; ++z) if (rows[z] == idx) q = z; }
        else { for (int z = 0; z < nc; ++z) if (cols[z] == idx) q = z; }
        long double Tr = P.g2n * P.norm[t2], Tf = FIXED_ATOL + FIXED_RTOL * std::fabs((long double)P.ref[t2]);
        if (q < 0) { visit(P.absd[t2], P.a[t2], P.c[t2], Tr, Tf, P.absd[t2]); continue; }
        T dn = (isrow ? rs_new[q] : cs_new[q]) - P.ref[t2];
        long double raw = std::fabs((long double)dn);
        if (!std::isfinite(dn)) { visit(INFINITY, INFINITY, INFINITY, Tr, Tf, raw); continue; }   // guarded: fires
        visit(raw, Tr > 0 ? raw / Tr : (raw > 0 ? INFINITY : 0), raw / Tf, Tr, Tf, raw);
      }
      r.amax = (double)amax; r.cmax = (double)cmax;
      r.bmax = (double)safe_ratio(absdmax, Mnew);
      r.cb = cbd < 0 ? 0.0 : (double)safe_ratio(cbd, Mnew);
      r.cb_rig = cbd_rig < 0 ? 0.0 : (double)safe_ratio(cbd_rig, Mnew);
      if (ung_rb) fl |= F_UNG_RB; if (ung_fx) fl |= F_UNG_FIXED; if (ung_c6) fl |= F_UNG_COMB6; if (ung_c3) fl |= F_UNG_COMB3;
      r.flags = fl;
      out.inj.push_back(r);
    }
  }
}

static std::vector<int> parse_precs(const std::string& s) {
  std::vector<int> v; size_t p = 0;
  while (p <= s.size()) { size_t q = s.find(',', p); std::string t = s.substr(p, q == std::string::npos ? std::string::npos : q - p);
    for (int k = 0; k < 4; ++k) if (t == PREC_NAME[k]) v.push_back(k);
    if (q == std::string::npos) break; p = q + 1; }
  return v;
}

int main(int argc, char** argv) {
  Opts o;
  for (int a = 1; a < argc; ++a) {
    std::string k = argv[a]; auto val = [&]() { return std::string(argv[++a]); };
    if (k == "--out") o.out = val(); else if (k == "--precisions") o.precs = parse_precs(val());
    else if (k == "--K") o.K = std::stoi(val()); else if (k == "--clean") o.clean = std::stoi(val());
    else if (k == "--threads") o.threads = std::stoi(val()); else if (k == "--fwd-probs") o.fwd_probs = std::stoi(val());
    else { std::fprintf(stderr, "unknown option %s\n", k.c_str()); return 2; }
  }
  if (o.threads <= 0) o.threads = std::max(1u, std::thread::hardware_concurrency());
  std::vector<int> cfgs;
  for (int p : o.precs) for (int si = 0; si < 3; ++si) for (int d = 0; d < 5; ++d) for (int s = 0; s < NSEEDS; ++s) cfgs.push_back(cfg_id(p, si, d, s));
  std::vector<CfgOut> outs(cfgs.size());
  std::atomic<size_t> next{0}; std::mutex mu; size_t done = 0; auto t0 = std::chrono::steady_clock::now();
  // largest problems first for load balance; output order is fixed afterwards
  std::vector<size_t> order(cfgs.size()); for (size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](size_t x, size_t y) { return cfg_of(cfgs[x]).n > cfg_of(cfgs[y]).n; });
  auto worker = [&]() {
    for (;;) {
      size_t w = next++; if (w >= order.size()) return;
      size_t idx = order[w]; int cfg = cfgs[idx];
      if (cfg_of(cfg).prec == P64) run_config<double>(cfg, o, outs[idx]); else run_config<float>(cfg, o, outs[idx]);
      std::lock_guard<std::mutex> lk(mu); ++done;
      if (done * 10 / cfgs.size() != (done - 1) * 10 / cfgs.size() || done == cfgs.size()) {
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::fprintf(stderr, "[campaign] %3zu%%  (%zu/%zu configurations, %.0f s)\n", done * 100 / cfgs.size(), done, cfgs.size(), el);
      }
    }
  };
  std::vector<std::thread> th; for (int t = 0; t < o.threads; ++t) th.emplace_back(worker); for (auto& x : th) x.join();
  std::string base = o.out + "/";
  FILE* fi = std::fopen((base + "injections.bin").c_str(), "wb");
  FILE* fc = std::fopen((base + "clean_checks.bin").c_str(), "wb");
  FILE* fp = std::fopen((base + "problems.csv").c_str(), "w");
  if (!fi || !fc || !fp) { std::fprintf(stderr, "cannot write to %s (does the directory exist?)\n", o.out.c_str()); return 1; }
  std::fprintf(fp, "cfg,prec,n,dist,seed,prob,maxC,frobC,fwd_err,ratio_med,ratio_max,n_big\n");
  size_t ni = 0, nc = 0;
  for (auto& co : outs) {
    std::fwrite(co.inj.data(), sizeof(InjRec), co.inj.size(), fi); ni += co.inj.size();
    std::fwrite(co.clean.data(), sizeof(CleanRec), co.clean.size(), fc); nc += co.clean.size();
    for (auto& p : co.probs) { Config c = cfg_of(p.cfg);
      std::fprintf(fp, "%d,%s,%d,%s,%d,%d,%.17Lg,%.17Lg,%.17g,%.17g,%.17g,%d\n", p.cfg, PREC_NAME[c.prec], c.n, DIST_NAME[c.dist], c.seed, p.prob, p.maxC, p.frobC, p.fwd_err, p.ratio_med, p.ratio_max, p.n_big); }
  }
  std::fclose(fi); std::fclose(fc); std::fclose(fp);
  FILE* fm = std::fopen((base + "campaign_meta.json").c_str(), "w");
  std::fprintf(fm, "{\"injections\": %zu, \"clean_checks\": %zu, \"K\": %d, \"clean_problems_per_config\": %d, \"precisions\": [", ni, nc, o.K, o.clean);
  for (size_t i = 0; i < o.precs.size(); ++i) std::fprintf(fm, "%s\"%s\"", i ? ", " : "", PREC_NAME[o.precs[i]]);
  std::fprintf(fm, "], \"injrec_bytes\": %zu, \"cleanrec_bytes\": %zu}\n", sizeof(InjRec), sizeof(CleanRec));
  std::fclose(fm);
  std::fprintf(stderr, "[campaign] wrote %zu injections, %zu fault-free checks\n", ni, nc);
  return 0;
}
