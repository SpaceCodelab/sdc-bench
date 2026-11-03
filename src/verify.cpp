// sdc-bench verifier (Algorithm 3): an independent re-implementation that rebuilds every
// sampled faulty product and all checksums from scratch, then classifies each discrepancy
// into the three phenomena that must be kept apart:
//   (1) genuine checksum blind spot: the corruption changes the output but leaves every
//       affected checksum discrepancy bit-identical to the fault-free one;
//   (2) floating-point ordering: a checker decision flips when the checksum side is
//       summed in a different, equally valid order (pairwise instead of sequential);
//   (3) implementation error: the verifier and the campaign disagree bit for bit.
// Only input generation (common.hpp) is shared with the campaign.
#include "common.hpp"
#include <thread>
#include <atomic>
#include <mutex>
#include <map>

using namespace sdc;

struct Counts {
  long verified = 0, mism_output = 0, mism_checks = 0;                 // (3)
  long changed = 0, cancel = 0, cancel_h6 = 0, cancel_h3 = 0;           // (1)
  long rbmiss_h6 = 0, h6 = 0;                                           // harmful at 1e-6 missed by RB
  long flip_rb = 0, flip_fixed = 0, flip_comb6 = 0, flip_comb3 = 0, flip_tol6 = 0, flip_tol3 = 0;   // (2)
  void add(const Counts& o) {
    verified += o.verified; mism_output += o.mism_output; mism_checks += o.mism_checks; changed += o.changed;
    cancel += o.cancel; cancel_h6 += o.cancel_h6; cancel_h3 += o.cancel_h3; rbmiss_h6 += o.rbmiss_h6; h6 += o.h6;
    flip_rb += o.flip_rb; flip_fixed += o.flip_fixed; flip_comb6 += o.flip_comb6; flip_comb3 += o.flip_comb3;
    flip_tol6 += o.flip_tol6; flip_tol3 += o.flip_tol3;
  }
};
struct Example { int cfg, fm, bit, bitclass; double r_max; long double dclean, dfault; };

template <class T> static T seq_sum(const T* x, int n, int stride) { T s = 0; for (int i = 0; i < n; ++i) s = s + x[(size_t)i * stride]; return s; }
template <class T> static T pair_sum(const T* x, int n, int stride) {
  if (n <= 2) { T s = 0; for (int i = 0; i < n; ++i) s = s + x[(size_t)i * stride]; return s; }
  int h = n / 2; T a = pair_sum(x, h, stride), b = pair_sum(x + (size_t)h * stride, n - h, stride); return a + b;
}

template <class T> struct Checks {   // all 2n checks of one output, computed from scratch
  std::vector<T> d, ref; std::vector<long double> norm;
};
template <class T> static void checks_from_scratch(const std::vector<T>& A, const std::vector<T>& B, const std::vector<T>& C, int n, bool pairwise, Checks<T>& K) {
  auto sum = [&](const T* x, int m, int st) { return pairwise ? pair_sum(x, m, st) : seq_sum(x, m, st); };
  std::vector<T> B1(n), A1(n), prod(n);
  for (int k = 0; k < n; ++k) B1[k] = sum(&B[(size_t)k * n], n, 1);
  for (int k = 0; k < n; ++k) A1[k] = sum(&A[k], n, n);
  K.d.resize(2 * n); K.ref.resize(2 * n); K.norm.resize(2 * n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) prod[k] = A[(size_t)i * n + k] * B1[k];
    K.ref[i] = sum(prod.data(), n, 1); K.d[i] = sum(&C[(size_t)i * n], n, 1) - K.ref[i];
  }
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) prod[k] = A1[k] * B[(size_t)k * n + j];
    K.ref[n + j] = sum(prod.data(), n, 1); K.d[n + j] = sum(&C[j], n, n) - K.ref[n + j];
  }
  std::vector<long double> nB1(n, 0), nA1(n, 0);
  for (int k = 0; k < n; ++k) for (int j = 0; j < n; ++j) nB1[k] += std::fabs((long double)B[(size_t)k * n + j]);
  for (int k = 0; k < n; ++k) for (int i = 0; i < n; ++i) nA1[k] += std::fabs((long double)A[(size_t)i * n + k]);
  for (int i = 0; i < n; ++i) { long double s = 0; for (int k = 0; k < n; ++k) s += std::fabs((long double)A[(size_t)i * n + k]) * nB1[k]; K.norm[i] = s; }
  for (int j = 0; j < n; ++j) { long double s = 0; for (int k = 0; k < n; ++k) s += nA1[k] * std::fabs((long double)B[(size_t)k * n + j]); K.norm[n + j] = s; }
}

static inline long double sr(long double a, long double b) { return b > 0 ? a / b : (a > 0 ? INFINITY : 0); }

template <class T> struct Stats { long double amax = 0, absdmax = 0, cmax = 0, cbd = -1, cbd_rig = -1; bool ung_rb = false, ung_fx = false, ung_c6 = false, ung_c3 = false; };
template <class T> static Stats<T> stats_of(const Checks<T>& K, int n, int prec, long double Mnaive) {
  Stats<T> s; long double u = unit_roundoff(prec), g = gamma_k(2.0L * n, u), kr = 2.0L * gamma_k(2.0L * n - 1, u) * (1.0L + u) / g;
  for (int t = 0; t < 2 * n; ++t) {
    long double raw = std::fabs((long double)K.d[t]), Tr = g * K.norm[t], Tf = FIXED_ATOL + FIXED_RTOL * std::fabs((long double)K.ref[t]);
    long double ad, a, c;
    if (!std::isfinite(K.d[t])) { ad = a = c = INFINITY; } else { ad = raw; a = Tr > 0 ? raw / Tr : (raw > 0 ? INFINITY : 0); c = raw / Tf; }
    s.amax = std::max(s.amax, a); s.absdmax = std::max(s.absdmax, ad); s.cmax = std::max(s.cmax, c);
    if (a > 1) s.cbd = std::max(s.cbd, ad);
    if (a > kr) s.cbd_rig = std::max(s.cbd_rig, ad);
    if (raw > Tr) s.ung_rb = true; if (raw > Tf) s.ung_fx = true;
    if (raw > Tr && raw > 1e-6L * Mnaive) s.ung_c6 = true; if (raw > Tr && raw > 1e-3L * Mnaive) s.ung_c3 = true;
  }
  return s;
}

template <class T> struct CfgData { std::vector<T> A, B, C; Checks<T> Kc; };
template <class T> static void prep(int cfg, CfgData<T>& D) {
  Config cf = cfg_of(cfg); int n = cf.n; std::vector<double> Ar, Br; gen_real(n, cf.dist, cf.seed, 0, Ar, Br);
  D.A.resize(Ar.size()); D.B.resize(Br.size());
  for (size_t t = 0; t < Ar.size(); ++t) {
    if constexpr (std::is_same<T, double>::value) { D.A[t] = Ar[t]; D.B[t] = Br[t]; }
    else { D.A[t] = round_in_f(Ar[t], cf.prec); D.B[t] = round_in_f(Br[t], cf.prec); }
  }
  D.C.assign((size_t)n * n, 0);
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) { T acc = 0; for (int k = 0; k < n; ++k) { T p = D.A[(size_t)i * n + k] * D.B[(size_t)k * n + j]; acc = acc + p; } D.C[(size_t)i * n + j] = acc; }
  checks_from_scratch(D.A, D.B, D.C, n, false, D.Kc);
}

static bool same(double x, double y) { return std::memcmp(&x, &y, 8) == 0; }

template <class T> static void verify_one(int cfg, const InjRec& r, const Site& s, const CfgData<T>& D, Counts& cnt, std::vector<Example>& ex) {
  Config cf = cfg_of(cfg); const int n = cf.n; const int fm = r.fm;
  // ---- full faulty product, every element recomputed
  std::vector<T> Cf((size_t)n * n);
  T ap = 0, bp = 0;
  if (fm == REUSE_A) ap = flip_in(D.A[(size_t)s.i * n + s.k], s.bit, cf.prec);
  if (fm == REUSE_B) bp = flip_in(D.B[(size_t)s.k * n + s.j], s.bit, cf.prec);
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
    bool hitA = fm == ACC && i == s.i && j == s.j;
    bool hitRA = fm == REUSE_A && i == s.i && j >= s.j && j < s.j + REUSE;
    bool hitRB = fm == REUSE_B && j == s.j && i >= s.i && i < s.i + REUSE;
    T acc = 0;
    for (int k = 0; k < n; ++k) {
      T a = (hitRA && k == s.k) ? ap : D.A[(size_t)i * n + k];
      T b = (hitRB && k == s.k) ? bp : D.B[(size_t)k * n + j];
      T p = a * b; acc = acc + p;
      if (hitA && k == s.k) acc = flip_acc(acc, s.bit);
    }
    Cf[(size_t)i * n + j] = acc;
  }
  // ---- harm and flags over the whole output
  long double maxC = 0, fr = 0; int nbig_clean = 0;
  for (auto v : D.C) { long double a = std::fabs((long double)v); maxC = std::max(maxC, a); fr += a * a; if (a >= RANGE_BOUND) nbig_clean++; }
  long double frob = std::sqrt(fr);
  bool changed = false, nan = false, inf = false, range = false; long double dmax = 0, dss = 0, cw = 0, M = 0, Mnaive = 0;
  for (size_t t = 0; t < Cf.size(); ++t) {
    T v = Cf[t], w = D.C[t];
    if (std::memcmp(&v, &w, sizeof(T)) != 0) changed = true;
    if (std::isnan(v)) nan = true; else if (std::isinf(v)) inf = true;
    if (!std::isfinite(v) || std::fabs((long double)v) >= RANGE_BOUND) range = true;
    long double av = std::fabs((long double)v); if (av > Mnaive) Mnaive = av;
    if (std::isfinite(v)) {
      M = std::max(M, av);
      long double dl = std::fabs((long double)v - (long double)w);
      if (dl > 0) { dmax = std::max(dmax, dl); dss += dl * dl; cw = std::max(cw, std::fabs((long double)w) > 0 ? dl / std::fabs((long double)w) : (long double)INFINITY); }
    }
  }
  bool nonfin = nan || inf;
  double rmax = nonfin ? INFINITY : (double)sr(dmax, maxC), rfro = nonfin ? INFINITY : (double)sr(std::sqrt(dss), frob), rcw = nonfin ? INFINITY : (double)cw;
  uint8_t fl = 0; if (changed) fl |= F_BITWISE; if (nan) fl |= F_NAN; if (inf) fl |= F_INF; if (range) fl |= F_RANGE;
  // ---- checks from scratch (sequential order, as specified)
  Checks<T> K; checks_from_scratch(D.A, D.B, Cf, n, false, K);
  Stats<T> st = stats_of(K, n, cf.prec, Mnaive);
  if (st.ung_rb) fl |= F_UNG_RB; if (st.ung_fx) fl |= F_UNG_FIXED; if (st.ung_c6) fl |= F_UNG_COMB6; if (st.ung_c3) fl |= F_UNG_COMB3;
  double amax = (double)st.amax, cmax = (double)st.cmax, bmax = (double)sr(st.absdmax, M);
  double cb = st.cbd < 0 ? 0.0 : (double)sr(st.cbd, M), cbr = st.cbd_rig < 0 ? 0.0 : (double)sr(st.cbd_rig, M);
  cnt.verified++;
  // (3) implementation error: bit-for-bit comparison with the campaign record
  bool out_ok = same(rmax, r.r_max) && same(rfro, r.r_fro) && same(rcw, r.r_cw) && ((fl & (F_BITWISE | F_NAN | F_INF | F_RANGE)) == (r.flags & (F_BITWISE | F_NAN | F_INF | F_RANGE)));
  bool chk_ok = same(amax, r.amax) && same(bmax, r.bmax) && same(cmax, r.cmax) && same(cb, r.cb) && same(cbr, r.cb_rig) && fl == r.flags;
  if (!out_ok) cnt.mism_output++; else if (!chk_ok) cnt.mism_checks++;
  // (1) genuine blind spot: output changed, every affected checksum discrepancy unchanged
  if (changed) {
    cnt.changed++;
    std::vector<int> aff;
    if (fm == ACC) { aff = {s.i, n + s.j}; }
    else if (fm == REUSE_A) { aff.push_back(s.i); for (int q = 0; q < REUSE; ++q) aff.push_back(n + s.j + q); }
    else { aff.push_back(n + s.j); for (int q = 0; q < REUSE; ++q) aff.push_back(s.i + q); }
    bool all_same = true; for (int t : aff) { T a = K.d[t], b = D.Kc.d[t]; if (std::memcmp(&a, &b, sizeof(T)) != 0) all_same = false; }
    if (all_same && !nonfin) {
      cnt.cancel++; if (rmax > 1e-6) cnt.cancel_h6++; if (rmax > 1e-3) cnt.cancel_h3++;
      if (ex.size() < 40) ex.push_back({cfg, fm, s.bit, r.bitclass, rmax, (long double)D.Kc.d[aff[0]], (long double)K.d[aff[0]]});
    }
  }
  bool rb_fire = nonfin || st.amax > 1;
  if (rmax > 1e-6) { cnt.h6++; if (!rb_fire) cnt.rbmiss_h6++; }
  // (2) ordering: same faulty output, checksum side summed pairwise
  Checks<T> Kp; checks_from_scratch(D.A, D.B, Cf, n, true, Kp);
  Stats<T> sp = stats_of(Kp, n, cf.prec, Mnaive);
  auto dec = [&](const Stats<T>& x, int which) -> bool {
    if (nonfin) return true;
    long double bm = sr(x.absdmax, M), cbm = x.cbd < 0 ? 0 : sr(x.cbd, M);
    switch (which) { case 0: return x.amax > 1; case 1: return x.cmax > 1; case 2: return cbm > 1e-6L / (1 + 1e-6L); case 3: return cbm > 1e-3L / (1 + 1e-3L); case 4: return bm > 1e-6L / (1 + 1e-6L); default: return bm > 1e-3L / (1 + 1e-3L); }
  };
  long* fl_[6] = {&cnt.flip_rb, &cnt.flip_fixed, &cnt.flip_comb6, &cnt.flip_comb3, &cnt.flip_tol6, &cnt.flip_tol3};
  for (int w = 0; w < 6; ++w) if (dec(st, w) != dec(sp, w)) (*fl_[w])++;
}

int main(int argc, char** argv) {
  std::string dir = "results"; int stride = 10, threads = 0;
  for (int a = 1; a < argc; ++a) { std::string k = argv[a];
    if (k == "--dir") dir = argv[++a]; else if (k == "--stride") stride = std::stoi(argv[++a]); else if (k == "--threads") threads = std::stoi(argv[++a]);
    else { std::fprintf(stderr, "unknown option %s\n", k.c_str()); return 2; } }
  if (threads <= 0) threads = std::max(1u, std::thread::hardware_concurrency());
  FILE* f = std::fopen((dir + "/injections.bin").c_str(), "rb"); if (!f) { std::fprintf(stderr, "no injections.bin in %s\n", dir.c_str()); return 1; }
  std::vector<InjRec> recs; { InjRec r; while (std::fread(&r, sizeof r, 1, f) == 1) recs.push_back(r); } std::fclose(f);
  // group records by configuration (they are stored contiguously)
  std::vector<std::pair<size_t, size_t>> groups; for (size_t i = 0; i < recs.size();) { size_t j = i; while (j < recs.size() && recs[j].cfg == recs[i].cfg) ++j; groups.push_back({i, j}); i = j; }
  // map (prec, fm, bitclass) -> counts
  std::map<std::tuple<int, int, int>, Counts> table; std::vector<Example> examples; std::mutex mu;
  std::atomic<size_t> next{0}; size_t done = 0;
  auto worker = [&]() {
    for (;;) {
      size_t g = next++; if (g >= groups.size()) return;
      auto [b, e] = groups[g]; int cfg = recs[b].cfg; Config cf = cfg_of(cfg);
      std::map<std::tuple<int, int, int>, Counts> loc; std::vector<Example> lex;
      auto run = [&](auto tag) {
        using T = decltype(tag); CfgData<T> D; prep<T>(cfg, D);
        for (int fm = 0; fm < 3; ++fm) {
          SiteStream ss(cf.n, cf.dist, cf.seed, fm, cf.prec);
          for (size_t i = b; i < e; ++i) {
            if (recs[i].fm != fm) continue;
            Site s = ss.next();
            if (i % (size_t)stride != 0) continue;
            verify_one<T>(cfg, recs[i], s, D, loc[{cf.prec, fm, recs[i].bitclass}], lex);
          }
        }
      };
      if (cf.prec == P64) run(double{}); else run(float{});
      std::lock_guard<std::mutex> lk(mu);
      for (auto& kv : loc) table[kv.first].add(kv.second);
      for (auto& x : lex) if (examples.size() < 60) examples.push_back(x);
      ++done; if (done * 10 / groups.size() != (done - 1) * 10 / groups.size() || done == groups.size()) std::fprintf(stderr, "[verify] %3zu%%  (%zu/%zu configurations)\n", done * 100 / groups.size(), done, groups.size());
    }
  };
  std::vector<std::thread> th; for (int t = 0; t < threads; ++t) th.emplace_back(worker); for (auto& x : th) x.join();
  FILE* o = std::fopen((dir + "/verify_summary.csv").c_str(), "w");
  std::fprintf(o, "prec,fm,bitclass,verified,mism_output,mism_checks,changed,cancel,cancel_h6,cancel_h3,h6,rbmiss_h6,flip_rb,flip_fixed,flip_comb6,flip_comb3,flip_tol6,flip_tol3\n");
  Counts tot;
  for (auto& kv : table) { auto [p, fm, bc] = kv.first; const Counts& c = kv.second; tot.add(c);
    std::fprintf(o, "%s,%s,%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n", PREC_NAME[p], FM_NAME[fm], bc, c.verified, c.mism_output, c.mism_checks, c.changed, c.cancel, c.cancel_h6, c.cancel_h3, c.h6, c.rbmiss_h6, c.flip_rb, c.flip_fixed, c.flip_comb6, c.flip_comb3, c.flip_tol6, c.flip_tol3); }
  std::fclose(o);
  FILE* x = std::fopen((dir + "/verify_examples.csv").c_str(), "w");
  std::fprintf(x, "cfg,prec,n,dist,seed,fm,bit,bitclass,r_max,d_clean,d_faulty\n");
  for (auto& e : examples) { Config c = cfg_of(e.cfg); std::fprintf(x, "%d,%s,%d,%s,%d,%s,%d,%d,%.17g,%.21Lg,%.21Lg\n", e.cfg, PREC_NAME[c.prec], c.n, DIST_NAME[c.dist], c.seed, FM_NAME[e.fm], e.bit, e.bitclass, e.r_max, e.dclean, e.dfault); }
  std::fclose(x);
  std::fprintf(stderr, "[verify] verified %ld injections: %ld output mismatches, %ld check mismatches (implementation errors); %ld blind-spot cancellations; ordering flips RB %ld, fixed %ld\n",
               tot.verified, tot.mism_output, tot.mism_checks, tot.cancel, tot.flip_rb, tot.flip_fixed);
  return (tot.mism_output + tot.mism_checks) ? 3 : 0;
}
