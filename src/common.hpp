// sdc-bench: shared definitions for the campaign (Algorithm 1) and the verifier (Algorithm 3).
// Only *input generation* lives here (RNG, matrices, fault-site draws, number formats).
// The arithmetic under test (GEMM, checksums, checks) is implemented separately in
// campaign.cpp and verify.cpp so that the verifier is an independent implementation.
//
// Build with: g++ -O2 -std=c++17 -ffp-contract=off -fno-fast-math
// (-ffp-contract=off forbids fused multiply-add, so every product and sum is rounded
// separately; the results are then bit-identical on any x86-64 machine.)
#pragma once
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <array>
#include <limits>
#include <algorithm>

namespace sdc {

// ------------------------------------------------------------------ campaign design
enum Prec : int { P64 = 0, P32 = 1, P16 = 2, PB16 = 3 };   // P16/PB16: FP16/BF16 inputs, FP32 accumulation
static const char* PREC_NAME[4] = {"fp64", "fp32", "fp16", "bf16"};
enum Dist : int { UNIF = 0, NORMAL = 1, LOGN = 2, ILLC = 3, NONNEG = 4 };
static const char* DIST_NAME[5] = {"uniform", "normal", "lognormal", "illcond", "nonneg"};
enum FM : int { ACC = 0, REUSE_A = 1, REUSE_B = 2 };
static const char* FM_NAME[3] = {"acc", "reuse_a", "reuse_b"};
static const int SIZES[3] = {64, 128, 256};
static const int NSEEDS = 3;
static const int REUSE = 8;            // outputs reached by one reused operand
static const double FIXED_ATOL = 1e-8, FIXED_RTOL = 1e-8;   // fixed-tolerance checker
static const double RANGE_BOUND = 1e6;                      // range invariant

struct Config { int prec, n, dist, seed; };
// Configuration id = ((prec*3 + size_index)*5 + dist)*3 + seed
inline int cfg_id(int prec, int si, int dist, int seed) { return ((prec * 3 + si) * 5 + dist) * NSEEDS + seed; }
inline Config cfg_of(int id) {
  Config c; c.seed = id % NSEEDS; id /= NSEEDS; c.dist = id % 5; id /= 5; int si = id % 3; c.prec = id / 3; c.n = SIZES[si]; return c;
}
inline int size_index(int n) { return n == 64 ? 0 : (n == 128 ? 1 : 2); }

// ------------------------------------------------------------------ portable RNG
// SplitMix64 for seed derivation; std-independent xoshiro256** for streams.
inline uint64_t splitmix64(uint64_t& x) {
  uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
inline uint64_t hash_seed(std::initializer_list<uint64_t> parts) {
  uint64_t s = 0x5DC0BE7C1234ABCDULL, h = 0;
  for (uint64_t p : parts) { s ^= p + 0x9E3779B97F4A7C15ULL + (s << 6) + (s >> 2); h = splitmix64(s); }
  return h;
}
struct Rng {
  uint64_t s[4];
  explicit Rng(uint64_t seed) { uint64_t x = seed; for (auto& v : s) v = splitmix64(x); }
  static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  uint64_t next() {
    uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
    return r;
  }
  double uniform01() { return (next() >> 11) * 0x1.0p-53; }            // [0,1)
  double uniform01_open() { return ((next() >> 11) + 1) * 0x1.0p-53; } // (0,1]
  uint64_t below(uint64_t m) {                                          // unbiased integer in [0,m)
    uint64_t lim = UINT64_MAX - (UINT64_MAX % m);
    for (;;) { uint64_t x = next(); if (x < lim) return x % m; }
  }
  double normal();                                                      // Marsaglia polar, defined below
};

// ------------------------------------------------------------------ deterministic elementary functions
// Built from +, -, *, / and exact operations (frexp, ldexp, floor, sqrt) only, so the
// generated inputs are bit-identical across C libraries (libm log/exp may differ in the
// last bit between glibc versions). Accuracy is about 1 ulp, ample for input generation.
inline double det_log(double x) {                       // x > 0, finite
  int e; double m = std::frexp(x, &e);                  // x = m * 2^e, m in [0.5, 1)
  if (m < 0.70710678118654752) { m *= 2.0; e -= 1; }    // m in [sqrt(1/2), sqrt(2))
  double s = (m - 1.0) / (m + 1.0), s2 = s * s, term = s, sum = 0.0;
  for (int k = 1; k <= 41; k += 2) { sum += term / k; term *= s2; }
  return 2.0 * sum + e * 0.69314718055994530942;
}
inline double det_exp(double x) {
  const double LN2_HI = 0.693147180369123816490, LN2_LO = 1.90821492927058770002e-10;
  double k = std::floor(x / 0.69314718055994530942 + 0.5);
  double r = (x - k * LN2_HI) - k * LN2_LO, term = 1.0, sum = 1.0;
  for (int i = 1; i <= 24; ++i) { term *= r / i; sum += term; }
  return std::ldexp(sum, (int)k);
}
inline double Rng::normal() {
  for (;;) {
    double v1 = 2.0 * uniform01() - 1.0, v2 = 2.0 * uniform01() - 1.0, s = v1 * v1 + v2 * v2;
    if (s > 0.0 && s < 1.0) return v1 * std::sqrt(-2.0 * det_log(s) / s);
  }
}

// ------------------------------------------------------------------ number formats
inline uint32_t f2u(float f) { uint32_t u; std::memcpy(&u, &f, 4); return u; }
inline float u2f(uint32_t u) { float f; std::memcpy(&f, &u, 4); return f; }
inline uint64_t d2u(double d) { uint64_t u; std::memcpy(&u, &d, 8); return u; }
inline double u2d(uint64_t u) { double d; std::memcpy(&d, &u, 8); return d; }

// IEEE binary16, round-to-nearest-even
inline uint16_t float_to_half(float f) {
  uint32_t x = f2u(f), sign = (x >> 16) & 0x8000u;
  uint32_t ax = x & 0x7FFFFFFFu;
  if (ax >= 0x7F800000u) return sign | (ax > 0x7F800000u ? 0x7E00u : 0x7C00u);   // NaN / Inf
  if (ax >= 0x477FF000u) return sign | 0x7C00u;                                   // overflow -> Inf
  if (ax < 0x38800000u) {                                                         // subnormal or zero
    if (ax < 0x33000000u) return sign;                                            // < 2^-25 -> 0
    uint32_t m = (ax & 0x7FFFFFu) | 0x800000u; int e = (int)(ax >> 23);
    int shift = 126 - e;                                                          // 14..24
    uint32_t r = m >> shift, rem = m & ((1u << shift) - 1), half = 1u << (shift - 1);
    if (rem > half || (rem == half && (r & 1u))) r++;
    return sign | r;
  }
  uint32_t r = ((ax >> 13) - (112u << 10)), rem = ax & 0x1FFFu;
  if (rem > 0x1000u || (rem == 0x1000u && (r & 1u))) r++;
  return sign | r;
}
inline float half_to_float(uint16_t h) {
  uint32_t sign = (uint32_t)(h & 0x8000u) << 16, e = (h >> 10) & 0x1F, m = h & 0x3FF;
  if (e == 0) { if (m == 0) return u2f(sign); float v = std::ldexp((float)m, -24); return sign ? -v : v; }
  if (e == 31) return u2f(sign | 0x7F800000u | (m << 13));
  return u2f(sign | ((e + 112u) << 23) | (m << 13));
}
// bfloat16, round-to-nearest-even (no NaN inputs occur)
inline uint16_t float_to_bf16(float f) { uint32_t x = f2u(f); x += 0x7FFFu + ((x >> 16) & 1u); return (uint16_t)(x >> 16); }
inline float bf16_to_float(uint16_t b) { return u2f((uint32_t)b << 16); }

// Width and class boundaries of the format in which a bit is flipped.
struct Fmt { int width, frac_bits, exp_bits; };
inline Fmt acc_fmt(int prec) { return prec == P64 ? Fmt{64, 52, 11} : Fmt{32, 23, 8}; }   // accumulator
inline Fmt in_fmt(int prec) {                                                               // operand register
  switch (prec) { case P64: return {64, 52, 11}; case P32: return {32, 23, 8}; case P16: return {16, 10, 5}; default: return {16, 7, 8}; }
}
// bit class: 0 = fraction, 1 = exponent, 2 = sign
inline int bit_class(const Fmt& f, int bit) { return bit < f.frac_bits ? 0 : (bit < f.frac_bits + f.exp_bits ? 1 : 2); }

// Flip one bit of a value held in the accumulator format.
inline double flip_acc(double v, int bit) { return u2d(d2u(v) ^ (1ULL << bit)); }
inline float flip_acc(float v, int bit) { return u2f(f2u(v) ^ (1u << bit)); }
// Flip one bit of an operand held in its *input* format (exact round trip for stored inputs).
inline double flip_in(double v, int bit, int) { return u2d(d2u(v) ^ (1ULL << bit)); }
inline float flip_in(float v, int bit, int prec) {
  if (prec == P32) return u2f(f2u(v) ^ (1u << bit));
  if (prec == P16) return half_to_float((uint16_t)(float_to_half(v) ^ (1u << bit)));
  return bf16_to_float((uint16_t)(float_to_bf16(v) ^ (1u << bit)));
}

// ------------------------------------------------------------------ input generation
// Real-valued draws depend on (n, dist, seed, problem) only, so every precision sees the
// same underlying matrices, rounded to its own input format.
//   uniform  : U[-1,1]
//   normal   : N(0,1)
//   lognormal: exp(N(0,1)) with a random sign (heavy-tailed)
//   illcond  : U[-1,1] with column k of A scaled by 10^(-6 k/(n-1)) (six decades)
//   nonneg   : U[0,1]
inline void gen_real(int n, int dist, int seed, int prob, std::vector<double>& A, std::vector<double>& B) {
  Rng g(hash_seed({0xA11CEULL, (uint64_t)n, (uint64_t)dist, (uint64_t)seed, (uint64_t)prob}));
  A.assign((size_t)n * n, 0.0); B.assign((size_t)n * n, 0.0);
  auto draw = [&](void) -> double {
    switch (dist) {
      case UNIF: return 2.0 * g.uniform01() - 1.0;
      case NORMAL: return g.normal();
      case LOGN: { double m = det_exp(g.normal()); return (g.next() >> 63) ? -m : m; }
      case ILLC: return 2.0 * g.uniform01() - 1.0;
      default: return g.uniform01();
    }
  };
  for (auto& v : A) v = draw();
  for (auto& v : B) v = draw();
  if (dist == ILLC)
    for (int i = 0; i < n; ++i) for (int k = 0; k < n; ++k) A[(size_t)i * n + k] *= det_exp(-6.0 * k / (n - 1) * 2.30258509299404568402);
}
inline double round_in(double v, int prec) { (void)prec; return v; }
inline float round_in_f(double v, int prec) {
  float f = (float)v;
  if (prec == P16) return half_to_float(float_to_half(f));
  if (prec == PB16) return bf16_to_float(float_to_bf16(f));
  return f;
}

// ------------------------------------------------------------------ fault-site schedule
// Site draws depend on (n, dist, seed, fm) only, so all precisions share fault sites;
// the bit is drawn from a separate, precision-specific stream.
struct Site { int i, j, k, bit; };
struct SiteStream {
  Rng site, bits; int n, fm, wacc, win;
  SiteStream(int n_, int dist, int seed, int fm_, int prec)
    : site(hash_seed({0x517EULL, (uint64_t)n_, (uint64_t)dist, (uint64_t)seed, (uint64_t)fm_})),
      bits(hash_seed({0xB175ULL, (uint64_t)n_, (uint64_t)dist, (uint64_t)seed, (uint64_t)fm_, (uint64_t)prec})),
      n(n_), fm(fm_), wacc(acc_fmt(prec).width), win(in_fmt(prec).width) {}
  Site next() {
    Site s;
    if (fm == ACC) { s.i = (int)site.below(n); s.j = (int)site.below(n); s.k = (int)site.below(n); s.bit = (int)bits.below(wacc); }
    else if (fm == REUSE_A) { s.i = (int)site.below(n); s.j = (int)site.below(n - REUSE + 1); s.k = (int)site.below(n); s.bit = (int)bits.below(win); }
    else { s.k = (int)site.below(n); s.j = (int)site.below(n); s.i = (int)site.below(n - REUSE + 1); s.bit = (int)bits.below(win); }
    return s;
  }
};

// ------------------------------------------------------------------ rounding constants
inline long double unit_roundoff(int prec) { return prec == P64 ? std::ldexp(1.0L, -53) : std::ldexp(1.0L, -24); }
inline long double gamma_k(long double k, long double u) { return k * u / (1.0L - k * u); }

// ------------------------------------------------------------------ binary record formats
#pragma pack(push, 1)
struct InjRec {               // one per injection (Algorithm 1 output)
  uint32_t cfg; uint32_t idx; uint8_t fm; uint8_t bit; uint8_t bitclass; uint8_t flags;
  double r_max, r_fro, r_cw;  // harm under max-norm, Frobenius, componentwise relative error
  double amax;                // max over all 2n checks of |d| / (gamma_2n * norm)         (rounding-bound)
  double bmax;                // max over all checks of |d| / max|C'|                         (tolerance term)
  double cmax;                // max over all checks of |d| / (atol + rtol*|ref|)             (fixed tolerance)
  double cb;                  // max of |d|/max|C'| over checks that also exceed the rounding bound   (combined)
  double cb_rig;              // same, with the rigorous bound 2*gamma_{2n-1}*(1+u)
};
struct CleanRec {             // one per fault-free check
  uint32_t cfg; uint16_t prob; uint16_t check;   // check < n: row; >= n: column
  double a, b, c;             // |d|/(gamma_2n*norm), |d|/max|C|, |d|/(atol+rtol*|ref|)
};
#pragma pack(pop)
enum Flag : uint8_t { F_BITWISE = 1, F_NAN = 2, F_INF = 4, F_RANGE = 8, F_UNG_RB = 16, F_UNG_FIXED = 32, F_UNG_COMB6 = 64, F_UNG_COMB3 = 128 };
// F_UNG_*: an *unguarded* comparison (|d| > T, IEEE semantics, no non-finite check) fired.

}  // namespace sdc
