// sdc-bench timing: cost of the detector separately from the cost of
// recovery, for the same scalar kernel and floating-point semantics as the campaign.
//   gemm        one scalar GEMM (i-j-k order, as in the campaign) = cost of one recovery
//   dup         two GEMMs plus a bitwise comparison of the outputs
//   abft_rb     checksum detector: references A(B1), (1^T A)B, row/column sums of C,
//               normalisers |A|(|B|1), (1^T|A|)|B|, and all 2n threshold comparisons
//   abft_comb   abft_rb plus max|C| and the tolerance comparison of the combined rule
// Output (stdout, CSV): prec,n,op,ns_median,ns_min,reps
// Run on an idle machine on mains power; each value is the median of 7 timed batches.
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>
#include <algorithm>
#include <functional>

static volatile double g_sink = 0;

template <class T> static void gemm(const std::vector<T>& A, const std::vector<T>& B, std::vector<T>& C, int n) {
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
    T acc = 0; for (int k = 0; k < n; ++k) { T p = A[(size_t)i * n + k] * B[(size_t)k * n + j]; acc = acc + p; }
    C[(size_t)i * n + j] = acc;
  }
}

template <class T> static int abft(const std::vector<T>& A, const std::vector<T>& B, const std::vector<T>& C, int n, bool comb, double tau) {
  std::vector<T> B1(n), A1(n); std::vector<double> nb(n), na(n);
  const double u = sizeof(T) == 8 ? std::ldexp(1.0, -53) : std::ldexp(1.0, -24);
  const double g = 2.0 * n * u / (1 - 2.0 * n * u);
  for (int k = 0; k < n; ++k) { T s = 0; double sa = 0; for (int j = 0; j < n; ++j) { T v = B[(size_t)k * n + j]; s = s + v; sa += std::fabs((double)v); } B1[k] = s; nb[k] = sa; }
  for (int k = 0; k < n; ++k) { A1[k] = 0; na[k] = 0; }
  for (int i = 0; i < n; ++i) for (int k = 0; k < n; ++k) { T v = A[(size_t)i * n + k]; A1[k] = A1[k] + v; na[k] += std::fabs((double)v); }
  double M = 0;
  if (comb) for (size_t t = 0; t < C.size(); ++t) { double v = std::fabs((double)C[t]); if (v > M) M = v; }
  int fires = 0;
  for (int i = 0; i < n; ++i) {
    T s = 0, rs = 0; double nr = 0;
    for (int k = 0; k < n; ++k) { T p = A[(size_t)i * n + k] * B1[k]; s = s + p; nr += std::fabs((double)A[(size_t)i * n + k]) * nb[k]; }
    for (int j = 0; j < n; ++j) rs = rs + C[(size_t)i * n + j];
    double d = std::fabs((double)(rs - s));
    bool f = !std::isfinite(d) || d > g * nr; if (comb) f = f && (!std::isfinite(d) || d > tau * M); fires += f;
  }
  std::vector<T> cs(n, 0);
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) cs[j] = cs[j] + C[(size_t)i * n + j];
  for (int j = 0; j < n; ++j) {
    T s = 0; double nc = 0;
    for (int k = 0; k < n; ++k) { T p = A1[k] * B[(size_t)k * n + j]; s = s + p; nc += na[k] * std::fabs((double)B[(size_t)k * n + j]); }
    double d = std::fabs((double)(cs[j] - s));
    bool f = !std::isfinite(d) || d > g * nc; if (comb) f = f && (!std::isfinite(d) || d > tau * M); fires += f;
  }
  return fires;
}

static double now_ns() { return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

static void measure(const char* prec, int n, const char* op, const std::function<void()>& f) {
  f();                                                   // warm-up
  double t0 = now_ns(); int reps = 0; do { f(); ++reps; } while (now_ns() - t0 < 2e8 && reps < 1000000);
  std::vector<double> batch;
  for (int b = 0; b < 7; ++b) { double s = now_ns(); for (int r = 0; r < reps; ++r) f(); batch.push_back((now_ns() - s) / reps); }
  std::sort(batch.begin(), batch.end());
  std::printf("%s,%d,%s,%.1f,%.1f,%d\n", prec, n, op, batch[3], batch[0], reps); std::fflush(stdout);
}

template <class T> static void run(const char* prec) {
  for (int n : {64, 128, 256, 512}) {
    std::vector<T> A((size_t)n * n), B((size_t)n * n), C((size_t)n * n), C2((size_t)n * n);
    uint64_t x = 88172645463325252ULL;
    auto rnd = [&]() { x ^= x << 13; x ^= x >> 7; x ^= x << 17; return (T)((x >> 11) * 0x1.0p-53 * 2.0 - 1.0); };
    for (auto& v : A) v = rnd();
    for (auto& v : B) v = rnd();
    gemm(A, B, C, n);
    measure(prec, n, "gemm", [&]() { gemm(A, B, C, n); g_sink = g_sink + (double)C[0]; });
    measure(prec, n, "dup", [&]() { gemm(A, B, C, n); gemm(A, B, C2, n); g_sink = g_sink + (std::memcmp(C.data(), C2.data(), C.size() * sizeof(T)) != 0); });
    measure(prec, n, "abft_rb", [&]() { g_sink = g_sink + abft(A, B, C, n, false, 1e-3); });
    measure(prec, n, "abft_comb", [&]() { g_sink = g_sink + abft(A, B, C, n, true, 1e-3); });
  }
}

int main() {
  std::printf("prec,n,op,ns_median,ns_min,reps\n");
  run<double>("fp64");
  run<float>("fp32");
  return 0;
}
