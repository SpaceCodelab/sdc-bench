# Changelog

## 1.0.0 (2026-11-03)

First public release.

- **Fault-injection campaign:** four precisions (FP64, FP32, and FP16 and BF16 with FP32 accumulation), three fault models and five input distributions.
- **Independent bit-level verifier** that separates implementation errors, checksum blind spots and summation-order effects.
- **Timing benchmark** for detector and recovery cost.
- **Analysis:**
  - coverage, false alarms (per check and per GEMM) and wasted recovery as functions of the tolerance;
  - in-sample and held-out calibration;
  - dependence-aware intervals;
  - the self-normalisation ablation;
  - the cost model.
- **Tolerance-aware checkers** use the threshold τ/(1+τ)·max|C′|.
- **Determinism regression test** and continuous integration.
