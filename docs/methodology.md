# Methodology

This page summarises the method that the code implements. The paper gives the full derivations and discussion.

## Kernel and workload

| | |
|---|---|
| Kernel | Scalar GEMM `C = AB`, i–j–k loop, sequential accumulation, no fused multiply–add |
| Sizes | n ∈ {64, 128, 256} |
| Precisions | FP64; FP32; FP16 and BF16 inputs with FP32 multiply and accumulate |
| Inputs | 5 distributions × 3 seeds: uniform U[−1,1]; standard normal; log-normal exp(N(0,1)) with a random sign; ill-conditioned (column k of A scaled by 10^(−6k/(n−1))); non-negative U[0,1] |
| Configurations | 4 precisions × 3 sizes × 5 distributions × 3 seeds = 180 |
| Fault-free controls | 50 GEMMs per configuration, giving 2n row and column checks each; problem 0 receives the injections |
| Injections | 2,000 per configuration and fault model (1,080,000 in the full run) |

**Exact emulation of FP16 and BF16.** A product of two FP16 values (at most 22 significant bits) or two BF16 values (at most 16) is exact in FP32. The software emulation therefore reproduces mixed-precision hardware semantics exactly.

**Same inputs in every precision.** The same real-valued matrices are rounded to each precision, so precisions are compared on identical inputs.

## Observation boundary and fault models

- **Where faults strike.** Faults corrupt values in flight in the arithmetic datapath. Operands in memory stay correct, so checksum references are formed from clean `A` and `B`.
- **Fault models.** Each injection flips one bit:
  - **ACC**: the running accumulator of one output element, after a random multiply–accumulate step. One element is affected.
  - **REUSE-A**: a register-held operand `A[i,k]` used for 8 consecutive output columns. Eight elements of one row are affected.
  - **REUSE-B**: the same for `B[k,j]`. Eight elements of one column are affected.
- **Bit selection.** Fault sites come from a seeded schedule shared by all precisions. The bit is drawn uniformly over the width of the format that holds the value: the accumulation format for ACC, the input format for the reuse models.

## Checksums and checkers

For every row `i` the discrepancy is `d_i = (C'1)_i − (A(B1))_i`; columns are symmetric, giving 2n checks per output.
- **Sums** are sequential.
- **References** are formed in the accumulation precision.
- **Normalisers** `ν_i = (|A|(|B|1))_i` are formed in 80-bit extended precision.

A checker fires if any of the 2n checks satisfies its condition. All checkers are guarded: a non-finite output fires every checker.

| Checker | Fires when |
|---|---|
| Duplication | `C' ≠ C` bitwise (idealised: one corrupted replica, reliable comparator) |
| Fixed tolerance | `|d_i| > 1e-8 + 1e-8·|ref_i|` |
| Rounding bound (RB) | `|d_i| > γ_{2n} ν_i`, where `γ_k = ku/(1−ku)` |
| Rigorous bound | `|d_i| > 2γ_{2n−1}(1+u) ν_i` (Proposition 1, references in working precision) |
| Range invariant | some element non-finite or `|c'_ij| ≥ 1e6` |
| Tolerance only | `|d_i| > τ/(1+τ) · M`, where `M = max|c'_ij|` over finite elements |
| Combined | both the RB condition and the tolerance-only condition hold for the same check |

**Why τ/(1+τ) (Proposition 2).** `M` is read from the faulty output, so a large corruption can raise its own threshold. Assume the corrupted elements lie in distinct rows or in distinct columns, which holds for all three fault models.
- **With threshold τ·M,** a harmful corruption can be missed whenever τ < r ≤ τ/(1−τ).
- **With τ/(1+τ)·M,** a miss requires r ≤ τ + (1+τ)η/max|C|, where η is the rounding-level part of the discrepancy.

The analysis reports the uncorrected rule as `tol_raw` and `comb_raw`, for comparison.

## Harm and metrics

The primary harm metric is the max-norm relative error, `r = max|C' − C| / max|C|`, where `C` is the fault-free output in the same precision. A non-finite output gives `r = ∞`. For a tolerance `τ`, an injection is **harmful** if `r > τ` and **benign** otherwise.

For each checker and tolerance the analysis reports:

| Metric | Definition |
|---|---|
| Coverage | harmful injections flagged / harmful injections |
| Wasted recovery | benign injections flagged / benign injections |
| False-alarm rate | fault-free checks that fire / fault-free checks; also per fault-free GEMM |

The Frobenius and componentwise relative errors are reported as alternative harm metrics.

## Calibration protocols

Each threshold family has one scale `s`. The scale chosen is the smallest with no false alarm on the calibration set, i.e. the largest calibration statistic.

| Protocol | Calibrate on | Score on |
|---|---|---|
| In sample | all 672,000 checks (2,250 GEMMs) of the precision | all 270,000 injections |
| Leave one seed out (3 folds) | 448,000 checks (1,500 GEMMs) | 224,000 checks (750 GEMMs), 90,000 injections |
| Leave one distribution out (5 folds) | 537,600 checks (1,800 GEMMs) | 134,400 checks (450 GEMMs), 54,000 injections |

Held-out data never influences `s`.

## Statistics

- **Why outcomes are dependent.** The checks of one GEMM share its matrices, and the injections of one configuration share one product.
- **Per-GEMM false alarms.** These carry binomial intervals: fault-free GEMMs are independent, and a zero count over N gives the 95% bound 3/N.
- **Coverage, wasted recovery and per-check false alarms** carry linearised cluster intervals:
  - clusters: configurations;
  - strata: size × distribution cells (stratified variant), or none (unstratified, conservative variant).
- **Wilson intervals** are also reported, for comparison.

## Algorithms and where they live

| Algorithm | Program | Role |
|---|---|---|
| 1. Injection campaign | `src/campaign.cpp` (`bin/campaign`) | Generates inputs, injects faults, and records every checker's normalised discrepancy |
| 2. Metrics at tolerance τ | `analysis/analyze.py` | All metrics, calibration, intervals, cost model; writes `report.md`, `summary.json`, `fig/` |
| 3. Independent verification | `src/verify.cpp` (`bin/verify`) | Rebuilds sampled injections from scratch. Separates implementation errors, checksum blind spots and summation-order effects |
| Timing | `src/bench_timing.cpp` (`bin/bench_timing`) | Cost of the GEMM, duplication, and the RB and combined detectors |

`src/common.hpp` contains only input generation, which the campaign and verifier share:
- random number generator and deterministic elementary functions;
- matrix generators;
- fault-site schedule;
- FP16/BF16 rounding;
- record layouts.

The arithmetic under test is implemented separately in the campaign and in the verifier.

**Incremental recomputation.** An injection recomputes only the outputs and checksums its fault reaches. It uses exactly the same sequence of rounded operations as a full recomputation, so each injection costs O(n) operations and still matches a full rebuild bit for bit.

**Re-thresholding without new injections.** The records store normalised discrepancies, not just decisions, so any new threshold can be evaluated without rerunning the campaign.
