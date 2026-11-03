# Results of the reference run

These are the headline numbers of the reference run of 2026-10-31. That run had 1,080,000 injections and 2,688,000 fault-free checks, and its records were bit-identical on two machines. Every value here, and many more, is in `summary.json` and `report.md`, which `analysis/analyze.py` produces from the archived records. Percentages are of the stated population.

## How much corruption is harmful

| | FP64 | FP32 | FP16 | BF16 |
|---|---|---|---|---|
| Changes the output bitwise | 95.65 | 91.41 | 98.28 | 98.57 |
| Harmful at τ = 10⁻⁶ | 37.71 | 66.50 | 87.44 | 88.70 |
| Harmful at τ = 10⁻³ | 21.48 | 35.12 | 45.29 | 55.42 |

## False alarms on fault-free checks

| Checker | FP64 | FP32 | FP16 | BF16 |
|---|---|---|---|---|
| Fixed tolerance 10⁻⁸ | 0.00 | 92.85 | 92.12 | 91.10 |
| Rounding bound | 0.00 | 0.00 | 0.00 | 0.00 |
| Tolerance only, τ = 10⁻⁶ | 0.00 | 46.53 | 42.14 | 37.92 |
| Tolerance only, τ = 10⁻³ | 0.00 | 0.00 | 0.00 | 0.00 |

A zero means none of 672,000 checks per precision fired. The 95% upper bound for a zero is 4.5 × 10⁻⁴ % per check, or 0.13% per fault-free GEMM.

## Coverage / wasted recovery at τ = 10⁻³

| Checker | FP64 | FP32 | FP16 | BF16 |
|---|---|---|---|---|
| Duplication (idealised) | 100.00 / 94.46 | 100.00 / 86.76 | 100.00 / 96.86 | 100.00 / 96.80 |
| Fixed tolerance 10⁻⁸ | 100.00 / 33.95 | 100.00 / 100.00 | 100.00 / 100.00 | 100.00 / 100.00 |
| Rounding bound | 100.00 / 58.39 | 89.10 / 4.29 | 87.99 / 8.33 | 88.58 / 10.47 |
| Rigorous bound | 100.00 / 56.45 | 81.94 / 2.19 | 80.71 / 4.35 | 81.07 / 5.60 |
| Tolerance only | 100.00 / 1.30 | 100.00 / 2.89 | 100.00 / 6.06 | 100.00 / 8.09 |
| Combined | 100.00 / 1.30 | 89.10 / 1.49 | 87.98 / 3.23 | 88.58 / 4.20 |
| Range invariant | 14.68 / 0.00 | 11.85 / 0.00 | 5.17 / 0.00 | 11.56 / 0.00 |

Coverage is the share of harmful corruptions that are flagged. Wasted recovery is the share of harmless corruptions that are flagged anyway.
The tolerance-only entries of 100.00 for FP32 and FP16 are rounded: each misses one corruption, and both lie within 0.09% of the tolerance.

For held-out calibration, per-GEMM false alarms, dependence-aware intervals, the self-normalisation ablation, timing and the cost model, see the corresponding sections of `report.md`.
