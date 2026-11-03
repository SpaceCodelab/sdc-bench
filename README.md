# sdc-bench

[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-informational.svg)](src)
[![Python 3](https://img.shields.io/badge/python-3.10%2B-informational.svg)](analysis)

**Tolerance-aware evaluation of silent data corruption (SDC) checking in matrix multiplication.**

sdc-bench is the experimental framework behind the paper *How Much Detection Is Enough? Tolerance-Aware Evaluation of Silent Data Corruption Checking in Numerical Kernels* (R. Singh). It injects single-bit faults into GEMM in four floating-point precisions and evaluates SDC checkers against the error an application can tolerate. The checkers are duplication, fixed-tolerance checksums, rounding-bound checksums, calibrated thresholds, and tolerance-aware checksums. For each one it reports three quantities:

- **coverage** of harmful corruptions;
- **false alarms** on fault-free runs, per check and per GEMM;
- **wasted recovery**: harmless corruptions that would still trigger recomputation.

The campaign is fully deterministic: the reference run reproduces bit for bit on different x86-64 machines and compilers.

## Highlights of the reference run

- **Bitwise change is not harm.** In FP64, 95.65% of injections change the output bitwise, but only 21.48% exceed a relative tolerance of 10⁻³.
- **Fixed tolerances do not port across precision.** A 10⁻⁸ tolerance raises no false alarm in FP64 and fires on 92.85% of fault-free FP32 checks.
- **Calibrated thresholds do not port across inputs.** Scaled on four input distributions, they raise an alarm on every fault-free GEMM of the held-out non-negative distribution.
- **Tolerance-aware checking.** At τ = 10⁻³, a checksum with threshold τ/(1+τ)·max|C′| misses 2 of 424,743 harmful corruptions across the four precisions. It raises no false alarm and wastes 1.30–8.09% of recoveries; duplication wastes 86.76–96.86%.

[`docs/results.md`](docs/results.md) has the full headline tables.

## Quick start

```bash
git clone <this repository> && cd sdc-bench
pip install -r analysis/requirements.txt
make test          # builds, runs a reduced campaign (~20 s), checks bit-identical output
make full          # full campaign: 1,080,000 injections, 2,688,000 fault-free checks
```

The results land in `results/`:

| File | Contents |
|---|---|
| `report.md` | Every table, in readable form |
| `summary.json` | Every number, with intervals |
| `fig/` | Data for every figure |

To re-analyse the archived records without rerunning anything, see [`docs/reproducing.md`](docs/reproducing.md).

## Repository layout

```
sdc-bench/
├── src/
│   ├── common.hpp          shared input generation: RNG, matrices, fault schedule, FP16/BF16, record layouts
│   ├── campaign.cpp        Algorithm 1: fault-injection campaign
│   ├── verify.cpp          Algorithm 3: independent bit-level verifier
│   └── bench_timing.cpp    detector and recovery cost
├── analysis/
│   ├── analyze.py          Algorithm 2: metrics, calibration, statistics, cost model, report
│   └── requirements.txt
├── scripts/
│   ├── run_all.sh          full pipeline: build → campaign → verify → timing → analysis
│   └── make_release.sh     packages the records of a run for archiving
├── tests/
│   ├── check_quick.sh      determinism regression test (make test)
│   └── quick_reference.sha256
├── docs/
│   ├── methodology.md      what is measured and how
│   ├── reproducing.md      step-by-step reproduction
│   ├── data-format.md      binary record layout
│   └── results.md          headline results of the reference run
├── data/README.md          where the archived records are and what they contain
├── .github/workflows/ci.yml
├── CITATION.cff  LICENSE  LICENSE-DATA  CHANGELOG.md  CONTRIBUTING.md
└── Makefile
```

## How it works

| Stage | Program | Input | Output |
|---|---|---|---|
| M1 Fault-injection campaign | `bin/campaign` | Matrices, precisions, fault models | Injection and fault-free check records |
| M2 Harm and checker evaluation | `bin/campaign`, `bin/verify` | Records | Harmful/benign label and every checker's decision, verified bit for bit |
| M3 Metrics, calibration, statistics | `analysis/analyze.py` | Records | Coverage, false alarms and wasted recovery vs. τ; held-out calibration; intervals; cost |

[`docs/methodology.md`](docs/methodology.md) covers the checker definitions, the fault models, the calibration protocols and the statistical treatment.

## Data

The records of the reference run (about 170 MB) are archived separately; [`data/README.md`](data/README.md) explains where. Records produced by a fresh run should match these SHA-256 digests:

```
6eaedfc4c1d7759139637dc9621e03b26ae58f71b3fd4133f66a9876b9a8b342  injections.bin
d18c4e1c10c58fec611bd98fc021b6b9dfd0af0a8d907ff6dc5cd1d1e77bf26c  clean_checks.bin
```

## Citation

If you use sdc-bench or its records, please cite the software using [`CITATION.cff`](CITATION.cff); GitHub's "Cite this repository" button reads it.

## License

- **Code:** [MIT](LICENSE).
- **Experimental records:** [CC BY 4.0](LICENSE-DATA).

## Contact

Ramanjeet Singh, Department of Computer Science and Engineering, Chandigarh University, Gharuan, Mohali, India. Email: ramanjeetsingh310@gmail.com
