# Reproducing the results

There are two ways to reproduce the results:
- **Re-analyse the archived records.** This takes minutes and needs no compilation.
- **Rerun the campaign from scratch.** On x86-64 Linux it produces the same records bit for bit.

## Requirements

- x86-64 Linux. Other platforms run correctly, but their records are not bit-identical to the reference.
- A C++17 compiler: GCC 13 and GCC 16 were tested.
- GNU make.
- Python 3.10 or newer, with NumPy (`pip install -r analysis/requirements.txt`).
- About 3 GB of RAM for the full analysis, and about 200 MB of disk for the full records.

## 1. Check the installation (about 20 seconds)

```bash
make test
```

This runs a reduced campaign, checks that its records match `tests/quick_reference.sha256` bit for bit, checks that the verifier finds no disagreement, and runs the analysis.

## 2a. Re-analyse the archived records

Download and unpack the records archive (see `data/README.md`), then:

```bash
cd sdc-bench-records-20261103 && sha256sum -c SHA256SUMS && cd -
python3 analysis/analyze.py sdc-bench-records-20261103
```

This writes three outputs into the records folder:
- `report.md`: every table, in readable form;
- `summary.json`: every number, with intervals;
- `fig/`: data for every figure.

## 2b. Rerun the full campaign

```bash
make full                        # or: THREADS=16 scripts/run_all.sh full
```

| Step | Program | Time on a 16-thread laptop |
|---|---|---|
| Campaign | `bin/campaign` | about 7 s |
| Verification (every 10th injection) | `bin/verify` | a few minutes |
| Timing | `bin/bench_timing` | a few minutes |
| Analysis | `analysis/analyze.py` | about 1.5 min |

The run should reproduce these digests:

```
6eaedfc4c1d7759139637dc9621e03b26ae58f71b3fd4133f66a9876b9a8b342  injections.bin
d18c4e1c10c58fec611bd98fc021b6b9dfd0af0a8d907ff6dc5cd1d1e77bf26c  clean_checks.bin
6bb8846b45d0b6d8fad7390ca166abaf4f9235b5c7df2a99f5d8470bf778f2f9  problems.csv
```

**Timing.** Measured times depend on the machine. Run the timing step on an idle machine on mains power. If possible, set the CPU frequency governor to `performance` (`sudo cpupower frequency-set -g performance`). The reported quantities are ratios between operations on the same matrices.

## 3. Archive a run

```bash
scripts/make_release.sh results
```

This writes `dist/sdc-bench-records-<date>.tar.gz`, containing:
- the records and their checksums;
- the metadata, logs, verifier output and timing;
- the data licence and the record-format description.
