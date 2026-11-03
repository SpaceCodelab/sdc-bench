# Data

Git does not hold the experimental records. They are about 170 MB of binary data, archived on Zenodo: https://doi.org/10.5281/zenodo.23911926.

The archive `sdc-bench-records-20261103.tar.gz` contains:

| File | Contents |
|---|---|
| `injections.bin` | 1,080,000 injection records, 76 bytes each |
| `clean_checks.bin` | 2,688,000 fault-free check records, 32 bytes each |
| `problems.csv` | One row per generated problem: output scale, fault-free forward error, normaliser ratio |
| `campaign_meta.json` | Campaign parameters |
| `verify_summary.csv`, `verify_examples.csv` | Output of the independent verifier |
| `timing.csv` | Timing medians (ns) |
| `system.txt`, `logs/` | Machine description and run logs |
| `SHA256SUMS` | Digests of the records |
| `LICENSE`, `DATA_FORMAT.md` | CC BY 4.0 licence, record layout |

The record layout is described in [`../docs/data-format.md`](../docs/data-format.md). To re-analyse the records, see [`../docs/reproducing.md`](../docs/reproducing.md).
