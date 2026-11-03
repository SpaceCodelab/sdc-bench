# Record formats

All files are little-endian, packed (no padding). They can be loaded with NumPy using the dtypes defined at the top of `analysis/analyze.py`.

## injections.bin: one 76-byte record per injection (1,080,000 in the reference run)

| Field | Type | Meaning |
|---|---|---|
| cfg | uint32 | Configuration id. prec = cfg // 45 (0 FP64, 1 FP32, 2 FP16, 3 BF16); size index = (cfg // 15) % 3 (64, 128, 256); dist = (cfg // 3) % 5 (uniform, normal, lognormal, illcond, nonneg); seed = cfg % 3 |
| idx | uint32 | Injection index within the configuration and fault model |
| fm | uint8 | Fault model: 0 accumulator, 1 reuse-A, 2 reuse-B |
| bit | uint8 | Flipped bit position within its format |
| bitclass | uint8 | 0 fraction, 1 exponent, 2 sign |
| flags | uint8 | 1 bitwise change, 2 NaN, 4 Inf, 8 range invariant fired, 16/32/64/128 unguarded RB / fixed / combined(1e-6) / combined(1e-3) fired |
| r_max, r_fro, r_cw | float64 | Harm: max-norm, Frobenius and componentwise relative error |
| amax | float64 | max over 2n checks of \|d\| / (gamma_2n · normaliser): RB fires iff > 1 |
| bmax | float64 | max over checks of \|d\| / max\|C'\|: tolerance-only fires iff > tau/(1+tau) |
| cmax | float64 | max over checks of \|d\| / (1e-8 + 1e-8·\|ref\|): fixed fires iff > 1 |
| cb | float64 | max of \|d\|/max\|C'\| over checks that also exceed the RB threshold: combined fires iff > tau/(1+tau) |
| cb_rig | float64 | Same with the rigorous bound 2·gamma_{2n-1}·(1+u) |

## clean_checks.bin: one 32-byte record per fault-free check (2,688,000)

| Field | Type | Meaning |
|---|---|---|
| cfg | uint32 | As above |
| prob | uint16 | Fault-free problem index, 0..49; injections were made on problem 0 |
| check | uint16 | < n: row check; >= n: column check |
| a, b, c | float64 | \|d\| / (gamma_2n · normaliser), \|d\| / max\|C\|, \|d\| / (1e-8 + 1e-8·\|ref\|) |

## Checksums of the reference run

```
6eaedfc4c1d7759139637dc9621e03b26ae58f71b3fd4133f66a9876b9a8b342  injections.bin
d18c4e1c10c58fec611bd98fc021b6b9dfd0af0a8d907ff6dc5cd1d1e77bf26c  clean_checks.bin
6bb8846b45d0b6d8fad7390ca166abaf4f9235b5c7df2a99f5d8470bf778f2f9  problems.csv
```
