#!/usr/bin/env python3
"""sdc-bench analysis (Algorithm 2): turns campaign records into every reported metric.

Usage: python3 analysis/analyze.py <results-dir>
Reads  : injections.bin, clean_checks.bin, problems.csv, campaign_meta.json,
         verify_summary.csv / verify_examples.csv (optional), timing.csv (optional)
Writes : report.md (human-readable), summary.json (all numbers), and figure data CSVs in fig/

Checker definitions (all comparisons guarded: any non-finite output fires every checker):
  dup      : output differs bitwise from the fault-free output (idealised duplication)
  fixed    : some check has |d| > 1e-8 + 1e-8*|ref|
  rb       : some check has |d| > gamma_2n * (|A|(|B|1))_i                 (classical threshold)
  rb_rig   : some check has |d| > 2*gamma_{2n-1}*(1+u) * (|A|(|B|1))_i     (rigorous bound)
  range    : some output element is non-finite or has magnitude >= 1e6
  tol(t)   : some check has |d| > t/(1+t) * max|C'|
  comb(t)  : some check has |d| above BOTH the rb threshold and t/(1+t) * max|C'|
  tol_raw, comb_raw: the same with t * max|C'| (ablation). With max|C'| taken from the faulty
             output, the raw rule can be masked by the corruption it should detect: it misses a
             harmful corruption only if t < r <= t/(1-t). The factor 1/(1+t) closes that window
             for every fault model whose corrupted elements lie in distinct rows or columns
             (see docs/methodology.md, Proposition 2).
  max|C'| is taken over the finite elements of the faulty output.
Harm: r = max|C'-C| / max|C| (primary); also Frobenius and componentwise relative error.
"""
import sys, os, json, csv, math
import numpy as np

D = sys.argv[1] if len(sys.argv) > 1 else "results"
FIG = os.path.join(D, "fig"); os.makedirs(FIG, exist_ok=True)
PREC = ["fp64", "fp32", "fp16", "bf16"]; DIST = ["uniform", "normal", "lognormal", "illcond", "nonneg"]
FMN = ["acc", "reuse_a", "reuse_b"]; BC = ["fraction", "exponent", "sign"]; SIZES = np.array([64, 128, 256])
U = np.array([2.0**-53, 2.0**-24, 2.0**-24, 2.0**-24])
F_BITWISE, F_NAN, F_INF, F_RANGE, F_UNG_RB, F_UNG_FIXED, F_UNG_C6, F_UNG_C3 = 1, 2, 4, 8, 16, 32, 64, 128

INJ = np.dtype([("cfg", "<u4"), ("idx", "<u4"), ("fm", "u1"), ("bit", "u1"), ("bitclass", "u1"), ("flags", "u1"),
                ("r_max", "<f8"), ("r_fro", "<f8"), ("r_cw", "<f8"), ("amax", "<f8"), ("bmax", "<f8"), ("cmax", "<f8"),
                ("cb", "<f8"), ("cb_rig", "<f8")])
CLN = np.dtype([("cfg", "<u4"), ("prob", "<u2"), ("check", "<u2"), ("a", "<f8"), ("b", "<f8"), ("c", "<f8")])
assert INJ.itemsize == 76 and CLN.itemsize == 32

def fields(cfg):
    cfg = cfg.astype(np.int64)
    return dict(seed=cfg % 3, dist=(cfg // 3) % 5, n=SIZES[(cfg // 15) % 3], prec=cfg // 45)

def gamma(k, u): return k * u / (1 - k * u)

def wilson(k, n, z=1.959963984540054):
    if n == 0: return [None, None, None]
    p = k / n; den = 1 + z * z / n; c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [100 * p, 100 * max(0.0, c - h), 100 * min(1.0, c + h)]

def rate(k, n):
    """Proportion with Wilson interval; for a zero count also the 95% upper bound 3/N (rule of three)."""
    w = wilson(int(k), int(n))
    out = {"k": int(k), "n": int(n), "pct": w[0], "ci": [w[1], w[2]]}
    if k == 0 and n > 0: out["upper95_rule_of_three_pct"] = 100 * 3.0 / n
    return out

def pct(x): return None if x is None else round(x, 2)
def fmt(r):
    if r is None or r.get("pct") is None: return "--"
    s = f"{r['pct']:.2f}"
    if r["k"] == 0 and "upper95_rule_of_three_pct" in r: s += f" (≤{r['upper95_rule_of_three_pct']:.2g})"
    return s

# --------------------------------------------------------------------------- load
inj = np.fromfile(os.path.join(D, "injections.bin"), dtype=INJ)
cln = np.fromfile(os.path.join(D, "clean_checks.bin"), dtype=CLN)
meta = json.load(open(os.path.join(D, "campaign_meta.json")))
fi, fc = fields(inj["cfg"]), fields(cln["cfg"])
nonfin = (inj["flags"] & (F_NAN | F_INF)) > 0
bitwise = (inj["flags"] & F_BITWISE) > 0
precs_present = sorted(set(fi["prec"].tolist()))
u_i, u_c = U[fi["prec"]], U[fc["prec"]]
g2n_i, g2n_c = gamma(2.0 * fi["n"], u_i), gamma(2.0 * fc["n"], u_c)
kr_i = 2 * gamma(2.0 * fi["n"] - 1, u_i) * (1 + u_i) / g2n_i
kr_c = 2 * gamma(2.0 * fc["n"] - 1, u_c) * (1 + u_c) / g2n_c

def teff(tau): return tau / (1.0 + tau)

def fire(name, tau=None, sel=None):
    """Firing decision of a checker on injections."""
    if name == "dup": f = bitwise.copy()
    elif name == "fixed": f = nonfin | (inj["cmax"] > 1)
    elif name == "rb": f = nonfin | (inj["amax"] > 1)
    elif name == "rb_rig": f = nonfin | (inj["amax"] > kr_i)
    elif name == "range": f = (inj["flags"] & F_RANGE) > 0
    elif name == "tol": f = nonfin | (inj["bmax"] > teff(tau))
    elif name == "comb": f = nonfin | (inj["cb"] > teff(tau))
    elif name == "comb_rig": f = nonfin | (inj["cb_rig"] > teff(tau))
    elif name == "tol_raw": f = nonfin | (inj["bmax"] > tau)
    elif name == "comb_raw": f = nonfin | (inj["cb"] > tau)
    else: raise ValueError(name)
    return f if sel is None else f[sel]

def clean_fire(name, tau=None):
    """Firing decision of a checker on fault-free checks."""
    if name == "fixed": return cln["c"] > 1
    if name == "rb": return cln["a"] > 1
    if name == "rb_rig": return cln["a"] > kr_c
    if name == "tol": return cln["b"] > teff(tau)
    if name == "comb": return (cln["a"] > 1) & (cln["b"] > teff(tau))
    if name == "comb_rig": return (cln["a"] > kr_c) & (cln["b"] > teff(tau))
    if name == "tol_raw": return cln["b"] > tau
    if name == "comb_raw": return (cln["a"] > 1) & (cln["b"] > tau)
    if name == "dup": return np.zeros(len(cln), bool)
    raise ValueError(name)

def harmful(tau, metric="r_max"): return inj[metric] > tau

def cov_wr(f, tau, sel, metric="r_max"):
    h = harmful(tau, metric)[sel]; ff = f[sel] if len(f) == len(inj) else f
    return rate(np.sum(ff & h), np.sum(h)), rate(np.sum(ff & ~h), np.sum(~h))

def fa(name, sel, tau=None):
    f = clean_fire(name, tau)[sel]; return rate(np.sum(f), np.sum(sel))

S = {"meta": meta, "definitions": __doc__}
ALL = np.ones(len(inj), bool); CALL = np.ones(len(cln), bool)
P = {p: fi["prec"] == p for p in precs_present}; PC = {p: fc["prec"] == p for p in precs_present}
P0 = cln["prob"] == 0     # the fault-free GEMM on which injections were made (one per configuration)
CHK = ["dup", "fixed", "rb", "rb_rig", "range", "tol", "comb", "comb_rig"]

# --------------------------------------------------------------------------- 1. headline per precision
head = {}
for p in precs_present:
    sp, sc = P[p], PC[p]; d = {"injections": int(sp.sum()), "fault_free_checks": int(sc.sum()), "fault_free_checks_prob0": int((sc & P0).sum())}
    d["bitwise"] = rate(np.sum(bitwise & sp), sp.sum()); d["nonfinite"] = int(np.sum(nonfin & sp))
    for t in (1e-6, 1e-3, 1e-1): d[f"harmful_{t:g}"] = rate(np.sum(harmful(t) & sp), sp.sum())
    for t in (1e-6, 1e-3):
        for ck in CHK:
            c, w = cov_wr(fire(ck, t), t, sp); d[f"{ck}_cov_{t:g}"] = c; d[f"{ck}_wr_{t:g}"] = w
    for ck in ["fixed", "rb", "rb_rig"]:
        d[f"{ck}_fa"] = fa(ck, sc); d[f"{ck}_fa_prob0"] = fa(ck, sc & P0)
    for t in (1e-6, 1e-3):
        for ck in ["tol", "comb", "comb_rig"]: d[f"{ck}_fa_{t:g}"] = fa(ck, sc, t)
    d["fixed_clean_c_p999"] = float(np.quantile(cln["c"][sc & P0], 0.999))
    head[PREC[p]] = d
S["headline"] = head

# --------------------------------------------------------------------------- 2. breakdowns
def breakdown(key, labels, sel_inj, sel_cln=None):
    out = {}
    for p in precs_present:
        for v, lab in enumerate(labels):
            si = P[p] & (sel_inj == v)
            if si.sum() == 0: continue
            d = {"injections": int(si.sum()), "bitwise": rate(np.sum(bitwise & si), si.sum())}
            for t in (1e-6, 1e-3): d[f"harmful_{t:g}"] = rate(np.sum(harmful(t) & si), si.sum())
            d["rb_cov_1e-06"], d["rb_wr_1e-06"] = cov_wr(fire("rb"), 1e-6, si)
            d["dup_cov_0.001"], d["dup_wr_0.001"] = cov_wr(fire("dup"), 1e-3, si)
            d["comb_cov_0.001"], d["comb_wr_0.001"] = cov_wr(fire("comb", 1e-3), 1e-3, si)
            d["tol_cov_0.001"], d["tol_wr_0.001"] = cov_wr(fire("tol", 1e-3), 1e-3, si)
            if sel_cln is not None:
                sc = PC[p] & (sel_cln == v)
                d["fixed_fa"] = fa("fixed", sc); d["rb_fa"] = fa("rb", sc); d["fixed_fa_prob0"] = fa("fixed", sc & P0)
                d["comb_fa_0.001"] = fa("comb", sc, 1e-3); d["tol_fa_0.001"] = fa("tol", sc, 1e-3)
            out[f"{PREC[p]}|{lab}"] = d
    return out
S["by_fault_model"] = breakdown("fm", FMN, inj["fm"])
S["by_distribution"] = breakdown("dist", DIST, fi["dist"], fc["dist"])
S["by_size"] = breakdown("n", ["64", "128", "256"], (fi["n"] == 128) * 1 + (fi["n"] == 256) * 2, (fc["n"] == 128) * 1 + (fc["n"] == 256) * 2)
S["by_bitclass"] = breakdown("bitclass", BC, inj["bitclass"])
S["by_seed"] = breakdown("seed", ["0", "1", "2"], fi["seed"], fc["seed"])

# --------------------------------------------------------------------------- 3. threshold families and calibration
def fam_stats(fam):
    """(injection statistic, clean statistic, n arrays): fires when statistic > s."""
    if fam == "tolrel": return inj["cmax"], cln["c"]
    shape = {"lin": lambda n: 2.0 * n, "sqrt": lambda n: np.sqrt(n), "const": lambda n: np.ones_like(n, float), "invsqrt": lambda n: 1 / np.sqrt(n)}[fam]
    si = inj["amax"] * g2n_i / u_i / shape(fi["n"]); sc = cln["a"] * g2n_c / u_c / shape(fc["n"])
    return si, sc
FAMS = ["lin", "sqrt", "const", "invsqrt", "tolrel"]

def eval_at(s_inj_thr, stat_i, sel_i, taus=(1e-6, 1e-3)):
    f = nonfin | (stat_i > s_inj_thr); out = {}
    for t in taus: out[f"cov_{t:g}"], out[f"wr_{t:g}"] = cov_wr(f, t, sel_i)
    return out

cal = {}
for p in precs_present:
    for fam in FAMS + ["pern"]:
        base = "lin" if fam == "pern" else fam
        st_i, st_c = fam_stats(base); res = {}
        # in-sample (all 50 problems per configuration, and problem 0 only)
        for proto, csel in (("insample_all", PC[p]), ("insample_prob0", PC[p] & P0)):
            if fam == "pern":
                thr_i = np.full(len(inj), np.inf); k = 0; nchk = 0
                for n in SIZES:
                    s_star = st_c[csel & (fc["n"] == n)].max(); thr_i[fi["n"] == n] = s_star
                r = eval_at(thr_i, st_i, P[p])
            else:
                s_star = st_c[csel].max(); r = eval_at(s_star, st_i, P[p]); r["s_star"] = float(s_star)
            res[proto] = r
        # held-out: leave one seed out, leave one distribution out
        for proto, gi, gc, vals in (("leave_one_seed_out", fi["seed"], fc["seed"], range(3)), ("leave_one_dist_out", fi["dist"], fc["dist"], range(5))):
            fa_k = fa_n = 0; hits = {t: [0, 0, 0, 0] for t in (1e-6, 1e-3)}; folds = []
            for v in vals:
                calc = PC[p] & (gc != v); held_c = PC[p] & (gc == v); held_i = P[p] & (gi == v)
                if fam == "pern":
                    thr_c = np.zeros(len(cln)); thr_i = np.zeros(len(inj))
                    for n in SIZES:
                        s_star = st_c[calc & (fc["n"] == n)].max(); thr_c[fc["n"] == n] = s_star; thr_i[fi["n"] == n] = s_star
                    fa_hits = np.sum(st_c[held_c] > thr_c[held_c]); f = nonfin | (st_i > thr_i)
                else:
                    s_star = st_c[calc].max(); fa_hits = np.sum(st_c[held_c] > s_star); f = nonfin | (st_i > s_star)
                fa_k += fa_hits; fa_n += held_c.sum()
                fold = {"held_out": int(v), "fa": rate(fa_hits, held_c.sum())}
                for t in (1e-6, 1e-3):
                    h = harmful(t) & held_i; b = ~harmful(t) & held_i
                    hits[t][0] += np.sum(f & h); hits[t][1] += h.sum(); hits[t][2] += np.sum(f & b); hits[t][3] += b.sum()
                    fold[f"cov_{t:g}"] = rate(np.sum(f & h), h.sum())
                folds.append(fold)
            res[proto] = {"fa": rate(fa_k, fa_n), "folds": folds}
            for t in (1e-6, 1e-3): res[proto][f"cov_{t:g}"] = rate(hits[t][0], hits[t][1]); res[proto][f"wr_{t:g}"] = rate(hits[t][2], hits[t][3])
        cal[f"{PREC[p]}|{fam}"] = res
S["calibration"] = cal

# frontier curves (figure data): false-alarm rate vs coverage at tau = 1e-6, per family
with open(os.path.join(FIG, "frontier.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["prec", "family", "s", "fa", "cov_1e-06", "cov_0.001"])
    for p in precs_present:
        for fam in FAMS:
            st_i, st_c = fam_stats(fam); c_sorted = np.sort(st_c[PC[p]]); nC = len(c_sorted)
            smax = c_sorted[-1]; grid = np.unique(np.concatenate([np.geomspace(smax * 1e-4, smax * 1e2, 181), [smax]]))
            for t_col in [None]:
                pass
            hs = {t: np.sort(st_i[P[p] & harmful(t) & ~nonfin]) for t in (1e-6, 1e-3)}
            hn = {t: np.sum(P[p] & harmful(t)) for t in (1e-6, 1e-3)}; nf = {t: np.sum(P[p] & harmful(t) & nonfin) for t in (1e-6, 1e-3)}
            for s in grid:
                fa_v = (nC - np.searchsorted(c_sorted, s, side="right")) / nC
                cv = [(nf[t] + len(hs[t]) - np.searchsorted(hs[t], s, side="right")) / hn[t] for t in (1e-6, 1e-3)]
                w.writerow([PREC[p], fam, f"{s:.6g}", f"{fa_v:.6g}", f"{100*cv[0]:.4f}", f"{100*cv[1]:.4f}"])

# --------------------------------------------------------------------------- 4. curves over tau (figure data)
TAUS = np.logspace(-16, 0, 65)
CHK_CURVE = CHK + ["tol_raw", "comb_raw"]
with open(os.path.join(FIG, "curves_tau.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["prec", "tau", "bitwise", "harmful"] + [f"cov_{c}" for c in CHK_CURVE] + [f"wr_{c}" for c in CHK_CURVE] + ["fa_tol", "fa_comb"])
    for p in precs_present:
        sp = P[p]; bw = 100 * np.mean(bitwise[sp])
        for t in TAUS:
            h = harmful(t)[sp]; row = [PREC[p], f"{t:.3g}", f"{bw:.4f}", f"{100*h.mean():.4f}"]
            covs, wrs = [], []
            for ck in CHK_CURVE:
                f = fire(ck, t)[sp]
                covs.append(f"{100*np.sum(f&h)/max(h.sum(),1):.4f}"); wrs.append(f"{100*np.sum(f&~h)/max((~h).sum(),1):.4f}")
            row += covs + wrs + [f"{100*np.mean(clean_fire('tol', t)[PC[p]]):.6f}", f"{100*np.mean(clean_fire('comb', t)[PC[p]]):.6f}"]
            w.writerow(row)

# --------------------------------------------------------------------------- 5. harm-metric sensitivity
sens = {}
for p in precs_present:
    sp = P[p]
    for t in (1e-6, 1e-3):
        hm = harmful(t, "r_max")[sp]; d = {}
        for m in ("r_max", "r_fro", "r_cw"):
            hx = harmful(t, m)[sp]
            d[m] = {"harmful": rate(hx.sum(), sp.sum()), "agree_with_max": rate(np.sum(hx == hm), sp.sum())}
            for ck in ("rb", "comb", "dup"):
                c, wv = cov_wr(fire(ck, t), t, sp, m); d[m][f"{ck}_cov"] = c; d[m][f"{ck}_wr"] = wv
        sens[f"{PREC[p]}|{t:g}"] = d
S["harm_metric_sensitivity"] = sens

# --------------------------------------------------------------------------- 6. non-finite outputs
nf = {}
for p in precs_present:
    sp = P[p] & nonfin; k = int(sp.sum()); fl = inj["flags"][sp]
    nf[PREC[p]] = {"count": k, "nan_only": int(np.sum((fl & F_NAN > 0) & (fl & F_INF == 0))), "inf_only": int(np.sum((fl & F_INF > 0) & (fl & F_NAN == 0))),
                   "both": int(np.sum((fl & F_NAN > 0) & (fl & F_INF > 0))),
                   "unguarded_rb_caught": rate(np.sum(fl & F_UNG_RB > 0), k), "unguarded_fixed_caught": rate(np.sum(fl & F_UNG_FIXED > 0), k),
                   "unguarded_comb_1e-6_caught": rate(np.sum(fl & F_UNG_C6 > 0), k), "unguarded_comb_1e-3_caught": rate(np.sum(fl & F_UNG_C3 > 0), k),
                   "guarded_caught": rate(k, k)}
S["nonfinite"] = nf

# --------------------------------------------------------------------------- 7. problems: scale ratio, forward error, distribution hypothesis
probs = list(csv.DictReader(open(os.path.join(D, "problems.csv"))))
cells = {}
for r in probs:
    key = (r["prec"], int(r["n"]), r["dist"]); c = cells.setdefault(key, {"ratio": [], "fwd": [], "big": 0, "zero": 0})
    c["ratio"].append(float(r["ratio_med"])); c["big"] += int(r["n_big"])
    if r["fwd_err"] != "nan": c["fwd"].append(float(r["fwd_err"]))
    if float(r["maxC"]) == 0.0: c["zero"] += 1
hyp = {}
def spearman(x, y):
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y)); return float(np.corrcoef(rx, ry)[0, 1])
for p in precs_present:
    rows = []
    for n in SIZES:
        for di, dn in enumerate(DIST):
            c = cells.get((PREC[p], int(n), dn)); si = P[p] & (fi["n"] == n) & (fi["dist"] == di)
            cv, _ = cov_wr(fire("rb"), 1e-6, si)
            rows.append({"n": int(n), "dist": dn, "norm_over_maxC_median": float(np.median(c["ratio"])), "fwd_err_median": float(np.median(c["fwd"])) if c["fwd"] else None,
                         "rb_cov_1e-06": cv, "problems_with_values_ge_1e6": c["big"], "problems_with_zero_output": c["zero"]})
    x = [r["norm_over_maxC_median"] / r["n"] for r in rows]; y = [r["rb_cov_1e-06"]["pct"] or 0 for r in rows]
    hyp[PREC[p]] = {"cells": rows, "spearman_ratio_per_n_vs_cov": spearman(np.array(x), np.array(y))}
S["distribution_hypothesis"] = hyp
# range invariant on fault-free GEMMs: fires if any element has magnitude >= 1e6
rng_fa = {}
for p in precs_present:
    rs = [r for r in probs if r["prec"] == PREC[p]]
    rng_fa[PREC[p]] = rate(sum(int(r["n_big"]) > 0 for r in rs), len(rs))
    for dn in DIST:
        rd = [r for r in rs if r["dist"] == dn]; rng_fa[f"{PREC[p]}|{dn}"] = rate(sum(int(r["n_big"]) > 0 for r in rd), len(rd))
S["range_invariant_fa_per_gemm"] = rng_fa

# --------------------------------------------------------------------------- 8. verification (Algorithm 3)
vpath = os.path.join(D, "verify_summary.csv")
if os.path.exists(vpath):
    V = list(csv.DictReader(open(vpath))); agg = {}
    for r in V:
        for key in ("ALL", r["prec"], f"{r['prec']}|{r['fm']}", f"{r['prec']}|bitclass{r['bitclass']}"):
            a = agg.setdefault(key, {}); [a.__setitem__(k, a.get(k, 0) + int(v)) for k, v in r.items() if k not in ("prec", "fm", "bitclass")]
    S["verification"] = agg

# --------------------------------------------------------------------------- 9. cost model (needs timing.csv from bench_timing)
tpath = os.path.join(D, "timing.csv")
if os.path.exists(tpath):
    T = list(csv.DictReader(open(tpath))); tm = {}
    for r in T: tm[(r["prec"], int(r["n"]), r["op"])] = float(r["ns_median"])
    cost = {}
    for p in ("fp64", "fp32"):
        pi = PREC.index(p)
        if pi not in precs_present: continue
        for n in SIZES:
            g = tm.get((p, int(n), "gemm"));
            if g is None: continue
            det = {"dup": tm[(p, int(n), "dup")] - g, "rb": tm[(p, int(n), "abft_rb")], "comb": tm[(p, int(n), "abft_comb")]}
            si = P[pi] & (fi["n"] == n); sc = PC[pi] & (fc["n"] == n)
            # per-GEMM false-alarm probability: any check of a fault-free GEMM fires
            grp = cln["cfg"][sc].astype(np.int64) * 100 + cln["prob"][sc]
            def gemm_fa(fl):
                ids = np.unique(grp); hit = np.unique(grp[fl[sc]]); return len(hit) / len(ids)
            t = 1e-3; h = harmful(t)[si]; Ph = h.mean(); entry = {"t_gemm_ns": g}
            for ck in ("dup", "rb", "comb"):
                f = fire(ck, t)[si]; covv = np.sum(f & h) / max(h.sum(), 1); wrv = np.sum(f & ~h) / max((~h).sum(), 1)
                fag = 0.0 if ck == "dup" else gemm_fa(clean_fire(ck, t))
                c_rec = g + det[ck]
                e = {"detector_overhead_pct": 100 * det[ck] / g, "gemm_false_alarm_prob": fag}
                for pf in (1e-9, 1e-6, 1e-3, 1e-1):
                    fire_p = (1 - pf) * fag + pf * (Ph * covv + (1 - Ph) * wrv); wasted_p = pf * (1 - Ph) * wrv + (1 - pf) * fag
                    e[f"p={pf:g}"] = {"expected_cost_rel_gemm": (g + det[ck] + fire_p * c_rec) / g, "wasted_recovery_cost_rel_gemm": wasted_p * c_rec / g}
                entry[ck] = e
            cost[f"{p}|{int(n)}"] = entry
    S["cost_model_tau_1e-3"] = cost

# --------------------------------------------------------------------------- 11. self-normalisation: raw tau*max|C'| vs corrected tau/(1+tau)*max|C'|
selfn = {}
for p in precs_present:
    sp = P[p]
    for t in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.3, 0.5, 1.0):
        h = harmful(t) & sp; b = ~harmful(t) & sp
        win = inj["r_max"] <= (t / (1 - t) if t < 1 else np.inf)
        d = {"harmful": int(h.sum())}
        for ck in ("tol_raw", "tol", "comb_raw", "comb"):
            f = fire(ck, t); miss = h & ~f
            d[ck] = {"cov": rate(np.sum(f & h), h.sum()), "wr": rate(np.sum(f & b), b.sum()), "miss": int(miss.sum()),
                     "miss_in_window": int(np.sum(miss & win)), "fa": fa(ck, PC[p], t)}
        selfn[f"{PREC[p]}|{t:g}"] = d
G = np.logspace(-3, 0, 61); gm = {"grid_points": len(G), "tau_min": 1e-3, "tau_max": 1.0}
for ck, stat in (("tol_raw", None), ("tol", None)):
    tot = 0; worst = 0.0; mincov = 100.0
    for p in precs_present:
        for t in G:
            h = harmful(t) & P[p]; f = fire(ck, t); miss = h & ~f; tot += int(miss.sum())
            mincov = min(mincov, 100 * np.sum(f & h) / h.sum())
            if miss.any(): worst = max(worst, float((inj["r_max"][miss] / t).max()))
    gm[ck] = {"misses_summed_over_grid_and_precisions": tot, "max_r_over_tau_among_misses": worst, "min_coverage_pct": mincov}
selfn["grid_1e-3_to_1"] = gm
m3 = harmful(1e-3) & ~fire("tol", 1e-3)
selfn["tol_misses_at_1e-3"] = [{"prec": PREC[int(fi["prec"][i])], "r_over_tau": float(inj["r_max"][i] / 1e-3)} for i in np.flatnonzero(m3)]
S["self_normalisation"] = selfn

# --------------------------------------------------------------------------- 12. false alarms with the fault-free GEMM as the unit
gid_c = cln["cfg"].astype(np.int64) * 100 + cln["prob"]
def gemm_level(fl, sel):
    """Number of fault-free GEMMs (in sel) with at least one firing check, and number of GEMMs."""
    ids = np.unique(gid_c[sel]); hit = np.unique(gid_c[sel & fl]); return len(hit), len(ids)
gfa = {}
for p in precs_present:
    d = {}
    for ck, t in (("fixed", None), ("rb", None), ("rb_rig", None), ("tol", 1e-6), ("tol", 1e-4), ("tol", 1e-3), ("comb", 1e-6), ("comb", 1e-3), ("comb_rig", 1e-3)):
        k, n = gemm_level(clean_fire(ck, t), PC[p]); d[ck if t is None else f"{ck}_{t:g}"] = rate(k, n)
    gfa[PREC[p]] = d
S["gemm_level_false_alarms"] = gfa

# held-out calibration, false alarms counted per fault-free GEMM
hgfa = {}
for p in precs_present:
    for fam in ("lin", "sqrt", "const", "invsqrt", "tolrel", "pern"):
        st_i, st_c = fam_stats("lin" if fam == "pern" else fam); res = {}
        for proto, gc, vals in (("leave_one_seed_out", fc["seed"], range(3)), ("leave_one_dist_out", fc["dist"], range(5))):
            k_tot = n_tot = 0; folds = []
            for v in vals:
                calc = PC[p] & (gc != v); held = PC[p] & (gc == v)
                if fam == "pern":
                    thr = np.zeros(len(cln))
                    for n_ in SIZES: thr[fc["n"] == n_] = st_c[calc & (fc["n"] == n_)].max()
                    fl = st_c > thr
                else:
                    fl = st_c > st_c[calc].max()
                k, n = gemm_level(fl, held); k_tot += k; n_tot += n; folds.append(rate(k, n))
            res[proto] = {"fa_gemm": rate(k_tot, n_tot), "folds": folds}
        hgfa[f"{PREC[p]}|{fam}"] = res
S["heldout_gemm_level_false_alarms"] = hgfa

# --------------------------------------------------------------------------- 13. dependence-aware intervals
# Injections of one configuration share the same matrices, and the 2n checks of one GEMM share them too, so
# outcomes are not independent. We treat each configuration (precision, size, distribution, seed) as a cluster
# and each size x distribution cell as a stratum, and use the Taylor-linearised variance of a ratio estimator:
#   z_c = (num_c - p*den_c) / sum(den);  var = sum_h k_h/(k_h-1) * sum_{c in h} (z_c - mean_h z)^2.
# "stratified": strata = 15 size x distribution cells, 3 seeds each (the design is fixed; matrices are random).
# "unstratified": all 45 configurations of a precision as one stratum (treats the design cells as random too).
def lin_ci(num_mask, den_mask, cfg_arr, sel):
    cfgs = np.unique(cfg_arr[sel]); idx = np.searchsorted(cfgs, cfg_arr)
    num = np.bincount(idx[sel & num_mask], minlength=len(cfgs)).astype(float)
    den = np.bincount(idx[sel & den_mask], minlength=len(cfgs)).astype(float)
    pt = num.sum() / den.sum(); z = (num - pt * den) / den.sum(); strata = (cfgs.astype(np.int64) // 3) % 15
    out = {"k": int(num.sum()), "n": int(den.sum()), "pct": 100 * pt}
    w = wilson(int(num.sum()), int(den.sum())); out["wilson_ci"] = [w[1], w[2]]
    for name, st in (("stratified", strata), ("unstratified", np.zeros_like(strata))):
        v = 0.0
        for h in np.unique(st):
            zh = z[st == h]; k = len(zh); v += k / (k - 1) * np.sum((zh - zh.mean()) ** 2)
        se = math.sqrt(v); out[f"{name}_ci"] = [100 * max(0.0, pt - 1.96 * se), 100 * min(1.0, pt + 1.96 * se)]
        binom = math.sqrt(pt * (1 - pt) / den.sum()) if 0 < pt < 1 else None
        out[f"{name}_design_effect"] = (se / binom) ** 2 if binom else None
    return out
boot = {}
for p in precs_present:
    sp = P[p]; sc = PC[p]; d = {}
    ALLI = np.ones(len(inj), bool); ALLC = np.ones(len(cln), bool)
    for t in (1e-6, 1e-3):
        d[f"harmful_{t:g}"] = lin_ci(harmful(t), ALLI, inj["cfg"], sp)
    for ck, t, kind in (("rb", 1e-6, "cov"), ("comb", 1e-3, "cov"), ("comb", 1e-3, "wr"), ("tol", 1e-3, "wr"), ("dup", 1e-3, "wr")):
        den = harmful(t) if kind == "cov" else ~harmful(t)
        d[f"{ck}_{kind}_{t:g}"] = lin_ci(fire(ck, t) & den, den, inj["cfg"], sp)
    for ck, t in (("fixed", None), ("tol", 1e-6)):
        d[f"{ck}_fa" + ("" if t is None else f"_{t:g}")] = lin_ci(clean_fire(ck, t), ALLC, cln["cfg"], sc)
    boot[PREC[p]] = d
S["dependence_aware_intervals"] = boot

# --------------------------------------------------------------------------- 14. calibration protocol sizes (per precision)
proto = {}
for name, g_i, g_c, vals in (("leave_one_seed_out", fi["seed"], fc["seed"], range(3)), ("leave_one_dist_out", fi["dist"], fc["dist"], range(5))):
    v = 0; calc = PC[0] & (g_c != v); held = PC[0] & (g_c == v); hi_ = P[0] & (g_i == v)
    proto[name] = {"folds": len(list(vals)), "calibration_checks": int(calc.sum()), "calibration_gemms": len(np.unique(gid_c[calc])),
                   "heldout_checks": int(held.sum()), "heldout_gemms": len(np.unique(gid_c[held])), "heldout_injections": int(hi_.sum())}
proto["insample"] = {"calibration_checks": int(PC[0].sum()), "calibration_gemms": len(np.unique(gid_c[PC[0]])), "scored_injections": int(P[0].sum())}
S["calibration_protocol_sizes"] = proto

# --------------------------------------------------------------------------- 15. timing per unit of work, and cost break-even
if os.path.exists(tpath):
    tw = {}
    for (pr, n, op), v in tm.items():
        tw[f"{pr}|{n}|{op}"] = {"ns": v, "ns_per_mac": v / n**3 if op in ("gemm", "dup") else None, "ns_per_output_element": v / n**2 if op.startswith("abft") else None}
    S["timing_per_unit"] = tw
    be = {}
    for key, e in S["cost_model_tau_1e-3"].items():
        pr, n = key.split("|"); g = tm[(pr, int(n), "gemm")]; pi = PREC.index(pr)
        si = P[pi] & (fi["n"] == int(n)); h = harmful(1e-3)[si]; Ph = h.mean()
        def pfire_rate(ck):
            f = fire(ck, 1e-3)[si]; return Ph * np.sum(f & h) / max(h.sum(), 1) + (1 - Ph) * np.sum(f & ~h) / max((~h).sum(), 1)
        out = {"P_h": Ph}
        for R in (1, 10, 100, 1000):
            # recovery cost R GEMMs; E = 1 + det + R * P(fire); break-even p between rb and comb
            d_rb = e["rb"]["detector_overhead_pct"] / 100; d_cb = e["comb"]["detector_overhead_pct"] / 100
            q_rb, q_cb = pfire_rate("rb"), pfire_rate("comb")
            out[f"R={R}"] = {"p_breakeven_rb_vs_comb": (d_cb - d_rb) / (R * (q_rb - q_cb)) if q_rb > q_cb else None,
                             "per_fault_recovery_rb": R * q_rb, "per_fault_recovery_comb": R * q_cb, "per_fault_recovery_dup": R * pfire_rate("dup")}
        # sensitivity: false-alarm probability per GEMM at its 95% upper bound (3/N GEMMs) instead of the observed zero
        ng = len(np.unique(gid_c[PC[pi] & (fc["n"] == int(n))])); out["gemm_fa_upper95"] = 3.0 / ng
        be[key] = out
    S["cost_breakeven_tau_1e-3"] = be

json.dump(S, open(os.path.join(D, "summary.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))

# --------------------------------------------------------------------------- report.md
L = []; A = L.append
A("# sdc-bench results\n")
A(f"Injections: {meta['injections']:,}; fault-free checks: {meta['clean_checks']:,} ({meta['clean_problems_per_config']} fault-free GEMMs per configuration); precisions: {', '.join(meta['precisions'])}.\n")
A("Percentages with 95% Wilson intervals in summary.json. A zero count is shown as `0.00 (≤x)`, where x is the 95% upper bound 3/N.\n")
A("## Headline by precision\n")
A("| Precision | bitwise | harmful 1e-6 | harmful 1e-3 | fixed FA | RB FA | RB cov 1e-6 | RB WR 1e-6 | comb cov 1e-3 | comb WR 1e-3 | comb FA 1e-3 | dup WR 1e-3 | non-finite |")
A("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for p, d in head.items():
    A(f"| {p} | {fmt(d['bitwise'])} | {fmt(d['harmful_1e-06'])} | {fmt(d['harmful_0.001'])} | {fmt(d['fixed_fa'])} | {fmt(d['rb_fa'])} | {fmt(d['rb_cov_1e-06'])} | {fmt(d['rb_wr_1e-06'])} | {fmt(d['comb_cov_0.001'])} | {fmt(d['comb_wr_0.001'])} | {fmt(d['comb_fa_0.001'])} | {fmt(d['dup_wr_0.001'])} | {d['nonfinite']} |")
A("\n## Trade-off at tau = 1e-3 (coverage / wasted recovery / fault-free false alarms)\n")
A("| Precision | Checker | Coverage | Wasted recovery | False alarms |\n|---|---|---|---|---|")
for p, d in head.items():
    for ck, fk in (("dup", None), ("rb", "rb_fa"), ("rb_rig", "rb_rig_fa"), ("tol", "tol_fa_0.001"), ("comb", "comb_fa_0.001"), ("comb_rig", "comb_rig_fa_0.001"), ("fixed", "fixed_fa"), ("range", None)):
        fa_txt = fmt(d[fk]) if fk else ("0 (by construction)" if ck == "dup" else fmt(rng_fa[p]) + " per GEMM")
        A(f"| {p} | {ck} | {fmt(d[ck + '_cov_0.001'])} | {fmt(d[ck + '_wr_0.001'])} | {fa_txt} |")
for title, key in (("fault model", "by_fault_model"), ("input distribution", "by_distribution"), ("matrix size", "by_size"), ("bit class", "by_bitclass"), ("seed", "by_seed")):
    A(f"\n## Breakdown by {title}\n")
    A("| Cell | n inj | bitwise | harmful 1e-6 | harmful 1e-3 | RB cov 1e-6 | comb cov 1e-3 | comb WR 1e-3 | tol cov 1e-3 | tol WR 1e-3 | fixed FA | RB FA |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for k, d in S[key].items():
        A(f"| {k} | {d['injections']} | {fmt(d['bitwise'])} | {fmt(d['harmful_1e-06'])} | {fmt(d['harmful_0.001'])} | {fmt(d['rb_cov_1e-06'])} | {fmt(d['comb_cov_0.001'])} | {fmt(d['comb_wr_0.001'])} | {fmt(d['tol_cov_0.001'])} | {fmt(d['tol_wr_0.001'])} | {fmt(d.get('fixed_fa'))} | {fmt(d.get('rb_fa'))} |")
A("\n## Threshold calibration: in-sample vs held-out\n")
A("Family `lin` is the classical rounding-bound shape (∝ 2n), `sqrt` ∝ √n, `const`, `invsqrt` ∝ n^-1/2, `pern` = separate scale per n, `tolrel` = scaled fixed tolerance. The zero-false-alarm scale is chosen on the calibration checks and scored on the held-out ones.\n")
A("| Precision|family | in-sample cov 1e-6 | LOSO held-out FA | LOSO cov 1e-6 | LODO held-out FA | LODO cov 1e-6 | LODO cov 1e-3 | LODO WR 1e-3 |")
A("|---|---|---|---|---|---|---|---|")
for k, r in cal.items():
    A(f"| {k} | {fmt(r['insample_all']['cov_1e-06'])} | {fmt(r['leave_one_seed_out']['fa'])} | {fmt(r['leave_one_seed_out']['cov_1e-06'])} | {fmt(r['leave_one_dist_out']['fa'])} | {fmt(r['leave_one_dist_out']['cov_1e-06'])} | {fmt(r['leave_one_dist_out']['cov_0.001'])} | {fmt(r['leave_one_dist_out']['wr_0.001'])} |")
A("\n## Harm-metric sensitivity (max-norm vs Frobenius vs componentwise)\n")
A("| Precision|tau | metric | harmful | agrees with max-norm | RB cov | comb cov | comb WR |\n|---|---|---|---|---|---|---|")
for k, d in sens.items():
    for m, x in d.items(): A(f"| {k} | {m} | {fmt(x['harmful'])} | {fmt(x['agree_with_max'])} | {fmt(x['rb_cov'])} | {fmt(x['comb_cov'])} | {fmt(x['comb_wr'])} |")
A("\n## Non-finite outputs\n")
A("| Precision | count | NaN only | Inf only | both | unguarded RB catches | unguarded fixed catches | unguarded comb(1e-3) catches |\n|---|---|---|---|---|---|---|---|")
for p, d in nf.items(): A(f"| {p} | {d['count']} | {d['nan_only']} | {d['inf_only']} | {d['both']} | {fmt(d['unguarded_rb_caught'])} | {fmt(d['unguarded_fixed_caught'])} | {fmt(d['unguarded_comb_1e-3_caught'])} |")
A("\n## Why coverage varies by distribution: normaliser / max|C| and forward error\n")
A("| Precision | n | distribution | median (|A|(|B|1))_i / max|C| | fault-free forward error | RB cov 1e-6 |\n|---|---|---|---|---|---|")
for p, h in hyp.items():
    for r in h["cells"]: A(f"| {p} | {r['n']} | {r['dist']} | {r['norm_over_maxC_median']:.4g} | {r['fwd_err_median']:.3g} | {fmt(r['rb_cov_1e-06'])} |")
    A(f"| {p} | | Spearman(ratio/n, coverage) | {h['spearman_ratio_per_n_vs_cov']:.3f} | | |")
if "verification" in S:
    A("\n## Independent verification (Algorithm 3)\n")
    A("Category (3) = implementation error (bit-level disagreement). Category (1) = genuine blind spot (output changed, affected checksums bit-identical to fault-free). Category (2) = decision flips when the checksum side is summed pairwise instead of sequentially.\n")
    A("| Group | verified | (3) output mismatch | (3) check mismatch | changed | (1) cancellations | ...harmful 1e-6 | ...harmful 1e-3 | harmful 1e-6 | RB misses | (2) RB flips | fixed flips | comb 1e-3 flips |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for k, a in sorted(S["verification"].items()):
        A(f"| {k} | {a['verified']} | {a['mism_output']} | {a['mism_checks']} | {a['changed']} | {a['cancel']} | {a['cancel_h6']} | {a['cancel_h3']} | {a['h6']} | {a['rbmiss_h6']} | {a['flip_rb']} | {a['flip_fixed']} | {a['flip_comb3']} |")
if "cost_model_tau_1e-3" in S:
    A("\n## Cost model at tau = 1e-3 (expected cost per GEMM relative to one unprotected GEMM)\n")
    A("| Precision|n | checker | detector overhead % | GEMM false-alarm prob | E[cost] p=1e-9 | p=1e-6 | p=1e-3 | p=1e-1 |\n|---|---|---|---|---|---|---|---|---|")
    for k, e in S["cost_model_tau_1e-3"].items():
        for ck in ("dup", "rb", "comb"):
            x = e[ck]; A(f"| {k} | {ck} | {x['detector_overhead_pct']:.2f} | {x['gemm_false_alarm_prob']:.3g} | " + " | ".join(f"{x[f'p={pf:g}']['expected_cost_rel_gemm']:.4f}" for pf in (1e-9, 1e-6, 1e-3, 1e-1)) + " |")
A("\n## Tolerance normalisation: raw tau*max|C'| vs corrected tau/(1+tau)*max|C'|\n")
A("Misses in window = harmful misses with tau < r <= tau/(1-tau), the only ones the raw rule can lose to self-normalisation.\n")
A("| Precision|tau | harmful | raw cov | raw misses (in window) | corrected cov | corrected misses (in window) | raw WR | corrected WR |\n|---|---|---|---|---|---|---|---|---|")
for k, d in S["self_normalisation"].items():
    if k.startswith("grid") or k.startswith("tol_misses"): continue
    A(f"| {k} | {d['harmful']} | {fmt(d['tol_raw']['cov'])} | {d['tol_raw']['miss']} ({d['tol_raw']['miss_in_window']}) | {fmt(d['tol']['cov'])} | {d['tol']['miss']} ({d['tol']['miss_in_window']}) | {fmt(d['tol_raw']['wr'])} | {fmt(d['tol']['wr'])} |")
A(f"\nGrid 1e-3..1: {json.dumps(S['self_normalisation']['grid_1e-3_to_1'])}\n")
A("\n## False alarms per fault-free GEMM (a GEMM counts once if any of its 2n checks fires)\n")
A("| Precision | " + " | ".join(next(iter(S["gemm_level_false_alarms"].values())).keys()) + " |")
A("|---|" + "---|" * len(next(iter(S["gemm_level_false_alarms"].values()))))
for p, d in S["gemm_level_false_alarms"].items(): A(f"| {p} | " + " | ".join(fmt(v) for v in d.values()) + " |")
A("\n## Held-out calibration, false alarms per fault-free GEMM\n")
A("| Precision|family | LOSO GEMMs with alarm | LODO GEMMs with alarm | LODO per fold (uniform, normal, lognormal, illcond, nonneg) |\n|---|---|---|---|")
for k, v in S["heldout_gemm_level_false_alarms"].items():
    A(f"| {k} | {fmt(v['leave_one_seed_out']['fa_gemm'])} | {fmt(v['leave_one_dist_out']['fa_gemm'])} | " + ", ".join(str(f['k']) for f in v['leave_one_dist_out']['folds']) + " |")
A("\n## Dependence-aware 95% intervals (clusters = configurations)\n")
A("| Precision | metric | estimate | Wilson (independence) | stratified by size x distribution | unstratified | design effect (stratified) |\n|---|---|---|---|---|---|---|")
for p, d in S["dependence_aware_intervals"].items():
    for k, v in d.items():
        de = v["stratified_design_effect"]
        A(f"| {p} | {k} | {v['pct']:.2f} | [{v['wilson_ci'][0]:.2f}, {v['wilson_ci'][1]:.2f}] | [{v['stratified_ci'][0]:.2f}, {v['stratified_ci'][1]:.2f}] | [{v['unstratified_ci'][0]:.2f}, {v['unstratified_ci'][1]:.2f}] | {'--' if de is None else f'{de:.1f}'} |")
A(f"\nCalibration protocol sizes (per precision): {json.dumps(S['calibration_protocol_sizes'])}\n")
if "cost_breakeven_tau_1e-3" in S:
    A("\n## Cost: break-even fault probability per GEMM between rounding bound and combined rule (tau = 1e-3)\n")
    A("Recovery costs R GEMMs. Below the break-even p the rounding bound is cheaper (smaller detector); above it the combined rule is cheaper (less wasted recovery).\n")
    A("| Precision|n | P_h | R=1 | R=10 | R=100 | R=1000 | GEMM false-alarm 95% upper bound |\n|---|---|---|---|---|---|---|---|")
    for k, v in S["cost_breakeven_tau_1e-3"].items():
        cell = lambda x: "none" if x is None or x > 1 else f"{x:.3g}"
        A(f"| {k} | {v['P_h']:.4f} | " + " | ".join(cell(v[f'R={R}']['p_breakeven_rb_vs_comb']) for R in (1, 10, 100, 1000)) + f" | {v['gemm_fa_upper95']:.3g} |")
    A("\n## Timing per unit of work\n")
    A("| Precision|n|op | ns | ns per multiply-accumulate | ns per output element |\n|---|---|---|---|")
    for k, v in sorted(S["timing_per_unit"].items()):
        mac = "" if v["ns_per_mac"] is None else f"{v['ns_per_mac']:.3f}"
        ele = "" if v["ns_per_output_element"] is None else f"{v['ns_per_output_element']:.2f}"
        A(f"| {k} | {v['ns']:.0f} | {mac} | {ele} |")
open(os.path.join(D, "report.md"), "w").write("\n".join(L) + "\n")
print(f"[analyze] wrote {os.path.join(D, 'report.md')}, summary.json and fig/*.csv")
