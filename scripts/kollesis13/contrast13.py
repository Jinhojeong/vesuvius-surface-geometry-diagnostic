#!/usr/bin/env python3
"""kollesis13, step 2: mass contrast on the census13 1.3-1.6 stratum.

Design = kollesis_lattice2, with the rescoped strata:
  TEST       = clusters with ratio in [1.3, 1.6)  (the band the earlier
               null never tested)
  BACKGROUND = clusters with ratio >= 1.6         (same coverage mask)
  stat = frac(test on fence) - frac(background on fence)

Arms:
  A) W = 160.0 fixed (manufacture prior), offset fitted on TEST
  B) W free 120-200 (fitted on TEST) - edge-drift check

Null: label permutation over converted clusters, fence held at the
observed fit (upward bias for the observed arm; a flat null despite the
bias is conservative). Positive p-values require the refit-per-permutation
null (refit_null13.py) before any claim.

Robustness: cuts to rings <= 38 and <= 16.
Sensitivity: injection - clusters per joint (fence at ~160 mm across the
measured arc) needed for arm A to detect at p < 0.05.

Reads results/kollesis13/census13_arc.csv.gz; writes contrast13_report.json.
"""
import gzip, csv, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import kollesis_lattice as K
from kollesis_lattice2 import fit_offset_only, on_frac

OUTDIR = os.path.join(REPO, "results", "kollesis13")
RNG = np.random.default_rng(20260809)
NPERM = 2000
TOL = 0.12

rows = []
with gzip.open(os.path.join(OUTDIR, "census13_arc.csv.gz"), "rt") as f:
    for r in csv.DictReader(f):
        rows.append((float(r["s_mm"]), float(r["ratio"]), int(r["ring"])))
s = np.array([r[0] for r in rows])
ratio = np.array([r[1] for r in rows])
ring = np.array([r[2] for r in rows])


def contrast(sv, is_test, W, o):
    return on_frac(sv[is_test], W, o) - on_frac(sv[~is_test], W, o)


def perm_null(sv, is_test, W, o, nperm=NPERM):
    n_test = int(is_test.sum())
    obs = contrast(sv, is_test, W, o)
    on = np.abs(((sv - o + W / 2) % W) - W / 2) <= TOL * W
    null = []
    idx = np.arange(len(sv))
    for _ in range(nperm):
        RNG.shuffle(idx)
        t = idx[:n_test]
        ft = on[t].mean()
        fb = on[idx[n_test:]].mean()
        null.append(ft - fb)
    null = np.array(null)
    p = float(np.mean(null >= obs))
    return obs, p, float(null.mean()), float(null.std())


def run_cut(name, keep):
    sv, isT = s[keep], (ratio[keep] < 1.6)
    nT, nB = int(isT.sum()), int((~isT).sum())
    out = {"n_test": nT, "n_background": nB}
    # arm A: W=160 fixed, offset on TEST
    oA, fA = fit_offset_only(sv[isT], 160.0)
    obsA, pA, muA, sdA = perm_null(sv, isT, 160.0, oA)
    out["W160"] = {"offset_mm": oA, "frac_test": on_frac(sv[isT], 160.0, oA),
                   "frac_background": on_frac(sv[~isT], 160.0, oA),
                   "contrast_pp": obsA * 100, "perm_p": pA,
                   "null_mean_pp": muA * 100, "null_sd_pp": sdA * 100}
    # arm B: W free on TEST
    W, o, sc = K.fence_fit_fast(sv[isT])
    obsB, pB, muB, sdB = perm_null(sv, isT, W, o)
    out["Wfree"] = {"W_mm": W, "offset_mm": o,
                    "edge_drift": bool(W <= 121 or W >= 199),
                    "frac_test": on_frac(sv[isT], W, o),
                    "frac_background": on_frac(sv[~isT], W, o),
                    "contrast_pp": obsB * 100, "perm_p": pB,
                    "null_mean_pp": muB * 100, "null_sd_pp": sdB * 100}
    print(f"[{name}] test {nT:,} bg {nB:,} | "
          f"W160: contrast {obsA*100:+.2f}pp p={pA:.3f} | "
          f"Wfree: W={W:.1f} contrast {obsB*100:+.2f}pp p={pB:.3f}",
          flush=True)
    return out


report = {"design": "mass contrast test(1.3-1.6) vs background(>=1.6), "
                    "census13 pass, label-permutation null at the "
                    f"observed-fitted fence (nperm={NPERM}); W free as "
                    "edge-drift check",
          "cuts": {}}
report["cuts"]["all"] = run_cut("all", np.ones(len(s), bool))
report["cuts"]["ring_le_38"] = run_cut("ring<=38", ring <= 38)
report["cuts"]["ring_le_16"] = run_cut("ring<=16", ring <= 16)

# ---- sensitivity by injection (all cut, W=160 arm) ----
isT = ratio < 1.6
sT = s[isT].copy()
smin, smax = s.min(), s.max()
joints = np.arange(smin + 80.0, smax, 160.0)
print(f"[sens] {len(joints)} synthetic joints in s range "
      f"[{smin:.0f},{smax:.0f}] mm", flush=True)
sens = {}
for per_joint in (5, 10, 15, 20, 30, 40, 60):
    hits = 0
    n_inj = int(per_joint * len(joints))
    if n_inj > len(sT):
        break
    for trial in range(20):
        rng = np.random.default_rng(1000 + trial)
        sm = sT.copy()
        pick = rng.choice(len(sm), size=n_inj, replace=False)
        jpos = joints[rng.integers(0, len(joints), size=n_inj)]
        sm[pick] = jpos + rng.normal(0, 0.12 * 160 / 2, size=n_inj)
        sv = np.concatenate([sm, s[~isT]])
        lab = np.concatenate([np.ones(len(sm), bool),
                              np.zeros((~isT).sum(), bool)])
        o, _ = fit_offset_only(sm, 160.0)
        _, p, _, _ = perm_null(sv, lab, 160.0, o, nperm=500)
        hits += (p < 0.05)
    sens[per_joint] = hits / 20
    print(f"  {per_joint}/joint ({n_inj:,} moved of {len(sT):,}): "
          f"power {hits/20:.2f}", flush=True)
    if hits / 20 >= 0.95:
        break
report["sensitivity_W160"] = {"joints_assumed": len(joints),
                              "power_by_clusters_per_joint": sens}

# ---- verdict ----
a = report["cuts"]["all"]
pos = (a["W160"]["perm_p"] < 0.05 and a["W160"]["contrast_pp"] > 0 and
       150.0 <= a["Wfree"]["W_mm"] <= 170.0 and a["Wfree"]["perm_p"] < 0.05)
report["verdict"] = "lattice_present" if pos else "no_lattice_detected"
json.dump(report, open(os.path.join(OUTDIR, "contrast13_report.json"), "w"),
          indent=1)
print(f"[verdict] {report['verdict']}", flush=True)
