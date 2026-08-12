#!/usr/bin/env python3
"""powered16_analyze.py - pooled 16-site run-level analysis of the extended
GrowPatch hazard experiment (8 archived sites + 8 new sites).

Method is randctl_analyze.py / powered_analyze.py verbatim (documented
there): unit of analysis is the RUN, per site x arm replicate means, paired
t across sites of the mean differences (df = n_sites - 1), scipy p, t-crit
CIs, MDE80 as in powered_analyze. NO ray-pooled CIs anywhere. score_ab.py
and powered_analyze.py are NOT modified; this script replicates their
documented method and must first reproduce the archived 8-site numbers
exactly (gate) before the 16-site pool is computed.

Inputs (read-only):
  /mnt/vesuvius/vcbuild/demo_out/powered/powered_scores.json    A/B, sites 0-7
  /mnt/vesuvius/vcbuild/demo_out/randctl/randctl_scores.json    C,   sites 0-7
  /mnt/vesuvius/vcbuild/demo_out/randctl/randctl_verdict.json   gate targets
  /mnt/vesuvius/vcbuild/demo_out/powered16/powered16_scores_new8.json
                                                                A/B/C, 8-15
Outputs:
  /mnt/vesuvius/vcbuild/demo_out/powered16/RESULTS.json
  /mnt/vesuvius/vcbuild/demo_out/powered16/RESULTS.md
"""
import json
import math
import os
import time
from collections import defaultdict

import numpy as np
from scipy import stats

PW = "/mnt/vesuvius/vcbuild/demo_out/powered/powered_scores.json"
RC = "/mnt/vesuvius/vcbuild/demo_out/randctl/randctl_scores.json"
RV = "/mnt/vesuvius/vcbuild/demo_out/randctl/randctl_verdict.json"
PW16 = "/mnt/vesuvius/vcbuild/demo_out/powered16"
NEW = os.path.join(PW16, "powered16_scores_new8.json")
NEW_SITES = "/mnt/vesuvius/hazard_zarr_smoke/demo_sites_new8.json"

EPS = [("area", "vx2"), ("cross", "pp"), ("dbl", "pp"), ("on", "pp"),
       ("per1000", "per-1000")]
EP_LABEL = {"area": "area_vx2", "cross": "crossing frac_quads_near (pp)",
            "dbl": "ray frac_double (pp of on-sheet rays)",
            "on": "on_sheet_rate (guardrail, pp)",
            "per1000": "n_double per 1000 sampled quads"}


def metric(r, key):
    """Identical to randctl_analyze.metric, plus per1000 from
    powered_analyze.endpoints (1000 * n_double / n_sampled)."""
    if key == "area":
        return r["area_vx2"]
    if key == "cross":
        return r["crossing"]["frac_quads_near"] * 100.0
    if key == "dbl":
        v = r["ray"].get("frac_double")
        return None if v is None else v * 100.0
    if key == "on":
        n = r["ray"].get("n_sampled")
        if not n:
            return None
        return 100.0 * r["ray"].get("n_on_sheet", 0) / n
    if key == "per1000":
        n = r["ray"].get("n_sampled")
        if not n:
            return None
        return 1000.0 * r["ray"].get("n_double", 0) / n
    raise KeyError(key)


def site_mean(runs, key):
    v = [metric(r, key) for r in runs]
    v = [x for x in v if x is not None]
    return (sum(v) / len(v)) if v else None


def paired(diffs):
    d = np.asarray(diffs, float)
    n = len(d)
    mean = float(d.mean())
    sd = float(d.std(ddof=1))
    se = sd / math.sqrt(n)
    t = mean / se if se > 0 else float("inf")
    p = float(2 * stats.t.sf(abs(t), n - 1))
    tcrit = float(stats.t.ppf(0.975, n - 1))
    mde80 = (tcrit + float(stats.t.ppf(0.80, n - 1))) * se
    return {"n_sites": n, "mean_diff": mean, "sd_diff": sd, "t": t,
            "df": n - 1, "p": p,
            "ci95": [mean - tcrit * se, mean + tcrit * se],
            "mde80": mde80,
            "sites_positive": int(sum(1 for x in d if x > 0))}


def contrast(M, sites, hi, lo, key):
    diffs, used = [], []
    for s in sites:
        a = site_mean(M[s][lo], key)
        b = site_mean(M[s][hi], key)
        if a is None or b is None:
            continue
        diffs.append(b - a)
        used.append(s)
    if len(diffs) < 2:
        return None
    out = paired(diffs)
    out["per_site_diff"] = {s: d for s, d in zip(used, diffs)}
    return out


# ------------------------------------------------------------------- load
pw = json.load(open(PW))
rc = json.load(open(RC))
M = defaultdict(lambda: defaultdict(list))
for r in pw["runs"] + rc["runs"]:
    M[r["site"]][r["arm"]].append(r)
old_sites = sorted(M.keys())
assert len(old_sites) == 8, f"expected 8 archived sites, got {len(old_sites)}"

# ------------------------------------------------------------------- GATE
# reproduce the archived 8-site contrasts exactly from the archived JSONs
print("=== GATE: reproduce archived 8-site numbers ===")
rv = json.load(open(RV))
gate = {"pass": True, "checks": []}


def close(a, b, tol=1e-9):
    if isinstance(a, list):
        return all(close(x, y, tol) for x, y in zip(a, b))
    if a is None or b is None:
        return a is b
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


for name_hi, name_lo in (("C", "A"), ("B", "C"), ("B", "A")):
    cname = f"{name_hi}-{name_lo}"
    for key, unit in EPS:
        if key == "per1000":
            continue  # not in the archived verdict
        got = contrast(M, old_sites, name_hi, name_lo, key)
        want = rv["contrasts"][cname][key]
        ok = all(close(got[k], want[k]) for k in
                 ("n_sites", "mean_diff", "sd_diff", "t", "p", "ci95"))
        gate["checks"].append({"contrast": cname, "endpoint": key, "ok": ok,
                               "got": {k: got[k] for k in
                                       ("mean_diff", "sd_diff", "t", "p", "ci95")},
                               "want": {k: want[k] for k in
                                        ("mean_diff", "sd_diff", "t", "p", "ci95")}})
        gate["pass"] &= ok
        print("  %-4s %-5s %s  got mean=%+.4f p=%.4f | want mean=%+.4f p=%.4f"
              % (cname, key, "OK  " if ok else "FAIL",
                 got["mean_diff"], got["p"], want["mean_diff"], want["p"]))

# headline spot checks against POWERED_VERDICT.md published rounding
ba_area = contrast(M, old_sites, "B", "A", "area")
ba_dbl = contrast(M, old_sites, "B", "A", "dbl")
spot = (round(ba_area["mean_diff"] / 1000, 1) == 234.0
        and round(ba_area["p"], 3) == 0.003
        and round(ba_dbl["mean_diff"], 2) == -2.65
        and round(ba_dbl["p"], 3) == 0.068)
gate["headline_spot_check"] = {
    "ok": spot,
    "area_mean_k": round(ba_area["mean_diff"] / 1000, 1), "area_p": round(ba_area["p"], 3),
    "dbl_mean_pp": round(ba_dbl["mean_diff"], 2), "dbl_p": round(ba_dbl["p"], 3),
    "want": "area +234.0k p=0.003; dbl -2.65pp p=0.068"}
gate["pass"] &= spot
print("  headline spot check:", "OK" if spot else "FAIL",
      gate["headline_spot_check"])
print("GATE", "PASS" if gate["pass"] else "FAIL")
if not gate["pass"]:
    json.dump(gate, open(os.path.join(PW16, "gate_FAIL.json"), "w"), indent=1)
    raise SystemExit("gate failed: archived 8-site numbers not reproduced")

# ------------------------------------------------------------- 16-site pool
if not os.path.exists(NEW):
    print("no powered16_scores_new8.json yet; gate-only run, exiting")
    raise SystemExit(0)

new = json.load(open(NEW))
for r in new["runs"]:
    M[r["site"]][r["arm"]].append(r)
sites16 = sorted(M.keys(), key=lambda s: int(s.split("_")[0][1:]))
print(f"\n=== 16-site pooled analysis ({len(sites16)} sites) ===")

contrasts = {}
for hi, lo in (("B", "A"), ("B", "C"), ("C", "A")):
    cname = f"{hi}-{lo}"
    contrasts[cname] = {}
    for key, unit in EPS:
        contrasts[cname][key] = contrast(M, sites16, hi, lo, key)

# cells: per site x arm mean/sd for the report
cells = {}
for s in sites16:
    cells[s] = {}
    for arm in ("A", "B", "C"):
        runs = M[s][arm]
        cells[s][arm] = {"n_runs": len(runs)}
        for key, _ in EPS:
            v = [metric(r, key) for r in runs]
            v = [x for x in v if x is not None]
            a = np.asarray(v, float)
            cells[s][arm][key] = {
                "n": len(v),
                "mean": float(a.mean()) if len(v) else None,
                "sd": float(a.std(ddof=1)) if len(v) > 1 else None}

# leave-one-site-out sensitivity on dbl (historical fragility: same-tile pair)
loo = {}
for cname, (hi, lo) in (("B-A", ("B", "A")), ("B-C", ("B", "C"))):
    rows = {}
    for drop in sites16:
        sub = [s for s in sites16 if s != drop]
        c = contrast(M, sub, hi, lo, "dbl")
        rows[drop] = {"mean_diff": c["mean_diff"], "p": c["p"],
                      "ci95": c["ci95"]}
    ps = [r["p"] for r in rows.values()]
    loo[cname] = {"per_dropped_site": rows,
                  "p_min": min(ps), "p_max": max(ps),
                  "all_below_0.05": all(p < 0.05 for p in ps)}

# block heterogeneity: archived 8 vs new 8 (the new block is all mid-scroll
# by construction; the archived block mostly was not, so pooled effects can
# be diluted by a real block difference rather than noise alone)
old8 = [s for s in sites16 if int(s.split("_")[0][1:]) < 8]
new8 = [s for s in sites16 if int(s.split("_")[0][1:]) >= 8]
blocks = {}
for bname, bsites in (("old8", old8), ("new8", new8)):
    blocks[bname] = {}
    for hi, lo in (("B", "A"), ("B", "C"), ("C", "A")):
        blocks[bname][f"{hi}-{lo}"] = {
            key: contrast(M, bsites, hi, lo, key) for key, _ in EPS}

guard = contrasts["B-A"]["on"]
guard_flag = guard["p"] < 0.05

# new-site selection record
new_sites_rec = json.load(open(NEW_SITES))
selection = {
    "rule": "select_sites16.py, deterministic, no RNG: candidate pool and "
            "thresholds identical to demo_sites.py (census8k clusters, "
            "ratio >= 3.0, n_sites >= 2, tile-interior margin 48, sort by "
            "(-n_sites, -ratio)); constraints: exclude the 7 tiles used by "
            "the original 8 sites, at most one new site per tile, slab cap "
            "2 across old+new, m7 L1 seed value >= 128 (feasibility); "
            "pass 1 picks top 8 among mid-scroll band candidates "
            "(0.25*zmax < gz < 0.8*zmax, zmax=11405), top-up pass unused "
            "(8/8 filled in pass 1); 0 candidates skipped for m7 < 128",
    "sites": new_sites_rec}

results = {
    "generated": time.strftime("%F %T"),
    "design": "16 sites (8 archived + 8 new, all-distinct tiles) x arms "
              "A (field only) / B (field + p1218_conf_v2_amp4.zarr) / "
              "C (random placebo) x 5 fresh reps per arm; "
              "unit = run; paired t across sites of replicate means, df=15",
    "arms": {"A": "direction field only",
             "B": "field + /mnt/vesuvius/hazard_zarr/p1218_conf_v2_amp4.zarr",
             "C": "sites 0-7: p1218_conf_rand_amp4.zarr (rng 20260806); "
                  "sites 8-15: p1218_conf_rand_amp4_ext8.zarr (rng 20260812); "
                  "both region-scoped random redraws, radius 288 L1, "
                  "amp_scale 4, matched count/sigma/amp/support"},
    "inputs": {"powered_scores": PW, "randctl_scores": RC,
               "new8_scores": NEW, "new_sites": NEW_SITES},
    "gate": gate,
    "n_runs": {"archived_AB": len(pw["runs"]), "archived_C": len(rc["runs"]),
               "new8_ABC": len(new["runs"]), "missing_new8": new["missing"]},
    "new_site_selection": selection,
    "cells": cells,
    "contrasts": contrasts,
    "block_contrasts": blocks,
    "leave_one_site_out_dbl": loo,
    "guardrail_on_sheet_confounded_BA": bool(guard_flag),
}
json.dump(results, open(os.path.join(PW16, "RESULTS.json"), "w"), indent=1)
print("wrote", os.path.join(PW16, "RESULTS.json"))

# --------------------------------------------------------------- RESULTS.md
def fmt(v, key):
    if v is None:
        return "n/a"
    if key == "area":
        return f"{v/1000:+.1f}k" if abs(v) < 1e7 else f"{v:.0f}"
    return f"{v:+.2f}"


L = ["# Powered16: GrowPatch hazard A/B/C pooled over 16 sites", ""]
L.append(f"Generated {time.strftime('%F %T')}. "
         "8 archived sites (powered + randctl blocks, 2026-08-06) + 8 new "
         "sites (this block, 2026-08-12), every site on its own census "
         "tile. Unit = run; paired t across the 16 site pairs of replicate "
         "means (df=15); scorer score_ab.score_run UNMODIFIED (2000 rays, "
         "near_r 6). Gate: this analyzer reproduced the archived 8-site "
         "contrasts exactly before pooling (see gate block in RESULTS.json).")
L.append("")
L.append("## Paired contrasts (16 sites)")
L.append("")
for cname in ("B-A", "B-C", "C-A"):
    L.append(f"### {cname}")
    L.append("")
    L.append("| endpoint | mean diff | sd | t(15) | p | 95% CI | MDE80 | sites>0 |")
    L.append("|---|---|---|---|---|---|---|---|")
    for key, unit in EPS:
        c = contrasts[cname][key]
        if c is None:
            L.append(f"| {EP_LABEL[key]} | n/a |")
            continue
        L.append("| %s | %s %s | %s | %.2f | %.4f | [%s, %s] | %s | %d/%d |" % (
            EP_LABEL[key], fmt(c["mean_diff"], key), unit,
            fmt(c["sd_diff"], key).lstrip("+"), c["t"], c["p"],
            fmt(c["ci95"][0], key), fmt(c["ci95"][1], key),
            fmt(c["mde80"], key).lstrip("+"),
            c["sites_positive"], c["n_sites"]))
    L.append("")

L.append("## Per-site replicate means (dbl pp / area k-vx2)")
L.append("")
L.append("| site | A dbl | B dbl | C dbl | A area | B area | C area |")
L.append("|---|---|---|---|---|---|---|")
for s in sites16:
    row = [s.split("_")[0]]
    for key in ("dbl", "area"):
        for arm in ("A", "B", "C"):
            m = cells[s][arm][key]["mean"]
            if m is None:
                row.append("n/a")
            elif key == "area":
                row.append(f"{m/1000:.0f}k")
            else:
                row.append(f"{m:.1f}")
    L.append("| " + " | ".join([row[0]] + row[1:4] + row[4:7]) + " |")
L.append("")

L.append("## Leave-one-site-out sensitivity (dbl)")
L.append("")
for cname in ("B-A", "B-C"):
    lo_ = loo[cname]
    L.append(f"### {cname}: p range [{lo_['p_min']:.4f}, {lo_['p_max']:.4f}]"
             f"{' (all < 0.05)' if lo_['all_below_0.05'] else ''}")
    L.append("")
    L.append("| dropped site | mean diff (pp) | p |")
    L.append("|---|---|---|")
    for s in sites16:
        r = lo_["per_dropped_site"][s]
        L.append(f"| {s.split('_')[0]} | {r['mean_diff']:+.2f} | {r['p']:.4f} |")
    L.append("")

L.append("## Block heterogeneity (archived 8 vs new 8, B-A)")
L.append("")
L.append("The new block is all mid-scroll by construction (the selection "
         "mid-band pass filled all 8 slots); the archived block mostly sat "
         "near the scroll ends. Same paired method within each block "
         "(df=7).")
L.append("")
L.append("| endpoint | old8 diff | old8 p | new8 diff | new8 p |")
L.append("|---|---|---|---|---|")
for key, unit in EPS:
    o = blocks["old8"]["B-A"][key]
    n = blocks["new8"]["B-A"][key]
    L.append("| %s | %s %s | %.4f | %s %s | %.4f |" % (
        EP_LABEL[key], fmt(o["mean_diff"], key), unit, o["p"],
        fmt(n["mean_diff"], key), unit, n["p"]))
L.append("")

L.append("## Guardrail")
L.append("")
if guard_flag:
    L.append(f"on_sheet_rate DIFFERS between B and A (p={guard['p']:.3f}); "
             "the frac_double denominator is confounded; lean on per1000.")
else:
    L.append(f"on_sheet_rate does not differ between B and A "
             f"(p={guard['p']:.3f}); the dbl denominator is arm-comparable.")
L.append("")
open(os.path.join(PW16, "RESULTS.md"), "w").write("\n".join(L))
print("wrote", os.path.join(PW16, "RESULTS.md"))
for cname in ("B-A", "B-C", "C-A"):
    for key, _ in EPS:
        c = contrasts[cname][key]
        if c:
            print("%-4s %-8s diff=%+12.3f sd=%10.3f t=%+6.2f p=%.4f" %
                  (cname, key, c["mean_diff"], c["sd_diff"], c["t"], c["p"]))
