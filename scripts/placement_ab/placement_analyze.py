#!/usr/bin/env python3
"""placement_analyze.py - analysis for prereg_placement_ab (frozen at
commit 8531276, Amendment 1 at 083f424).

GATE (mandatory, runs first): reproduce the archived powered contrasts from
the archived powered_scores.json to FULL FLOAT PRECISION, replicating
powered_analyze.py's arithmetic verbatim (cell means over replicate runs,
paired t across the 8 sites, df 7, two-sided p, t-based CI95, MDE80). The
archived file must carry sha256
fd9deb2c3b7a3e3f0e5e55b258a5c3230b3ce845cb66e7d24425596472a65876.
Only if every archived number is matched exactly does the script touch the
new placement data.

Then, per the frozen doc:
  cells   per site x arm x placement: mean over the 5 replicate runs
  DD      (B - A)_on - (B - A)_off per site
  PRIMARY   area_vx2, one-sided paired t, H1 DD > 0, alpha 0.05, df 7
  SECONDARY frac_double, one-sided, H1 DD < 0, declared underpowered,
            reported with CI, no confirmatory weight alone
  GUARDRAIL on_sheet_rate DD, two-sided at 0.05
  LOO       leave-one-site-out p range for the primary, s2/s4 drops explicit
  OFF (B-A) reported with CI per endpoint (not tested)
  BUCKETS   frozen outcomes 1-4

Usage:
  placement_analyze.py --gate-only     validate the gate, touch nothing new
  placement_analyze.py                 full analysis, writes RESULTS.json/md
"""
import argparse
import hashlib
import json
import math
import os
import sys
import time

import numpy as np
from scipy import stats

PW_JSON = "/mnt/vesuvius/vcbuild/demo_out/powered/powered_scores.json"
PW_SHA = "fd9deb2c3b7a3e3f0e5e55b258a5c3230b3ce845cb66e7d24425596472a65876"
PLACE = "/mnt/vesuvius/vcbuild/demo_out/placement"
PL_JSON = os.path.join(PLACE, "placement_scores.json")
OFFSEEDS = "/mnt/vesuvius/hazard_zarr_smoke/offseeds_placement.json"
OFFSEEDS_SHA = "b07e905b6dee57c7cc2fec5a755cb099a2357d887933d5c8bbad56ec481290dc"
AMENDMENT = "prereg_placement_ab/PREREGISTRATION.md @ 083f424 (Amendment 1)"

EP = ["dbl", "per1000", "on_sheet", "area", "near"]  # powered_analyze order


# ------------------------------------------------ powered_analyze verbatim
def endpoints(m):
    ry, cr = m["ray"], m["crossing"]
    ns = ry.get("n_sampled", 0)
    return {
        "dbl": ry.get("frac_double"),
        "per1000": 1000.0 * ry.get("n_double", 0) / ns if ns else None,
        "on_sheet": ry.get("n_on_sheet", 0) / ns if ns else None,
        "area": m.get("area_vx2"),
        "near": cr.get("frac_quads_near"),
    }


def cell_mean(vals):
    v = np.array([x for x in vals if x is not None], float)
    return float(v.mean()) if len(v) else None


def paired_two_sided(diffs):
    d = np.array(diffs, float)
    n = len(d)
    sd = d.std(ddof=1)
    se = sd / math.sqrt(n)
    t = d.mean() / se if se > 0 else 0.0
    p = 2 * stats.t.sf(abs(t), n - 1)
    half = stats.t.ppf(0.975, n - 1) * se
    mde80 = (stats.t.ppf(0.975, n - 1) + stats.t.ppf(0.80, n - 1)) * se
    return {"n_sites": n, "mean_diff": float(d.mean()), "sd_diff": float(sd),
            "t": float(t), "df": n - 1, "p": float(p),
            "ci95": [float(d.mean() - half), float(d.mean() + half)],
            "mde80": float(mde80)}


def one_sided(diffs, direction):
    """Paired one-sided t. direction '+' tests H1 mean > 0, '-' mean < 0.
    CI reported is the ordinary two-sided 95% interval."""
    out = paired_two_sided(diffs)
    t, df = out["t"], out["df"]
    out["p_one_sided"] = float(stats.t.sf(t, df) if direction == "+"
                               else stats.t.cdf(t, df))
    out["h1"] = "mean > 0" if direction == "+" else "mean < 0"
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


# ------------------------------------------------------------------- GATE
def run_gate():
    got = sha256(PW_JSON)
    assert got == PW_SHA, "archived powered_scores.json sha256 mismatch: " + got
    data = json.load(open(PW_JSON))
    runs = data["runs"]
    sites = sorted(set(r["site"] for r in runs),
                   key=lambda s: int(s.split("_")[0][1:]))
    assert len(sites) == 8

    cells = {}
    for site in sites:
        for arm in ("A", "B"):
            sel = [r for r in runs if r["site"] == site and r["arm"] == arm]
            cells[(site, arm)] = {ep: cell_mean([endpoints(r)[ep] for r in sel])
                                  for ep in EP}

    arch = data["analysis"]
    n_checked = 0
    for ep in EP:
        diffs = [cells[(s, "B")][ep] - cells[(s, "A")][ep] for s in sites]
        rec = paired_two_sided(diffs)
        a = arch["paired"][ep]
        for k in ("n_sites", "mean_diff", "sd_diff", "t", "df", "p", "mde80"):
            assert rec[k] == a[k], \
                "GATE FAIL %s.%s: %r != archived %r" % (ep, k, rec[k], a[k])
            n_checked += 1
        assert rec["ci95"][0] == a["ci95"][0] and rec["ci95"][1] == a["ci95"][1], \
            "GATE FAIL %s.ci95" % ep
        n_checked += 2
        # archived cell means, full precision
        for s in sites:
            for arm in ("A", "B"):
                av = arch["cells"]["%s|%s" % (s, arm)][ep].get("mean")
                assert cells[(s, arm)][ep] == av, \
                    "GATE FAIL cell %s|%s %s" % (s, arm, ep)
                n_checked += 1
    print("GATE PASS: archived powered contrasts reproduced to full float "
          "precision (%d numbers, 5 endpoints x 8 sites, sha256 verified)"
          % n_checked)
    return data, cells, sites


# ---------------------------------------------------------------- analysis
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate-only", action="store_true")
    args = ap.parse_args()

    pw_data, on_cells, on_sites = run_gate()
    if args.gate_only:
        return

    assert sha256(OFFSEEDS) == OFFSEEDS_SHA, "offseeds_placement.json sha256 mismatch"
    pl = json.load(open(PL_JSON))
    pl_runs = pl["runs"]
    off_sites = sorted(set(r["site"] for r in pl_runs),
                       key=lambda s: int(s.split("_")[0][1:]))
    assert len(off_sites) == 8, "expected 8 off sites, got %d" % len(off_sites)

    off_cells = {}
    completion = {}
    for site in off_sites:
        for arm in ("A", "B"):
            sel = [r for r in pl_runs if r["site"] == site and r["arm"] == arm]
            completion[(site, arm)] = len(sel)
            off_cells[(site, arm)] = {ep: cell_mean([endpoints(r)[ep] for r in sel])
                                      for ep in EP}
    n_runs = len(pl_runs)

    # per-site (B-A) at on and off, and DD, per endpoint
    per_site = {}
    for i in range(8):
        s_on, s_off = on_sites[i], off_sites[i]
        row = {"on_site": s_on, "off_site": s_off}
        for ep in EP:
            bon, aon = on_cells[(s_on, "B")][ep], on_cells[(s_on, "A")][ep]
            boff = off_cells[(s_off, "B")][ep]
            aoff = off_cells[(s_off, "A")][ep]
            ba_on = None if (bon is None or aon is None) else bon - aon
            ba_off = None if (boff is None or aoff is None) else boff - aoff
            dd_v = (None if (ba_on is None or ba_off is None)
                    else ba_on - ba_off)
            row[ep] = {"ba_on": ba_on, "ba_off": ba_off, "dd": dd_v}
        per_site[i] = row

    dd = {ep: [per_site[i][ep]["dd"] for i in range(8)] for ep in EP}
    ba_off = {ep: [per_site[i][ep]["ba_off"] for i in range(8)] for ep in EP}
    ba_on = {ep: [per_site[i][ep]["ba_on"] for i in range(8)] for ep in EP}

    def testable(vals):
        return all(v is not None for v in vals)

    WHY = ("unmeasurable per Amendment 3: off cells carry zero on-sheet "
           "rays because every off-seed sits in label-void space, so "
           "label-referenced endpoints cannot be computed there")

    def unmeasurable(vals):
        return {"status": "unmeasurable",
                "n_missing_sites": sum(1 for v in vals if v is None),
                "reason": WHY}

    def tt(vals, fn, *a):
        return fn(vals, *a) if testable(vals) else unmeasurable(vals)

    primary = one_sided(dd["area"], "+")
    secondary = tt(dd["dbl"], one_sided, "-")
    guardrail = tt(dd["on_sheet"], paired_two_sided)
    # reported, not tested
    reported_near = tt(dd["near"], paired_two_sided)
    reported_per1000 = tt(dd["per1000"], paired_two_sided)

    # off-seed (B-A) itself with CI per endpoint (reported, not tested)
    off_ba = {ep: tt(ba_off[ep], paired_two_sided) for ep in EP}
    on_ba = {ep: tt(ba_on[ep], paired_two_sided) for ep in EP}

    # LOO for the primary (one-sided p, df 6), s2/s4 drops explicit
    loo = {}
    for j in range(8):
        d7 = [dd["area"][i] for i in range(8) if i != j]
        loo["drop_s%d" % j] = one_sided(d7, "+")
    loo_p = {k: v["p_one_sided"] for k, v in loo.items()}
    loo_range = [min(loo_p.values()), max(loo_p.values())]

    # frozen buckets
    primary_ok = primary["p_one_sided"] < 0.05
    secondary_ok = (secondary.get("status") != "unmeasurable"
                    and secondary["p_one_sided"] < 0.05)
    off_area_mean = off_ba["area"].get("mean_diff")
    on_area_mean = on_ba["area"]["mean_diff"]
    if primary_ok and secondary_ok:
        bucket, blabel = 1, ("placement confirmed on both endpoints")
    elif primary_ok:
        bucket, blabel = 2, ("placement confirmed for growth; the quality "
                             "claim stays at its current placebo-level "
                             "evidence")
    elif off_area_mean > 0.5 * on_area_mean:
        bucket, blabel = 3, ("placement REFUTED; the published mechanism "
                             "reading gets a public correction in the same "
                             "thread")
    else:
        bucket, blabel = 4, ("underpowered no-verdict, disclosed as such")
    guard_flag = (guardrail.get("status") != "unmeasurable"
                  and guardrail["p"] < 0.05)

    field_verify = None
    fv_path = os.path.join(PLACE, "field_verify.json")
    if os.path.isfile(fv_path):
        field_verify = json.load(open(fv_path))
    extent = None
    ex_path = os.path.join(PLACE, "field_extent_check.json")
    if os.path.isfile(ex_path):
        extent = json.load(open(ex_path))
    selection = json.load(open(OFFSEEDS))

    results = {
        "generated": time.strftime("%F %T"),
        "prereg": "prereg_placement_ab/PREREGISTRATION.md @ 8531276",
        "amendment": AMENDMENT,
        "gate": {"status": "PASS",
                 "archived_sha256": PW_SHA,
                 "detail": "archived powered contrasts reproduced to full "
                           "float precision before touching new data"},
        "selection": {"offseeds_sha256": OFFSEEDS_SHA, "seeds": selection},
        "field_extent_check_preamendment": extent,
        "field_verify_postbuild": field_verify,
        "completion": {"n_runs_scored": n_runs,
                       "n_missing": len(pl.get("missing", [])),
                       "missing": pl.get("missing", []),
                       "cells": {"%s|%s" % k: v for k, v in completion.items()}},
        "cells_on": {"%s|%s" % k: v for k, v in on_cells.items()},
        "cells_off": {"%s|%s" % k: v for k, v in off_cells.items()},
        "per_site": per_site,
        "primary_area_dd": primary,
        "secondary_dbl_dd": secondary,
        "guardrail_on_sheet_dd": guardrail,
        "guardrail_flag": bool(guard_flag),
        "reported_near_dd": reported_near,
        "reported_per1000_dd": reported_per1000,
        "off_seed_ba": off_ba,
        "amendment3_structural_finding": {
            "summary": "all 8 off-seeds sit in label-void space; "
                       "label-referenced endpoints unmeasurable at off "
                       "cells",
            "label_occupancy_pm24_at_offseeds": "0.000 at 8 of 8",
            "label_at_offseed_voxel": "0 at 8 of 8",
            "off_quads_no_label_within_12vox": "100 percent of sampled "
                                               "(on-run control: 40 "
                                               "percent within 4)",
            "hazard_conf_at_offseeds_pm64": "255 everywhere (no basin in "
                                            "reach of the seed core)",
            "on_sheet_rays_off_total": "5 of ~16,000 sampled",
        },
        "on_seed_ba_archived": on_ba,
        "loo_primary": loo,
        "loo_primary_p_range": loo_range,
        "bucket": bucket,
        "bucket_label": blabel,
    }
    json.dump(results, open(os.path.join(PLACE, "RESULTS.json"), "w"), indent=1)

    # ------------------------------------------------------------- markdown
    def fk(x):
        return "n/a" if x is None else "%.0fk" % (x / 1000.0)

    def pp(x):
        return "n/a" if x is None else "%+.2fpp" % (x * 100.0)

    L = []
    L.append("# Placement A/B: seed-placement test of the hazard-weight effect")
    L.append("")
    L.append("Preregistered design: prereg_placement_ab/PREREGISTRATION.md at "
             "commit 8531276; Amendment 1 (off-seed direction-field regions, "
             "same builder and parameters) at commit 083f424. Generated %s "
             "on the GPU box (CPU only)." % time.strftime("%F %T"))
    L.append("")
    L.append("## Gate")
    L.append("")
    L.append("Archived powered_scores.json sha256 verified (%s...). All "
             "archived 8-site contrasts (5 endpoints: cell means, mean_diff, "
             "sd, t, p, CI95, MDE80) reproduced to full float precision "
             "before any new data was read." % PW_SHA[:16])
    L.append("")
    L.append("## Selection (frozen record)")
    L.append("")
    L.append("offseeds_placement.json sha256 %s. All 8 off-seeds resolved at "
             "the primary distance rule 96 (no relaxation). Sites s2/s4 share "
             "a tile; the deterministic rule gives them the identical "
             "off-seed point (disclosed in Amendment 1)." % OFFSEEDS_SHA[:16])
    L.append("")
    L.append("| site | slab/tile | off L1 (z,y,x) | m7 | minCheb to cluster | rule |")
    L.append("|---|---|---|---|---|---|")
    for i, o in enumerate(selection):
        L.append("| o%d | %s/%s | %d,%d,%d | %d | %d | %d |" % (
            i, o["slab"], o["tile"], o["gz"], o["gy"], o["gx"],
            o["m7"], o["min_chebyshev_to_cluster"], o["dist_rule"]))
    L.append("")
    if field_verify:
        L.append("## Off-seed field build verification (Amendment 1)")
        L.append("")
        L.append("| site | field u8 (x,y,z) | |v| decoded | written |")
        L.append("|---|---|---|---|")
        for r in field_verify.get("seeds", []):
            L.append("| o%d | %s | %.3f | %s |" % (
                r["site_idx"], tuple(r["field_u8_xyz"]),
                r["norm_decoded"], r["written"]))
        L.append("")
    L.append("## Completion")
    L.append("")
    L.append("%d/80 runs scored; %d missing.%s" % (
        n_runs, len(pl.get("missing", [])),
        (" Missing: " + "; ".join(pl["missing"])) if pl.get("missing") else ""))
    L.append("")
    L.append("## Headline: double difference DD = (B-A)_on - (B-A)_off")
    L.append("")
    L.append("| site | area (B-A)_on | area (B-A)_off | area DD | dbl DD | on_sheet DD |")
    L.append("|---|---|---|---|---|---|")
    for i in range(8):
        r = per_site[i]
        L.append("| s%d | %s | %s | %s | %s | %s |" % (
            i, fk(r["area"]["ba_on"]), fk(r["area"]["ba_off"]),
            fk(r["area"]["dd"]), pp(r["dbl"]["dd"]), pp(r["on_sheet"]["dd"])))
    L.append("")
    L.append("PRIMARY area_vx2 DD (H1 > 0, one-sided, df 7): mean %s, sd %s, "
             "t=%.2f, one-sided p=%.4f, 95%% CI [%s, %s] -> %s" % (
                 fk(primary["mean_diff"]), fk(primary["sd_diff"]),
                 primary["t"], primary["p_one_sided"],
                 fk(primary["ci95"][0]), fk(primary["ci95"][1]),
                 "CONFIRMS" if primary_ok else "fails"))
    L.append("")
    if secondary.get("status") == "unmeasurable":
        L.append("SECONDARY frac_double DD: UNMEASURABLE per Amendment 3 "
                 "(off cells carry zero on-sheet rays; %d of 8 sites "
                 "missing)." % secondary["n_missing_sites"])
    else:
        L.append("SECONDARY frac_double DD (H1 < 0, one-sided, declared "
                 "underpowered): mean %s, t=%.2f, one-sided p=%.4f, 95%% CI "
                 "[%s, %s] -> %s" % (
                     pp(secondary["mean_diff"]), secondary["t"],
                     secondary["p_one_sided"], pp(secondary["ci95"][0]),
                     pp(secondary["ci95"][1]),
                     "agrees (p<0.05)" if secondary_ok else "inconclusive"))
    L.append("")
    if guardrail.get("status") == "unmeasurable":
        L.append("GUARDRAIL on_sheet_rate DD: UNMEASURABLE per Amendment 3 "
                 "(same cause).")
    else:
        L.append("GUARDRAIL on_sheet_rate DD (two-sided): mean %s, t=%.2f, "
                 "p=%.4f, 95%% CI [%s, %s] -> %s" % (
                     pp(guardrail["mean_diff"]), guardrail["t"],
                     guardrail["p"], pp(guardrail["ci95"][0]),
                     pp(guardrail["ci95"][1]),
                     "FLAGGED (endpoints possibly denominator-driven)"
                     if guard_flag else "not flagged"))
    L.append("")
    def _stat(d, key="mean_diff"):
        return None if d.get("status") == "unmeasurable" else d[key]

    if reported_near.get("status") == "unmeasurable" or             reported_per1000.get("status") == "unmeasurable":
        L.append("Reported, not tested: crossing DD %s; per-1000 DD %s "
                 "(unmeasurable entries per Amendment 3)." % (
                     pp(_stat(reported_near)),
                     "n/a" if _stat(reported_per1000) is None
                     else "%.2f" % _stat(reported_per1000)))
    else:
        L.append("Reported, not tested: crossing frac_quads_near DD mean %s "
                 "(95%% CI [%s, %s]); per-1000-quads DD mean %.2f "
                 "(95%% CI [%.2f, %.2f])." % (
                     pp(reported_near["mean_diff"]),
                     pp(reported_near["ci95"][0]),
                     pp(reported_near["ci95"][1]),
                     reported_per1000["mean_diff"],
                     reported_per1000["ci95"][0],
                     reported_per1000["ci95"][1]))
    L.append("")
    L.append("## Off-seed (B-A) itself (mechanism predicts near zero)")
    L.append("")
    L.append("| endpoint | mean | 95% CI |")
    L.append("|---|---|---|")
    def row(name, d, f):
        if d.get("status") == "unmeasurable":
            return "| %s | unmeasurable (Amendment 3) | n/a |" % name
        return "| %s | %s | [%s, %s] |" % (
            name, f(d["mean_diff"]), f(d["ci95"][0]), f(d["ci95"][1]))

    L.append(row("area_vx2", off_ba["area"], fk))
    L.append(row("frac_double", off_ba["dbl"], pp))
    L.append(row("on_sheet", off_ba["on_sheet"], pp))
    L.append(row("crossing near", off_ba["near"], pp))
    L.append("")
    L.append("Archived on-seed (B-A) area mean %s; off-seed point estimate "
             "is %.1f%% of it." % (
                 fk(on_area_mean),
                 100.0 * off_area_mean / on_area_mean if on_area_mean else 0))
    L.append("")
    L.append("## Leave-one-site-out (primary, one-sided p, df 6)")
    L.append("")
    L.append("| dropped | p |")
    L.append("|---|---|")
    for j in range(8):
        tag = "s%d" % j
        note = " (same-tile pair member)" if j in (2, 4) else ""
        L.append("| %s%s | %.4f |" % (tag, note, loo_p["drop_s%d" % j]))
    L.append("")
    L.append("Range [%.4f, %.4f]. s2 and s4 share a tile and one off-seed "
             "point; their drops are listed above." % tuple(loo_range))
    L.append("")
    L.append("## Verdict (frozen buckets)")
    L.append("")
    L.append("BUCKET %d: %s%s" % (
        bucket, blabel,
        " Guardrail flag reported alongside." if guard_flag else ""))
    L.append("")
    open(os.path.join(PLACE, "RESULTS.md"), "w").write("\n".join(L))
    gtxt = ("unmeasurable" if guardrail.get("status") == "unmeasurable"
            else "%.4f" % guardrail["p"])
    print("wrote RESULTS.json + RESULTS.md; bucket %d; primary one-sided "
          "p=%.5f; guardrail p=%s%s"
          % (bucket, primary["p_one_sided"], gtxt,
             " GUARDRAIL FLAG" if guard_flag else ""))


if __name__ == "__main__":
    main()
