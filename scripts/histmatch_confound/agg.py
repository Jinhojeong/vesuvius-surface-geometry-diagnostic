"""Aggregate the four-arm histogram-matching x normalization 2x2.

Reads the jsonl written by run_4arm.py and emits confound.json plus a printed
table. Nothing here re-runs the model; every number is a function of the shipped
jsonl, control_his60.jsonl (reproduction check) and the plans constants.

Definitions used throughout:
  A = original volume, plans CTNormalization
  B = original volume, instance z-score      (the villa wrapper path)
  C = histogram-matched copy, plans CTNormalization
  D = histogram-matched copy, instance z-score

  wrapper-path recovery      dBD = D - B   (what a bench run through the wrapper
                                            would report as "recovery from matching")
  correctly-normalised level dAC = C - A   (what the model actually does with a
                                            level-matched input)
  confound                   dBD - dAC     paired per volume, then median

Medians are per-volume medians (median over volumes of the per-volume statistic),
which is TAUIL-Abd-Elilah's convention in surface_bench_m7.json and the one every
earlier artifact in this line uses.
"""
from __future__ import annotations
import argparse, json, statistics as st
from pathlib import Path

import numpy as np
from scipy import stats

D = Path("/mnt/vesuvius/experiments/histmatch_confound")
SS = Path("/mnt/vesuvius/experiments/shell_split")
ARMS = ["A_orig_ctplans", "B_orig_instzs", "C_match_ctplans", "D_match_instzs"]
AFFINE_ARMS = ["E_affine_ctplans", "F_affine_instzs"]
ALL_ARMS = ARMS + AFFINE_ARMS
METRICS = ["recall", "precision", "ppf"]


def load(*paths: Path) -> dict:
    """Merge one or more jsonl files keyed (sample, arm). Rows with status
    'degenerate' (no labelled sheet in the scored crop) carry recall None and are
    dropped by every downstream median, but they stay in the file."""
    out: dict = {}
    for p in paths:
        if not p.exists():
            continue
        for ln in p.read_text().splitlines():
            if not ln.strip():
                continue
            r = json.loads(ln)
            out.setdefault(r["sample"], {})[r["arm"]] = r
    return out


def med(v):
    return st.median(v) if v else None


def paired(byv: dict, arm1: str, arm2: str, metric: str):
    """arm2 - arm1, per volume, only where both arms scored ok."""
    out = []
    for s, d in sorted(byv.items()):
        if arm1 in d and arm2 in d and d[arm1].get(metric) is not None \
                and d[arm2].get(metric) is not None:
            out.append((s, d[arm2][metric] - d[arm1][metric]))
    return out


def wilc(v):
    if len(v) < 6 or all(x == 0 for x in v):
        return None
    return float(stats.wilcoxon(v).pvalue)


def arm_table(byv: dict, arms: list[str]) -> dict:
    t = {}
    for arm in arms:
        vals = {m: [d[arm][m] for d in byv.values()
                    if arm in d and d[arm].get(m) is not None] for m in METRICS}
        if not vals["recall"]:
            continue
        t[arm] = {"n": len(vals["recall"]),
                  **{f"median_{m}": round(med(vals[m]), 4) for m in METRICS},
                  **{f"q1_{m}": round(float(np.percentile(vals[m], 25)), 4) for m in METRICS},
                  **{f"q3_{m}": round(float(np.percentile(vals[m], 75)), 4) for m in METRICS}}
    return t


def delta_block(byv: dict, a1: str, a2: str) -> dict:
    out = {}
    for m in METRICS:
        pr = paired(byv, a1, a2, m)
        v = [x for _, x in pr]
        if not v:
            continue
        out[m] = {"n": len(v), "median": round(med(v), 4),
                  "q1": round(float(np.percentile(v, 25)), 4),
                  "q3": round(float(np.percentile(v, 75)), 4),
                  "n_positive": int(sum(x > 0 for x in v)),
                  "wilcoxon_p": wilc(v)}
    return out


def confound_block(byv: dict) -> dict:
    """(D - B) - (C - A), paired within volume."""
    out = {}
    for m in METRICS:
        bd = dict(paired(byv, "B_orig_instzs", "D_match_instzs", m))
        ac = dict(paired(byv, "A_orig_ctplans", "C_match_ctplans", m))
        common = sorted(set(bd) & set(ac))
        v = [bd[s] - ac[s] for s in common]
        if not v:
            continue
        out[m] = {"n": len(v),
                  "median_dBD": round(med([bd[s] for s in common]), 4),
                  "median_dAC": round(med([ac[s] for s in common]), 4),
                  "median_confound_paired": round(med(v), 4),
                  "confound_of_medians": round(med([bd[s] for s in common])
                                               - med([ac[s] for s in common]), 4),
                  "q1": round(float(np.percentile(v, 25)), 4),
                  "q3": round(float(np.percentile(v, 75)), 4),
                  "n_positive": int(sum(x > 0 for x in v)),
                  "wilcoxon_p": wilc(v)}
    return out


def sigma_block(byv: dict, plans_mean: float, plans_std: float) -> dict:
    o = [d[a]["orig_blk_mean"] for d in byv.values() for a in ARMS if a in d][:0]
    om = [next(iter(d.values()))["orig_blk_mean"] for d in byv.values()]
    osd = [next(iter(d.values()))["orig_blk_std"] for d in byv.values()]
    mm = [r["match_blk_mean"] for d in byv.values() for r in d.values()
          if "match_blk_mean" in r]
    msd = [r["match_blk_std"] for d in byv.values() for r in d.values()
           if "match_blk_std" in r]
    # sigma = (block mean - plans mean) / plans std, i.e. how far the volume's own
    # centre sits from the intensity fingerprint the plans normalization assumes.
    # scale = block std / plans std, the second half of the same discrepancy: instance
    # z-score divides by the block's own std, plans divides by 47.7438.
    out = {"n": len(om),
           "orig_blk_mean_median": round(med(om), 2),
           "orig_blk_std_median": round(med(osd), 2),
           "orig_sigma_of_median": round((med(om) - plans_mean) / plans_std, 3),
           "orig_sigma_per_volume_median": round(
               med([(x - plans_mean) / plans_std for x in om]), 3),
           "orig_scale_ratio_median": round(med([x / plans_std for x in osd]), 3)}
    if mm:
        out.update({"match_blk_mean_median": round(med(mm), 2),
                    "match_blk_std_median": round(med(msd), 2),
                    "match_sigma_of_median": round((med(mm) - plans_mean) / plans_std, 3),
                    "match_sigma_per_volume_median": round(
                        med([(x - plans_mean) / plans_std for x in mm]), 3),
                    "match_scale_ratio_median": round(med([x / plans_std for x in msd]), 3)})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(D / "confound.json"))
    a = ap.parse_args()

    plans = json.loads(Path("/mnt/vesuvius/models/surface_m7_nnunet/plans.json").read_text())
    ip = plans["foreground_intensity_properties_per_channel"]["0"]
    pm, ps = float(ip["mean"]), float(ip["std"])

    loc = load(D / "loc60_4arm.jsonl", D / "loc60_affine.jsonl")
    abl16 = json.loads((D / "abl16.json").read_text())
    loc16 = {s: v for s, v in loc.items() if s in abl16}
    nonloc = load(D / "nonloc60_2arm.jsonl")
    patch16 = load(D / "patch16_4arm.jsonl")

    res: dict = {"plans": {"mean": pm, "std": ps,
                           "percentile_00_5": ip["percentile_00_5"],
                           "percentile_99_5": ip["percentile_99_5"]},
                 "arm_definitions": {
                     "A_orig_ctplans": "original volume, plans CTNormalization",
                     "B_orig_instzs": "original volume, instance z-score (villa wrapper path)",
                     "C_match_ctplans": "histogram-matched copy, plans CTNormalization",
                     "D_match_instzs": "histogram-matched copy, instance z-score"}}

    for tag, byv, arms in (("located16", loc16, ALL_ARMS),
                           ("located60", loc, ALL_ARMS),
                           ("nonlocated60", nonloc, ARMS[:2]),
                           ("located16_patchnorm", patch16, ARMS)):
        if not byv:
            continue
        blk = {"arms": arm_table(byv, arms),
               "sigma": sigma_block(byv, pm, ps),
               "samples": sorted(byv)}
        if len(arms) >= 4:
            blk["deltas"] = {
                "dAC_plans_matched_response": delta_block(byv, "A_orig_ctplans", "C_match_ctplans"),
                "dBD_wrapper_recovery": delta_block(byv, "B_orig_instzs", "D_match_instzs"),
                "dAB_norm_penalty_orig": delta_block(byv, "A_orig_ctplans", "B_orig_instzs"),
                "dCD_norm_penalty_matched": delta_block(byv, "C_match_ctplans", "D_match_instzs"),
                "dAE_plans_pure_level_scale": delta_block(byv, "A_orig_ctplans", "E_affine_ctplans"),
                "dEC_plans_lut_nonlinearity": delta_block(byv, "E_affine_ctplans", "C_match_ctplans"),
                "dBF_affine_invariance_selftest": delta_block(byv, "B_orig_instzs", "F_affine_instzs")}
            blk["confound"] = confound_block(byv)
            lv = [r["lut_distinct_levels_on_block"] for d in byv.values()
                  for r in d.values() if "lut_distinct_levels_on_block" in r]
            bl = [r["block_distinct_levels"] for d in byv.values()
                  for r in d.values() if "block_distinct_levels" in r]
            if lv:
                blk["lut_compression"] = {
                    "block_distinct_levels_median": med(bl),
                    "lut_distinct_levels_on_block_median": med(lv),
                    "levels_collapsed_median": med(bl) - med(lv),
                    "surviving_fraction_median": round(med(lv) / med(bl), 4)}
        res[tag] = blk

    # ---- the gap, both ends, in both normalization worlds ----
    if nonloc and loc:
        g = {}
        for arm, world in (("A_orig_ctplans", "plans_CTNormalization"),
                           ("B_orig_instzs", "instance_zscore_wrapper_path")):
            def col(byv, a):  # degenerate rows carry recall None and drop out here
                return [d[a]["recall"] for d in byv.values()
                        if a in d and d[a]["recall"] is not None]
            lv, nv = col(loc, arm), col(nonloc, arm)
            dv, cv = col(loc, "D_match_instzs"), col(loc, "C_match_ctplans")
            matched = med(dv if arm == "B_orig_instzs" else cv)
            gap = med(nv) - med(lv)
            g[world] = {"located_n": len(lv), "nonlocated_n": len(nv),
                        "located_median_recall": round(med(lv), 4),
                        "nonlocated_median_recall": round(med(nv), 4),
                        "gap": round(gap, 4),
                        "matched_located_median_recall": round(matched, 4),
                        "recovery": round(matched - med(lv), 4),
                        "recovered_share_of_gap": round((matched - med(lv)) / gap, 4)
                        if abs(gap) > 1e-9 else None}
        res["gap"] = g

    # ---- reproduction checks ----
    chk = {}
    old = {}
    for ln in (SS / "control_his60.jsonl").read_text().splitlines():
        if ln.strip():
            r = json.loads(ln)
            d = r["desc"]
            tp = d["n_pred_scored"] - d["n_fp"]
            old[r["sample"]] = {"recall": tp / d["n_sheet"], "precision": tp / d["n_pred_scored"],
                                "ppf": d["n_pred_scored"] / d["n_scored"]}
    dev = {m: [abs(loc[s]["A_orig_ctplans"][m] - old[s][m]) for s in loc
               if s in old and "A_orig_ctplans" in loc[s]] for m in METRICS}
    chk["armA_vs_control_his60"] = {"n": len(dev["recall"]),
                                    **{f"max_abs_dev_{m}": float(max(dev[m])) for m in METRICS}}
    na = json.loads((SS / "norm_ablation16.json").read_text())
    if isinstance(na, dict):  # the file gained a {legend, rows} wrapper on 08-07
        na = na["rows"]
    nad = {r["sample"]: r for r in na}
    for arm, col in (("A_orig_ctplans", "ct_plans"), ("B_orig_instzs", "instance_zscore")):
        dv = [max(abs(round(loc16[s][arm][m], 4) - nad[s][col][i])
                  for i, m in enumerate(METRICS))
              for s in loc16 if s in nad and arm in loc16[s]]
        chk[f"{arm}_vs_norm_ablation16"] = {"n": len(dv), "max_abs_dev": round(max(dv), 5)}
    if patch16:
        dv = [abs(patch16[s][a]["recall"] - loc16[s][a]["recall"])
              for s in patch16 for a in patch16[s]
              if s in loc16 and a in loc16[s] and "instzs" in a]
        if dv:
            chk["instance_zscore_patch_vs_volume"] = {
                "n": len(dv), "max_abs_recall_dev": round(max(dv), 5),
                "median_abs_recall_dev": round(med(dv), 5)}
    res["reproduction_checks"] = chk

    Path(a.out).write_text(json.dumps(res, indent=1))

    # per-volume table, so every cohort median has a visible row behind it
    cols = ["cohort", "sample", "orig_blk_mean", "orig_blk_std", "match_blk_mean",
            "match_blk_std"]
    for arm in ALL_ARMS:
        cols += [f"{arm}_{m}" for m in METRICS]
    cols += ["dAC_recall", "dBD_recall", "confound_recall"]
    lines = [",".join(cols)]
    for tag, byv in (("located60", loc), ("nonlocated60", nonloc),
                     ("located16_patchnorm", patch16)):
        for s in sorted(byv):
            d = byv[s]
            any_r = next(iter(d.values()))
            row = [tag, s, any_r.get("orig_blk_mean"), any_r.get("orig_blk_std"),
                   any_r.get("match_blk_mean"), any_r.get("match_blk_std")]
            for arm in ALL_ARMS:
                row += [d[arm][m] if arm in d else None for m in METRICS]
            def g(arm, m="recall"):
                return d[arm][m] if arm in d else None
            dac = (g("C_match_ctplans") - g("A_orig_ctplans")
                   if g("C_match_ctplans") is not None and g("A_orig_ctplans") is not None
                   else None)
            dbd = (g("D_match_instzs") - g("B_orig_instzs")
                   if g("D_match_instzs") is not None and g("B_orig_instzs") is not None
                   else None)
            row += [dac, dbd, (dbd - dac) if dac is not None and dbd is not None else None]
            lines.append(",".join("" if x is None else
                                  (f"{x:.6f}" if isinstance(x, float) else str(x))
                                  for x in row))
    Path(a.out).with_name("per_volume.csv").write_text("\n".join(lines) + "\n")

    # ---- print ----
    for tag in ("located16", "located60", "nonlocated60", "located16_patchnorm"):
        if tag not in res:
            continue
        print(f"\n### {tag} ###")
        for arm, v in res[tag]["arms"].items():
            print("  %-16s n=%-3d recall %.4f  precision %.4f  ppf %.4f"
                  % (arm, v["n"], v["median_recall"], v["median_precision"], v["median_ppf"]))
        s = res[tag]["sigma"]
        print("  sigma: orig blk-mean median %.2f (%+.3f sigma)%s"
              % (s["orig_blk_mean_median"], s["orig_sigma_of_median"],
                 "  matched %.2f (%+.3f sigma)" % (s["match_blk_mean_median"],
                                                   s["match_sigma_of_median"])
                 if "match_blk_mean_median" in s else ""))
        if "lut_compression" in res[tag]:
            lc = res[tag]["lut_compression"]
            print("  lut: %g of %g distinct grey levels on the block survive the map "
                  "(%.1f%%)" % (lc["lut_distinct_levels_on_block_median"],
                                lc["block_distinct_levels_median"],
                                100 * lc["surviving_fraction_median"]))
        if "confound" in res[tag]:
            for m in METRICS:
                c = res[tag]["confound"].get(m)
                if not c:
                    continue
                print("  %-9s dBD %+.4f  dAC %+.4f  confound %+.4f (paired med, "
                      "n=%d, %d>0, p=%s)"
                      % (m, c["median_dBD"], c["median_dAC"],
                         c["median_confound_paired"], c["n"], c["n_positive"],
                         "%.3g" % c["wilcoxon_p"] if c["wilcoxon_p"] is not None else "na"))
            for k, v in res[tag]["deltas"].items():
                if "recall" in v:
                    r = v["recall"]
                    print("    %-32s recall %+.4f (n=%d, %d>0, p=%s)"
                          % (k, r["median"], r["n"], r["n_positive"],
                             "%.3g" % r["wilcoxon_p"] if r["wilcoxon_p"] is not None else "na"))
    if "gap" in res:
        print("\n### gap, both ends ###")
        for w, v in res["gap"].items():
            print("  %-32s located %.4f  nonlocated %.4f  gap %+.4f  matched %.4f  "
                  "recovery %+.4f  share %s"
                  % (w, v["located_median_recall"], v["nonlocated_median_recall"],
                     v["gap"], v["matched_located_median_recall"], v["recovery"],
                     "%.3f" % v["recovered_share_of_gap"]
                     if v["recovered_share_of_gap"] is not None else "na"))
    print("\n### reproduction checks ###")
    print(json.dumps(res["reproduction_checks"], indent=1))
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
