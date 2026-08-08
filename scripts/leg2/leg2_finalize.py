"""Leg-2 range check, readout stage.

Everything here comes off the retained 2000-bin histograms, so no GPU and no
re-threshold. The question is whether `surface_recto_3dunet` carries signal on
TAUIL's 892 at his geometry, and the estimator that answers it without depending
on a threshold choice is AUC. AUC 0.5 is chance. Dice and precision are also
reported, at the fixed grid and at a per-volume oracle threshold, because a model
can have real AUC and still be useless at every usable operating point.

Leg 1 (m7) is scored the same way from its own retained histograms on the same
volumes, so the comparison is head to head rather than against a remembered number.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path

import numpy as np

LEG2 = Path("/mnt/vesuvius/experiments/leg2/rows.jsonl")
LEG1_ROWS = Path("/mnt/vesuvius/experiments/shell892/rows.jsonl")
LEG1_HIST = Path("/mnt/vesuvius/experiments/shell892/pred_histograms.jsonl")
OUT = Path("/mnt/vesuvius/experiments/leg2")
GRID = ("0.05", "0.10", "0.20", "0.30", "0.50", "0.70")


def loadl(p: Path):
    for ln in p.read_text().splitlines():
        if ln.strip():
            yield json.loads(ln)


def auc_from_hist(h_sheet, h_scored):
    """AUC of P(surface) separating labelled sheet from scored non-sheet, computed
    from binned counts with ties handled by the mid-rank convention. Returns None
    when either class is empty."""
    pos = np.asarray(h_sheet, dtype=np.float64)
    neg = np.asarray(h_scored, dtype=np.float64) - pos
    neg = np.maximum(neg, 0.0)           # binning cannot make sheet exceed scored
    P, N = pos.sum(), neg.sum()
    if P <= 0 or N <= 0:
        return None
    # negatives strictly below each bin, plus half the ties in the bin
    below = np.concatenate([[0.0], np.cumsum(neg)[:-1]])
    return float((pos * (below + 0.5 * neg)).sum() / (P * N))


def curve(h_sheet, h_scored):
    """Recall, precision and dice at every bin edge, from the tail sums."""
    pos = np.asarray(h_sheet, dtype=np.float64)
    tot = np.asarray(h_scored, dtype=np.float64)
    P = pos.sum()
    tail_p = np.concatenate([np.cumsum(pos[::-1])[::-1], [0.0]])   # TP at each edge
    tail_t = np.concatenate([np.cumsum(tot[::-1])[::-1], [0.0]])   # predicted at edge
    if P <= 0:
        return None
    rec = tail_p / P
    with np.errstate(invalid="ignore", divide="ignore"):
        prec = np.where(tail_t > 0, tail_p / np.maximum(tail_t, 1e-9), np.nan)
        dice = np.where((tail_t + P) > 0, 2.0 * tail_p / (tail_t + P), np.nan)
    return rec, prec, dice, tail_p, tail_t, P


def per_volume(h_sheet, h_scored, nbin: int):
    c = curve(h_sheet, h_scored)
    if c is None:
        return None
    rec, prec, dice, tail_p, tail_t, P = c
    i = int(np.nanargmax(dice))
    # precision ceiling: best precision at any threshold that still recovers a
    # non-trivial slice of the sheet, so a single lucky voxel cannot set it
    ok = rec >= 0.05
    pc = float(np.nanmax(prec[ok])) if ok.any() else None
    return {
        "auc": auc_from_hist(h_sheet, h_scored),
        "dice_best": float(dice[i]),
        "dice_best_threshold": round(i / nbin, 5),
        "recall_at_dice_best": float(rec[i]),
        "precision_at_dice_best": float(prec[i]) if not math.isnan(prec[i]) else None,
        "precision_ceiling_at_recall_05": pc,
    }


def pooled(rows, key_sheet="hist_sheet", key_scored="hist_scored"):
    """Pooled dice: aggregate counts across volumes, then one dice curve. This is the
    July estimator, kept so 0.23 is comparable to whatever comes out here."""
    ns = None
    for r in rows:
        if r.get(key_sheet) is None:
            continue
        s = np.asarray(r[key_sheet], dtype=np.float64)
        t = np.asarray(r[key_scored], dtype=np.float64)
        ns = (s, t) if ns is None else (ns[0] + s, ns[1] + t)
    if ns is None:
        return None
    c = curve(ns[0], ns[1])
    rec, prec, dice, tail_p, tail_t, P = c
    i = int(np.nanargmax(dice))
    nbin = len(ns[0])
    return {
        "auc": auc_from_hist(ns[0], ns[1]),
        "pooled_dice_best": float(dice[i]),
        "pooled_dice_best_threshold": round(i / nbin, 5),
        "pooled_recall_at_best": float(rec[i]),
        "pooled_precision_at_best": float(prec[i]),
        "pooled_dice_at_0.50": float(dice[int(0.50 * nbin)]),
        "pooled_base_rate": float(P / ns[1].sum()),
    }


def budget_point(h_sheet, h_scored, budget: float = 0.12):
    """TAUIL's matched-budget operating point: threshold each volume so that at most
    `budget` of its scored voxels are predicted positive, then read precision lift
    against that volume's own base rate. His registered floor for leg 2 was a lift of
    2.0x, and lift 1.0x is chance by construction, so this is the estimator his floor
    was written in and the one his 20-volume result is quoted in."""
    pos = np.asarray(h_sheet, dtype=np.float64)
    tot = np.asarray(h_scored, dtype=np.float64)
    n_s, n_sh = tot.sum(), pos.sum()
    if not n_s or not n_sh:
        return None
    tail_t = np.concatenate([np.cumsum(tot[::-1])[::-1], [0.0]])
    tail_p = np.concatenate([np.cumsum(pos[::-1])[::-1], [0.0]])
    i = int(np.argmax(tail_t / n_s <= budget))
    npred, ntp = tail_t[i], tail_p[i]
    if not npred:
        return None
    prec, base = ntp / npred, n_sh / n_s
    return {
        "budget": budget,
        "budget_threshold": round(i / len(tot), 5),
        "budget_pred_positive_fraction": float(npred / n_s),
        "budget_recall": float(ntp / n_sh),
        "budget_precision": float(prec),
        "budget_precision_lift": float(prec / base),
    }


def sweep_lift(row: dict, t: str):
    """Precision, dice and lift at a fixed threshold, from the exact counts rather
    than the histogram."""
    sw = (row.get("threshold_sweep") or {}).get(t)
    if not sw or not row.get("n_sheet") or not row.get("n_scored"):
        return None
    npred, ntp = sw["n_pred_scored"], sw["n_tp"]
    n_sh, n_s = row["n_sheet"], row["n_scored"]
    base = n_sh / n_s
    out = {"recall": ntp / n_sh, "pred_positive_fraction": npred / n_s,
           "dice": (2.0 * ntp / (npred + n_sh)) if (npred + n_sh) else None}
    out["precision"] = (ntp / npred) if npred else None
    out["precision_lift"] = (out["precision"] / base) if (npred and base) else None
    return out


def med(v):
    v = [x for x in v if x is not None]
    return round(float(np.median(v)), 5) if v else None


def block(rows, tag):
    ok = [r for r in rows if r.get("status") == "ok"]
    per = [r["leg2"] for r in ok if r.get("leg2")]
    b = {
        "tag": tag,
        "n_rows": len(rows),
        "n_ok": len(ok),
        "median_denominator": "n_ok",
        "median_auc": med([x["auc"] for x in per]),
        "frac_auc_above_0.5": (round(sum(1 for x in per if (x["auc"] or 0) > 0.5)
                                    / len(per), 4) if per else None),
        "median_dice_best": med([x["dice_best"] for x in per]),
        "median_dice_best_threshold": med([x["dice_best_threshold"] for x in per]),
        "median_precision_ceiling": med([x["precision_ceiling_at_recall_05"]
                                         for x in per]),
        "median_base_rate": med([r.get("base_rate") for r in ok]),
        "median_p_sheet_minus_p_scored": med(
            [(r["p_mean_sheet"] - r["p_mean_scored"]) for r in ok
             if r.get("p_mean_sheet") is not None]),
        "pooled": pooled(ok),
    }
    # his registered estimator: matched-budget precision lift against a 2.0x floor
    bp = [budget_point(r["hist_sheet"], r["hist_scored"]) for r in ok
          if r.get("hist_sheet")]
    bp = [x for x in bp if x]
    b["matched_budget_0.12"] = {
        "n": len(bp),
        "median_budget_recall": med([x["budget_recall"] for x in bp]),
        "median_budget_precision": med([x["budget_precision"] for x in bp]),
        "median_budget_precision_lift": med([x["budget_precision_lift"] for x in bp]),
        "frac_lift_at_or_above_2.0": (round(sum(
            1 for x in bp if x["budget_precision_lift"] >= 2.0) / len(bp), 4)
            if bp else None),
        "registered_floor": "his preregistered leg-2 floor was lift >= 2.0x at "
                            "matched budget; 1.0x is chance by construction",
    }
    lifts = [sweep_lift(r, "0.20") for r in ok]
    lifts = [x for x in lifts if x]
    b["at_threshold_0.20"] = {
        "n": len(lifts),
        "median_recall": med([x["recall"] for x in lifts]),
        "median_precision": med([x["precision"] for x in lifts]),
        "median_precision_lift": med([x["precision_lift"] for x in lifts]),
        "median_pred_positive_fraction": med(
            [x["pred_positive_fraction"] for x in lifts]),
    }
    # dice at the fixed grid, from the exact counts rather than the histogram
    grid = {}
    for t in GRID:
        ds = []
        for r in ok:
            sw = (r.get("threshold_sweep") or {}).get(t)
            if not sw:
                continue
            tp, npred, nsh = sw["n_tp"], sw["n_pred_scored"], r["n_sheet"]
            if npred + nsh:
                ds.append(2.0 * tp / (npred + nsh))
        grid[t] = med(ds)
    b["median_dice_by_threshold"] = grid
    return b


def main() -> None:
    rows = list(loadl(LEG2))
    nbin = 2000
    for r in rows:
        if r.get("hist_sheet") and r.get("n_sheet"):
            r["leg2"] = per_volume(r["hist_sheet"], r["hist_scored"], nbin)
        else:
            r["leg2"] = None

    # ---- leg 1 on the same volumes, from its own retained histograms ----
    l1h = {h["sample"]: h for h in loadl(LEG1_HIST)}
    l1r = {r["sample"]: r for r in loadl(LEG1_ROWS)}
    leg1 = []
    for nm, h in l1h.items():
        pv = per_volume(h["hist_sheet"], h["hist_scored"], len(h["hist_scored"]))
        if pv is None:
            continue
        r = l1r.get(nm, {})
        leg1.append({"sample": nm, "status": r.get("status"),
                     "population": r.get("population"),
                     "n_sheet": r.get("n_sheet"), "base_rate": r.get("base_rate"),
                     "hist_sheet": h["hist_sheet"], "hist_scored": h["hist_scored"],
                     "leg2": pv})

    loc = [r for r in rows if r.get("population") == "located"]
    non = [r for r in rows if r.get("population") == "nonlocated"]
    l1_ok = [r for r in leg1 if r.get("status") == "ok"]

    # paired comparison on volumes where both legs are scorable
    by2 = {r["sample"]: r for r in rows if r.get("leg2")}
    pairs = [(by2[r["sample"]]["leg2"], r["leg2"]) for r in l1_ok
             if r["sample"] in by2]
    d_auc = [b["auc"] - a["auc"] for a, b in pairs
             if a["auc"] is not None and b["auc"] is not None]
    d_dice = [b["dice_best"] - a["dice_best"] for a, b in pairs]

    summary = {
        "what_this_is": (
            "Leg-2 range check. surface_recto_3dunet over the 892 public volumes at "
            "TAUIL's geometry, under the checkpoint's own declared z-score "
            "normalisation, one 256^3 forward pass per volume. Scored on the inner "
            "128^3 with label class 2 excluded, the same voxels leg 1 used."),
        "provenance": {
            "checkpoint": "/mnt/vesuvius/models/surface_recto_3dunet/"
                          "checkpoint_inference_ready.pth",
            "normalization_scheme": "zscore, intensity_properties empty, so CT "
                                    "normalisation is not defined for this checkpoint",
            "patch": [256, 256, 256],
            "geometry": "centred 256^3 block, TRIM=64, inner 128^3, class 2 excluded",
            "nbin": nbin,
            "note": "geometry cross-checked against leg 1 row by row on n_sheet and "
                    "n_scored; see gates.label_side_vs_leg1",
        },
        "populations": {
            "all": block(rows, "all 892"),
            "located": block(loc, "located population"),
            "nonlocated": block(non, "non-located population"),
        },
        "leg1_m7_same_estimators": {
            "tag": "m7, from its own retained histograms, same volumes",
            "n": len(l1_ok),
            "median_auc": med([r["leg2"]["auc"] for r in l1_ok]),
            "median_dice_best": med([r["leg2"]["dice_best"] for r in l1_ok]),
            "median_precision_ceiling": med(
                [r["leg2"]["precision_ceiling_at_recall_05"] for r in l1_ok]),
            "pooled": pooled(l1_ok),
        },
        "paired_leg1_minus_leg2": {
            "what": "m7 minus surface_recto_3dunet, per volume, positive means m7 "
                    "is higher",
            "n_pairs": len(pairs),
            "median_delta_auc": med(d_auc),
            "median_delta_dice_best": med(d_dice),
            "n_leg1_auc_higher": sum(1 for x in d_auc if x > 0),
            "n_leg1_dice_higher": sum(1 for x in d_dice if x > 0),
        },
    }

    # ---- gates ----
    mismatch = []
    for r in rows:
        a = l1r.get(r["sample"])
        if a and (a.get("n_sheet") != r.get("n_sheet")
                  or a.get("n_scored") != r.get("n_scored")):
            mismatch.append(r["sample"])
    summary["gates"] = {
        "label_side_vs_leg1": {
            "n_checked": len(rows),
            "fields": ["n_sheet", "n_scored"],
            "n_mismatch": len(mismatch),
            "mismatches": mismatch[:10],
            "what": "the two legs must see the same voxels, or no per-volume "
                    "comparison between them is meaningful",
        },
        "sweep_vs_histogram": None,
    }
    # the fixed-grid counts and the histogram must agree at 0.50
    worst = 0.0
    for r in rows:
        if not r.get("hist_sheet"):
            continue
        sw = (r.get("threshold_sweep") or {}).get("0.50")
        if not sw:
            continue
        h = np.asarray(r["hist_scored"], dtype=np.float64)
        from_hist = float(h[int(0.50 * nbin):].sum())
        if sw["n_pred_scored"]:
            worst = max(worst, abs(from_hist - sw["n_pred_scored"])
                        / sw["n_pred_scored"])
    summary["gates"]["sweep_vs_histogram"] = {
        "what": "predicted-positive count at 0.50 from the histogram against the "
                "exact count, per volume",
        "max_rel_delta": round(worst, 8),
    }

    slim = [{k: v for k, v in r.items()
             if k not in ("hist_scored", "hist_sheet")} for r in rows]
    (OUT / "rows_slim.jsonl").write_text(
        "".join(json.dumps(r) + "\n" for r in slim))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: summary[k] for k in
                      ("leg1_m7_same_estimators", "paired_leg1_minus_leg2", "gates")},
                     indent=1))
    for k in ("all", "located", "nonlocated"):
        b = summary["populations"][k]
        print(f"\n{k}: n_ok={b['n_ok']} auc={b['median_auc']} "
              f"frac_auc>0.5={b['frac_auc_above_0.5']} "
              f"dice_best={b['median_dice_best']} "
              f"prec_ceiling={b['median_precision_ceiling']} "
              f"base={b['median_base_rate']} sep={b['median_p_sheet_minus_p_scored']}")
        print("   pooled:", json.dumps(b["pooled"]))
        print("   dice by threshold:", json.dumps(b["median_dice_by_threshold"]))


if __name__ == "__main__":
    main()
