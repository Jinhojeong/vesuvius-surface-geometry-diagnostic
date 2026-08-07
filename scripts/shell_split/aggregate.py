"""Turn per-volume rows into the statistics his comment quotes, plus the ones it doesn't.

Conventions are named explicitly because the two obvious ones disagree:

  per-volume median  -- compute the statistic inside each volume, then take the median
                        across volumes. This is what m7_margin_fp.py does for the shell
                        enrichments and what reproduces his 0.7554 obs/pred exactly.
  pooled             -- sum the counts across volumes first, then form one ratio. This
                        weights big-sheet volumes more.

obs/pred: ordinary least squares of log(enrichment) on shell index k for k in 2..5,
extrapolated to k = 1, then observed shell 1 divided by that. Fitting his shipped rows
this way returns 0.7554, which is the number in his comment, so this is his convention.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.stats import wilcoxon, pearsonr, spearmanr


def loadl(p):
    return [json.loads(l) for l in Path(p).read_text().splitlines() if l.strip()]


def pred_shell1(sh: dict) -> float | None:
    ks = [k for k in range(2, 6)
          if sh.get(f"shell_{k}") is not None and sh[f"shell_{k}"] > 0]
    if len(ks) < 2:
        return None
    x = np.array(ks, dtype=float)
    y = np.log([sh[f"shell_{k}"] for k in ks])
    b, a = np.polyfit(x, y, 1)
    return float(np.exp(a + b))


def med(v):
    return round(float(np.median(v)), 4) if len(v) else None


def summarize(rows: list[dict], tag: str) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    out = {"tag": tag, "n_rows": len(rows), "n_ok": len(ok)}

    shells = {k: [r["shells"][f"shell_{k}"] for r in ok
                  if r.get("shells", {}).get(f"shell_{k}") is not None] for k in range(1, 6)}
    out["shell_enrichment_median"] = {f"shell_{k}": med(shells[k]) for k in range(1, 6)}
    out["shell_enrichment_n"] = {f"shell_{k}": len(shells[k]) for k in range(1, 6)}

    op, s1 = [], []
    for r in ok:
        p = pred_shell1(r.get("shells", {}))
        o = r.get("shells", {}).get("shell_1")
        if p and o is not None:
            op.append(o / p); s1.append(o)
    op = np.array(op)
    out["obs_over_pred"] = {
        "n": len(op), "median": med(op), "mean": round(float(op.mean()), 4) if len(op) else None,
        "q10": round(float(np.quantile(op, .1)), 4) if len(op) else None,
        "q90": round(float(np.quantile(op, .9)), 4) if len(op) else None,
        "frac_above_1": round(float(np.mean(op > 1)), 4) if len(op) else None,
    }
    if len(op) >= 6:
        w = wilcoxon(op - 1.0)
        out["obs_over_pred"]["wilcoxon_stat"] = float(w.statistic)
        out["obs_over_pred"]["wilcoxon_p"] = float(w.pvalue)

    en = [r["enrichment"] for r in ok]
    out["margin_enrichment"] = {
        "median": med(en), "q10": round(float(np.quantile(en, .1)), 3) if en else None,
        "q90": round(float(np.quantile(en, .9)), 3) if en else None,
        "frac_above_1": round(float(np.mean(np.array(en) > 1)), 4) if en else None,
        "median_fp_share_in_margin": med([r["fp_share_in_margin"] for r in ok]),
        "median_margin_share_of_nonsheet": med([r["margin_share_of_nonsheet"] for r in ok]),
    }

    # FP load
    load_scored = [r["desc"]["n_fp"] / r["desc"]["n_scored"] for r in ok]
    load_ns = [r["desc"]["n_fp"] / r["desc"]["n_nonsheet_scored"] for r in ok]
    tp = [(r["desc"]["n_pred_scored"] - r["desc"]["n_fp"]) for r in ok]
    out["fp_load"] = {
        "median_fp_frac_of_scored": med(load_scored),
        "median_fp_frac_of_nonsheet_scored": med(load_ns),
        "pooled_fp_frac_of_scored": round(sum(r["desc"]["n_fp"] for r in ok)
                                          / sum(r["desc"]["n_scored"] for r in ok), 5),
        "pooled_fp_frac_of_nonsheet_scored": round(sum(r["desc"]["n_fp"] for r in ok)
                                                   / sum(r["desc"]["n_nonsheet_scored"] for r in ok), 5),
        "median_recall": med([t / max(1, r["desc"]["n_sheet"]) for t, r in zip(tp, ok)]),
        "median_precision": med([t / max(1, r["desc"]["n_pred_scored"]) for t, r in zip(tp, ok)]),
        "median_base_rate": med([r["desc"]["n_sheet"] / r["desc"]["n_scored"] for r in ok]),
    }

    # descriptive: FP mass share by shell, and beyond-k volume share
    per_vol = []
    for r in ok:
        d = r["desc"]
        if d["n_fp"] == 0:
            continue
        per_vol.append([f / d["n_fp"] for f in d["shell_fp"]])
    pv = np.array(per_vol) if per_vol else np.zeros((0, 5))
    tot_fp = sum(r["desc"]["n_fp"] for r in ok)
    pooled = [sum(r["desc"]["shell_fp"][i] for r in ok) / tot_fp for i in range(5)] if tot_fp else []
    out["fp_mass_share_by_shell"] = {
        "per_volume_median_pct": [round(float(np.median(pv[:, i])) * 100, 2) for i in range(5)] if len(pv) else None,
        "pooled_pct": [round(p * 100, 2) for p in pooled],
        "per_volume_median_pct_beyond5": round(float(np.median(1 - pv.sum(axis=1))) * 100, 2) if len(pv) else None,
        "pooled_pct_beyond5": round((1 - sum(pooled)) * 100, 2) if pooled else None,
    }
    bt = {}
    for k in (2, 3, 4):
        pvv = [r["desc"][f"beyond{k}_vox"] / r["desc"]["n_nonsheet_scored"] for r in ok]
        pf = [r["desc"][f"beyond{k}_fp"] / r["desc"]["n_fp"] for r in ok if r["desc"]["n_fp"]]
        bt[f"beyond{k}"] = {
            "vol_share_per_volume_median_pct": round(float(np.median(pvv)) * 100, 2),
            "vol_share_pooled_pct": round(sum(r["desc"][f"beyond{k}_vox"] for r in ok)
                                          / sum(r["desc"]["n_nonsheet_scored"] for r in ok) * 100, 2),
            "fp_share_per_volume_median_pct": round(float(np.median(pf)) * 100, 3),
            "fp_share_pooled_pct": round(sum(r["desc"][f"beyond{k}_fp"] for r in ok) / tot_fp * 100, 3),
        }
    out["beyond_k_shares"] = bt
    out["shell_vol_share"] = {
        "per_volume_median_pct": [round(float(np.median([r["desc"]["shell_vox"][i]
                                                         / r["desc"]["n_nonsheet_scored"] for r in ok])) * 100, 2)
                                  for i in range(5)],
        "pooled_pct": [round(sum(r["desc"]["shell_vox"][i] for r in ok)
                             / sum(r["desc"]["n_nonsheet_scored"] for r in ok) * 100, 2) for i in range(5)],
    }
    return out


def compare(ours: list[dict], his_json: str) -> dict:
    his = {r["sample"]: r for r in json.loads(Path(his_json).read_text())["rows"]
           if r.get("status") == "ok"}
    o = {r["sample"]: r for r in ours if r.get("status") == "ok"}
    common = sorted(set(his) & set(o))
    res = {"n_common": len(common)}
    for k in range(1, 6):
        a = [his[s]["shells"].get(f"shell_{k}") for s in common]
        b = [o[s]["shells"].get(f"shell_{k}") for s in common]
        pair = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
        if len(pair) < 4:
            continue
        x = np.array([p[0] for p in pair]); y = np.array([p[1] for p in pair])
        res[f"shell_{k}"] = {
            "n": len(pair), "his_median": round(float(np.median(x)), 3),
            "ours_median": round(float(np.median(y)), 3),
            "ratio_of_medians": round(float(np.median(y) / np.median(x)), 3),
            "pearson_r": round(float(pearsonr(x, y)[0]), 4),
            "spearman_r": round(float(spearmanr(x, y)[0]), 4),
            "median_abs_rel_err": round(float(np.median(np.abs(y - x) / np.abs(x))), 4),
        }
    x = np.array([his[s]["enrichment"] for s in common])
    y = np.array([o[s]["enrichment"] for s in common])
    res["margin_enrichment"] = {
        "his_median": round(float(np.median(x)), 3), "ours_median": round(float(np.median(y)), 3),
        "pearson_r": round(float(pearsonr(x, y)[0]), 4),
        "spearman_r": round(float(spearmanr(x, y)[0]), 4),
        "median_abs_rel_err": round(float(np.median(np.abs(y - x) / np.abs(x))), 4)}
    # his implied FP mass share by shell, from his enrichments x our (identical) shell volumes
    imp, oursh = [], []
    for s in common:
        d = o[s]["desc"]
        sv = [v / d["n_nonsheet_scored"] for v in d["shell_vox"]]
        hs = [his[s]["shells"].get(f"shell_{k}") for k in range(1, 6)]
        if any(h is None for h in hs):
            continue
        imp.append([h * v for h, v in zip(hs, sv)])
        oursh.append([f / d["n_fp"] for f in d["shell_fp"]])
    imp, oursh = np.array(imp), np.array(oursh)
    res["fp_mass_share_by_shell_pct"] = {
        "his_implied_per_volume_median": [round(float(np.median(imp[:, i])) * 100, 2) for i in range(5)],
        "ours_per_volume_median": [round(float(np.median(oursh[:, i])) * 100, 2) for i in range(5)],
        "note": ("his implied share = his shipped enrichment_k x the shell-k volume share "
                 "measured off the labels, which are byte-identical between the two runs"),
    }
    return res


if __name__ == "__main__":
    rows = loadl(sys.argv[1])
    tag = sys.argv[2]
    res = summarize(rows, tag)
    if len(sys.argv) > 3 and sys.argv[3] != "-":
        res["vs_his"] = compare(rows, sys.argv[3])
    print(json.dumps(res, indent=1))
    if len(sys.argv) > 4:
        Path(sys.argv[4]).write_text(json.dumps({"summary": res, "rows": rows}, indent=1))
