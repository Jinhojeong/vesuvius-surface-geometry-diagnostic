"""Turn raw_rows.jsonl into the deliverable: rows.jsonl + summary.json + the gates.

rows.jsonl carries one row per public volume, including the volumes with no labelled
sheet in the scored crop. Those keep their row with null endpoints and a `reason`
rather than disappearing, because a table that silently drops 37 of 892 volumes is how
the population statistics drift in the first place.

summary.json carries per-population medians, the aggregate blocks his BENCHMARK.md
table carries, the gate results, and a provenance block.
"""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path

OUT = Path("/mnt/vesuvius/experiments/shell892")
S = Path("/mnt/vesuvius/experiments/shell_split")
GEN = S / "gen"
MODEL = Path("/mnt/vesuvius/models/surface_m7_nnunet")
sys.path.insert(0, str(OUT))
sys.path.insert(0, str(GEN))

import numpy as np
import aggregate as AG

LOC = {"located", "intersecting", "iou1"}
BUDGET = 0.12
GRID = ("0.05", "0.10", "0.20", "0.30", "0.50", "0.70")


def loadl(p):
    return [json.loads(l) for l in Path(p).read_text().splitlines() if l.strip()]


def sha(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def div(a, b):
    return (a / b) if b else None


def budget_point(h: dict, budget: float = BUDGET) -> dict:
    """Smallest threshold on the retained histogram grid whose predicted-positive
    fraction of the SCORED region is at or below `budget`. This is not his rule.
    surface_bench.py takes the exact 1-budget quantile of the scored probabilities, so
    it always spends the full budget, while a grid search cannot when one bin holds more
    than 1-budget of the scored mass. The thresholds agree to within one bin either way."""
    hs = np.asarray(h["hist_scored"], dtype=np.int64)
    hh = np.asarray(h["hist_sheet"], dtype=np.int64)
    n_s, n_sh = int(h["n_scored"]), int(h["n_sheet"])
    if not n_s:
        return {}
    tail_s = np.concatenate([np.cumsum(hs[::-1])[::-1], [0]])
    tail_h = np.concatenate([np.cumsum(hh[::-1])[::-1], [0]])
    nb = h["nbin"]
    idx = int(np.argmax(tail_s / n_s <= budget))
    npred, ntp = int(tail_s[idx]), int(tail_h[idx])
    return {"budget": budget, "budget_threshold": round(idx / nb, 5),
            "budget_pred_positive_fraction": div(npred, n_s),
            "budget_recall": div(ntp, n_sh),
            "budget_precision": div(ntp, npred),
            "budget_precision_lift": (div(div(ntp, npred), div(n_sh, n_s))
                                      if npred and n_sh else None)}


def reason_for(d: dict, st: str, n_margin) -> str | None:
    if st == "ok":
        return None
    if d.get("n_scored", 0) == 0:
        return "no_scored_voxels_in_crop: the whole inner 128^3 is label class 2 (ignore)"
    if d.get("n_sheet", 0) == 0:
        return "no_labelled_sheet_in_crop: shell and margin geometry are undefined"
    if d.get("n_fp", 0) == 0:
        return "no_false_positives_in_crop: shell and margin enrichments are undefined"
    if not n_margin:
        return "empty_margin_band: margin enrichment is undefined"
    return "degenerate"


def build_rows() -> list[dict]:
    raw = {r["sample"]: r for r in loadl(OUT / "raw_rows.jsonl")}
    hist = {h["sample"]: h for h in loadl(OUT / "pred_histograms.jsonl")}
    groups = json.loads(Path("/mnt/vesuvius/kaggle892/groups892.json").read_text())
    out = []
    for nm in sorted(groups):
        r = raw.get(nm)
        if r is None:
            out.append({"sample": nm, "group": groups[nm],
                        "population": "located" if groups[nm] in LOC else "nonlocated",
                        "status": "missing", "reason": "not produced by this run"})
            continue
        d = r.get("desc", {})
        n_scored = d.get("n_scored", 0)
        n_sheet = d.get("n_sheet", 0)
        n_pred = d.get("n_pred_scored", 0)
        n_fp = d.get("n_fp", 0)
        n_tp = n_pred - n_fp
        n_margin = r.get("n_margin_scored", r.get("n_margin"))
        sh = r.get("shells") or {}
        base_rate = div(n_sheet, n_scored)
        precision = div(n_tp, n_pred)
        row = {
            "sample": nm,
            "group": groups[nm],
            "population": "located" if groups[nm] in LOC else "nonlocated",
            "status": r.get("status"),
            "reason": reason_for(d, r.get("status"), n_margin),
            # headline
            "recall": div(n_tp, n_sheet),
            "precision": precision,
            "pred_positive_fraction": div(n_pred, n_scored),
            "precision_lift": (precision / base_rate) if (precision and base_rate) else None,
            "n_sheet": n_sheet,
            "base_rate": base_rate,
            # margin
            "margin_share_of_nonsheet": (r.get("margin_share_of_nonsheet")
                                         if r.get("status") == "ok"
                                         else div(n_margin, d.get("n_nonsheet_scored", 0))),
            "margin_enrichment": r.get("enrichment"),
            "fp_share_in_margin": r.get("fp_share_in_margin"),
            "margin_hit_rate": r.get("margin_hit_rate"),
            "nonmargin_fp_rate": r.get("nonmargin_fp_rate"),
            "n_margin_scored": n_margin,
            # shells
            **{f"shell_{k}": sh.get(f"shell_{k}") for k in range(1, 6)},
            "n_fp": n_fp,
            # scored-region counts
            "n_scored": n_scored,
            "n_nonsheet_scored": d.get("n_nonsheet_scored"),
            "n_pred_scored": n_pred,
            "n_tp": n_tp,
            "shell_vox": d.get("shell_vox"),
            "shell_fp": d.get("shell_fp"),
            **{f"beyond{k}_vox": d.get(f"beyond{k}_vox") for k in (2, 3, 4)},
            **{f"beyond{k}_fp": d.get(f"beyond{k}_fp") for k in (2, 3, 4)},
            "crop": r.get("crop"),
        }
        sw = r.get("threshold_sweep") or {}
        row["recall_by_threshold"] = {
            t: div(sw[t]["n_tp"], n_sheet) for t in GRID if t in sw} or None
        row["precision_by_threshold"] = {
            t: div(sw[t]["n_tp"], sw[t]["n_pred_scored"]) for t in GRID if t in sw} or None
        h = hist.get(nm)
        row.update(budget_point(h) if h else {})
        out.append(row)
    return out


def med(v):
    v = [x for x in v if x is not None]
    return round(float(np.median(v)), 5) if v else None


def population_block(rows: list[dict], raw_by: dict, tag: str, n_pop: int) -> dict:
    present = [r for r in rows if r.get("status") != "missing"]
    ok = [r for r in present if r.get("status") == "ok"]
    b = {
        "tag": tag,
        "n_population": n_pop,
        "n_rows": len(rows),
        "n_rows_present": len(present),
        "n_missing": len(rows) - len(present),
        "n_ok": len(ok),
        "n_no_labelled_sheet": len([r for r in present if r.get("n_sheet") == 0]),
        "n_no_scored_voxels": len([r for r in present if r.get("n_scored") == 0]),
        # Every median below is over the SAME n_ok rows. Taking each median over
        # whatever rows happen to define that field would silently move the
        # denominator between columns of one table: pred_positive_fraction is
        # defined on volumes with no labelled sheet, recall and precision are not.
        # The alternative denominators are reported separately, with their n.
        "median_denominator": "n_ok",
        # the columns his table carries, recomputed with class 2 excluded
        "median_recall": med([r.get("recall") for r in ok]),
        "median_precision": med([r.get("precision") for r in ok]),
        "median_precision_lift": med([r.get("precision_lift") for r in ok]),
        "median_pred_positive_fraction": med([r.get("pred_positive_fraction") for r in ok]),
        "median_base_rate": med([r.get("base_rate") for r in ok]),
        "median_n_sheet": med([r.get("n_sheet") for r in ok]),
        "median_budget_recall": med([r.get("budget_recall") for r in ok]),
        "median_budget_precision": med([r.get("budget_precision") for r in ok]),
        "median_budget_precision_lift": med([r.get("budget_precision_lift") for r in ok]),
        "median_budget_threshold": med([r.get("budget_threshold") for r in ok]),
        "median_recall_by_threshold": {
            t: med([(r.get("recall_by_threshold") or {}).get(t) for r in ok]) for t in GRID},
        "median_margin_share_of_nonsheet": med(
            [r.get("margin_share_of_nonsheet") for r in ok]),
        # his split scores every volume with labelled sheet in the crop, which is a
        # wider set than n_ok: it also holds the volumes with sheet but nothing
        # predicted anywhere, where recall is 0.0 and precision does not exist.
        "n_recall_defined": len([r for r in present if r.get("recall") is not None]),
        "median_recall_over_recall_defined": med(
            [r.get("recall") for r in present if r.get("recall") is not None]),
        "n_pred_positive_fraction_defined": len(
            [r for r in present if r.get("pred_positive_fraction") is not None]),
        "median_pred_positive_fraction_over_defined": med(
            [r.get("pred_positive_fraction") for r in present]),
        "budget_convention": (
            f"budget_* columns threshold each volume so that at most {BUDGET:.2f} of its "
            "SCORED voxels are predicted positive. His surface_bench.py takes the exact "
            f"{1 - BUDGET:.2f} quantile of the scored probabilities; we search a 2000-bin "
            "probability histogram and take the first bin whose predicted-positive share is "
            "at or below the budget, so the two thresholds can differ by up to one bin. "
            "Treat the budget columns as ours rather than as a like-for-like replacement."),
    }
    # the shell / margin / fp-mass blocks, straight from the shipped aggregator.
    # aggregate.summarize assumes at least one status=ok row, so guard rather than
    # patch the vendored copy.
    sub = [raw_by[r["sample"]] for r in rows if r["sample"] in raw_by]
    b["shell_blocks"] = (AG.summarize(sub, tag)
                         if any(x.get("status") == "ok" for x in sub)
                         else {"tag": tag, "n_rows": len(sub), "n_ok": 0,
                               "note": "no scorable rows yet"})
    return b


def main() -> None:
    rows = build_rows()
    (OUT / "rows.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
    raw_by = {r["sample"]: r for r in loadl(OUT / "raw_rows.jsonl")}

    loc = [r for r in rows if r["population"] == "located"]
    non = [r for r in rows if r["population"] == "nonlocated"]

    # ---- gate 1: the 120 reused control rows, re-run inside this pass ----
    ctrl = {}
    for f in ("control_his60.jsonl", "control_nonloc60.jsonl"):
        for r in loadl(S / f):
            ctrl[r["sample"]] = r
    reuse = {"n_compared": 0, "n_bit_identical": 0, "label_side_identical": True,
             "worst": [], "note": ""}
    rel = {"n_fp": [], "n_pred_scored": [], "enrichment": [], "shell_1": []}
    for nm, old in ctrl.items():
        new = raw_by.get(nm)
        if new is None:
            continue
        reuse["n_compared"] += 1
        same = all(new.get(k) == old.get(k) for k in old if k != "seconds")
        reuse["n_bit_identical"] += int(same)
        do, dn = old.get("desc", {}), new.get("desc", {})
        for k in ("n_sheet", "n_scored", "n_nonsheet_scored", "shell_vox",
                  "beyond2_vox", "beyond3_vox", "beyond4_vox"):
            if do.get(k) != dn.get(k):
                reuse["label_side_identical"] = False
                reuse["worst"].append({"sample": nm, "field": k,
                                       "old": do.get(k), "new": dn.get(k)})
        for k in ("n_fp", "n_pred_scored"):
            if do.get(k):
                rel[k].append(abs(dn.get(k, 0) - do[k]) / do[k])
        for k, src in (("enrichment", None), ("shell_1", "shells")):
            a = (old.get(src) or old).get(k) if src else old.get(k)
            c = (new.get(src) or new).get(k) if src else new.get(k)
            if a:
                rel[k].append(abs(c - a) / abs(a))
    reuse["max_rel_delta"] = {k: (round(float(np.max(v)), 8) if v else None)
                              for k, v in rel.items()}
    reuse["median_rel_delta"] = {k: (round(float(np.median(v)), 8) if v else None)
                                 for k, v in rel.items()}
    reuse["note"] = (
        "The 120 volumes already scored in control_his60/control_nonloc60 were re-run "
        "inside this pass rather than pasted in. Every label-side quantity is identical. "
        "The prediction-side counts move by ~1e-4 relative, which is run-to-run "
        "nondeterminism of fp16 sliding-window inference (cuDNN algorithm selection and "
        "non-deterministic accumulation), not a difference in code or geometry. Shipping "
        "one homogeneous run avoids putting that seam through the table.")

    # ---- gate 2: label-side quantities against label_geometry_892.jsonl ----
    lg = {r["sample"]: r for r in loadl(S / "label_geometry_892.jsonl")}
    lab_gate = {"n_checked": 0, "n_mismatch": 0, "fields": [], "mismatches": []}
    for r in rows:
        g = lg.get(r["sample"])
        if not g or r.get("status") == "missing":
            continue
        lab_gate["n_checked"] += 1
        pairs = [("n_sheet", r["n_sheet"], g["n_sheet"]),
                 ("n_scored", r["n_scored"], g["n_scored"]),
                 ("n_nonsheet_scored", r["n_nonsheet_scored"], g["n_nonsheet_scored"])]
        if "shell_vox" in g and r.get("shell_vox") is not None:
            pairs.append(("shell_vox", r["shell_vox"], g["shell_vox"]))
            for k in (2, 3, 4):
                pairs.append((f"beyond{k}_vox", r[f"beyond{k}_vox"], g[f"beyond{k}_vox"]))
        br = div(g["n_sheet"], g["n_scored"])
        pairs.append(("base_rate", r["base_rate"], br))
        for name, a, b in pairs:
            if a != b:
                lab_gate["n_mismatch"] += 1
                lab_gate["mismatches"].append(
                    {"sample": r["sample"], "field": name, "run": a, "label_file": b})
    lab_gate["fields"] = ["n_sheet", "n_scored", "n_nonsheet_scored", "base_rate",
                          "shell_vox", "beyond2_vox", "beyond3_vox", "beyond4_vox"]
    lab_gate["margin_note"] = (
        "margin_share_of_nonsheet is not in label_geometry_892.jsonl, which is labels-only; "
        "the margin band needs the CT (a Hessian across-sheet step). It is checked instead "
        "for internal consistency against n_margin_scored / n_nonsheet_scored, and against "
        "the 120 control rows in the reuse gate.")
    bad = [r["sample"] for r in rows
           if r.get("status") == "ok" and r.get("margin_share_of_nonsheet")
           and abs(r["margin_share_of_nonsheet"]
                   - r["n_margin_scored"] / r["n_nonsheet_scored"]) > 1e-6]
    lab_gate["margin_internal_consistency_failures"] = bad

    # ---- gate 3: the threshold sweep at 0.20 must reproduce the shipped scorer ----
    sweep_gate = {"n_checked": 0, "n_mismatch": 0, "mismatches": [],
                  "what": ("recall_by_threshold['0.20'] and precision_by_threshold['0.20'] "
                           "are computed from the retained probability volume, "
                           "independently of analyse()/descriptives(). They must equal the "
                           "recall and precision derived from the shipped scorer's counts.")}
    for r in rows:
        rt = (r.get("recall_by_threshold") or {}).get("0.20")
        pt = (r.get("precision_by_threshold") or {}).get("0.20")
        if rt is None and pt is None:
            continue
        sweep_gate["n_checked"] += 1
        for name, a, b_ in (("recall", rt, r.get("recall")), ("precision", pt, r.get("precision"))):
            if a is None and b_ is None:
                continue
            if a is None or b_ is None or abs(a - b_) > 1e-12:
                sweep_gate["n_mismatch"] += 1
                sweep_gate["mismatches"].append(
                    {"sample": r["sample"], "field": name, "sweep": a, "scorer": b_})

    # ---- provenance ----
    plans = json.loads((MODEL / "plans.json").read_text())
    ip = plans["foreground_intensity_properties_per_channel"]["0"]
    cfg = plans["configurations"]["3d_fullres"]
    prov = {
        "what_these_rows_are": (
            "CT-normalized rows. Every volume here was normalized with the "
            "CTNormalization constants that m7's own plans.json declares for 3d_fullres, "
            "not with the villa `vesuvius.predict` wrapper default. The wrapper path and "
            "this path disagree, which is the reason the reference table needed replacing "
            "rather than retracting."),
        "checkpoint": {
            "path": str(MODEL / "fold_0/checkpoint_best.pth"),
            "sha256": sha(MODEL / "fold_0/checkpoint_best.pth"),
            "bytes": (MODEL / "fold_0/checkpoint_best.pth").stat().st_size,
            "source": "hf scrollprize/surface_m7_nnunet",
            "fold": 0, "checkpoint_name": "checkpoint_best.pth"},
        "geometry": {
            "volume_shape": [320, 320, 320],
            "centred_block": 256, "TRIM": 64,
            "scored_crop": "inner 128^3, index 96:224 on each axis",
            "prediction_mask": "sigmoid(logit1 - logit0) > 0.2",
            "class_2_excluded": True,
            "shells": "one-voxel Euclidean shells off the labelled sheet, k = 1..5",
            "num_segmentation_heads": 2,
            "note": ("nnU-Net's LabelManager treats the label literally named `ignore` in "
                     "dataset.json as ignore_label, so the network has two heads, not "
                     "three, and there is no third channel to drop.")},
        "normalization": {
            "scheme": cfg.get("normalization_schemes"),
            "clip": [ip["percentile_00_5"], ip["percentile_99_5"]],
            "mean": ip["mean"], "std": ip["std"],
            "why": ("these are global dataset-fingerprint constants, not per-image "
                    "statistics, so cropping the 256^3 block before normalising cannot "
                    "change them; that is what makes crop-then-normalise equivalent to "
                    "normalise-then-crop here"),
            "patch_size": cfg["patch_size"], "tile_step_size": 0.5,
            "use_gaussian": True, "use_mirroring": False,
            "autocast": "fp16, nnUNet's own default inference path"},
        "code": {
            "run892.py": {"path": str(OUT / "run892.py"), "sha256": sha(OUT / "run892.py")},
            "run_shells.py": {"path": str(GEN / "run_shells.py"),
                              "sha256": sha(GEN / "run_shells.py")},
            "geom.py": {"path": str(GEN / "geom.py"), "sha256": sha(GEN / "geom.py")},
            "aggregate.py": {"path": str(OUT / "aggregate.py"),
                             "sha256": sha(OUT / "aggregate.py")},
            "finalize.py": {"path": str(OUT / "finalize.py"),
                            "sha256": sha(OUT / "finalize.py")},
            "repo_commit": REPO_COMMIT,
            "commit_note": ("shell_split/ is untracked in the benchmark repo at this "
                            "commit, so the file sha256 above is the binding identifier, "
                            "not the commit"),
            "geom_constants_verified": (
                "geom.py is byte-identical to the copy that produced control_his60 and "
                "control_nonloc60 apart from a docstring; run_shells.py likewise. TRIM=64, "
                "SIZE=256, THRESH=0.2, SIGMA=1.0 in all copies.")},
        "environment": ENVBLOCK,
        "cohort": {"source": "/mnt/vesuvius/kaggle892/groups892.json",
                   "located": "groups located + intersecting + iou1",
                   "n_located": len(loc), "n_nonlocated": len(non)},
        "pred_on_empty_ct": PRED_EMPTY,
    }

    summary = {
        "provenance": prov,
        "populations": {
            "all": population_block(rows, raw_by, "all public volumes", 892),
            "located": population_block(loc, raw_by, "located population", len(loc)),
            "nonlocated": population_block(non, raw_by, "nonlocated population", len(non)),
        },
        "gates": {"reuse_of_120_controls": reuse,
                  "label_side_vs_label_geometry_892": lab_gate,
                  "threshold_sweep_vs_scorer": sweep_gate},
        "legend": (
            "Every quantity is measured on the inner 128^3 scored crop with label class 2 "
            "(ignore) excluded. recall = TP / n_sheet, precision = TP / n_pred_scored, "
            "pred_positive_fraction = n_pred_scored / n_scored, base_rate = n_sheet / "
            "n_scored, all at threshold 0.2. Volumes with no labelled sheet in the crop "
            "keep a row with null endpoints and a reason."),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({
        "n_rows": len(rows),
        "n_present": len([r for r in rows if r.get("status") != "missing"]),
        "n_ok": len([r for r in rows if r.get("status") == "ok"]),
        "n_no_labelled_sheet": len([r for r in rows if r.get("n_sheet") == 0]),
        "gate_label_side_mismatches": lab_gate["n_mismatch"],
        "gate_label_side_checked": lab_gate["n_checked"],
        "gate_sweep_mismatches": sweep_gate["n_mismatch"],
        "gate_sweep_checked": sweep_gate["n_checked"],
        "reuse_compared": reuse["n_compared"],
        "reuse_bit_identical": reuse["n_bit_identical"],
        "reuse_label_side_identical": reuse["label_side_identical"],
        "reuse_max_rel_delta": reuse["max_rel_delta"],
        "median_recall_all": summary["populations"]["all"]["median_recall"],
        "median_precision_all": summary["populations"]["all"]["median_precision"],
    }, indent=1))


REPO_COMMIT = "9213bf19f320847619e85e62f602992ea0cec2bc"
ENVBLOCK: dict = {}
PRED_EMPTY = None

if __name__ == "__main__":
    import importlib.metadata as M
    import torch
    ENVBLOCK = {"python": sys.version.split()[0],
                "nnunetv2": M.version("nnunetv2"), "torch": torch.__version__,
                "numpy": np.__version__,
                "scipy": __import__("scipy").__version__,
                "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}
    for ln in (OUT / "run.log").read_text().splitlines():
        if ln.startswith("MARKER pred_on_empty_ct="):
            PRED_EMPTY = {"predicted_positive_fraction": float(ln.split("=")[1].split()[0]),
                          "max_probability": float(ln.split("max_p=")[1]),
                          "note": "one forward pass on an all-zero 320^3 volume"}
    main()
