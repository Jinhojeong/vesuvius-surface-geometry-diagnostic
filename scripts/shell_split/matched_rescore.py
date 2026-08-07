"""Adversarial re-scoring. Is our advantage a threshold artifact or a ranking difference?

Our path and his published path disagree on recall and precision over the same 60 volumes.
One cheap explanation is that the two runs simply sit at different points on the same
curve, in which case matching the operating point should make the gap close. This scores
OUR probability field four ways against the same labels and reports the paired differences.

  (a) his own endpoints, threshold 0.2 and the matched-budget point at 0.12
  (b) the threshold that reproduces HIS per-volume predicted-positive fraction exactly,
      then read recall off it
  (c) the threshold that reproduces HIS per-volume recall exactly, then read precision
      off it

Raw per-shell false-positive counts are carried along so the shell obs/pred question can
be recomputed from counts rather than from three-decimal enrichments.

What it needs
  --cache  a directory of cached probability maps, one {sample}.npy per volume, float16,
           already cropped to the inner 128^3. run_shells.py writes these when given
           --cache, so this script needs no GPU and no checkpoint. It is deterministic and
           reruns byte for byte.
  --labels the label tif directory used for the run
  --his    his published results/surface_bench_m7.json from
           github.com/TAUIL-Abd-Elilah/vesuvius-repro. His rows are copied in as the
           target operating points; his stack is not rerun here.

    python matched_rescore.py --cache cache_his60 --labels /mnt/vesuvius/kaggle892/labels \
        --his surface_bench_m7.json --out matched_rescore.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

import numpy as np
import tifffile
from scipy.ndimage import distance_transform_edt

TRIM, SIZE, THRESH, BUDGET = 64, 256, 0.2, 0.12


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--his", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cache, lbdir = Path(a.cache), Path(a.labels)
    his = {r["sample"]: r for r in json.loads(Path(a.his).read_text())["rows"]
           if r.get("status") == "ok"}

    rows = []
    for nm in sorted(his):
        f = cache / f"{nm}.npy"
        if not f.exists():
            continue
        p = np.load(f).astype(np.float32)
        lab = np.asarray(tifffile.imread(str(lbdir / f"{nm}.tif")))
        off = (lab.shape[0] - SIZE) // 2
        lo, hi = off + TRIM, off + SIZE - TRIM
        labc = lab[lo:hi, lo:hi, lo:hi]
        assert labc.shape == p.shape, (labc.shape, p.shape)
        sheet = labc == 1
        scored = labc != 2
        n_scored = float(scored.sum())
        n_sheet = float(sheet.sum())
        ps = p[scored]
        sh_s = sheet[scored]

        def at(t):
            pred = ps > t
            tp = float((pred & sh_s).sum())
            npp = float(pred.sum())
            return {"recall": tp / n_sheet, "precision": tp / max(npp, 1.0),
                    "ppf": npp / n_scored, "thr": float(t)}

        r = {"sample": nm, "n_sheet": int(n_sheet), "base_rate": n_sheet / n_scored,
             "his_n_sheet": his[nm]["n_sheet"], "his_base_rate": his[nm]["base_rate"]}
        v = at(THRESH)
        r.update({f"ours_{k}": x for k, x in v.items() if k != "thr"})
        r.update({"his_recall": his[nm]["recall"], "his_precision": his[nm]["precision"],
                  "his_ppf": his[nm]["pred_positive_fraction"]})

        # (a) his matched-budget endpoint, recomputed on our field
        b = at(float(np.quantile(ps, 1.0 - BUDGET)))
        r.update({"ours_budget_recall": b["recall"], "ours_budget_precision": b["precision"],
                  "ours_budget_thr": b["thr"], "ours_budget_ppf": b["ppf"],
                  "his_budget_recall": his[nm]["budget_recall"],
                  "his_budget_precision": his[nm]["budget_precision"],
                  "his_budget_thr": his[nm]["budget_threshold"]})

        # (b) match his predicted-positive fraction exactly
        m = at(float(np.quantile(ps, 1.0 - his[nm]["pred_positive_fraction"])))
        r.update({"ours_at_his_ppf_recall": m["recall"],
                  "ours_at_his_ppf_precision": m["precision"],
                  "ours_at_his_ppf_ppf": m["ppf"]})

        # (c) match his recall exactly, then compare precision. The threshold is his recall
        # quantile of our probabilities taken over sheet voxels only.
        p_sheet = np.sort(p[sheet])
        k = int(np.floor((1.0 - his[nm]["recall"]) * len(p_sheet)))
        c = at(float(p_sheet[min(max(k, 0), len(p_sheet) - 1)]))
        r.update({"ours_at_his_recall_recall": c["recall"],
                  "ours_at_his_recall_precision": c["precision"],
                  "ours_at_his_recall_ppf": c["ppf"]})

        d = distance_transform_edt(~sheet)
        ns = scored & ~sheet
        fp = (p > THRESH) & ns
        r["shell_fp"] = [int((fp & ((d > j - 1) & (d <= j) & ns)).sum()) for j in range(1, 6)]
        r["shell_vox"] = [int((((d > j - 1) & (d <= j)) & ns).sum()) for j in range(1, 6)]
        r["n_fp"] = int(fp.sum())
        r["n_ns"] = int(ns.sum())
        rows.append(r)
        print("MARK", nm, flush=True)

    Path(a.out).write_text(json.dumps(rows, indent=1))
    print("wrote", len(rows))


if __name__ == "__main__":
    main()
