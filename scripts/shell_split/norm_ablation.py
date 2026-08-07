"""Normalization ablation on our own path, everything else held fixed.

The question is narrow. `vesuvius.predict` defaults to `--normalization instance_zscore`
and never reads the scheme out of an nnU-Net checkpoint's plans.json, so a checkpoint whose
plans declare CTNormalization gets per-volume z-scoring instead. This runs the same
published checkpoint over the same centred crop with the same threshold and the same
sliding window, and changes only the input normalization, so any difference in the
endpoints is the normalization and nothing else.

Two presets.

  --arms pair  the two schemes that matter, plans CTNormalization against instance
               z-score. This is what norm_ablation16.json holds.
  --arms all   adds instance z-score over non-zero voxels only, instance min-max, and the
               raw uint8 block with no normalization at all, as a spread check on how much
               of the effect is z-scoring specifically. This is what norm_ablation.json
               holds, on a smaller cohort because it costs five forward passes a volume.

Both presets record the centred block's own mean and standard deviation next to the plans
constants, because how far a volume sits from the training fingerprint is what governs how
much the two schemes can differ.

What it needs
  a CUDA device and the published checkpoint at run_shells.MODEL, loadable by nnUNetv2
  the volumes and labels at run_shells.IM and run_shells.LB
  --his, his published results/surface_bench_m7.json from
         github.com/TAUIL-Abd-Elilah/vesuvius-repro. It fixes the cohort, since the rows
         are taken in sorted order off the volumes he scored, and it carries his own
         endpoints so each row can be read against them. Nothing in this script runs his
         stack; his numbers are copied in, not reproduced.

    python norm_ablation.py --his surface_bench_m7.json --arms pair  --limit 16 \
        --out norm_ablation16.json
    python norm_ablation.py --his surface_bench_m7.json --arms all   --limit 4  \
        --out norm_ablation.json

The sliding window runs under fp16 autocast, which is nnU-Net's own inference default, so
this is not bit-reproducible. Reruns move individual endpoints by about 1e-4. Every median
reported off these files is unchanged at four decimals across reruns, but do not expect a
byte-identical file. matched_rescore.py, which scores cached probabilities, is exact.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

import numpy as np
import tifffile
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geom import TRIM, SIZE, THRESH
from run_shells import build_predictor, ct_normalize, IM, LB

LEGEND = ("every arm is [recall, precision, predicted_positive_fraction], measured over the "
          "scored crop at threshold %.1f with class 2 excluded. 'his' is copied from his "
          "published rows and is not rerun here. 'ct_plans' is the CTNormalization the "
          "checkpoint's own plans.json declares; 'instance_zscore' is what "
          "vesuvius.predict applies by default." % THRESH)

ARMS = {"pair": ["ct_plans", "instance_zscore"],
        "all": ["ct_plans", "instance_zscore", "instance_zscore_nonzero",
                "instance_minmax", "raw_uint8_none"]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--his", required=True, help="his published surface_bench_m7.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--arms", choices=sorted(ARMS), default="pair")
    ap.add_argument("--limit", type=int, default=16)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()

    pred_obj, ip = build_predictor(a.device)
    his = {r["sample"]: r for r in json.loads(Path(a.his).read_text())["rows"]
           if r.get("status") == "ok"}
    names = sorted(his)[:a.limit]

    def forward(x: np.ndarray) -> np.ndarray:
        t = torch.from_numpy(np.ascontiguousarray(x, dtype=np.float32))[None]
        with torch.no_grad():
            return pred_obj.predict_sliding_window_return_logits(t).float().cpu().numpy()

    rows = []
    for nm in names:
        ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
        lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
        off = (ct.shape[0] - SIZE) // 2
        lo, hi = off + TRIM, off + SIZE - TRIM
        blk = ct[off:off + SIZE, off:off + SIZE, off:off + SIZE].astype(np.float32)
        labc = lab[lo:hi, lo:hi, lo:hi]
        sheet, scored = labc == 1, labc != 2

        def endpoints(lg: np.ndarray) -> list:
            s = slice(TRIM, SIZE - TRIM)
            p = 1.0 / (1.0 + np.exp(-(lg[1, s, s, s] - lg[0, s, s, s])))
            m = p > THRESH
            tp = float((m & sheet).sum())
            npp = float((m & scored).sum())
            return [round(tp / float(sheet.sum()), 4), round(tp / max(npp, 1.0), 4),
                    round(npp / float(scored.sum()), 4)]

        mu, sd = float(blk.mean()), float(blk.std())
        nz = blk[blk > 0]
        arm_input = {
            "ct_plans": lambda: ct_normalize(blk, ip),
            "instance_zscore": lambda: (blk - mu) / max(sd, 1e-8),
            "instance_zscore_nonzero": lambda: (blk - float(nz.mean())) / max(float(nz.std()), 1e-8),
            "instance_minmax": lambda: (blk - blk.min()) / max(float(blk.max() - blk.min()), 1e-8),
            "raw_uint8_none": lambda: blk.copy(),
        }
        r = {"sample": nm, "blk_mean": round(mu, 2), "blk_std": round(sd, 2),
             "plans_mean": round(float(ip["mean"]), 2), "plans_std": round(float(ip["std"]), 2),
             "his": [round(his[nm]["recall"], 4), round(his[nm]["precision"], 4),
                     round(his[nm]["pred_positive_fraction"], 4)]}
        for arm in ARMS[a.arms]:
            r[arm] = endpoints(forward(arm_input[arm]()))
        rows.append(r)
        print(json.dumps(r), flush=True)

    Path(a.out).write_text(json.dumps({"legend": LEGEND, "rows": rows}, indent=1))

    g = lambda k, i: np.array([r[k][i] for r in rows], float)
    print("\n=== medians over %d volumes (recall / precision / pred-pos-frac) ===" % len(rows))
    for k in ["his"] + ARMS[a.arms]:
        print(" %-24s %.4f  %.4f  %.4f" % (k, np.median(g(k, 0)), np.median(g(k, 1)),
                                           np.median(g(k, 2))))
    for k in ARMS[a.arms]:
        print(" |%s - his| median abs: recall %.4f  prec %.4f  ppf %.4f"
              % (k, np.median(abs(g(k, 0) - g("his", 0))),
                 np.median(abs(g(k, 1) - g("his", 1))),
                 np.median(abs(g(k, 2) - g("his", 2)))))


if __name__ == "__main__":
    main()
