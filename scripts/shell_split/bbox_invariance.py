"""Is OUR readout window sensitive to how much context the sliding window is given?

Run the same checkpoint over a restricted box and over the whole volume, then compare
the two on the region they share. This isolates one alternative explanation for our
divergence from a published run, namely that our restricted box is what makes our
numbers look the way they do.
What came back, over n=3 volumes. The aggregate endpoints move by at most 0.01 between the
two paths, recall and precision alike, while the thresholded masks agree only at Dice 0.923
to 0.942. So the summary statistics this study reports are stable to how much context the
sliding window is given, and the mask itself is patch-grid sensitive at the several-percent
Dice level. Read that as a bound on the aggregates, not as a claim that the two paths emit
the same voxels.
"""
import json, sys, time
from pathlib import Path
import numpy as np, tifffile, torch
sys.path.insert(0, str(Path(__file__).resolve().parent))
from geom import distance_profile, TRIM, SIZE, THRESH
from run_shells import build_predictor, ct_normalize, IM, LB

pred_obj, ip = build_predictor("cuda")
names = json.loads(Path(sys.argv[1]).read_text())
rows = []
for nm in names:
    ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
    lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
    off = (ct.shape[0] - SIZE) // 2
    lo, hi = off + TRIM, off + SIZE - TRIM

    def prob(block, o):
        x = torch.from_numpy(ct_normalize(block, ip))[None]
        with torch.no_grad():
            lg = pred_obj.predict_sliding_window_return_logits(x).float().cpu().numpy()
        a, b = lo - o, hi - o
        return 1.0 / (1.0 + np.exp(-(lg[1, a:b, a:b, a:b] - lg[0, a:b, a:b, a:b])))

    t0 = time.time()
    p_box = prob(ct[off:off + SIZE, off:off + SIZE, off:off + SIZE], off)
    p_full = prob(ct, 0)
    m1, m2 = p_box > THRESH, p_full > THRESH
    dice = 2 * (m1 & m2).sum() / max(1, m1.sum() + m2.sum())

    labc = lab[lo:hi, lo:hi, lo:hi]
    sheet, scored = labc == 1, labc != 2
    def rp(m):
        tp = float((m & sheet).sum())
        return (round(tp / max(1.0, float(sheet.sum())), 4),
                round(tp / max(1.0, float((m & scored).sum())), 4))
    row = {"sample": nm,
           "dice_box_vs_full": round(float(dice), 5),
           "max_abs_prob_diff": round(float(np.abs(p_box - p_full).max()), 5),
           "mean_abs_prob_diff": round(float(np.abs(p_box - p_full).mean()), 6),
           "recall_prec_box": rp(m1), "recall_prec_full": rp(m2),
           "base_rate": round(float(sheet.sum()) / float(scored.sum()), 5),
           "shells_box": distance_profile(labc, m1),
           "shells_full": distance_profile(labc, m2),
           "sec": round(time.time() - t0, 1)}
    rows.append(row)
    print(json.dumps(row), flush=True)
Path("/mnt/vesuvius/experiments/shell_split/bbox_invariance.json").write_text(
    json.dumps(rows, indent=1))
