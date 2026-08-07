"""Does our hand-fed normalisation match nnUNet's own preprocessing?

predict_single_npy_array runs DefaultPreprocessor end to end: crop-to-nonzero,
CTNormalization from plans, resample to target spacing, then the same sliding window,
then resample back and softmax. If p1/(p0+p1) off that path matches sigmoid(l1-l0) off
our path, then nothing in our preprocessing is home-made.
"""
import json, sys
from pathlib import Path
import numpy as np, tifffile, torch
sys.path.insert(0, str(Path(__file__).resolve().parent))
from geom import TRIM, SIZE, THRESH
from run_shells import build_predictor, ct_normalize, IM, LB

pred_obj, ip = build_predictor("cuda")
rows = []
for nm in json.loads(Path(sys.argv[1]).read_text()):
    ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
    lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
    off = (ct.shape[0] - SIZE) // 2
    lo, hi = off + TRIM, off + SIZE - TRIM

    # ours, full volume
    x = torch.from_numpy(ct_normalize(ct, ip))[None]
    with torch.no_grad():
        lg = pred_obj.predict_sliding_window_return_logits(x).float().cpu().numpy()
    p_ours = 1.0 / (1.0 + np.exp(-(lg[1, lo:hi, lo:hi, lo:hi] - lg[0, lo:hi, lo:hi, lo:hi])))

    # nnUNet's own preprocessing path
    seg, probs = pred_obj.predict_single_npy_array(
        ct[None].astype(np.float32), {"spacing": [1.0, 1.0, 1.0]},
        None, None, True)
    p0 = probs[0, lo:hi, lo:hi, lo:hi].astype(np.float64)
    p1 = probs[1, lo:hi, lo:hi, lo:hi].astype(np.float64)
    p_nn = p1 / np.maximum(p0 + p1, 1e-12)

    labc = lab[lo:hi, lo:hi, lo:hi]
    sheet, scored = labc == 1, labc != 2
    def rp(m):
        tp = float((m & sheet).sum())
        return (round(tp / max(1.0, float(sheet.sum())), 4),
                round(tp / max(1.0, float((m & scored).sum())), 4))
    m1, m2 = p_ours > THRESH, p_nn > THRESH
    row = {"sample": nm,
           "dice_ourprep_vs_nnunetprep": round(float(2 * (m1 & m2).sum() / max(1, m1.sum() + m2.sum())), 6),
           "max_abs_prob_diff": round(float(np.abs(p_ours - p_nn).max()), 6),
           "mean_abs_prob_diff": round(float(np.abs(p_ours - p_nn).mean()), 8),
           "recall_prec_ours": rp(m1), "recall_prec_nnunet": rp(m2)}
    rows.append(row)
    print(json.dumps(row), flush=True)
Path("/mnt/vesuvius/experiments/shell_split/preproc_crosscheck.json").write_text(json.dumps(rows, indent=1))
