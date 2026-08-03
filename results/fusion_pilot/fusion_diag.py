"""Secondary decomposition for the fusion readout (exploratory, labeled so).

The preregistered fused/separated verdict counts a broken predicted run as
separated regardless of where it breaks. This addendum splits non-fused
neighbour sites into two kinds. A separation is a predicted gap when the
break lies strictly between the two sheets' true material, and a miss when
the prediction already fails on the site's own sheet material. Also reports
detection strength at site centers per cell.
"""
import glob
import json
import os

import numpy as np
from scipy import ndimage as ndi

SPAN, STEP = 12.0, 0.5
K = int(2 * SPAN / STEP) + 1
OFFS = np.linspace(-SPAN, SPAN, K)
CTR = K // 2
N_SAMPLE = 20000
THRESHOLDS = (0.4, 0.5, 0.6)


def diag(path):
    with np.load(path) as d:
        prob = d["prob"].astype(np.float32)
        tid = d["turn_id"].astype(np.int16)
    mask = tid > 0
    sm = ndi.gaussian_filter(
        ndi.distance_transform_edt(~mask).astype(np.float32)
        - ndi.distance_transform_edt(mask).astype(np.float32), 1.0)
    g = np.stack(np.gradient(sm), 0).astype(np.float32)
    surf = mask & ~ndi.binary_erosion(mask)
    pts = np.argwhere(surf)
    rng = np.random.default_rng(1218)
    pts = pts[rng.choice(len(pts), size=min(N_SAMPLE, len(pts)),
                         replace=False)]
    nv = g[:, pts[:, 0], pts[:, 1], pts[:, 2]].T
    nv /= np.linalg.norm(nv, axis=1, keepdims=True) + 1e-6
    cc = pts[:, None, :].astype(np.float32) + OFFS[None, :, None] * nv[:, None, :]
    flat = cc.reshape(-1, 3).T
    rt = ndi.map_coordinates(tid.astype(np.float32), flat, order=0,
                             mode="constant").reshape(len(pts), K).astype(np.int16)
    rp = ndi.map_coordinates(prob, flat, order=1,
                             mode="constant").reshape(len(pts), K)

    ctr_probs = []
    n_nb = 0
    res = {t: {"fused": 0, "sep_gap": 0, "sep_miss": 0, "sep_nbmiss": 0}
           for t in THRESHOLDS}
    for i in range(len(pts)):
        own = rt[i, CTR]
        if own <= 0:
            continue
        other = (rt[i] > 0) & (rt[i] != own)
        if not other.any():
            continue
        n_nb += 1
        ctr_probs.append(float(rp[i, CTR]))
        idx = np.where(other)[0]
        kn = idx[np.argmin(np.abs(idx - CTR))]
        a, b = (CTR, kn) if kn > CTR else (kn, CTR)
        seg_t = rt[i, a:b + 1]
        seg_p = rp[i, a:b + 1]
        own_part = seg_t == own
        gap_part = seg_t == 0
        nb_part = (seg_t > 0) & (seg_t != own)
        for t in THRESHOLDS:
            pred = seg_p >= t
            if pred.all():
                res[t]["fused"] += 1
            elif (~pred & own_part).any() or rp[i, CTR] < t:
                res[t]["sep_miss"] += 1
            elif (~pred & gap_part).any():
                res[t]["sep_gap"] += 1
            else:
                res[t]["sep_nbmiss"] += 1

    cp = np.array(ctr_probs)
    out = {"n_neighbour_sites": n_nb,
           "center_prob_quartiles": [round(float(np.percentile(cp, q)), 3)
                                     for q in (25, 50, 75)]}
    for t in THRESHOLDS:
        r = res[t]
        n = max(n_nb, 1)
        det = int((cp >= t).sum())
        out[f"t{t}"] = {
            "fused_pct": round(100.0 * r["fused"] / n, 2),
            "sep_true_gap_kept_pct": round(100.0 * r["sep_gap"] / n, 2),
            "sep_neighbour_miss_pct": round(100.0 * r["sep_nbmiss"] / n, 2),
            "sep_own_miss_pct": round(100.0 * r["sep_miss"] / n, 2),
            "detected_at_center_pct": round(100.0 * det / n, 2),
            "fused_given_detected_pct": round(
                100.0 * r["fused"] / max(det, 1), 2),
        }
    return out


def main():
    files = sorted(glob.glob("/mnt/vesuvius/fusion_pilot/*.npz"))
    out = {}
    for f in files:
        name = os.path.basename(f)[:-4]
        out[name] = diag(f)
        print(name, json.dumps(out[name]), flush=True)
    json.dump(out, open("/mnt/vesuvius/fusion_pilot/fusion_diag.json", "w"),
              indent=1)
    print("DIAG_DONE", flush=True)


main()
