"""Preregistered fusion readout for the 8-cell physical-arm pilot.

Implements exactly the definition posted in villa#191 comment 5162354182
before the data existed. At surface sites where turn_id puts a neighbouring
sheet within the ray span, threshold the probability and walk the same
along-normal ray as the PHerc 1218 work, across the cell's true air gap. The
site counts as fused when predicted material runs unbroken from one sheet
into the other across that gap, and as separated when the prediction keeps a
gap there. Per cell, fusion rate over those neighbour sites and false-split
rate over single-sheet control sites, at thresholds 0.4 / 0.5 / 0.6.
"""
import glob
import json
import os
import sys

import numpy as np
from scipy import ndimage as ndi

SPAN, STEP = 12.0, 0.5
K = int(2 * SPAN / STEP) + 1
OFFS = np.linspace(-SPAN, SPAN, K)
CTR = K // 2
N_SAMPLE = 20000
THRESHOLDS = (0.4, 0.5, 0.6)


def readout(path):
    with np.load(path) as d:
        prob = d["prob"].astype(np.float32)
        gt = d["gt_surface"].astype(bool)
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
    ray_tid = ndi.map_coordinates(tid.astype(np.float32), flat, order=0,
                                  mode="constant").reshape(len(pts), K)
    ray_tid = ray_tid.astype(np.int16)
    ray_prob = ndi.map_coordinates(prob, flat, order=1,
                                   mode="constant").reshape(len(pts), K)

    n_nb = 0
    fused = {t: 0 for t in THRESHOLDS}
    n_ct = 0
    fsplit = {t: 0 for t in THRESHOLDS}
    for i in range(len(pts)):
        own = ray_tid[i, CTR]
        if own <= 0:
            continue
        row = ray_tid[i]
        other = (row > 0) & (row != own)
        if other.any():
            # neighbour site: nearest other-sheet sample along the ray
            idx = np.where(other)[0]
            kn = idx[np.argmin(np.abs(idx - CTR))]
            a, b = (CTR, kn) if kn > CTR else (kn, CTR)
            n_nb += 1
            for t in THRESHOLDS:
                if (ray_prob[i, a:b + 1] >= t).all():
                    fused[t] += 1
        else:
            # control site: the true contiguous own-sheet run through p0
            a = CTR
            while a > 0 and row[a - 1] == own:
                a -= 1
            b = CTR
            while b < K - 1 and row[b + 1] == own:
                b += 1
            n_ct += 1
            for t in THRESHOLDS:
                if not (ray_prob[i, a:b + 1] >= t).all():
                    fsplit[t] += 1

    return {
        "n_neighbour_sites": n_nb,
        "fusion_rate": {str(t): round(100.0 * fused[t] / max(n_nb, 1), 2)
                        for t in THRESHOLDS},
        "n_control_sites": n_ct,
        "false_split_rate": {str(t): round(100.0 * fsplit[t] / max(n_ct, 1), 2)
                             for t in THRESHOLDS},
    }


def main():
    files = sorted(glob.glob("/mnt/vesuvius/fusion_pilot/*.npz"))
    print(f"{len(files)} cells", flush=True)
    out = {}
    for f in files:
        name = os.path.basename(f)[:-4]
        r = readout(f)
        out[name] = r
        print(name, json.dumps(r), flush=True)
    json.dump(out, open("/mnt/vesuvius/fusion_pilot/fusion_readout.json", "w"),
              indent=1)
    print("DONE", flush=True)


main()
