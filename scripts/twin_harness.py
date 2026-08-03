"""Finite-thickness, instance-labelled build of Diego-dcv's synthetic twin.

The shipped make_volume paints each crushed turn as a zero-thickness curve
(half_t computed, never used) and its --fuse maps the welded sector onto the
target ellipse coincidentally, so no voxel twin of a two-sheet stack exists.
This harness reuses the repo's geometry (build_twin, ellipse tables) and
paints what our instruments need:

  intensity.npy  uint8 (z,y,x), papyrus 90, air 0
  labels.npy     int16 turn ids (1-based), the true sheet instances
  weld_truth.json  the welded sector (turns, angle range) for recall checks

Weld model: the source turn's sector is placed at contact distance (2*half_t
centre separation) OUTSIDE the target turn's ellipse, so the two sheets touch
with no air gap and stay distinct instances in the labels. That is the
double-thickness single-mask column our flag is built for, with true GT.

Usage (on the box):
  python twin_harness.py --columns 24 --voxel-um 20 --z-mm 2 \
      --fuse 20,21,150,210 --out /mnt/vesuvius/twin/run1
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, "/mnt/vesuvius/twin/vesuvius-topological-grid/scripts")
import synthetic_scroll_twin as T


def paint(args):
    truth, meta = T.build_twin(args.columns)
    g = T.G
    vox = args.voxel_um / 1000.0            # mm per voxel
    half_t = g["sheet_um"] / 2000.0         # half sheet thickness, mm
    p = g["pitch_um"] / 1000.0
    a_out, b_out = meta["section_w_mm"] / 2, meta["section_h_mm"] / 2
    pad = 1.0
    nx = int((2 * a_out + 2 * pad) / vox)
    ny = int((2 * b_out + 2 * pad) / vox)
    nz = max(int(args.z_mm / vox), 8)
    n_turns = int(np.ceil(meta["n_turns"]))
    print(f"grid ({nz},{ny},{nx}) vox={args.voxel_um}um turns={n_turns} "
          f"sheet={g['sheet_um']}um gap={g['pitch_um']-g['sheet_um']}um",
          flush=True)

    fuse = None
    if args.fuse:
        ti, tj, a0, a1 = (float(v) for v in args.fuse.split(","))
        fuse = (int(ti), int(tj), a0, a1)

    # dense parametric samples per turn so curve raster has no gaps
    t_hi = np.linspace(0, 2 * np.pi, 20000)
    ang = np.degrees(t_hi) % 360

    labels2d = np.zeros((ny, nx), np.int16)
    dist_best = np.full((ny, nx), 1e9, np.float32)
    half_vox = half_t / vox

    for t in range(n_turns):
        r_t = g["r0_mm"] + (t + 0.5) * p
        a, b = T.ellipse_axes_for_perimeter(2 * np.pi * r_t, g["ratio"])
        xs, ys = a * np.cos(t_hi), b * np.sin(t_hi)
        if fuse and fuse[0] <= t <= fuse[1] and t != fuse[1]:
            # weld: stack this turn at contact OUTSIDE the target ellipse,
            # each source turn at its own multiple of the sheet thickness
            # (centre separation 2*half_t per layer), keeping identity
            a2, b2 = T.ellipse_axes_for_perimeter(
                2 * np.pi * (g["r0_mm"] + (fuse[1] + 0.5) * p), g["ratio"])
            rr = np.hypot(a2 * np.cos(t_hi), b2 * np.sin(t_hi))
            sc = (rr + 2 * half_t * (fuse[1] - t)) / rr
            sel = (ang >= fuse[2]) & (ang <= fuse[3])
            xs = np.where(sel, a2 * np.cos(t_hi) * sc, xs)
            ys = np.where(sel, b2 * np.sin(t_hi) * sc, ys)
        ix = np.clip(((xs + a_out + pad) / vox), 0, nx - 1).astype(int)
        iy = np.clip(((ys + b_out + pad) / vox), 0, ny - 1).astype(int)
        curve = np.zeros((ny, nx), bool)
        curve[iy, ix] = True
        d = ndi.distance_transform_edt(~curve).astype(np.float32)
        ring = d <= half_vox
        take = ring & (d < dist_best)
        labels2d[take] = t + 1
        dist_best[take] = d[take]
        if (t + 1) % 10 == 0:
            print(f"  turn {t+1}/{n_turns}", flush=True)

    labels = np.repeat(labels2d[None, :, :], nz, axis=0)
    intensity = np.where(labels > 0, 90, 0).astype(np.uint8)

    os.makedirs(args.out, exist_ok=True)
    np.save(os.path.join(args.out, "labels.npy"), labels)
    np.save(os.path.join(args.out, "intensity.npy"), intensity)
    rec = {
        "columns": args.columns, "voxel_um": args.voxel_um, "z_mm": args.z_mm,
        "grid": [int(nz), int(ny), int(nx)], "n_turns": n_turns,
        "sheet_um": g["sheet_um"], "pitch_um": g["pitch_um"],
        "fuse": list(fuse) if fuse else None,
        "half_t_vox": round(half_vox, 3),
    }
    json.dump(rec, open(os.path.join(args.out, "weld_truth.json"), "w"),
              indent=1)
    n_inst = len(np.unique(labels2d[labels2d > 0]))
    fg = int((labels2d > 0).sum())
    print(f"done: {n_inst} instances, fg {fg:,} px/slice, "
          f"saved to {args.out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--columns", type=int, default=24)
    ap.add_argument("--voxel-um", type=float, default=20.0)
    ap.add_argument("--z-mm", type=float, default=2.0)
    ap.add_argument("--fuse", default=None,
                    help="ti,tj,ang0,ang1: weld turns ti..tj-1 onto tj over the angle range")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    paint(args)


main()
