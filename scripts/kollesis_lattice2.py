#!/usr/bin/env python3
"""kollesis_lattice2.py - mass-based contrast test for the kollesis lattice.

Why v2. The v1 column stratification saturated: 262k join-band clusters fill
essentially every 5 mm arc bin, so an on-fence fraction over occupied bins
degenerates to coverage geometry and cannot see ~27 joins even if present,
and its redraw null compares different objects (bin representatives vs raw
cluster draws). Diego's twin fence worked on 47 sparse detector flags; our
census flags all fusion, a different regime.

V2 design. Compare cluster MASS at the fence between the join-band stratum
(thickness ratio 1.6-2.6, his detect() band) and the background stratum
(everything else), which shares the same coverage mask by construction:

    stat = frac(band clusters on fence) - frac(background clusters on fence)

with the fence (W, offset) fitted on the band stratum. Null = permutation of
the band/background labels over all converted clusters at the fitted fence
(the fence is selected on the observed band, which biases the observed stat
UPWARD; a null result despite that bias is conservative evidence of absence).
If the label-permutation comes out positive, a full refit-per-permutation
null is required before any claim.

Also reports the same contrast at W fixed to 160.0 mm (Diego's manufacture
prior, offset fitted), which removes the W-selection freedom entirely, and
the per-cluster converted table is written out so the whole thing can be
rerun from the artifact.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os

import numpy as np

import kollesis_lattice as K


def fit_offset_only(svals, W, tol=0.12):
    offs = np.arange(0.0, W, 2.0)
    r = np.abs(((np.asarray(svals, float)[:, None] - offs[None, :] + W / 2)
                % W) - W / 2)
    frac = (r <= tol * W).mean(axis=0)
    i = int(np.argmax(frac))
    return float(offs[i]), float(frac[i])


def on_frac(svals, W, o, tol=0.12):
    r = np.abs(((np.asarray(svals, float) - o + W / 2) % W) - W / 2)
    return float((r <= tol * W).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crossings",
                    default="/mnt/vesuvius/kollesis/positions_merged.csv.gz")
    ap.add_argument("--origins",
                    default="/mnt/vesuvius/kollesis/origins_merged.csv")
    ap.add_argument("--census-dir", default="/mnt/vesuvius/census8k")
    ap.add_argument("--out", default="/mnt/vesuvius/kollesis")
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260808)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    origins = K.load_origins(args.origins)
    rings, zs, ths = K.load_rings(args.crossings)
    tabs = K.ring_tables(rings, zs, ths)
    rows = K.scan_census(args.census_dir)
    conv, drops = K.convert(rows, origins, rings, tabs, zs, ths)
    print(f"[convert] {len(conv):,} clusters (drops {drops})", flush=True)

    with gzip.open(os.path.join(args.out, "clusters_arc.csv.gz"), "wt",
                   newline="") as f:
        f.write("s_mm,z,ratio,th,n_sites,ring,dr_um,gz,gy,gx\n")
        for c in conv:
            f.write(f"{c['s']:.2f},{c['z']},{c['ratio']},{c['th']},"
                    f"{c['n_sites']},{c['ring']},{c['dr_um']:.1f},"
                    f"{c['gz']},{c['gy']},{c['gx']}\n")

    s_all = np.array([c["s"] for c in conv])
    band = np.array([K.JOIN_BAND[0] <= c["ratio"] <= K.JOIN_BAND[1]
                     for c in conv])
    sA, sB = s_all[band], s_all[~band]
    print(f"[strata] band {len(sA):,} background {len(sB):,}", flush=True)

    report = {"design": "mass contrast band-vs-background at the fence; "
                        "label-permutation null (anti-conservative for the "
                        "fitted-W arm, see docstring)",
              "n_band": int(len(sA)), "n_background": int(len(sB)),
              "drops": drops, "arms": {}}

    for arm, (W, o, note) in {
        "fitted": (*K.fence_fit_fast(sA)[:2],
                   "W and offset fitted on the band stratum"),
        "prior160": (160.0, None, "W fixed at 160.0 mm (manufacture prior), "
                                  "offset fitted on the band stratum"),
    }.items():
        if o is None:
            o, _ = fit_offset_only(sA, W)
        fA, fB = on_frac(sA, W, o), on_frac(sB, W, o)
        obs = fA - fB
        lab = band.copy()
        null = np.empty(args.nperm)
        for i in range(args.nperm):
            rng.shuffle(lab)
            null[i] = on_frac(s_all[lab], W, o) - on_frac(s_all[~lab], W, o)
        p = float((null >= obs).mean())
        report["arms"][arm] = {
            "W_mm": W, "offset_mm": o, "note": note,
            "frac_band": fA, "frac_background": fB,
            "contrast_pp": 100 * obs, "perm_p": p,
            "null_mean_pp": 100 * float(null.mean()),
            "null_sd_pp": 100 * float(null.std()),
            "null_q95_pp": 100 * float(np.quantile(null, 0.95)),
        }
        print(f"[{arm}] W={W:.1f} o={o:.1f} band {fA:.4f} bg {fB:.4f} "
              f"contrast {100 * obs:+.3f}pp perm_p={p:.4f} "
              f"(null sd {100 * null.std():.3f}pp)", flush=True)

    json.dump(report, open(os.path.join(args.out, "lattice2_report.json"),
                           "w"), indent=1)
    print("[done] lattice2_report.json", flush=True)


if __name__ == "__main__":
    main()
