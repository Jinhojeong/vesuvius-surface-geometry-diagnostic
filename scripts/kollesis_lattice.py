#!/usr/bin/env python3
"""kollesis_lattice.py - test Diego's kollesis-lattice hypothesis on our census.

Diego-dcv (villa #191, comment 5200998841): kollesis joins are benign
double-thickness by manufacture and repeat every kollema width (~160 mm)
along the arc, so flagged columns that sit on a ~160 mm lattice in arc
coordinates are candidates for a benign must-cross third state in the
hazard weight volume.

Pipeline (all inputs generated on the box, no external downloads):
  1. crossings: pitch_qa_ray_positions.csv.gz (iyando's pitch_qa.py
     --positions on the stitched p1218_full run; z, theta_deg, k,
     r_l1_vox, r_um; 6-deg rays, 8 slices per slab)
  2. origins: slice_origins.csv (per-slice mask centroid, the exact frame
     pitch_qa casts rays from)
  3. census: census8k tile jsons (cluster centroids + thickness ratio)
  4. rings per (z, theta): crossings sorted by descending r (ring 0 =
     outermost), ring perimeters via Diego's arc_fraction (chord sum)
  5. arc position s per cluster: anchored at the outermost resolved ring,
     s = sum(per_0..per_{j-1}) + frac_j(theta) * per_j, only where rings
     0..j are complete at that slice (exclusions counted)
  6. candidates: thickness ratio in the join band (1.6-2.6x, Diego's
     detect() band) AND z-coherent (a manufacture seam runs the full
     height; local fusion does not)
  7. blind fence_fit (Diego's mode 13, verbatim grid: W 120-200 step 0.5,
     offset step 2 mm, tol 0.12*W, ties prefer larger W)
  8. nulls: (A) 200 redraws of the candidate count from the non-band
     cluster arc pool (inherits coverage), full refit each; (B) offset
     permutation at the fitted W.

Attribution: fence_fit/on_fence/arc_fraction ported from Diego-dcv's
vesuvius-topological-grid (kollesis_detector.py, displacement_field.py).
A vectorized fence_fit is used for the nulls after being asserted equal
to the verbatim implementation on the observed data.
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import json
import math
import os
import re
from collections import defaultdict

import numpy as np

UM_PER_VOX = 17.28  # L1, pitch_qa convention
TILE_RE = re.compile(r"tile_y(\d+)_x(\d+)\.json$")
JOIN_BAND = (1.6, 2.6)
Z_LO, Z_HI = 1000, 11000


# ---------------------------------------------------- Diego's code, verbatim

def fence_fit(svals, wmin=120.0, wmax=200.0, tol=0.12):
    svals = np.asarray(svals, float)
    best = None
    for W in np.arange(wmin, wmax + 0.25, 0.5):
        for o in np.arange(0.0, W, 2.0):
            r = np.abs(((svals - o + W / 2) % W) - W / 2)
            score = (round(float(np.mean(r <= tol * W)), 3),
                     -float(r.mean()), W)
            if best is None or score > best[0]:
                best = (score, W, o)
    return best[1], best[2], best[0]


def on_fence(s, W, o, tol=0.12):
    return abs(((s - o + W / 2) % W) - W / 2) <= tol * W


def arc_fraction(theta_deg, r):
    th = np.radians(np.asarray(theta_deg, float))
    r = np.asarray(r, float)
    thc = np.concatenate([th, th[:1]])
    rc = np.concatenate([r, r[:1]])
    x, y = rc * np.cos(thc), rc * np.sin(thc)
    seg = np.hypot(np.diff(x), np.diff(y))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    return s[:-1] / s[-1], s[-1]


# ------------------------------------------------- vectorized fit (for nulls)

def fence_fit_fast(svals, wmin=120.0, wmax=200.0, tol=0.12):
    svals = np.asarray(svals, float)
    best = None
    for W in np.arange(wmin, wmax + 0.25, 0.5):
        offs = np.arange(0.0, W, 2.0)
        r = np.abs(((svals[:, None] - offs[None, :] + W / 2) % W) - W / 2)
        frac = np.round((r <= tol * W).mean(axis=0), 3)
        mres = r.mean(axis=0)
        i = int(np.lexsort((-mres, frac))[-1])
        score = (float(frac[i]), -float(mres[i]), W)
        if best is None or score > best[0]:
            best = (score, W, float(offs[i]))
    return best[1], best[2], best[0]


# --------------------------------------------------------------- input loads

def load_origins(path):
    org = {}
    for row in csv.DictReader(open(path)):
        org[int(float(row["z"]))] = (float(row["cy"]), float(row["cx"]))
    return org


def load_rings(path):
    """(z, theta) -> descending-sorted crossing radii in um."""
    d = defaultdict(list)
    with gzip.open(path, "rt", newline="") as f:
        for row in csv.DictReader(f):
            z = int(float(row["z"]))
            if not (Z_LO <= z <= Z_HI):
                continue
            d[(z, float(row["theta_deg"]))].append(float(row["r_um"]))
    rings = {k: sorted(v, reverse=True) for k, v in d.items()}
    zs = sorted({k[0] for k in rings})
    ths = sorted({k[1] for k in rings})
    return rings, zs, ths


def ring_tables(rings, zs, ths):
    """Per z: complete-ring count J(z), per-ring perimeter and cumulative
    outer arc, per-ring arc-fraction interpolator inputs."""
    tabs = {}
    for z in zs:
        depth = min((len(rings.get((z, t), [])) for t in ths), default=0)
        if depth == 0:
            continue
        R = np.array([[rings[(z, t)][j] for t in ths] for j in range(depth)])
        pers, fracs = [], []
        ok = True
        for j in range(depth):
            r = R[j]
            if not np.isfinite(r).all() or (r <= 0).any():
                depth = j
                ok = j > 0
                break
            f, per = arc_fraction(ths, r)
            pers.append(per / 1000.0)  # um -> mm
            fracs.append(f)
        if not ok or not pers:
            continue
        cum = np.concatenate([[0.0], np.cumsum(pers)])  # outer anchor
        tabs[z] = {"depth": len(pers), "pers": pers, "cum": cum,
                   "fracs": fracs, "R": R[:len(pers)]}
    return tabs


def scan_census(census_dir):
    rows = []
    for p in sorted(glob.glob(os.path.join(census_dir, "z*", "tile_*.json"))):
        slab = os.path.basename(os.path.dirname(p))
        m = TILE_RE.search(p)
        if not m or not slab.startswith("z"):
            continue
        z0, y0, x0 = int(slab[1:]), int(m.group(1)), int(m.group(2))
        rec = json.load(open(p))
        if rec.get("thin"):
            continue
        for c in rec.get("clusters", []):
            rows.append((z0 + c["z"], y0 + c["y"], x0 + c["x"],
                         c["ratio"], c["th"], int(c.get("n_sites", 1))))
    return rows


# ------------------------------------------------------------ arc conversion

def convert(rows, origins, rings, tabs, zs, ths, max_dr_um=260.0):
    zarr_ = np.array(zs)
    tharr = np.array(ths)
    out = []
    drops = defaultdict(int)
    for gz, gy, gx, ratio, th, ns in rows:
        if not (Z_LO <= gz <= Z_HI):
            drops["z_range"] += 1
            continue
        z = int(zarr_[np.argmin(np.abs(zarr_ - gz))])
        if abs(z - gz) > 20 or z not in tabs or z not in origins:
            drops["no_slice"] += 1
            continue
        cy, cx = origins[z]
        r_um = math.hypot(gy - cy, gx - cx) * UM_PER_VOX
        theta = math.degrees(math.atan2(gy - cy, gx - cx)) % 360.0
        it = int(np.argmin(np.minimum(np.abs(tharr - theta),
                                      360 - np.abs(tharr - theta))))
        t = ths[it]
        cr = rings.get((z, t))
        if not cr:
            drops["no_ray"] += 1
            continue
        tab = tabs[z]
        j = int(np.argmin([abs(r_um - c) for c in cr]))
        if abs(r_um - cr[j]) > max_dr_um:
            drops["dr_too_far"] += 1
            continue
        if j >= tab["depth"]:
            drops["ring_unresolved"] += 1
            continue
        frac = float(np.interp(theta, np.array(ths), tab["fracs"][j],
                               period=360.0))
        s_mm = float(tab["cum"][j] + frac * tab["pers"][j])
        out.append({"s": s_mm, "z": z, "ratio": ratio, "th": th,
                    "n_sites": ns, "ring": j, "dr_um": abs(r_um - cr[j]),
                    "gz": gz, "gy": gy, "gx": gx})
    return out, dict(drops)


# ------------------------------------------------------------- stratification

def columns(cands, bin_mm=5.0, min_slices=6, min_span=50.0):
    """z-coherent columns: s-bin with wide z support."""
    bins = defaultdict(list)
    for c in cands:
        bins[int(c["s"] // bin_mm)].append(c)
    reps, members = [], []
    for b, cs in sorted(bins.items()):
        zsl = sorted({c["z"] for c in cs})
        if len(zsl) >= min_slices and (zsl[-1] - zsl[0]) * UM_PER_VOX / 1000.0 >= min_span:
            reps.append(float(np.median([c["s"] for c in cs])))
            members.append({"s_bin_mm": b * bin_mm, "n_clusters": len(cs),
                            "n_slices": len(zsl),
                            "z_span_mm": round((zsl[-1] - zsl[0]) * UM_PER_VOX / 1000.0, 1)})
    return reps, members


# --------------------------------------------------------------------- nulls

def null_redraw(pool, n, rng, draws=200):
    scores = []
    pool = np.asarray(pool, float)
    for _ in range(draws):
        s = rng.choice(pool, size=n, replace=False) if len(pool) >= n else \
            rng.choice(pool, size=n, replace=True)
        _, _, sc = fence_fit_fast(s)
        scores.append(sc[0])
    return scores


def null_offsets(svals, W, rng, draws=2000, tol=0.12):
    svals = np.asarray(svals, float)
    offs = rng.uniform(0, W, size=draws)
    r = np.abs(((svals[:, None] - offs[None, :] + W / 2) % W) - W / 2)
    return (r <= tol * W).mean(axis=0).tolist()


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crossings",
                    default="/mnt/vesuvius/p1218_full/pitch_qa_ray_positions.csv.gz")
    ap.add_argument("--origins",
                    default="/mnt/vesuvius/p1218_full/slice_origins.csv")
    ap.add_argument("--census-dir", default="/mnt/vesuvius/census8k")
    ap.add_argument("--out", default="/mnt/vesuvius/kollesis")
    ap.add_argument("--seed", type=int, default=20260807)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    origins = load_origins(args.origins)
    rings, zs, ths = load_rings(args.crossings)
    print(f"[crossings] {sum(len(v) for v in rings.values()):,} crossings, "
          f"{len(zs)} slices, {len(ths)} rays", flush=True)
    tabs = ring_tables(rings, zs, ths)
    depths = [t["depth"] for t in tabs.values()]
    print(f"[rings] {len(tabs)} slices with complete rings, depth "
          f"median {np.median(depths):.0f} range {min(depths)}-{max(depths)}",
          flush=True)

    rows = scan_census(args.census_dir)
    conv, drops = convert(rows, origins, rings, tabs, zs, ths)
    print(f"[convert] {len(rows):,} clusters -> {len(conv):,} converted; "
          f"drops {drops}", flush=True)

    cands = [c for c in conv if JOIN_BAND[0] <= c["ratio"] <= JOIN_BAND[1]]
    bg = [c["s"] for c in conv if not (JOIN_BAND[0] <= c["ratio"] <= JOIN_BAND[1])]
    print(f"[strata] join-band {len(cands):,}, background {len(bg):,}",
          flush=True)

    reps, cols = columns(cands)
    print(f"[columns] {len(reps)} z-coherent columns", flush=True)

    report = {"inputs": {"crossings": args.crossings, "n_slices": len(zs),
                         "n_rays": len(ths), "ring_depth_median": float(np.median(depths))},
              "convert": {"total": len(rows), "converted": len(conv),
                          "drops": drops},
              "strata": {"join_band": len(cands), "background": len(bg),
                         "columns": cols}}

    fit_target = reps if len(reps) >= 8 else [c["s"] for c in cands]
    report["fit_on"] = "columns" if len(reps) >= 8 else "all_band_clusters"
    if len(fit_target) < 8:
        report["verdict"] = "insufficient candidates"
        json.dump(report, open(os.path.join(args.out, "lattice_report.json"), "w"), indent=1)
        print("[verdict] insufficient candidates", flush=True)
        return

    W, o, sc = fence_fit(fit_target)
    Wf, of, scf = fence_fit_fast(fit_target)
    assert (W, o) == (Wf, of) and abs(sc[0] - scf[0]) < 1e-9, \
        f"fast fit mismatch: {(W, o, sc)} vs {(Wf, of, scf)}"
    obs = sc[0]
    print(f"[fence] W={W:.1f}mm offset={o:.1f} on-fence frac={obs:.3f} "
          f"(n={len(fit_target)})", flush=True)

    nA = null_redraw(bg, len(fit_target), rng)
    pA = float(np.mean([s >= obs for s in nA]))
    nB = null_offsets(fit_target, W, rng)
    pB = float(np.mean([s >= obs for s in nB]))
    print(f"[nulls] redraw p={pA:.3f} (null mean {np.mean(nA):.3f}); "
          f"offset-perm p={pB:.3f} (null mean {np.mean(nB):.3f})", flush=True)

    report["fence"] = {"W_mm": W, "offset_mm": o, "on_fence_frac": obs,
                       "n_fit": len(fit_target)}
    report["null_redraw"] = {"draws": len(nA), "p": pA,
                             "mean": float(np.mean(nA)),
                             "q95": float(np.quantile(nA, 0.95))}
    report["null_offsets"] = {"draws": len(nB), "p": pB,
                              "mean": float(np.mean(nB))}
    in_range = 150.0 <= W <= 170.0
    positive = in_range and pA < 0.05 and pB < 0.05
    report["verdict"] = ("lattice_present" if positive else
                        "no_lattice_detected")
    json.dump(report, open(os.path.join(args.out, "lattice_report.json"), "w"),
              indent=1)
    print(f"[verdict] {report['verdict']} (W in 150-170: {in_range})",
          flush=True)


if __name__ == "__main__":
    main()
