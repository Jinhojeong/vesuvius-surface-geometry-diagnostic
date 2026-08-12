#!/usr/bin/env python3
"""Scroll-wide label coverage of the m7-predicted sheet surface, PHerc1218.

For each label tile under /mnt/vesuvius/p1218_repair_v3/blocks_repaired,
take the fixed center window (tile-local z 64:192, y 192:320, x 192:320),
fetch the matching m7 L1 crop (label grid == m7 prediction level 1, zero
offset), and measure what fraction of predicted-surface voxels carry an
instance label (strict and within Chebyshev 4).

Outputs:
  /mnt/vesuvius/labelcov1218/coverage_tiles.csv    per-tile rows
  /mnt/vesuvius/labelcov1218/coverage_summary.json aggregate + anchors
  /mnt/vesuvius/labelcov1218/coverage.log          progress log (stdout)

Deterministic: fixed window, fixed thresholds (pred = m7 >= 128,
lab = labels > 0, near = Chebyshev distance <= 4 via maximum_filter size 9),
tiles processed in sorted path order.
"""

import csv
import glob
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import maximum_filter

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

M7_URL = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
          "representations/predictions/surfaces/"
          "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")
CACHE_DIR = "/mnt/vesuvius/hazard_zarr_smoke/m7L1_cache"
BLOCKS_GLOB = "/mnt/vesuvius/p1218_repair_v3/blocks_repaired/z*/tile_*.npz"
OUT_DIR = "/mnt/vesuvius/labelcov1218"
DEMO_SITES = "/mnt/vesuvius/hazard_zarr_smoke/demo_sites.json"
OFFSEEDS = "/mnt/vesuvius/hazard_zarr_smoke/offseeds_placement_v2.json"

WIN_LO = (64, 192, 192)   # tile-local window lo (z, y, x)
WIN_HI = (192, 320, 320)  # tile-local window hi
PRED_TH = 128
NEAR_CHEB = 4
MIN_PRED_FOR_DIST = 1000
LOW_COV = 0.1
N_Z_BANDS = 10
MAX_FAIL_FRAC = 0.10

TILE_SHAPE = (256, 512, 512)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def window_metrics(lab_bool, m7_crop):
    """Counts inside one window. lab_bool and m7_crop have identical shape."""
    pred = m7_crop >= PRED_TH
    near = maximum_filter(lab_bool.astype(np.uint8),
                          size=2 * NEAR_CHEB + 1) > 0
    n_pred = int(pred.sum())
    n_lab = int(lab_bool.sum())
    n_and = int((pred & lab_bool).sum())
    n_near = int((pred & near).sum())
    return n_pred, n_lab, n_and, n_near


def fetch_crop(lvl, lo, hi):
    """read_crop with one extra attempt on failure (reader retries inside)."""
    try:
        return lvl.read_crop(lo, hi)
    except Exception as e:
        log(f"    fetch retry after: {e}")
        time.sleep(5.0)
        return lvl.read_crop(lo, hi)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    lvl = RemoteZarrLevel(M7_URL, 1, cache_dir=CACHE_DIR)
    log(f"m7 L1 shape={lvl.shape} chunks={lvl.chunks} dtype={lvl.dtype}")

    tiles = sorted(glob.glob(BLOCKS_GLOB))
    if len(sys.argv) > 1:  # smoke-run limit
        tiles = tiles[:int(sys.argv[1])]
    log(f"{len(tiles)} tiles")

    csv_path = os.path.join(OUT_DIR, "coverage_tiles.csv")
    fields = ["path", "slab", "tile", "z0", "y0", "x0",
              "win_lo_z", "win_lo_y", "win_lo_x",
              "win_hi_z", "win_hi_y", "win_hi_x", "clamped", "n_vox",
              "n_pred", "n_lab", "n_pred_and_lab", "n_pred_near_lab",
              "coverage_strict", "coverage_near4", "reverse_strict",
              "status"]
    rows = []
    n_fail = 0
    t0 = time.time()

    with open(csv_path, "w", newline="") as cf:
        w = csv.DictWriter(cf, fieldnames=fields)
        w.writeheader()
        for i, path in enumerate(tiles):
            rel = path.split("blocks_repaired/")[1]
            slab, tile = rel.split("/")
            tile = tile[:-4]
            row = {"path": rel, "slab": slab, "tile": tile}
            try:
                with np.load(path) as f:
                    z0 = int(f["z0"]); y0 = int(f["y0"]); x0 = int(f["x0"])
                    labels = f["labels"]
            except Exception as e:
                log(f"  {rel}: npz load failed: {e}")
                row["status"] = "npz_failed"
                w.writerow(row); cf.flush()
                rows.append(row)
                n_fail += 1
                continue
            origin = (z0, y0, x0)
            row.update(z0=z0, y0=y0, x0=x0)
            glo = tuple(o + l for o, l in zip(origin, WIN_LO))
            ghi = tuple(o + h for o, h in zip(origin, WIN_HI))
            # clamp against the m7/label grid extent
            ghi_c = tuple(min(h, s) for h, s in zip(ghi, lvl.shape))
            clamped = ghi_c != ghi
            if any(l >= h for l, h in zip(glo, ghi_c)):
                row["status"] = "empty_after_clamp"
                w.writerow(row); cf.flush()
                rows.append(row)
                continue
            row.update(win_lo_z=glo[0], win_lo_y=glo[1], win_lo_x=glo[2],
                       win_hi_z=ghi_c[0], win_hi_y=ghi_c[1],
                       win_hi_x=ghi_c[2], clamped=int(clamped))
            lz, ly, lx = WIN_LO
            dz, dy, dx = (h - l for l, h in zip(glo, ghi_c))
            lab = labels[lz:lz + dz, ly:ly + dy, lx:lx + dx] > 0
            del labels
            try:
                m7 = fetch_crop(lvl, glo, ghi_c)
            except Exception as e:
                log(f"  {rel}: fetch_failed: {e}")
                row["status"] = "fetch_failed"
                w.writerow(row); cf.flush()
                rows.append(row)
                n_fail += 1
                if n_fail > MAX_FAIL_FRAC * len(tiles):
                    log(f"ABORT: {n_fail} failures exceed "
                        f"{MAX_FAIL_FRAC:.0%} of {len(tiles)} tiles")
                    summary = {"status": "aborted_failures",
                               "n_tiles": len(tiles), "n_done": i + 1,
                               "n_fail": n_fail}
                    with open(os.path.join(OUT_DIR,
                                           "coverage_summary.json"),
                              "w") as sf:
                        json.dump(summary, sf, indent=1)
                    sys.exit(2)
                continue
            n_pred, n_lab, n_and, n_near = window_metrics(lab, m7)
            row.update(
                n_vox=int(dz * dy * dx),
                n_pred=n_pred, n_lab=n_lab,
                n_pred_and_lab=n_and, n_pred_near_lab=n_near,
                coverage_strict=(round(n_and / n_pred, 6)
                                 if n_pred else None),
                coverage_near4=(round(n_near / n_pred, 6)
                                if n_pred else None),
                reverse_strict=(round(n_and / n_lab, 6)
                                if n_lab else None),
                status="ok")
            w.writerow(row); cf.flush()
            rows.append(row)
            if (i + 1) % 25 == 0 or i + 1 == len(tiles):
                el = time.time() - t0
                log(f"  {i + 1}/{len(tiles)} tiles, {n_fail} failed, "
                    f"{el / 60:.1f} min elapsed, "
                    f"{el / (i + 1):.2f} s/tile")

    ok = [r for r in rows if r.get("status") == "ok"]
    log(f"main pass done: {len(ok)} ok, {n_fail} failed, "
        f"{len(rows) - len(ok) - n_fail} other")

    # ---- pooled global numbers ------------------------------------------
    tot_pred = sum(r["n_pred"] for r in ok)
    tot_lab = sum(r["n_lab"] for r in ok)
    tot_and = sum(r["n_pred_and_lab"] for r in ok)
    tot_near = sum(r["n_pred_near_lab"] for r in ok)
    tot_vox = sum(r["n_vox"] for r in ok)
    pooled = {
        "n_windows": len(ok),
        "n_vox": tot_vox,
        "n_pred": tot_pred,
        "n_lab": tot_lab,
        "n_pred_and_lab": tot_and,
        "n_pred_near_lab": tot_near,
        "P_label_given_pred": tot_and / tot_pred if tot_pred else None,
        "P_label_near4_given_pred": tot_near / tot_pred if tot_pred else None,
        "P_pred_given_label": tot_and / tot_lab if tot_lab else None,
    }

    # ---- tile distribution (n_pred >= MIN_PRED_FOR_DIST) ------------------
    dist_tiles = [r for r in ok if r["n_pred"] >= MIN_PRED_FOR_DIST]
    covs = np.array([r["coverage_near4"] for r in dist_tiles])
    deciles = (np.percentile(covs, np.arange(0, 101, 10)).round(4).tolist()
               if len(covs) else [])
    dist = {
        "min_pred": MIN_PRED_FOR_DIST,
        "n_tiles": len(dist_tiles),
        "coverage_near4_deciles_p0_to_p100": deciles,
        "n_below_0.1": int((covs < 0.1).sum()),
        "n_0.1_to_0.5": int(((covs >= 0.1) & (covs < 0.5)).sum()),
        "n_0.5_to_0.9": int(((covs >= 0.5) & (covs <= 0.9)).sum()),
        "n_above_0.9": int((covs > 0.9).sum()),
    }

    # ---- z-band structure of low-coverage tiles ---------------------------
    zc = np.array([(r["win_lo_z"] + r["win_hi_z"]) / 2 for r in dist_tiles])
    low = covs < LOW_COV
    zmin, zmax = (float(zc.min()), float(zc.max())) if len(zc) else (0, 1)
    edges = np.linspace(zmin, zmax + 1e-9, N_Z_BANDS + 1)
    z_bands = []
    for b in range(N_Z_BANDS):
        m = (zc >= edges[b]) & (zc < edges[b + 1])
        z_bands.append({
            "z_lo": round(float(edges[b]), 1),
            "z_hi": round(float(edges[b + 1]), 1),
            "n_tiles": int(m.sum()),
            "n_low_cov": int((m & low).sum()),
            "frac_low": (round(float((m & low).sum() / m.sum()), 4)
                         if m.sum() else None),
            "median_cov_near4": (round(float(np.median(covs[m])), 4)
                                 if m.sum() else None),
        })

    # ---- radial structure (axis approximated by volume xy center) ---------
    cy = lvl.shape[1] / 2.0
    cx = lvl.shape[2] / 2.0
    yc = np.array([(r["win_lo_y"] + r["win_hi_y"]) / 2 for r in dist_tiles])
    xc = np.array([(r["win_lo_x"] + r["win_hi_x"]) / 2 for r in dist_tiles])
    rr = np.sqrt((yc - cy) ** 2 + (xc - cx) ** 2)
    r_bands = []
    if len(rr):
        redges = np.percentile(rr, [0, 33.34, 66.67, 100.0])
        redges[-1] += 1e-9
        for b in range(3):
            m = (rr >= redges[b]) & (rr < redges[b + 1])
            r_bands.append({
                "r_lo": round(float(redges[b]), 1),
                "r_hi": round(float(redges[b + 1]), 1),
                "n_tiles": int(m.sum()),
                "n_low_cov": int((m & low).sum()),
                "frac_low": (round(float((m & low).sum() / m.sum()), 4)
                             if m.sum() else None),
            })

    # ---- sanity anchors ----------------------------------------------------
    # (a) center-window coverage of the tiles named in the two site files
    # (b) 128-cubes centered on each on-seed / off-seed global point,
    #     labels read from the named tile (window intersected w/ tile extent)
    with open(DEMO_SITES) as f:
        onseeds = json.load(f)
    with open(OFFSEEDS) as f:
        offseeds = json.load(f)

    by_key = {(r["slab"], r["tile"]): r for r in rows}
    tile_origin = {}

    def anchor_entry(site, kind):
        key = (site["slab"], site["tile"])
        g = (site["gz"], site["gy"], site["gx"])
        ent = {"kind": kind, "slab": site["slab"], "tile": site["tile"],
               "g_L1": list(g)}
        tr = by_key.get(key)
        ent["tile_center_window_coverage_near4"] = (
            tr.get("coverage_near4") if tr else None)
        ent["tile_center_window_n_pred"] = tr.get("n_pred") if tr else None
        # point-centered 128-cube, intersected with tile extent and grid
        if key not in tile_origin:
            p = (f"/mnt/vesuvius/p1218_repair_v3/blocks_repaired/"
                 f"{site['slab']}/{site['tile']}.npz")
            with np.load(p) as f:
                tile_origin[key] = ((int(f["z0"]), int(f["y0"]),
                                     int(f["x0"])), f["labels"] > 0)
        origin, labbool = tile_origin[key]
        lo = [max(g[a] - 64, origin[a], 0) for a in range(3)]
        hi = [min(g[a] + 64, origin[a] + TILE_SHAPE[a], lvl.shape[a])
              for a in range(3)]
        if any(l >= h for l, h in zip(lo, hi)):
            ent["point_window"] = None
            return ent
        lab = labbool[lo[0] - origin[0]:hi[0] - origin[0],
                      lo[1] - origin[1]:hi[1] - origin[1],
                      lo[2] - origin[2]:hi[2] - origin[2]]
        try:
            m7 = fetch_crop(lvl, tuple(lo), tuple(hi))
        except Exception as e:
            ent["point_window"] = f"fetch_failed: {e}"
            return ent
        n_pred, n_lab, n_and, n_near = window_metrics(lab, m7)
        ent["point_window"] = {
            "lo": lo, "hi": hi,
            "n_pred": n_pred, "n_lab": n_lab,
            "coverage_strict": (round(n_and / n_pred, 6)
                                if n_pred else None),
            "coverage_near4": (round(n_near / n_pred, 6)
                               if n_pred else None),
        }
        return ent

    anchors = []
    log("anchor windows (on-seed sites)")
    for s in onseeds:
        anchors.append(anchor_entry(s, "on_seed"))
    log("anchor windows (off-seed points)")
    for s in offseeds:
        anchors.append(anchor_entry(s, "off_seed"))
    # free anchor label arrays
    tile_origin.clear()

    summary = {
        "status": "ok",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "method": {
            "m7_url": M7_URL,
            "m7_level": 1,
            "grid_shape": list(lvl.shape),
            "window_tile_local": {"z": [64, 192], "y": [192, 320],
                                  "x": [192, 320]},
            "pred_threshold": PRED_TH,
            "near_chebyshev": NEAR_CHEB,
            "note": ("One fixed 128^3 center window per tile, clamped at "
                     "the grid boundary; not a whole-volume integral. "
                     "Near-window-edge labels outside the window are not "
                     "seen by the dilation."),
        },
        "n_tiles_total": len(tiles),
        "n_tiles_ok": len(ok),
        "n_tiles_failed": n_fail,
        "n_tiles_other": len(rows) - len(ok) - n_fail,
        "n_tiles_clamped": sum(1 for r in ok if r.get("clamped")),
        "pooled": pooled,
        "tile_distribution": dist,
        "z_bands_low_coverage": z_bands,
        "radial_bands_low_coverage_axis_approx_volume_center": r_bands,
        "anchors": anchors,
    }
    spath = os.path.join(OUT_DIR, "coverage_summary.json")
    with open(spath, "w") as sf:
        json.dump(summary, sf, indent=1)
    log(f"summary written to {spath}")
    log(f"total runtime {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
