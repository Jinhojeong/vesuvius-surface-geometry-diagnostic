#!/usr/bin/env python3
"""CT-based partial discriminant for unlabeled m7 surface, PHerc1218.

Question (villa #193): the m7-predicted surface is only ~45% covered by
instance labels. Is the unlabeled remainder real papyrus the labels missed,
or m7 over-prediction? CT intensity/texture over real papyrus differs from
air / non-sheet fill, so measure CT statistics at m7-positive voxels in
low-coverage (VOID) windows against fully-labeled (CONTROL) windows.

Sampling (deterministic, declared):
  pool  = coverage_tiles.csv rows with status ok, clamped == 0,
          n_pred >= 100000
  VOID  = pool with coverage_near4 < 0.1   -> 60 windows
  CTRL  = pool with coverage_near4 > 0.98  -> 30 windows
  Stratification: sort candidates by (z-band decile, radial third, path),
  take every k-th (k = n_candidates // n_target), first n_target.
  z-band deciles over pool window-center z range; radial terciles over the
  pool distribution of in-plane distance from the grid center (1898.5,
  1898.5).

Per window: the coverage run's 128^3 window is reconstructed from the csv
win_lo/win_hi columns (L1 grid coords, label grid == m7 L1 == CT L1, zero
offset). A 64^3 cube is placed inside the window so it lies in a single CT
chunk (CT chunks are 128^3 u1 uncompressed, so this bounds the fetch at one
2 MB chunk per window): per axis, with w = win_lo and b = next multiple of
128 above w, cube_lo = w+32 if b == w else (b-64 if b-w >= 64 else b).

Statistics per window over m7-positive voxels (m7 >= 128) in the cube:
  mean/std CT intensity, 256-bin histogram, and a sheet-texture proxy:
  fraction of positive voxels within 4 steps (along the most-radial in-plane
  axis, and separately along z) of a significant intensity extremum
  (alternation of |diff| >= 3 sign). CONTROL windows additionally split
  positives into labeled vs unlabeled. m7-negative voxel stats are kept as
  an air/background reference.

Outputs under /mnt/vesuvius/voidct1218/:
  windows.csv        per-window stats
  window_hists.npz   per-window 256-bin histograms
  summary.json       selection rule, CT URL, sanity check, references
  run.log            progress log

Resource discipline: 3 worker threads, 0.25 s stagger between fetch starts,
CT cache under this directory (never the m7 cache dir for CT), run under
nice. CPU-only.
"""

import csv
import json
import math
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from scipy.ndimage import maximum_filter

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

CT_URL = ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
          "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
M7_URL = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
          "representations/predictions/surfaces/"
          "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")
LEVEL = 1
COV_CSV = "/mnt/vesuvius/labelcov1218/coverage_tiles.csv"
BLOCKS_DIR = "/mnt/vesuvius/p1218_repair_v3/blocks_repaired"
OUT_DIR = "/mnt/vesuvius/voidct1218"
CT_CACHE = os.path.join(OUT_DIR, "ct_cache")
M7_CACHE = "/mnt/vesuvius/hazard_zarr_smoke/m7L1_cache"

N_VOID = 60
N_CTRL = 30
MIN_PRED = 100000
VOID_COV = 0.1
CTRL_COV = 0.98
PRED_TH = 128
CUBE = 64
GRID_CENTER = (1898.5, 1898.5)  # (y, x) of the 3797^2 in-plane grid
N_Z_BANDS = 10
TEX_AMP = 3       # min |diff| for a significant intensity step
TEX_NEAR = 4      # extremum within this many voxels along the profile axis
N_WORKERS = 3
STAGGER_S = 0.25

_log_lock = threading.Lock()
_stagger_lock = threading.Lock()
_last_start = [0.0]


def log(msg):
    with _log_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def stagger():
    with _stagger_lock:
        now = time.time()
        wait = _last_start[0] + STAGGER_S - now
        if wait > 0:
            time.sleep(wait)
        _last_start[0] = time.time()


def dir_bytes(d):
    total = 0
    for root, _, files in os.walk(d):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


# ---- selection -----------------------------------------------------------

def load_pool():
    rows = list(csv.DictReader(open(COV_CSV)))
    pool = []
    for r in rows:
        if r["status"] != "ok" or r["clamped"] != "0":
            continue
        if int(r["n_pred"]) < MIN_PRED:
            continue
        r = dict(r)
        r["zc"] = (int(r["win_lo_z"]) + int(r["win_hi_z"])) / 2
        yc = (int(r["win_lo_y"]) + int(r["win_hi_y"])) / 2
        xc = (int(r["win_lo_x"]) + int(r["win_hi_x"])) / 2
        r["rc"] = math.hypot(yc - GRID_CENTER[0], xc - GRID_CENTER[1])
        pool.append(r)
    return pool


def stratify(cands, pool, n_target):
    zs = [r["zc"] for r in pool]
    zmin, zmax = min(zs), max(zs)
    rs = sorted(r["rc"] for r in pool)
    t1, t2 = rs[len(rs) // 3], rs[2 * len(rs) // 3]

    def zband(r):
        return min(N_Z_BANDS - 1,
                   int((r["zc"] - zmin) / (zmax - zmin + 1e-9) * N_Z_BANDS))

    def rthird(r):
        return 0 if r["rc"] < t1 else (1 if r["rc"] < t2 else 2)

    for r in cands:
        r["zband"] = zband(r)
        r["rthird"] = rthird(r)
    cands = sorted(cands, key=lambda r: (r["zband"], r["rthird"], r["path"]))
    if len(cands) <= n_target:
        return cands, (zmin, zmax, t1, t2)
    k = len(cands) // n_target
    return [cands[i * k] for i in range(n_target)], (zmin, zmax, t1, t2)


# ---- per-window ----------------------------------------------------------

def cube_lo_axis(w):
    """64-range within [w, w+128) contained in a single 128-chunk."""
    b = ((w + 127) // 128) * 128
    if b == w:
        return w + 32
    if b - w >= CUBE:
        return b - CUBE
    return b


def tex_frac(ct, pos, axis):
    """Fraction of pos voxels within TEX_NEAR steps (along axis) of a
    significant intensity extremum (sign alternation of |diff|>=TEX_AMP)."""
    d = np.diff(ct.astype(np.int16), axis=axis)
    sgn = np.zeros_like(d, dtype=np.int8)
    sgn[d >= TEX_AMP] = 1
    sgn[d <= -TEX_AMP] = -1
    sl_a = [slice(None)] * 3
    sl_b = [slice(None)] * 3
    sl_a[axis] = slice(None, -1)
    sl_b[axis] = slice(1, None)
    alt = (sgn[tuple(sl_a)] * sgn[tuple(sl_b)] == -1)
    altp = np.zeros(ct.shape, bool)
    sl_c = [slice(None)] * 3
    sl_c[axis] = slice(1, 1 + alt.shape[axis])
    altp[tuple(sl_c)] = alt
    size = [1, 1, 1]
    size[axis] = 2 * TEX_NEAR + 1
    near = maximum_filter(altp, size=tuple(size))
    n = int(pos.sum())
    return float(near[pos].sum() / n) if n else None


def vox_stats(vals):
    if vals.size == 0:
        return None, None, np.zeros(256, dtype=np.int64)
    return (float(vals.mean()), float(vals.std()),
            np.bincount(vals, minlength=256).astype(np.int64))


def process_window(r, ct_lvl, m7_lvl, group):
    stagger()
    lo = [cube_lo_axis(int(r[f"win_lo_{a}"])) for a in "zyx"]
    hi = [c + CUBE for c in lo]
    ct = ct_lvl.read_crop(tuple(lo), tuple(hi))
    m7 = m7_lvl.read_crop(tuple(lo), tuple(hi))
    with np.load(os.path.join(BLOCKS_DIR, r["path"])) as f:
        org = (int(f["z0"]), int(f["y0"]), int(f["x0"]))
        lab_t = f["labels"]
        sl = tuple(slice(l - o, h - o) for l, h, o in zip(lo, hi, org))
        lab = lab_t[sl] > 0
        del lab_t

    pos = m7 >= PRED_TH
    n_pos = int(pos.sum())
    pos_mean, pos_std, pos_hist = vox_stats(ct[pos])
    neg_mean, neg_std, neg_hist = vox_stats(ct[~pos])

    # profile axis: the in-plane axis most aligned with the local radial dir
    yc = (int(r["win_lo_y"]) + int(r["win_hi_y"])) / 2 - GRID_CENTER[0]
    xc = (int(r["win_lo_x"]) + int(r["win_hi_x"])) / 2 - GRID_CENTER[1]
    rad_axis = 1 if abs(yc) >= abs(xc) else 2

    out = {
        "path": r["path"], "group": group,
        "zband": r["zband"], "rthird": r["rthird"],
        "coverage_near4": float(r["coverage_near4"]),
        "n_pred_window": int(r["n_pred"]),
        "cube_lo_z": lo[0], "cube_lo_y": lo[1], "cube_lo_x": lo[2],
        "n_pos": n_pos,
        "n_pos_lab": int((pos & lab).sum()),
        "pos_mean": pos_mean, "pos_std": pos_std,
        "neg_mean": neg_mean, "neg_std": neg_std,
        "frac_zero_all": float((ct == 0).mean()),
        "tex_rad": tex_frac(ct, pos, rad_axis),
        "tex_z": tex_frac(ct, pos, 0),
        "rad_axis": {1: "y", 2: "x"}[rad_axis],
    }
    hists = {"pos": pos_hist, "neg": neg_hist}
    if group == "control":
        lm, ls, lh = vox_stats(ct[pos & lab])
        um, us, uh = vox_stats(ct[pos & ~lab])
        out.update(lab_mean=lm, lab_std=ls, n_pos_unlab=n_pos - out["n_pos_lab"],
                   unlab_mean=um, unlab_std=us)
        hists["lab"] = lh
        hists["unlab"] = uh
    return out, hists


# ---- main ---------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(CT_CACHE, exist_ok=True)
    t0 = time.time()
    ct_bytes0 = dir_bytes(CT_CACHE)
    m7_bytes0 = dir_bytes(M7_CACHE)

    pool = load_pool()
    voids = [r for r in pool if float(r["coverage_near4"]) < VOID_COV]
    ctrls = [r for r in pool if float(r["coverage_near4"]) > CTRL_COV]
    sel_void, strata = stratify(voids, pool, N_VOID)
    sel_ctrl, _ = stratify(ctrls, pool, N_CTRL)
    log(f"pool={len(pool)} void_cands={len(voids)} ctrl_cands={len(ctrls)} "
        f"selected {len(sel_void)} void + {len(sel_ctrl)} control")

    ct_lvl = RemoteZarrLevel(CT_URL, LEVEL, cache_dir=CT_CACHE)
    m7_lvl = RemoteZarrLevel(M7_URL, LEVEL, cache_dir=M7_CACHE)
    log(f"CT L{LEVEL} shape={ct_lvl.shape} chunks={ct_lvl.chunks} "
        f"dtype={ct_lvl.dtype} compressor={ct_lvl.codec}")
    if ct_lvl.shape != m7_lvl.shape:
        log(f"ABORT: CT L{LEVEL} shape {ct_lvl.shape} != m7 L{LEVEL} "
            f"shape {m7_lvl.shape}")
        sys.exit(2)

    # ---- sanity check: one cube from a known-labeled control window ------
    sc = sel_ctrl[0]
    lo = tuple(cube_lo_axis(int(sc[f"win_lo_{a}"])) for a in "zyx")
    hi = tuple(c + CUBE for c in lo)
    arr = ct_lvl.read_crop(lo, hi)
    u = np.unique(arr)
    pcts = np.percentile(arr, [1, 10, 25, 50, 75, 90, 99]).tolist()
    sanity = {
        "window": sc["path"], "cube_lo": list(lo),
        "n_unique_values": int(u.size), "std": float(arr.std()),
        "mean": float(arr.mean()),
        "percentiles_1_10_25_50_75_90_99": pcts,
    }
    log(f"sanity cube {sc['path']} lo={lo}: mean={arr.mean():.1f} "
        f"std={arr.std():.1f} unique={u.size} pcts={pcts}")
    if u.size < 30 or arr.std() < 5.0:
        log("ABORT: sanity check failed (near-constant or degenerate cube)")
        with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
            json.dump({"status": "sanity_failed", "ct_url": CT_URL,
                       "sanity": sanity}, f, indent=1)
        sys.exit(2)
    sanity["verdict"] = "ok"

    # ---- sweep -----------------------------------------------------------
    jobs = [(r, "void") for r in sel_void] + [(r, "control") for r in sel_ctrl]
    results = []
    hist_store = {}
    n_fail = 0
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(process_window, r, ct_lvl, m7_lvl, g): (r, g)
                for r, g in jobs}
        for i, fut in enumerate(as_completed(futs)):
            r, g = futs[fut]
            try:
                out, hists = fut.result()
            except Exception as e:
                log(f"  FAIL {g} {r['path']}: {e}")
                n_fail += 1
                continue
            results.append(out)
            key = out["path"].replace("/", "__").replace(".npz", "")
            for name, h in hists.items():
                hist_store[f"{key}__{name}"] = h
            if (i + 1) % 15 == 0 or i + 1 == len(jobs):
                log(f"  {i + 1}/{len(jobs)} windows, {n_fail} failed, "
                    f"{(time.time() - t0) / 60:.1f} min")

    results.sort(key=lambda o: (o["group"], o["path"]))
    fields = ["path", "group", "zband", "rthird", "coverage_near4",
              "n_pred_window", "cube_lo_z", "cube_lo_y", "cube_lo_x",
              "n_pos", "n_pos_lab", "n_pos_unlab", "pos_mean", "pos_std",
              "lab_mean", "lab_std", "unlab_mean", "unlab_std",
              "neg_mean", "neg_std", "frac_zero_all", "tex_rad", "tex_z",
              "rad_axis"]
    with open(os.path.join(OUT_DIR, "windows.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for o in results:
            w.writerow(o)
    np.savez_compressed(os.path.join(OUT_DIR, "window_hists.npz"),
                        **hist_store)

    ct_dl = dir_bytes(CT_CACHE) - ct_bytes0
    m7_dl = dir_bytes(M7_CACHE) - m7_bytes0
    summary = {
        "status": "ok" if n_fail == 0 else f"ok_with_{n_fail}_failures",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "ct_url": CT_URL,
        "m7_url": M7_URL,
        "level": LEVEL,
        "sanity_check": sanity,
        "selection_rule": {
            "pool": f"coverage_tiles.csv status==ok, clamped==0, "
                    f"n_pred>={MIN_PRED}",
            "void": f"coverage_near4 < {VOID_COV}, target {N_VOID}",
            "control": f"coverage_near4 > {CTRL_COV}, target {N_CTRL}",
            "stratification": "sort by (z-band decile, radial third, path), "
                              "take every k-th (k = n_cands // n_target)",
            "z_range": [strata[0], strata[1]],
            "radial_tercile_edges": [strata[2], strata[3]],
            "cube_rule": "64-cube inside the 128-window placed in a single "
                         "128 CT chunk: per axis lo = w+32 if aligned else "
                         "(b-64 if b-w>=64 else b), b = next 128 multiple",
            "n_void_candidates": len(voids),
            "n_ctrl_candidates": len(ctrls),
            "selected_void": [r["path"] for r in sel_void],
            "selected_control": [r["path"] for r in sel_ctrl],
        },
        "n_windows_done": len(results),
        "n_failed": n_fail,
        "download_bytes": {"ct": ct_dl, "m7": m7_dl},
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    log(f"done: {len(results)} windows, {n_fail} failed, "
        f"CT dl {ct_dl / 1e6:.1f} MB, m7 dl {m7_dl / 1e6:.1f} MB, "
        f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
