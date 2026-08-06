"""Double-thickness ray flag over the 892 public label volumes, split by
TAUIL's recall groups (derived from our overlap_report.json, the same file
his overlap_recall_split.json is sha-pinned to).

Same instrument as twin_measure.py / the PHerc 1218 work: sheet mask,
EDT-SDF (sigma 1.0) normals, SPAN 12 / STEP 0.5 rays, contiguous run through
p0, per-volume median thickness, flag ratio > 1.6.

Label format (verified empirically): semantic classes, NOT instance ids.
0 = labeled background (dark, image mean ~74), 1 = sheet (bright, ~102),
2 = ignore/unlabeled (dark, ~72, one giant component, 40-75% of volume).
Sheet mask = (lab == 1); the >=2-ids criterion is unavailable, so
multi_id_rate is null. sheet_fraction and ignore_fraction are recorded
(sheet_fraction joins TAUIL's label_sheet_fraction confound).

Per volume: {sample, group, n_valid, median_thick, flag_rate, multi_id_rate,
sheet_fraction, ignore_fraction}. Resumable via flag_rows.jsonl. Aggregate:
per-group medians + Mann-Whitney located-vs-nonlocated on flag_rate.
"""
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy import stats

BASE = "/mnt/vesuvius/kaggle892"
ROWS = f"{BASE}/flag_rows.jsonl"
OUT = f"{BASE}/flag_split.json"
N_SAMPLE = 2000
SPAN, STEP = 12.0, 0.5
K = int(2 * SPAN / STEP) + 1
OFFS = np.linspace(-SPAN, SPAN, K)
CTR = K // 2
WORKERS = 3

GROUPS = json.load(open(f"{BASE}/groups892.json"))
URL = ("https://huggingface.co/buckets/scrollprize/datasets/resolve/"
       "surfaces/kaggle/labels/{s}.tif")


def valid_tif(path):
    try:
        with open(path, "rb") as f:
            return f.read(4) in (b"II*\x00", b"MM\x00*")
    except OSError:
        return False


def ensure_labels(samples):
    """Purge HTML-poisoned files and re-fetch until every label validates."""
    for attempt in range(6):
        bad = [s for s in samples
               if not valid_tif(f"{BASE}/labels/{s}.tif")]
        if not bad:
            return []
        print(f"refetch attempt {attempt}: {len(bad)} bad/missing labels",
              flush=True)
        for s in bad:
            p = f"{BASE}/labels/{s}.tif"
            try:
                os.remove(p)
            except OSError:
                pass
            subprocess.run(
                ["curl", "-fsS", "-L", "--retry", "5",
                 "--retry-all-errors", "--connect-timeout", "30",
                 "--max-time", "600", "-o", p, URL.format(s=s)])
            if not valid_tif(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        time.sleep(30)
    return [s for s in samples
            if not valid_tif(f"{BASE}/labels/{s}.tif")]


def measure(sample):
    path = f"{BASE}/labels/{sample}.tif"
    try:
        lab = tifffile.imread(path)
    except Exception as e:
        return {"sample": sample, "group": GROUPS[sample], "error": str(e)}
    mask = lab == 1
    sheet_frac = float(mask.mean())
    ignore_frac = float((lab == 2).mean())
    empty = {"sample": sample, "group": GROUPS[sample], "n_valid": 0,
             "median_thick": None, "flag_rate": None, "multi_id_rate": None,
             "sheet_fraction": sheet_frac, "ignore_fraction": ignore_frac}
    if not mask.any():
        return empty
    sm = ndi.gaussian_filter(
        ndi.distance_transform_edt(~mask).astype(np.float32)
        - ndi.distance_transform_edt(mask).astype(np.float32), 1.0)
    g = np.stack(np.gradient(sm), 0).astype(np.float32)

    surf = mask & ~ndi.binary_erosion(mask)
    pts = np.argwhere(surf)
    if len(pts) == 0:
        return empty
    rng = np.random.default_rng(1218)
    pts = pts[rng.choice(len(pts), size=min(N_SAMPLE, len(pts)),
                         replace=False)]
    nv = g[:, pts[:, 0], pts[:, 1], pts[:, 2]].T
    nv /= np.linalg.norm(nv, axis=1, keepdims=True) + 1e-6
    cc = (pts[:, None, :].astype(np.float32)
          + OFFS[None, :, None] * nv[:, None, :])
    onmask = ndi.map_coordinates(
        mask.astype(np.float32), cc.reshape(-1, 3).T, order=0,
        mode="constant").reshape(len(pts), K) > 0.5

    th_l = []
    for i in range(len(pts)):
        if not onmask[i, CTR]:
            continue
        a = CTR
        while a > 0 and onmask[i, a - 1]:
            a -= 1
        b = CTR
        while b < K - 1 and onmask[i, b + 1]:
            b += 1
        th_l.append((b - a) * STEP)

    if not th_l:
        return empty
    th = np.array(th_l)
    med = float(np.median(th))
    flag = th > 1.6 * med
    return {"sample": sample, "group": GROUPS[sample],
            "n_valid": int(len(th)),
            "median_thick": med,
            "flag_rate": float(flag.mean()),
            "multi_id_rate": None,
            "sheet_fraction": sheet_frac,
            "ignore_fraction": ignore_frac}


def main():
    done = set()
    if os.path.exists(ROWS):
        with open(ROWS) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    if "error" not in r:
                        done.add(r["sample"])
                except Exception:
                    pass
    todo = sorted(s for s in GROUPS if s not in done)
    missing = ensure_labels(todo)
    if missing:
        print(f"WARNING: {len(missing)} label files unrecoverable, e.g. "
              f"{missing[:5]}", flush=True)
        todo = [s for s in todo if s not in missing]
    print(f"{len(GROUPS)} samples, {len(done)} done, {len(todo)} to run",
          flush=True)

    with open(ROWS, "a") as fout, \
            ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(measure, s): s for s in todo}
        for k, fut in enumerate(as_completed(futs)):
            row = fut.result()
            fout.write(json.dumps(row) + "\n")
            fout.flush()
            if (k + 1) % 25 == 0:
                print(f"{k + 1}/{len(todo)}", flush=True)

    rows = [json.loads(l) for l in open(ROWS)]
    rows = {r["sample"]: r for r in rows}  # dedupe, last wins
    rows = list(rows.values())
    valid = [r for r in rows if r.get("n_valid")]
    print(f"aggregating over {len(valid)} valid volumes "
          f"({len(rows) - len(valid)} empty/error)", flush=True)

    def agg(rs):
        if not rs:
            return {"n": 0}
        return {
            "n": len(rs),
            "median_thick_med": float(np.median(
                [r["median_thick"] for r in rs])),
            "flag_rate_med": float(np.median([r["flag_rate"] for r in rs])),
            "flag_rate_mean": float(np.mean([r["flag_rate"] for r in rs])),
            "sheet_fraction_med": float(np.median(
                [r["sheet_fraction"] for r in rs])),
            "ignore_fraction_med": float(np.median(
                [r["ignore_fraction"] for r in rs])),
        }

    by = {}
    for gname in ("iou1", "intersecting", "located", "nonlocated"):
        by[gname] = agg([r for r in valid if r["group"] == gname])
    located_all = [r for r in valid if r["group"] != "nonlocated"]
    inter_all = [r for r in valid if r["group"] in ("iou1", "intersecting")]
    nonloc = [r for r in valid if r["group"] == "nonlocated"]
    by["located_all_189"] = agg(located_all)
    by["intersecting_all_122"] = agg(inter_all)

    fr_loc = [r["flag_rate"] for r in located_all]
    fr_non = [r["flag_rate"] for r in nonloc]
    mwu = stats.mannwhitneyu(fr_loc, fr_non, alternative="two-sided")
    mwu_g = stats.mannwhitneyu(fr_loc, fr_non, alternative="greater")
    mt_loc = [r["median_thick"] for r in located_all]
    mt_non = [r["median_thick"] for r in nonloc]
    mwu_mt = stats.mannwhitneyu(mt_loc, mt_non, alternative="two-sided")

    out = {
        "instrument": {
            "n_sample_sites": N_SAMPLE, "span": SPAN, "step": STEP,
            "flag_ratio": 1.6, "seed": 1218,
            "label_format": "uint8 semantic classes 320^3: 0=background, "
                            "1=sheet, 2=ignore/unlabeled; mask=(lab==1); "
                            "no instance ids so multi_id_rate is null",
            "groups_source": "overlap_report.json (sha-matched to TAUIL "
                             "overlap_recall_split.json source)",
        },
        "n_volumes": len(rows), "n_valid_volumes": len(valid),
        "groups": by,
        "tests": {
            "flag_rate_located_vs_nonlocated": {
                "U": float(mwu.statistic), "p_two_sided": float(mwu.pvalue),
                "p_greater_located": float(mwu_g.pvalue),
                "n_located": len(fr_loc), "n_nonlocated": len(fr_non),
            },
            "median_thick_located_vs_nonlocated": {
                "U": float(mwu_mt.statistic),
                "p_two_sided": float(mwu_mt.pvalue),
            },
        },
        "rows": sorted(rows, key=lambda r: r["sample"]),
    }
    json.dump(out, open(OUT, "w"), indent=1)
    print("wrote", OUT, flush=True)
    print(json.dumps({k: v for k, v in out["groups"].items()}, indent=1))
    print(json.dumps(out["tests"], indent=1))


if __name__ == "__main__":
    main()
