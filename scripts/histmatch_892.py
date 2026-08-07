"""Histogram matching for the sixth check: map each located volume's intensity
distribution onto the pooled non-located reference, so the recall bench can
rerun on intensity-normalized copies.

Deterministic, no rng. Reference = pooled uint8 histogram of the 703
non-located volumes (stride 2). Per located volume: own full-volume histogram
-> monotone LUT v -> smallest r with refCDF[r] >= ownCDF[v]. --write-dir emits
matched tifs; by default the script only writes luts + a post-match check
(the four intensity endpoints recomputed on the matched located volumes),
which is the self-test that the mapping lands the located population on the
reference regime."""
import argparse
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import tifffile

BASE = "/mnt/vesuvius/kaggle892"
GROUPS = json.load(open(f"{BASE}/groups892.json"))
LOCATED = sorted(s for s, g in GROUPS.items()
                 if g in ("located", "intersecting", "iou1"))
NONLOC = sorted(s for s, g in GROUPS.items() if g == "nonlocated")


def hist_of(sample):
    im = tifffile.imread(f"{BASE}/images/{sample}.tif")[::2, ::2, ::2]
    return np.bincount(im.ravel(), minlength=256).astype(np.int64)


def match_and_check(args):
    sample, ref_cdf, write_dir = args
    im = tifffile.imread(f"{BASE}/images/{sample}.tif")
    h = np.bincount(im.ravel(), minlength=256).astype(np.float64)
    own_cdf = np.cumsum(h) / h.sum()
    lut = np.searchsorted(ref_cdf, own_cdf, side="left").clip(0, 255).astype(np.uint8)
    matched = lut[im]
    if write_dir:
        tifffile.imwrite(f"{write_dir}/{sample}.tif", matched)
    lab = tifffile.imread(f"{BASE}/labels/{sample}.tif")
    sub = matched[::2, ::2, ::2]
    lsub = lab[::2, ::2, ::2]
    bg = sub[lsub == 0]
    sh = sub[lsub == 1]
    return {"sample": sample, "lut": lut.tolist(),
            "bg_median": float(np.median(bg)) if bg.size else None,
            "sheet_median": float(np.median(sh)) if sh.size else None,
            "iqr_all": float(np.subtract(*np.percentile(sub, [75, 25])))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-dir", default=None,
                    help="emit matched tifs here (default: luts + check only)")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        hists = list(ex.map(hist_of, NONLOC, chunksize=16))
    ref = np.sum(hists, axis=0).astype(np.float64)
    ref_cdf = np.cumsum(ref) / ref.sum()

    jobs = [(s, ref_cdf, args.write_dir) for s in LOCATED]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(match_and_check, jobs, chunksize=4))

    for r in rows:
        r["contrast"] = (r["sheet_median"] - r["bg_median"]
                         if r["sheet_median"] is not None and r["bg_median"] is not None
                         else None)
    out = {"n_located_matched": len(rows),
           "reference": "pooled histogram of the 703 non-located volumes, stride 2",
           "post_match_located_medians": {
               k: float(np.median([r[k] for r in rows if r[k] is not None]))
               for k in ("bg_median", "sheet_median", "contrast", "iqr_all")},
           "luts": {r["sample"]: r["lut"] for r in rows},
           "per_volume": [{k: r[k] for k in
                           ("sample", "bg_median", "sheet_median", "contrast", "iqr_all")}
                          for r in rows]}
    json.dump(out, open(f"{BASE}/histmatch_check.json", "w"))
    print("post-match located medians:", out["post_match_located_medians"])
    print(f"wrote histmatch_check.json ({len(rows)} volumes)")


if __name__ == "__main__":
    main()
