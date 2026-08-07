"""Diego's fifth-door check: per-volume intensity statistics, located vs
non-located, on the 892 image+label pairs. Prespecified endpoints, in order:
background median (class 0), sheet median (class 1), contrast (sheet minus
background), and overall IQR. Mann-Whitney two-sided per endpoint. Groups
come from the same groups892.json the flag run used."""
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import tifffile
from scipy.stats import mannwhitneyu

BASE = "/mnt/vesuvius/kaggle892"
GROUPS = json.load(open(f"{BASE}/groups892.json"))
STRIDE = 2  # every 2nd voxel per axis = 4.1M voxels per volume


def one(sample):
    im = tifffile.imread(f"{BASE}/images/{sample}.tif")[::STRIDE, ::STRIDE, ::STRIDE]
    lab = tifffile.imread(f"{BASE}/labels/{sample}.tif")[::STRIDE, ::STRIDE, ::STRIDE]
    bg = im[lab == 0]
    sh = im[lab == 1]
    r = {"sample": sample,
         "bg_median": float(np.median(bg)) if bg.size else None,
         "sheet_median": float(np.median(sh)) if sh.size else None,
         "iqr_all": float(np.subtract(*np.percentile(im, [75, 25]))),
         "n_sheet": int(sh.size), "n_bg": int(bg.size)}
    r["contrast"] = (r["sheet_median"] - r["bg_median"]
                     if r["sheet_median"] is not None and r["bg_median"] is not None
                     else None)
    return r

samples = sorted(GROUPS.keys())
with ProcessPoolExecutor(max_workers=8) as ex:
    rows = list(ex.map(one, samples, chunksize=8))
for r in rows:
    r["group"] = GROUPS[r["sample"]]
    r["located"] = GROUPS[r["sample"]] in ("located", "intersecting", "iou1")

out = {"stride": STRIDE, "n": len(rows),
       "endpoints_prespecified": ["bg_median", "sheet_median", "contrast", "iqr_all"],
       "groups": {}, "tests": {}, "rows": rows}
loc = [r for r in rows if r["located"]]
non = [r for r in rows if not r["located"]]
out["n_located"], out["n_nonlocated"] = len(loc), len(non)
for key in out["endpoints_prespecified"]:
    a = [r[key] for r in loc if r[key] is not None]
    b = [r[key] for r in non if r[key] is not None]
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    out["groups"][key] = {"located_median": float(np.median(a)),
                          "nonlocated_median": float(np.median(b)),
                          "located_iqr": [float(np.percentile(a, 25)), float(np.percentile(a, 75))],
                          "nonlocated_iqr": [float(np.percentile(b, 25)), float(np.percentile(b, 75))]}
    out["tests"][key] = {"U": float(U), "p_two_sided": float(p),
                         "n_located": len(a), "n_nonlocated": len(b)}
    print(f"{key:14s} loc {np.median(a):7.2f} non {np.median(b):7.2f} p={p:.3e}")
json.dump(out, open(f"{BASE}/intensity_split.json", "w"), indent=1)
print("wrote intensity_split.json")
