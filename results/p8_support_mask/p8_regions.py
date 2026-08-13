"""P8 day1 step 3: pick phantom-adjacent regions with real surface nearby.

Rule (deterministic): from voidct1218/windows_class.csv keep class
masked_empty; for each, look at the m7 L1 neighbourhood and require BOTH
(a) the window itself is phantom-dense and (b) a CT-supported, labelled
neighbourhood exists within one L1 chunk of the window, so a seed placed on
real surface can grow into the phantom zone. Rank by phantom voxels and take
the top N on distinct z-slabs.
"""
import csv
import json
import sys

import numpy as np

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

M7 = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
      "representations/predictions/surfaces/"
      "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")
CT = ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
      "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4

m7 = RemoteZarrLevel(M7, 1, cache_dir="/mnt/vesuvius/hazard_zarr_smoke/m7L1_cache")
ct = RemoteZarrLevel(CT, 1, cache_dir="/mnt/vesuvius/p8_sprint/ctL1_cache")

rows = list(csv.DictReader(open("/mnt/vesuvius/voidct1218/windows_class.csv")))
cand = [r for r in rows if r.get("cls","") == "masked_empty"]
print("masked_empty windows:", len(cand))
if not cand:
    print("cols:", list(rows[0].keys()))
    sys.exit(1)


def cube_lo(r):
    return (int(float(r["cube_lo_z"])), int(float(r["cube_lo_y"])),
            int(float(r["cube_lo_x"])))


out = []
for r in cand:
    z, y, x = cube_lo(r)
    # 64-cube of the window plus a 128-vox skirt, at L1
    lo = (max(0, z - 64), max(0, y - 64), max(0, x - 64))
    hi = (z + 128, y + 128, x + 128)
    try:
        p = m7.read_crop(lo, hi) >= 128
        c = ct.read_crop(lo, hi) > 0
    except Exception as e:
        continue
    phantom = int((p & ~c).sum())
    supported = int((p & c).sum())
    if phantom < 20000 or supported < 5000:
        continue
    out.append({"win": r.get("path","?"),
                "cube_lo_l1": [z, y, x], "phantom": phantom,
                "supported": supported,
                "slab": z // 1000})

out.sort(key=lambda d: -d["phantom"])
picked, slabs = [], set()
for d in out:
    if d["slab"] in slabs:
        continue
    picked.append(d)
    slabs.add(d["slab"])
    if len(picked) >= N:
        break

print("candidates with both phantom and supported mass:", len(out))
for d in picked:
    print(" PICK", d["win"], d["cube_lo_l1"], "phantom", d["phantom"],
          "supported", d["supported"])
json.dump({"rule": "masked_empty, phantom>=20000 and supported>=5000 in "
                   "window+128 skirt at L1, rank by phantom, distinct z-slab",
           "picked": picked, "n_candidates": len(out)},
          open("/mnt/vesuvius/p8_sprint/regions.json", "w"), indent=1)
