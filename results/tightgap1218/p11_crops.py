"""P11 step 3: extract the crops, per PREREGISTRATION.md (frozen 9d5d71dbaf45ab85).

Rules applied here, unchanged from the frozen document:
  4. 128^3 at level 0, centred on the site, wholly inside the volume, CT present
     in all eight octants.
  5. census row order within band, skip if overlap with an accepted cube is more
     than a quarter of its volume, target 60 per band, realised count reported.

Usage: p11_crops.py <per_band_target>
"""
import json
import os
import sys

import numpy as np
import zarr

OUT = "/mnt/vesuvius/tightgap1218"
CT = ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
      "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
BLOCKS = "/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
N = 128
TARGET = int(sys.argv[1]) if len(sys.argv) > 1 else 60

rows = json.load(open(OUT + "/sites_gaps.json"))
ct = zarr.open(CT, mode="r")
ct0 = ct["0"] if hasattr(ct, "keys") and "0" in ct else ct
Z, Y, X = ct0.shape
print("CT", ct0.shape, ct0.dtype, "| sites", len(rows), flush=True)

# slab/tile origin: names encode level-0 offsets, e.g. z4480 / tile_y1344_x1792
def origin(slab, tile):
    z = int(slab[1:])
    parts = tile.split("_")
    y = int(parts[1][1:])
    x = int(parts[2][1:])
    return z, y, x

accepted, per_band, skipped = [], {}, {"oob": 0, "octant": 0, "overlap": 0, "read": 0}
os.makedirs(OUT + "/crops", exist_ok=True)

for r in rows:
    b = r["band"]
    if per_band.get(b, 0) >= TARGET:
        continue
    oz, oy, ox = origin(r["slab"], r["tile"])
    cz, cy, cx = oz + r["z"], oy + r["y"], ox + r["x"]
    z0, y0, x0 = cz - N // 2, cy - N // 2, cx - N // 2
    if z0 < 0 or y0 < 0 or x0 < 0 or z0 + N > Z or y0 + N > Y or x0 + N > X:
        skipped["oob"] += 1
        continue
    ov = False
    for a in accepted:
        d = [abs(z0 - a[0]), abs(y0 - a[1]), abs(x0 - a[2])]
        inter = 1.0
        for k in range(3):
            inter *= max(0, N - d[k])
        if inter > 0.25 * N ** 3:
            ov = True
            break
    if ov:
        skipped["overlap"] += 1
        continue
    try:
        vol = np.asarray(ct0[z0:z0 + N, y0:y0 + N, x0:x0 + N])
    except Exception:
        skipped["read"] += 1
        continue
    h = N // 2
    octs = [vol[:h, :h, :h], vol[:h, :h, h:], vol[:h, h:, :h], vol[:h, h:, h:],
            vol[h:, :h, :h], vol[h:, :h, h:], vol[h:, h:, :h], vol[h:, h:, h:]]
    if any((o > 0).mean() < 0.01 for o in octs):
        skipped["octant"] += 1
        continue
    tag = "%s_%s_%d_%d_%d" % (r["slab"], r["tile"], r["z"], r["y"], r["x"])
    np.savez_compressed("%s/crops/%s.npz" % (OUT, tag), intensity=vol,
                        gap=np.float32(r["gap"]), band=b,
                        site=np.array([cz, cy, cx]), A=r["A"], B=r["B"])
    accepted.append((z0, y0, x0))
    per_band[b] = per_band.get(b, 0) + 1
    if len(accepted) % 10 == 0:
        print("  accepted %d | %s" % (len(accepted), per_band), flush=True)

out = dict(prereg="9d5d71dbaf45ab85", target_per_band=TARGET,
           realised=per_band, accepted=len(accepted), skipped=skipped)
json.dump(out, open(OUT + "/crops_summary.json", "w"), indent=1)
print(json.dumps(out, indent=1))
