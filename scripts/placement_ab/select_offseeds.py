#!/usr/bin/env python3
"""select_offseeds.py - deterministic OFF-cluster seeds for prereg_placement_ab.

For each of the 8 original demo sites (demo_sites.json), pick one seed in the
SAME tile, far from every census cluster centroid in that tile, on the m7
predicted sheet. Rule frozen in prereg_placement_ab/PREREGISTRATION.md:

  - tile-local interior 48 <= z <= 208, 48 <= y <= 464, 48 <= x <= 464
  - stride-16 lattice over that interior
  - Chebyshev distance >= 96 (L1 vox) from EVERY census cluster centroid
    listed for the tile in census8k (relax 96 -> 64 -> 48 if empty, recorded)
  - order survivors by md5("slab/tile:off:z:y:x"), take the first with
    m7 prediction L1 value >= 128 (same reader as select_sites16.py)

Runs on the box. Outputs:
  /mnt/vesuvius/hazard_zarr_smoke/offseeds_placement.json
  /mnt/vesuvius/vcbuild/demo_out/placement/off_sites.txt ("i slab tile gz gy gx")
"""
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

SMOKE = "/mnt/vesuvius/hazard_zarr_smoke"
M7_URL = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
          "representations/predictions/surfaces/"
          "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")

sites = json.load(open(SMOKE + "/demo_sites.json"))
lvl = RemoteZarrLevel(M7_URL, 1, cache_dir=SMOKE + "/m7L1_cache")

out = []
for si, s in enumerate(sites):
    slab, tile = s["slab"], s["tile"]
    cj = json.load(open(f"/mnt/vesuvius/census8k/{slab}/{tile}.json"))
    cents = [(c["z"], c["y"], c["x"]) for c in cj["clusters"]]
    npz = np.load(f"/mnt/vesuvius/p1218_repair_v3/blocks_repaired/{slab}/{tile}.npz")
    z0, y0, x0 = int(npz["z0"]), int(npz["y0"]), int(npz["x0"])
    npz.close()

    chosen = None
    for dist in (96, 64, 48):
        cand = []
        for z in range(48, 209, 16):
            for y in range(48, 465, 16):
                for x in range(48, 465, 16):
                    if all(max(abs(z - cz), abs(y - cy), abs(x - cx)) >= dist
                           for cz, cy, cx in cents):
                        cand.append((z, y, x))
        cand.sort(key=lambda p: hashlib.md5(
            f"{slab}/{tile}:off:{p[0]}:{p[1]}:{p[2]}".encode()).hexdigest())
        n_m7_checked = 0
        for z, y, x in cand:
            gz, gy, gx = z0 + z, y0 + y, x0 + x
            v = int(lvl.read_crop((gz, gy, gx),
                                  (gz + 1, gy + 1, gx + 1))[0, 0, 0])
            n_m7_checked += 1
            if v >= 128:
                chosen = {"slab": slab, "tile": tile,
                          "z": z, "y": y, "x": x,
                          "gz": gz, "gy": gy, "gx": gx,
                          "seed_L0": [gz * 2, gy * 2, gx * 2],
                          "m7": v, "dist_rule": dist,
                          "n_candidates": len(cand),
                          "n_m7_checked": n_m7_checked,
                          "min_chebyshev_to_cluster": int(min(
                              max(abs(z - cz), abs(y - cy), abs(x - cx))
                              for cz, cy, cx in cents))}
                break
        if chosen:
            break
    assert chosen is not None, f"no off-seed found for {slab}/{tile}"
    if chosen["dist_rule"] != 96:
        print(f"RELAXED to {chosen['dist_rule']} for {slab}/{tile}")
    out.append(chosen)
    print(f"o{si} {slab}/{tile} L1 ({chosen['gz']},{chosen['gy']},{chosen['gx']})"
          f" m7 {chosen['m7']} minCheb {chosen['min_chebyshev_to_cluster']}"
          f" rule {chosen['dist_rule']} cand {chosen['n_candidates']}")

json.dump(out, open(SMOKE + "/offseeds_placement.json", "w"), indent=1)
os.makedirs("/mnt/vesuvius/vcbuild/demo_out/placement", exist_ok=True)
with open("/mnt/vesuvius/vcbuild/demo_out/placement/off_sites.txt", "w") as f:
    for i, p in enumerate(out):
        f.write(f"{i} {p['slab']} {p['tile']} {p['gz']} {p['gy']} {p['gx']}\n")
print(f"\nwrote {len(out)} off-seeds")
