#!/usr/bin/env python3
"""apply_offseed_nudge.py - Amendment 2 of prereg_placement_ab.

For each off-seed whose voxel reads unwritten in the built direction field,
move it to the field-written voxel at minimal Chebyshev distance from the
frozen point, preferring m7 >= 128 within the minimal shell, ties broken by
md5("nudge:z:y:x"). Re-verify the 96-voxel cluster-distance rule at the
nudged point; abort if violated.

Outputs:
  /mnt/vesuvius/hazard_zarr_smoke/offseeds_placement_v2.json
  /mnt/vesuvius/vcbuild/demo_out/placement/off_sites.txt   (REWRITTEN)
  /mnt/vesuvius/vcbuild/demo_out/placement/nudge_record.json
"""
import hashlib
import json
import sys

import numpy as np
import zarr

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

SMOKE = "/mnt/vesuvius/hazard_zarr_smoke"
M7_URL = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
          "representations/predictions/surfaces/"
          "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")

g = zarr.open(SMOKE + "/m7_normals_L1.zarr", mode="r")
AX, AY, AZ = g["x"]["1"], g["y"]["1"], g["z"]["1"]
lvl = RemoteZarrLevel(M7_URL, 1, cache_dir=SMOKE + "/m7L1_cache")

seeds = json.load(open(SMOKE + "/offseeds_placement.json"))


def written_mask(cz, cy, cx, R):
    sl = (slice(max(0, cz - R), cz + R + 1),
          slice(max(0, cy - R), cy + R + 1),
          slice(max(0, cx - R), cx + R + 1))
    wx = np.asarray(AX[sl]).astype(np.float32)
    wy = np.asarray(AY[sl]).astype(np.float32)
    wz = np.asarray(AZ[sl]).astype(np.float32)
    v = np.stack([(wz - 128) / 127, (wy - 128) / 127, (wx - 128) / 127], -1)
    norm = np.sqrt((v ** 2).sum(-1))
    origin = (max(0, cz - R), max(0, cy - R), max(0, cx - R))
    return norm > 0.5, origin


def is_written(cz, cy, cx):
    m, o = written_mask(cz, cy, cx, 0)
    return bool(m[cz - o[0], cy - o[1], cx - o[2]])


def m7_at(p):
    return int(lvl.read_crop(p, (p[0] + 1, p[1] + 1, p[2] + 1))[0, 0, 0])


record = []
out = []
for si, s in enumerate(seeds):
    cz, cy, cx = s["gz"], s["gy"], s["gx"]
    if is_written(cz, cy, cx):
        out.append(dict(s))
        record.append({"site": si, "nudged": False})
        continue
    R = 48
    mask, origin = written_mask(cz, cy, cx, R)
    assert mask.any(), f"o{si}: no written voxel within {R}"
    c0 = np.array([cz - origin[0], cy - origin[1], cx - origin[2]])
    idx = np.argwhere(mask)
    cheb = np.abs(idx - c0).max(1)
    best = int(cheb.min())
    shell = [tuple(int(v) for v in p) for p in idx[cheb == best]]

    def absco(p):
        return (p[0] + origin[0], p[1] + origin[1], p[2] + origin[2])

    def md5key(p):
        A = absco(p)
        return hashlib.md5(f"nudge:{A[0]}:{A[1]}:{A[2]}".encode()).hexdigest()

    shell.sort(key=md5key)
    pick = None
    for p in shell:
        if m7_at(absco(p)) >= 128:
            pick = p
            break
    if pick is None:
        pick = shell[0]
    A = absco(pick)

    # re-verify cluster distance rule at the nudged point
    slab, tile = s["slab"], s["tile"]
    cj = json.load(open(f"/mnt/vesuvius/census8k/{slab}/{tile}.json"))
    npz = np.load(
        f"/mnt/vesuvius/p1218_repair_v3/blocks_repaired/{slab}/{tile}.npz")
    z0, y0, x0 = int(npz["z0"]), int(npz["y0"]), int(npz["x0"])
    npz.close()
    mind = min(max(abs(A[0] - (z0 + c["z"])), abs(A[1] - (y0 + c["y"])),
                   abs(A[2] - (x0 + c["x"]))) for c in cj["clusters"])
    assert mind >= 96, f"o{si}: nudged point violates distance rule ({mind})"

    ns = dict(s)
    ns["gz"], ns["gy"], ns["gx"] = A
    ns["z"], ns["y"], ns["x"] = A[0] - z0, A[1] - y0, A[2] - x0
    ns["seed_L0"] = [A[0] * 2, A[1] * 2, A[2] * 2]
    ns["m7"] = m7_at(A)
    ns["nudged_from"] = [cz, cy, cx]
    ns["nudge_cheb"] = best
    ns["min_chebyshev_to_cluster"] = int(mind)
    out.append(ns)
    record.append({"site": si, "nudged": True, "from": [cz, cy, cx],
                   "to": list(A), "cheb": best, "m7_at_pick": ns["m7"],
                   "min_dist_after": int(mind),
                   "shell_ties": len(shell)})
    print(f"o{si} nudged {cz},{cy},{cx} -> {A[0]},{A[1]},{A[2]} "
          f"(cheb {best}, m7 {ns['m7']}, minDist {mind})")

json.dump(out, open(SMOKE + "/offseeds_placement_v2.json", "w"), indent=1)
json.dump(record, open(
    "/mnt/vesuvius/vcbuild/demo_out/placement/nudge_record.json", "w"),
    indent=1)
with open("/mnt/vesuvius/vcbuild/demo_out/placement/off_sites.txt", "w") as f:
    for i, p in enumerate(out):
        f.write(f"{i} {p['slab']} {p['tile']} {p['gz']} {p['gy']} {p['gx']}\n")
print("nudge applied:", sum(1 for r in record if r["nudged"]), "of", len(out))
