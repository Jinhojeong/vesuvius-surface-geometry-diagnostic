#!/usr/bin/env python3
"""select_sites16.py - pick 8 NEW seed sites to extend the powered A/B block.

Replicates demo_sites.py verbatim (candidate pool, thresholds, sort order),
then applies the powered16 constraints:
  - exclude every (slab, tile) already used by the original 8 sites
  - at most ONE new site per tile (the original s2/s4 same-tile pair was a
    sensitivity problem)
  - slab cap 2 TOTAL across old+new (demo_sites.py capped 2 per slab within
    its own selection; here the seen-counter is initialised with the
    original sites' slabs so the spread rule holds over all 16)
  - feasibility check: m7 prediction L1 value at the seed voxel >= 128
    (all 8 original sites had 255; a seed off the predicted sheet cannot
    grow). Skips are recorded.
Selection is fully deterministic: no RNG, sort key (-n_sites, -ratio),
then the same mid-scroll filter (0.25*zmax < gz < 0.8*zmax, fallback to
picks[:8] when fewer than 6 mid picks) as demo_sites.py.

Outputs:
  /mnt/vesuvius/hazard_zarr_smoke/demo_sites16.json   (orig 8 verbatim + new 8)
  /mnt/vesuvius/hazard_zarr_smoke/demo_sites_new8.json (new 8 only)
  /mnt/vesuvius/vcbuild/demo_out/powered16/new_sites.txt ("i slab tile gz gy gx",
      i = 8..15, same format as ab2/sites.txt)
"""
import glob
import json
import sys

import numpy as np

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

SMOKE = "/mnt/vesuvius/hazard_zarr_smoke"
M7_URL = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
          "representations/predictions/surfaces/"
          "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")

orig = json.load(open(SMOKE + "/demo_sites.json"))
used_tiles = {(s["slab"], s["tile"]) for s in orig}
slab_seen = {}
for s in orig:
    slab_seen[s["slab"]] = slab_seen.get(s["slab"], 0) + 1
print("original tiles excluded:", sorted(used_tiles))
print("original slab counts:", slab_seen)

# ---- candidate pool: identical to demo_sites.py ------------------------
CAND = []
n_no_npz = 0
for f in glob.glob("/mnt/vesuvius/census8k/*/*.json"):
    r = json.load(open(f))
    if r.get("thin"):
        continue
    slab, tile = r["tile"].split("/")
    npz = f"/mnt/vesuvius/p1218_repair_v3/blocks_repaired/{slab}/{tile}.npz"
    z0 = y0 = x0 = None
    for c in r["clusters"]:
        if c["ratio"] < 3.0 or c["n_sites"] < 2:
            continue
        if z0 is None:
            try:
                d = np.load(npz)
            except FileNotFoundError:
                n_no_npz += 1
                break
            z0, y0, x0 = int(d["z0"]), int(d["y0"]), int(d["x0"])
            d.close()
        # interior of the tile only (margin 48) so the column is solvable
        if not (48 <= c["z"] <= 208 and 48 <= c["y"] <= 464 and 48 <= c["x"] <= 464):
            continue
        CAND.append({
            "slab": slab, "tile": tile,
            "gz": z0 + c["z"], "gy": y0 + c["y"], "gx": x0 + c["x"],
            "ratio": c["ratio"], "th": c["th"], "n_sites": c["n_sites"],
        })

CAND.sort(key=lambda c: (-c["n_sites"], -c["ratio"]))
zmax = max(c["gz"] for c in CAND) if CAND else 0
print(f"candidates {len(CAND):,} (tiles without npz skipped: {n_no_npz}), "
      f"zmax {zmax}")

# ---- constrained pick --------------------------------------------------
# demo_sites.py prefers the mid-scroll band (0.25*zmax < gz < 0.8*zmax) and
# needs 8 sites; here that preference is applied directly: pass 1 picks the
# top 8 (same sort order) among MID-BAND candidates under the tile/slab
# constraints, pass 2 tops up from non-mid candidates only if pass 1 finds
# fewer than 8. Deterministic, no RNG.
lvl = RemoteZarrLevel(M7_URL, 1, cache_dir=SMOKE + "/m7L1_cache")
picks, tile_seen, skipped_m7 = [], set(), []


def eligible(c):
    key = (c["slab"], c["tile"])
    if key in used_tiles or key in tile_seen:
        return False
    if slab_seen.get(c["slab"], 0) >= 2:
        return False
    v = int(lvl.read_crop((c["gz"], c["gy"], c["gx"]),
                          (c["gz"] + 1, c["gy"] + 1, c["gx"] + 1))[0, 0, 0])
    if v < 128:
        skipped_m7.append({"slab": c["slab"], "tile": c["tile"],
                           "gzyx": [c["gz"], c["gy"], c["gx"]], "m7": v})
        return False
    c["m7_center"] = v
    return True


def take(c):
    picks.append(c)
    tile_seen.add((c["slab"], c["tile"]))
    slab_seen[c["slab"]] = slab_seen.get(c["slab"], 0) + 1


in_mid = lambda c: 0.25 * zmax < c["gz"] < 0.8 * zmax
n_mid_pass = 0
for c in CAND:
    if len(picks) >= 8:
        break
    if in_mid(c) and eligible(c):
        take(c)
        n_mid_pass += 1
if len(picks) < 8:
    for c in CAND:
        if len(picks) >= 8:
            break
        if not in_mid(c) and eligible(c):
            take(c)

out = picks
print(f"picked {len(picks)} ({n_mid_pass} mid-scroll pass-1, "
      f"{len(picks) - n_mid_pass} top-up)")
print(f"m7 seed-value skips: {len(skipped_m7)}")
for s in skipped_m7:
    print("  SKIP", s)

for p in out:
    p["seed_L0"] = [p["gz"] * 2, p["gy"] * 2, p["gx"] * 2]

sites16 = orig + out
json.dump(sites16, open(SMOKE + "/demo_sites16.json", "w"), indent=1)
json.dump(out, open(SMOKE + "/demo_sites_new8.json", "w"), indent=1)
import os
os.makedirs("/mnt/vesuvius/vcbuild/demo_out/powered16", exist_ok=True)
with open("/mnt/vesuvius/vcbuild/demo_out/powered16/new_sites.txt", "w") as f:
    for i, p in enumerate(out):
        f.write(f"{i + 8} {p['slab']} {p['tile']} {p['gz']} {p['gy']} {p['gx']}\n")

print(f"\nnew sites (indices 8..{7 + len(out)}):")
for i, p in enumerate(out):
    print(f"  s{i + 8} {p['slab']}/{p['tile']} L1 ({p['gz']},{p['gy']},{p['gx']})"
          f" ratio {p['ratio']} th {p['th']} n_sites {p['n_sites']}"
          f" m7 {p['m7_center']}")
assert len(out) == 8, f"expected 8 new sites, got {len(out)}"
tiles_all = [(s["slab"], s["tile"]) for s in sites16]
assert len(set(tiles_all[8:])) == 8, "new sites share a tile"
assert not (set(tiles_all[:8]) & set(tiles_all[8:])), "new site reuses an original tile"
print("\nconstraint check PASS: 8 new sites, all tiles distinct from each "
      "other and from the original 8")
