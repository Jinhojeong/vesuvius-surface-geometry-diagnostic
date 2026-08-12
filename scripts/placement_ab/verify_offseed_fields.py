#!/usr/bin/env python3
"""verify_offseed_fields.py - Amendment 1 post-build control check: every
off-seed voxel must read back a WRITTEN unit vector from the direction
field (same check the original seeds pass). Writes
demo_out/placement/field_verify.json; run_placement.sh gates on
"all_written": true. Read-only w.r.t. the field."""
import json
import math
import os
import time

import numpy as np
import zarr

SMOKE = "/mnt/vesuvius/hazard_zarr_smoke"
FIELD = SMOKE + "/m7_normals_L1.zarr"
PLACE = "/mnt/vesuvius/vcbuild/demo_out/placement"

import os
_V2 = SMOKE + "/offseeds_placement_v2.json"
offs = json.load(open(_V2 if os.path.exists(_V2) else
                      SMOKE + "/offseeds_placement.json"))
arrs = {n: zarr.open_array(os.path.join(FIELD, n, "1"), mode="r") for n in "xyz"}

rows, all_written = [], True
for i, o in enumerate(offs):
    z, y, x = o["gz"], o["gy"], o["gx"]
    v = [int(arrs[n][z, y, x]) for n in "xyz"]
    written = v != [128, 128, 128]
    norm = math.sqrt(sum(((c - 128) / 127.0) ** 2 for c in v))
    unit = written and 0.9 < norm < 1.1
    # neighborhood written fraction (+-16), same probe as the extent check
    h = 16
    sl = (slice(z - h, z + h + 1), slice(y - h, y + h + 1),
          slice(x - h, x + h + 1))
    vx = arrs["x"][sl].astype(np.int16)
    vy = arrs["y"][sl].astype(np.int16)
    vz = arrs["z"][sl].astype(np.int16)
    w = (vx != 128) | (vy != 128) | (vz != 128)
    rows.append({"site_idx": i, "off_L1_zyx": [z, y, x],
                 "field_u8_xyz": v, "written": bool(written),
                 "norm_decoded": round(norm, 4), "unit_vector": bool(unit),
                 "written_frac_pm16": round(float(w.mean()), 4)})
    all_written = all_written and unit
    print("o%d %s u8=%s |v|=%.3f written=%s frac16=%.3f"
          % (i, (z, y, x), tuple(v), norm, written, float(w.mean())))

rec = {"generated": time.strftime("%F %T"),
       "amendment": "prereg_placement_ab @ 083f424 Amendment 1",
       "check": "off-seed voxel reads back a written unit vector",
       "all_written": bool(all_written), "seeds": rows}
json.dump(rec, open(os.path.join(PLACE, "field_verify.json"), "w"), indent=1)
print("all_written:", all_written)
if not all_written:
    raise SystemExit("BAND MISS: at least one off-seed voxel unwritten after "
                     "its region build; STOP and report per instructions")
