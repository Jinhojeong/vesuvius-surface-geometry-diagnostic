"""P8 day1 step 3-4: build before/support zarrs for ONE region, run the
PR #1156 masking, and run the label oracle check.

argv: <region_index> <halfsize_L0>
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import zarr

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

M7 = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
      "representations/predictions/surfaces/"
      "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")
CT = ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
      "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
ROOT = Path("/mnt/vesuvius/p8_sprint")
IDX = int(sys.argv[1])
R = int(sys.argv[2]) if len(sys.argv) > 2 else 640

reg = json.load(open(ROOT / "regions.json"))["picked"][IDX]
z1, y1, x1 = reg["cube_lo_l1"]                     # L1 window corner
# centre of the region at L0 = (corner + 32) * 2
CZ, CY, CX = (z1 + 32) * 2, (y1 + 32) * 2, (x1 + 32) * 2
tag = "r%d" % IDX
print("region", tag, reg["win"], "centre L0", (CZ, CY, CX))

m7 = RemoteZarrLevel(M7, 0, cache_dir=str(ROOT / "m7L0_cache"))
ct = RemoteZarrLevel(CT, 0, cache_dir=str(ROOT / "ctL0_cache"))
SHAPE = tuple(int(v) for v in m7.shape)
CH = tuple(int(v) for v in m7.chunks)

lo = [max(0, c - R) for c in (CZ, CY, CX)]
hi = [min(SHAPE[i], (CZ, CY, CX)[i] + R) for i in range(3)]
for i in range(3):
    lo[i] -= lo[i] % CH[i]
    hi[i] = min(SHAPE[i], ((hi[i] + CH[i] - 1) // CH[i]) * CH[i])
print("region L0 bounds", lo, hi)


def make(path, levels_attr=True):
    p = ROOT / path
    root = zarr.open_group(str(p), mode="w")
    a = root.create_dataset("0", shape=SHAPE, chunks=CH, dtype="u1",
                            fill_value=0, write_empty_chunks=False,
                            compressor=zarr.Blosc(cname="zstd", clevel=1))
    if levels_attr:
        root.attrs["multiscales"] = [{
            "version": "0.4",
            "axes": [{"name": n, "type": "space"} for n in ("z", "y", "x")],
            "datasets": [{"path": "0", "coordinateTransformations":
                          [{"type": "scale", "scale": [1.0, 1.0, 1.0]}]}]}]
    return p, a


bpath, barr = make("%s_before.zarr" % tag)
spath, sarr = make("%s_support.zarr" % tag)

n_pred_chunks = n_ct_chunks = 0
for zz in range(lo[0], hi[0], CH[0]):
    for yy in range(lo[1], hi[1], CH[1]):
        for xx in range(lo[2], hi[2], CH[2]):
            e = (min(zz + CH[0], SHAPE[0]), min(yy + CH[1], SHAPE[1]),
                 min(xx + CH[2], SHAPE[2]))
            p = m7.read_crop((zz, yy, xx), e)
            if p.any():
                barr[zz:e[0], yy:e[1], xx:e[2]] = p
                n_pred_chunks += 1
            c = ct.read_crop((zz, yy, xx), e)
            if c.any():
                sarr[zz:e[0], yy:e[1], xx:e[2]] = c
                n_ct_chunks += 1
print("chunks written: pred", n_pred_chunks, "ct", n_ct_chunks, flush=True)

# ---- run the PR's CLI on the pair
apath = ROOT / ("%s_after.zarr" % tag)
cmd = ["/home/jinhojeong/Vesuvius/.venv/bin/python3", "-c",
       "import sys; sys.path.insert(0, '/home/jinhojeong/Vesuvius/"
       "villa-pr1156/vesuvius/src');"
       "from vesuvius.models.run.mask_predictions import "
       "mask_finalized_predictions as f;"
       "import json; print(json.dumps(f(%r, %r, %r, verbose=False)))"
       % (str(bpath / "0"), str(apath), str(spath / "0"))]
print("running PR mask CLI ...", flush=True)
res = subprocess.run(cmd, capture_output=True, text=True)
print("rc", res.returncode)
print(res.stdout[-2000:])
print(res.stderr[-2000:])
stats = json.loads(res.stdout.strip().splitlines()[-1]) if res.returncode == 0 else None
json.dump({"region": reg, "tag": tag, "bounds_l0": [lo, hi],
           "chunks_pred": n_pred_chunks, "chunks_ct": n_ct_chunks,
           "mask_stats": stats},
          open(ROOT / ("%s_repair.json" % tag), "w"), indent=1)
print("saved", ROOT / ("%s_repair.json" % tag))
