"""P8 day1 control: can the tracer consume a LOCAL sparse copy of m7 L0?

Builds a full-shape, level-0-only local OME-Zarr whose chunks are populated
only inside one region around an existing demo seed, then the caller runs
the tracer against it and compares with the archived S3 run.
"""
import json
import sys
from pathlib import Path

import numpy as np
import zarr

M7 = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
      "representations/predictions/surfaces/"
      "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")
OUT = Path(sys.argv[1])            # local zarr path
CZ, CY, CX = (int(v) for v in sys.argv[2:5])   # seed in L0 zyx
R = int(sys.argv[5]) if len(sys.argv) > 5 else 640   # half-size in L0 vox

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

lvl0 = RemoteZarrLevel(M7, 0, cache_dir="/mnt/vesuvius/p8_sprint/m7L0_cache")
SHAPE = tuple(int(v) for v in lvl0.shape)
CH = tuple(int(v) for v in lvl0.chunks)
print("m7 L0 shape", SHAPE, "chunks", CH)

OUT.mkdir(parents=True, exist_ok=True)
root = zarr.open_group(str(OUT), mode="w")
a = root.create_dataset("0", shape=SHAPE, chunks=CH, dtype="u1",
                        fill_value=0, write_empty_chunks=False,
                        compressor=zarr.Blosc(cname="zstd", clevel=1))
# OME-NGFF multiscales with a single level, mirroring the source convention
root.attrs["multiscales"] = [{
    "version": "0.4",
    "axes": [{"name": n, "type": "space"} for n in ("z", "y", "x")],
    "datasets": [{"path": "0",
                  "coordinateTransformations": [
                      {"type": "scale", "scale": [1.0, 1.0, 1.0]}]}],
}]

z0, z1 = max(0, CZ - R), min(SHAPE[0], CZ + R)
y0, y1 = max(0, CY - R), min(SHAPE[1], CY + R)
x0, x1 = max(0, CX - R), min(SHAPE[2], CX + R)
# snap to chunk boundaries so partial chunks are never half-written
z0 -= z0 % CH[0]; y0 -= y0 % CH[1]; x0 -= x0 % CH[2]
z1 = min(SHAPE[0], ((z1 + CH[0] - 1) // CH[0]) * CH[0])
y1 = min(SHAPE[1], ((y1 + CH[1] - 1) // CH[1]) * CH[1])
x1 = min(SHAPE[2], ((x1 + CH[2] - 1) // CH[2]) * CH[2])
print("region L0 z", z0, z1, "y", y0, y1, "x", x0, x1)

nz = 0
for zz in range(z0, z1, CH[0]):
    for yy in range(y0, y1, CH[1]):
        for xx in range(x0, x1, CH[2]):
            ze, ye, xe = (min(zz + CH[0], SHAPE[0]), min(yy + CH[1], SHAPE[1]),
                          min(xx + CH[2], SHAPE[2]))
            blk = lvl0.read_crop((zz, yy, xx), (ze, ye, xe))
            if blk.any():
                a[zz:ze, yy:ye, xx:xe] = blk
                nz += 1
print("chunks written (nonzero):", nz)
json.dump({"shape": list(SHAPE), "chunks": list(CH),
           "region": [z0, z1, y0, y1, x0, x1], "chunks_written": nz,
           "seed_l0": [CZ, CY, CX]},
          open(OUT.parent / (OUT.name + ".region.json"), "w"), indent=1)
print("done")
