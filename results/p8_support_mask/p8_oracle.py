"""Label oracle: did the repair remove any LABELED voxel? Expect 0."""
import glob, json, sys
import numpy as np, zarr

ROOT = "/mnt/vesuvius/p8_sprint"
tag = sys.argv[1]
rep = json.load(open("%s/%s_repair.json" % (ROOT, tag)))
lo, hi = rep["bounds_l0"]
b = zarr.open("%s/%s_before.zarr/0" % (ROOT, tag), mode="r")
a = zarr.open("%s/%s_after.zarr" % (ROOT, tag), mode="r")
if hasattr(a, "keys") and "0" in a: a = a["0"]

l1lo = [v // 2 for v in lo]; l1hi = [v // 2 for v in hi]
tot_lab = tot_lab_pos = removed_lab = removed_all = 0
tiles = 0
for p in glob.glob("/mnt/vesuvius/p1218_repair_v3/blocks_repaired/z*/tile_*.npz"):
    import re
    slab = p.split("/")[-2]; m = re.search(r"tile_y(\d+)_x(\d+)", p)
    z0 = int(slab[1:]); y0 = int(m.group(1)); x0 = int(m.group(2))
    if not (l1lo[0] < z0 + 256 and z0 < l1hi[0] and l1lo[1] < y0 + 512
            and y0 < l1hi[1] and l1lo[2] < x0 + 512 and x0 < l1hi[2]):
        continue
    with np.load(p) as d: lab = d["labels"]
    zs, ze = max(z0, l1lo[0]), min(z0 + lab.shape[0], l1hi[0])
    ys, ye = max(y0, l1lo[1]), min(y0 + lab.shape[1], l1hi[1])
    xs, xe = max(x0, l1lo[2]), min(x0 + lab.shape[2], l1hi[2])
    if zs >= ze or ys >= ye or xs >= xe: continue
    sub = lab[zs-z0:ze-z0, ys-y0:ye-y0, xs-x0:xe-x0] > 0
    if not sub.any(): continue
    tiles += 1
    # L1 voxel -> L0 block of 2^3; sample the L0 corner voxel (2*idx)
    bb = b[2*zs:2*ze:2, 2*ys:2*ye:2, 2*xs:2*xe:2]
    aa = a[2*zs:2*ze:2, 2*ys:2*ye:2, 2*xs:2*xe:2]
    n = min(sub.shape[0], bb.shape[0]), min(sub.shape[1], bb.shape[1]), min(sub.shape[2], bb.shape[2])
    sub = sub[:n[0], :n[1], :n[2]]; bb = bb[:n[0], :n[1], :n[2]]; aa = aa[:n[0], :n[1], :n[2]]
    posb = bb >= 128; posa = aa >= 128
    tot_lab += int(sub.sum()); tot_lab_pos += int((sub & posb).sum())
    removed_lab += int((sub & posb & ~posa).sum())
    removed_all += int((posb & ~posa).sum())
out = {"tag": tag, "tiles_overlapping": tiles, "labeled_voxels_L1": tot_lab,
       "labeled_and_m7pos_before": tot_lab_pos,
       "labeled_positives_REMOVED": removed_lab,
       "all_positives_removed_in_sampled_grid": removed_all}
print(json.dumps(out, indent=1))
json.dump(out, open("%s/%s_oracle.json" % (ROOT, tag), "w"), indent=1)
