# What the seam agreement field in export_summary.json actually covers

`export_summary.json` carries `"overlap_id_agreement": "4249/5126 = 0.8289"`.
That key name is broader than the measurement behind it, so this note states
the scope and points at the number that supersedes it.

## Scope of the 82.9 percent

Reading `export_points.py` back, the block that computes it does three things
the key name does not say. It takes the first fifteen y-adjacent block pairs it
encounters and stops there. It compares only the y overlap, so no z seam and no
x seam is ever looked at. And it subsamples one voxel in sixteen on each axis.
The 5,126 comparisons are therefore drawn from one corner of the tree along one
axis, and they are not a sample of the seam population.

## The figure that supersedes it

Diego-dcv ran the full-resolution pass over the same tree and measured id
agreement on **all 106,629,002 overlap voxels** lying on the crossing table's
323 planes. That comes out at **88.8 percent**, so roughly one seam voxel in
nine carries a merge fault rather than one in six.

The gap is not sampling error. A 5,126-voxel sample sits inside about plus or
minus one point at 95 percent confidence, and 88.8 is far outside that, which
is what a narrow non-representative measurement looks like when a wide one
replaces it.

## What was done about it

The Kaggle card for `jhjeong0815/pherc1218-label-points` now carries the 88.8
percent figure with this scope note beside it. The correction was posted on
PR #1 on 2026-08-27. `export_summary.json` is left as the script emitted it,
since it is a run artifact rather than a claim, and this note is the place the
scope is stated.
