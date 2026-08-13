# PHerc1203 thickness zone map

`zones.npz` carries five grids of shape 21 by 36 by 36 (cells of 96 cubed
voxels, 1.8 mm) on the exact grid of the physical-truth labels, origin_l1
[3936, 0, 0], label shape [2016, 3456, 3456]. See `provenance.json` for
input and output hashes, probe constants and the deterministic sampling
rule. The map was regenerated on 2026-08-12 after a read-path bug in the
original release; the retraction and redelivery are in villa issue 191 and
the fixed reader is `scripts/thickness_map_1203.py` in this repo.

One grid-padding note, from villa issue 1407. The label array is 3456 in y
and x but the underlying lo volume ends at 3422; the final 34 voxels in
both axes are chunk padding, silently zero rather than marked invalid, and
labeled content in this release actually ends before y 2646 and x 2896. In
this map that means the last cell row and column (index 35 in y and x)
contain no material by construction. Measured on the shipped file, all
1,491 cells with y or x index 35 carry zero runs and zero valid share, so
per-cell readers should treat index 35 as padding, not as a material
deficit. The coordinate transform for the label grid is in issue 1407.
