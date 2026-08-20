# Amendment 2 to the tight-contact validation set preregistration

Amends PREREGISTRATION.md, sha256 prefix `9d5d71dbaf45ab85`. AMENDMENT_1
(`a07c2f86f39b3bd8`) stands unchanged. Written 2026-08-20, after versions 1 and
2 were published and before any corrected crop was extracted.

## What was wrong

Rule 4 of the frozen document says the crop is "128^3 at level 0, centred on the
site". The site coordinates come from the repaired label block names, and those
blocks are on the CT's level-1 grid, extent 11,624 by 3,797 by 3,797, which
matches the level-1 array exactly. `p11_crops.py` and `p11_control2.py` took
those level-1 coordinates and indexed the level-0 array, shape 23,247 by 7,593
by 7,593. Every published crop therefore carries CT from about half the true
offset. It is a real region of PHerc1218 and not the region its own labels
describe.

The clause "at level 0" was written on a wrong premise about which grid the
sites live on. The faithful reading of rule 4 is that the crop is 128 cubed on
the grid the sites and labels are given in, and that is what this amendment
fixes it to.

Evidence, reproducible with `p11_gridcheck.py`. Correlation between the crop's
CT and its own label mask is 0.11 to 0.19 for the shipped arrays and 0.46 to
0.65 for level 1 read at the same indices. Reading level 0 at doubled indices
over a 256 cube and mean-pooling back to 128 reproduces the level-1 values
exactly, which fixes level 1 as the correct region rather than merely the
better-correlated one. Zero fraction is 0.0000 to 0.0625 in the shipped arrays
against 0.35 to 0.42 in the correct CT.

## What changes

One line. The CT array read for both arms becomes level 1, and the volume bounds
used by the out-of-bounds test become the level-1 shape. The crop is still 128
cubed in label voxels, so nothing about the label content of a crop changes.

## What does not change, stated so that nothing here looks tuned

The per-band target of 60. The eight-octant rule and its 1 percent threshold.
The overlap exclusion at a quarter of cube volume. Census row order within band.
The band edges. The control arm rule from AMENDMENT_1 and its seed. The gap
measurement, which never touched the CT and is not re-run.

No parameter is being chosen. The only edit is forced by the defect.

## Disclosure, because I looked before freezing this

A feasibility probe was run before this document was written, 150 sites sampled
evenly through each band's census order, applying the unchanged rule against the
correct CT. It is disclosed here rather than left out because the reader should
know what I had seen. On the two bands that had finished when this was frozen,
the octant rule rejected 1 of 150 and 0 of 150, against 409 of 429 on the
displaced CT. So the corrected run is expected to accept far more sites,
particularly in the tightest band, which shipped 14 crops of a targeted 60.

Seeing that changed no rule, because there is nothing here to change. It does
change what I expect, and the honest form of that is to write the expectation
down now rather than to present it as a finding later.

## What this predicts, frozen here

The 0 to 2 band should reach or approach its 60 target rather than stopping at
14, and the realised counts should be closer to 60 across all five bands. The
`ct_empty_frac` distribution should move up substantially, since the correct CT
at these sites is 35 to 42 percent zero against nearly zero in the shipped
crops. Membership will differ from version 2 by more than a handful of files,
so this is a new version rather than a re-upload.

## What would make this amendment wrong

If the corrected extraction accepts fewer sites than the displaced one, the
premise that the octant rule was mostly rejecting on the wrong region is wrong
and the results file says so. If crops from the corrected run still show CT to
label correlation near 0.1, the grid diagnosis is wrong and the whole repair is
withdrawn rather than shipped. If the realised band counts come out identical to
version 2, the change did not take effect and the run is void.

## What stays published

Versions 1 and 2 stay on Kaggle rather than being deleted, marked as defective
in the card, so that anyone who pulled them can tell what they have. The gap
population, all 49,295 sites with their band histogram, is unaffected and is not
re-derived, since it was measured along label normals in the label volume.
