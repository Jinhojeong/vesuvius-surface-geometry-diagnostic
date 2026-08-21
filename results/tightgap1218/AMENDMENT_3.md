# Amendment 3 to the tight-contact validation set preregistration

Amends PREREGISTRATION.md, sha256 prefix `9d5d71dbaf45ab85`. AMENDMENT_1
(`a07c2f86f39b3bd8`) and AMENDMENT_2 (`f9d57773c715cb65`) stand unchanged.
Written 2026-08-21, after version 3 was published and before any version 4 crop
was extracted.

## What is wrong

The repaired label blocks are 256 by 512 by 512 arrays laid on strides of 224 in
z and 448 in y and x, so neighbouring blocks physically overlap by 32 voxels in z
and 64 in y and x. The stitching step walks blocks in sorted order and writes

    out[window][sub > 0] = sub[sub > 0] + base

with a different `base` per block. So wherever two blocks cover the same voxel,
the later block's write wins and stamps that voxel with the later block's base.
The two ids the split separated, A_id and B_id, are defined through the base of
the block the site falls in. Any of their voxels lying in an overlap region can
therefore be relabelled to a different base and disappear from the pair.

Measured over all 300 version-3 contact crops, with the site block's own
contribution as the reference:

| | |
| --- | --- |
| crops containing at least one doubly covered voxel | 262 of 300, 87.3 percent |
| median share of a crop that is doubly covered | 37.5 percent |
| contributing blocks per crop | median 4, maximum 8 |
| share of the site block's pair voxels lost to the overwrite | median 0.0 percent, mean 10.4 percent |
| crops losing more than 10 percent of the pair | 76 |
| crops losing more than half the pair | 26 |
| crops where the site block holds both ids but the stitched crop lost one entirely | 3 |

So the defect is real and it degrades the shipped instance labels, but it is not
the explanation for the 52 crops that carry neither id. Only 3 of those are
caused by it. The rest are crops whose pair does not lie inside the 128 cube
at all. An earlier reading that attributed all of them to this mechanism was
checked against the data and does not hold, and this document records the
measured number rather than the larger one.

This defect is present in versions 1, 2 and 3 alike, since all three used the
same stitching step. It is independent of the pyramid-level defect AMENDMENT_2
fixed.

## What changes

The stitching step gains a write order and a priority. The block containing the
site is written first. Every other contributing block then writes only where
nothing has been written yet, that is with the mask `(sub > 0) & (out == 0)`.
This makes the site block's ids authoritative inside the crop, so the pair can
no longer be overwritten, and it makes the result independent of block sort
order in the overlap regions.

## What does not change

Every selection rule. The per-band target of 60, the eight-octant rule and its 1
percent threshold, the overlap exclusion at a quarter of cube volume, census row
order within band, the band edges, the control-arm rule and seed from
AMENDMENT_1, and the level-1 CT read from AMENDMENT_2. The crops selected are
the same crops. Only the contents of the `instance` and `surface` arrays change,
and only inside overlap regions.

The gap values do not change and are not recomputed. They were measured inside a
single label block before any stitching, so this defect never touched them.

## What this predicts, frozen here

The both-instances count should rise from 245 of 300 by about 3, since 3 crops
lose an id entirely to the overwrite today. It should not rise by 55, and a much
larger jump would mean I have misunderstood the mechanism. The pair-voxel loss
measured the same way should fall to zero for every crop. The `surface` array,
being `instance > 0`, changes only where an overlap region had a labelled voxel
in one block and not the other, so it should barely move.

## What would make this amendment wrong

If the corrected stitch changes any gap value, the change has reached something
it should not have and the run is void. If the both-instances count moves by
much more than 3, the mechanism described above is not the one operating. If any
crop's pair-voxel loss is still non-zero after the fix, the fix does not do what
this document says.

## Disclosure

The defect was found by an adversarial audit of the pre-publication state, run
because two earlier errors had already been found the same day. The audit's own
estimate of the impact was several times larger than the numbers above; those
numbers are my own recount over all 300 crops and are the ones this document is
frozen against.

Versions 1 to 3 stay on Kaggle rather than being deleted, marked defective, so
that anyone who pulled them can tell what they have.
