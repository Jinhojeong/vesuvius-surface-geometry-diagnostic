# Do the repaired PHerc1218 labels sit on the CT ridge?

Executed against PREREGISTRATION.md (`ad186ffbf127216b`). Measured on the 209
contact crops of the published tight-contact set that carry both split
instances, per sheet, with the normal oriented outward from the scroll axis so
that a negative offset means the CT ridge sits inward of the label-run centre.

## Answer

No systematic displacement is detectable.

| corridor | crops | median offset | 95 percent interval | share negative | median absolute |
| --- | --- | --- | --- | --- | --- |
| ±3 vox | 209 | 0.00 | −0.25 to +0.12 | 49.8 percent | 0.63 |
| ±4 vox | 209 | 0.00 | −0.25 to +0.12 | 49.8 percent | 0.75 |
| ±8 vox | 209 | −0.25 | −0.62 to +0.12 | 51.2 percent | 1.25 |

The interval contains zero at every corridor and the sign split is even to
within half a percent at two of the three. **The preregistration named this
outcome in advance**, its third failure condition reads "if the sign
distribution is near even, there is no systematic relation and the correct
report is that". So this is the frozen reading rather than a post-hoc one.

## What the spread means, which is not nothing

A zero median is not zero error. Half of the per-crop medians sit more than
0.63 voxels from the ridge at the tightest corridor, and that spread is the sum
of real annotation scatter and estimator noise, which this design cannot
separate. Anyone scoring a model against these labels inherits that scatter as
a floor on how precisely the comparison can resolve surface placement. What
they do not inherit is a systematic bias in one direction.

## An implementation error, caught by the preregistered failure conditions

The first pass anchored the labelled run at each crop's centre voxel and
required that voxel to carry one of the two split ids. It carries a different
instance 24 times in 40, so 81 of 209 crops were discarded, the discard rate
passed the preregistered one-third threshold, and the offsets drifted with the
corridor. That looked like the estimator failing on the data. It was the
anchoring being wrong, which a check of the label at the crop centre settled
before anything was written up. The frozen text says "the labelled run through
that site", and a site is a point on a labelled sheet, so points are now
sampled on each of the two instances directly. With that fixed, no crop is
discarded for a run reason and the corridor sweep behaves.

## Context, not a comparison

This question was prompted by a correction on villa#193, where a signed
placement summary on Scroll1A was withdrawn for unanchored normals and reissued
as a median of −0.5150 voxels with 186 of 189 cubes negative. That is a
different scroll, a different label source and a different sampling of points.
The measurement here does not reproduce that relation on these labels, and it
is not evidence against it either.

## Limits

These are repaired automatic labels, not hand annotation, so a null says the
repair and the CT agree on placement to within the scatter, not that either is
correct. The crops come from contact sites, so this samples the regime where
sheets are close rather than the scroll as a whole. The corridor bounds the
reportable offset by construction, which is why the sweep is reported rather
than a single number. One scroll.
