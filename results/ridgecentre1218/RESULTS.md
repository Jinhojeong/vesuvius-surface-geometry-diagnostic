# Do the repaired PHerc1218 labels sit on the CT ridge?

Executed against PREREGISTRATION.md (`ad186ffbf127216b`). Measured on the 209
contact crops of the published tight-contact set that carry both split
instances, per sheet, with the normal oriented away from an assumed winding
centre so that a negative offset means the CT ridge sits inward of the
label-run centre. That orientation rule turned out to rest on a constant I had
not checked, which is the subject of the sensitivity section below. It does not
change the answer, and the reason it cannot is worth more than the answer.

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

## The other two failure conditions, with their numbers

The preregistration set three ways for this to be wrong and only the third is
quoted above, so here are the other two. Sites are discarded when the corridor
leaves the crop, when the labelled run has no interior, or when the profile is
flat, and that came to 556 of 5,408 sites, a rate of 10.3 percent against the
one-third threshold. The median absolute offset does grow with the corridor,
0.63 then 0.75 then 1.25 voxels, but sub-proportionally, falling from 0.21 to
0.19 to 0.16 of the half-width. Proportional growth would have meant the number
was made by the corridor. Growth slower than the corridor is what bounded noise
around a true zero looks like.

## Sensitivity, including a constant I should have checked first

The orientation rule needs a winding centre and I supplied one, half of 7,593
in each of y and x, without verifying that against the volume. It appears
nowhere else in this project, and the crops' own coordinates do not support it.
Every site sits at y below 3,754 and x below 3,738, so a point at 3,796 is off
the edge of the sampled region rather than at its centre. Fitting an axis from
the crops' own normals, treating each instance-to-instance direction as radial
and solving for the point that best explains them, gives y 2,380 and x 2,390.

The answer holds under all three rules, and under one that needs no axis at
all:

| orientation rule | median | share negative |
| --- | --- | --- |
| the constant I used | 0.00 | 49.8 percent |
| axis fitted from the normals | 0.00 | 46.4 percent |
| centroid of the sites | 0.00 | 46.4 percent |

The axis-free test is the one that settles it. Before any orientation is
applied, the roughly 24 point measurements inside a single crop do not agree
with each other on a sign. In 162 of 209 crops between 35 and 65 percent of the
points are negative, no crop is outside 15 to 85 percent, and the median crop
sits at 47.6 percent. **There is no signed displacement at the crop level for an
orientation rule to orient.** A wrong axis can only mislabel which way a real
effect points. It cannot manufacture a null out of one, and here the effect is
absent before the axis is consulted at all.

One more sensitivity. The walk that finds a labelled run's extent stops at 20
voxels, which 173 of the runs reach, mostly where the instance-to-instance
normal runs shallow across its own sheet. Dropping those leaves the median at
0.00 with 45.0 percent negative and the median absolute at 0.75.

`p14_sensitivity.py` reproduces this section.

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

The orientation rule is sound in construction but its winding centre was
guessed rather than measured, and the section above is the repair. Anyone
extending this should fit the axis rather than assume it. These are repaired
automatic labels, not hand annotation, so a null says the
repair and the CT agree on placement to within the scatter, not that either is
correct. The crops come from contact sites, so this samples the regime where
sheets are close rather than the scroll as a whole. The corridor bounds the
reportable offset by construction, which is why the sweep is reported rather
than a single number. One scroll.
