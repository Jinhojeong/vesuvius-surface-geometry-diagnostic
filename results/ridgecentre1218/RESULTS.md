# Do the repaired PHerc1218 labels sit on the CT ridge?

Executed against PREREGISTRATION.md (`ad186ffbf127216b`). Measured on the 209
contact crops of the published tight-contact set that carry both split
instances, per sheet, with the normal oriented away from an assumed winding
centre so that a negative offset means the CT ridge sits inward of the
label-run centre.

**Re-run on the version-5 membership, 2026-08-22.** The table below is the
209-crop run on an earlier release. The current answer is at the end of this
file, measured on all 300 version-5 crops with correctly stitched labels, and it
agrees.

**This file was rewritten on 2026-08-20 after the first version was withdrawn.**
The first version read the CT out of the shipped crop arrays, which turned out
to hold CT from a different physical region than their own labels. The answer
below is the same, but the first version had no right to it, and the section on
the defect is the part worth reading.

## Answer

No systematic displacement is detectable.

| corridor | crops | median offset | 95 percent interval | share negative | median absolute |
| --- | --- | --- | --- | --- | --- |
| ±3 vox | 209 | 0.00 | 0.00 to 0.00 | 43.1 percent | 0.50 |
| ±4 vox | 209 | 0.00 | −0.12 to +0.12 | 44.0 percent | 0.50 |
| ±8 vox | 209 | 0.00 | −0.25 to +0.25 | 48.3 percent | 1.13 |

The median is zero at every corridor and the interval contains zero at every
corridor, tightly enough at ±3 that both bootstrap bounds land on zero. **The
preregistration named this outcome in advance**, its third failure condition
reads "if the sign distribution is near even, there is no systematic relation
and the correct report is that".

## The defect that made the first version worthless

The repaired instance labels live on the CT's level-1 grid. The repair blocks
run to z0 11,368 with 256-voxel blocks, an extent of 11,624, which is the
level-1 z extent exactly. In y and x the block grid covers 4,096, which spans
level 1's 3,797 and falls far short of level 0's 7,593.
`p11_crops.py`, which built the tight-contact crops, took site coordinates from
those level-1 block names and then read the CT out of the level-0 array, shape
23,247 by 7,593 by 7,593. So every crop shipped CT from roughly half the true
offset, a real region of the scroll but not the one its labels describe.

Two things fix it, and a third that I first reached for does not.

The block arithmetic. The repair blocks run to z0 11,368 with 256-voxel blocks,
an extent of 11,624, which is the
level-1 z extent exactly. In y and x the block grid covers 4,096, which spans
level 1's 3,797 and falls far short of level 0's 7,593. No site can reach the midpoint of a 7,593-wide axis.

The supersample identity. Reading level 0 at doubled indices over a 256 cube
and mean-pooling back to 128 reproduces the level-1 values exactly, so level 1
is the physically matching region rather than merely a better-correlated one.
The shipped arrays are not equal to either.

What does not work is correlating a crop's CT with its own label mask, which is
what I used first. That contrast depends almost entirely on how much empty
space the crop contains, because in a fully dense region the CT barely
separates sheet from gap by absolute intensity. Across crops it runs 0.64 at 39
percent emptiness down to 0.05 at zero emptiness. Quoting a range from a handful
of crops as though it described the set was wrong, and the corrected form is the
shift test below.

The shift test is the one that settles alignment. Reading the CT at the same
site but offset by 20 or 40 voxels scores worse than reading it at the site
itself on 30 of 30 sampled crops. The mean CT under the label exceeds the mean
under the background on 97.0 percent of all 300 crops, by a median of 6.1 grey
levels. A
misaligned CT would not care where it was read.

**A null measured against an unrelated array is not a null, it is an artefact.**
The first version's headline evidence was that the roughly 24 point
measurements inside a single crop did not agree on a sign. That is exactly what
two unrelated arrays produce, so the strongest-sounding part of the argument was
the clearest symptom of the bug. Re-running the same estimator with only the CT
source changed keeps the within-crop disagreement, 173 of 209 crops between 35
and 65 percent negative against 162 before, but now it sits on a CT that
actually tracks the labels.

`p11_gridcheck.py` and `v3_diag.py` reproduce the alignment evidence, `p14_correctct.py` is the
published estimator with one line changed.

## What the spread means, which is not nothing

A zero median is not zero error. Half of the per-crop medians sit more than 0.50
voxels from the ridge at the tightest corridor, and that spread is the sum of
real annotation scatter and estimator noise, which this design cannot separate.
Anyone scoring a model against these labels inherits that scatter as a floor on
how precisely the comparison can resolve surface placement. What they do not
inherit is a systematic bias in one direction. The corrected figure is tighter
than the withdrawn one, 0.50 against 0.63 voxels, which is what a CT that
actually tracks its labels should give.

## The other two failure conditions, with their numbers

Sites are discarded when the corridor leaves the crop, when the labelled run has
no interior, or when the profile is flat. On the correct CT that is 351 of 5,276
sites, a rate of 6.7 percent against the preregistered one-third threshold, and
the flat-profile discard falls to zero because a CT that contains the sheet is
never flat across it. The median absolute offset grows with the corridor, 0.50
then 0.50 then 1.13 voxels, but sub-proportionally, falling from 0.17 to 0.13 to
0.14 of the half-width. Proportional growth would have meant the number was made
by the corridor.

## Sensitivity, including a constant I should have checked first

The orientation rule needs a winding centre and I supplied one, half of 7,593
in each of y and x, without verifying it. In hindsight that constant was the
first visible symptom of the grid bug, since 7,593 is the level-0 width and the
sites are level-1 coordinates, so no site could ever reach its midpoint.
Fitting an axis from the crops' own normals gives y 2,380 and x 2,390.

The answer held under all three rules on the withdrawn data, at 49.8, 46.4 and
46.4 percent negative, and the axis-free within-crop test behaves the same way
on the corrected CT. A wrong axis can only mislabel which way a real effect
points. It cannot manufacture a null out of one.

## An implementation error, caught by the preregistered failure conditions

The first pass anchored the labelled run at each crop's centre voxel and
required that voxel to carry one of the two split ids. It carries a different
instance 24 times in 40, so 81 of 209 crops were discarded, the discard rate
passed the preregistered one-third threshold, and the offsets drifted with the
corridor. That looked like the estimator failing on the data. It was the
anchoring being wrong. The frozen text says "the labelled run through that
site", and a site is a point on a labelled sheet, so points are sampled on each
of the two instances directly.

## Context, not a comparison

This question was prompted by a correction on villa#193, where a signed
placement summary on Scroll1A was withdrawn for unanchored normals and reissued
as a median of −0.5150 voxels with 186 of 189 cubes negative. That is a
different scroll, a different label source and a different sampling of points.
The measurement here does not reproduce that relation on these labels, and it
is not evidence against it either.

## Limits

The 209 crops were selected by a rule whose CT test ran on the displaced volume,
so the sample is correctly measured but not the sample the preregistration
intended. That is a property of the published set rather than of this
measurement, and it is being fixed in the set rather than here. These are
repaired automatic labels, not hand annotation, so a null says the repair and
the CT agree on placement to within the scatter, not that either is correct. The
crops come from contact sites, so this samples the regime where sheets are close
rather than the scroll as a whole. The corridor bounds the reportable offset by
construction, which is why the sweep is reported rather than a single number.
One scroll.

## Re-run on version 5, 2026-08-22

Per AMENDMENT_1 (`953017184b93ba59`), frozen before this ran. The earlier answer
sat on 209 crops from a release whose labels had lost part of the split pair to
an overlap overwrite, so it had to be redone once the labels were fixed. This run
uses all 300 version-5 contact crops, every one of which carries both split ids.

The winding-centre constant is retired. The primary is oriented by an axis fitted
from the crops' own normals, at y 1,767 and x 1,776. That is close to the level-1
volume centre of 1,898.5, which is what a scroll roughly centred in its volume
should give, and it is a long way from the 3,796.5 the first run used.

| corridor | crops | median | 95 percent interval | share negative | median absolute |
| --- | --- | --- | --- | --- | --- |
| ±3 vox | 300 | −0.06 | −0.25 to 0.00 | 50.0 percent | 0.50 |
| ±4 vox | 300 | −0.13 | −0.25 to 0.00 | 50.3 percent | 0.63 |
| ±8 vox | 300 | −0.25 | −0.38 to 0.00 | 53.3 percent | 1.19 |

**The reading is unchanged. There is no systematic displacement.** The
preregistration's third failure condition fires, the sign split is 50.0, 50.3 and
53.3 percent negative, and it says in advance that a near-even split means no
systematic relation. The axis-free test agrees, 248 of 300 crops have between 35
and 65 percent of their own points negative and the median crop sits at exactly
50 percent. The other two conditions do not fire. The discard rate is 296 of
7,417 sites, 4.0 percent against the one-third threshold, and the median absolute
offset falls as a share of the corridor, 0.167 then 0.156 then 0.148.

**One thing did move and it should be said plainly.** On 209 crops the medians
were 0.00 at every corridor with zero interior to the interval. On 300 they are
−0.06, −0.13 and −0.25, and zero sits at the upper bound of all three intervals
rather than inside them. Every point estimate is on the same side. That is a
hint of a small inward bias that the earlier run did not show, and it is smaller
than the sampling step of the profile at the tightest corridor. It does not
survive the sign test, so it is not being claimed. It is being recorded so that
anyone who repeats this on more data knows where to look.

Every AMENDMENT_1 prediction held. Medians within a quarter voxel of zero,
intervals containing zero, sign split inside 40 to 60 percent, and median
absolute offsets within 0.2 voxels of the earlier figures. The amendment's added
condition also passes, the fitted axis and the volume centre differ by at most
0.06 voxels against interval widths of 0.25 and 0.38, so the answer is not an
artefact of the axis choice.

The scatter is what a consumer inherits, and it is essentially unchanged. Half
the per-crop medians sit more than 0.50 voxels from the ridge at the tightest
corridor. `p14_v5.py` reproduces this section.

