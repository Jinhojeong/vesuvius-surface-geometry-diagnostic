# Preregistration: do my repaired PHerc1218 labels sit on the CT ridge?

Frozen 2026-08-20, before any offset was computed. Written because a signed
quantity is only meaningful if its sign convention is fixed in advance, which
is exactly the failure that made this question worth asking.

## Why

The tight-contact validation set ships instance labels as the thing models are
scored against, and another contributor is pinning it as a held-out set. If
those labels sit systematically off the CT ridge, every model scored against
them inherits that displacement, and the set's card should say so with a
number rather than the qualitative limit it carries today.

The question was prompted by TAUIL-Abd-Elilah's correction on villa#193, where
a signed placement summary on Scroll1A was withdrawn because Hessian normals
were never sign-anchored, then reissued oriented as a median of −0.5150 voxels
with 186 of 189 cubes negative. That is a different scroll, a different label
source and a different sampling of points, so it is context here and not a
result to confirm or deny. What is measured below is my own labels.

## Data

The 254 contact crops of the published tight-contact set, version 2, each
carrying the level-0 CT and the repaired v2.0 instance labels on the same grid.
Only crops whose two split instances are both present are used, since the
measurement needs a labelled run with an identifiable centre; that is 209 of
254 by the set's own manifest field.

## Sign convention, fixed here

The normal at a site is taken from the two instance centroids, oriented from
instance A to instance B, which is a construction with a definite direction and
not an eigenvector. The reference direction is outward from the scroll axis:
the volume's centre in the crop's own z slice defines an inward vector, and the
normal is flipped when needed so that positive always means outward from the
scroll centre. The offset is `ridge − centre` along that oriented normal, so
**negative means the CT ridge sits inward of the label-run centre** and
positive means it sits outward. Both terms are measured along the same oriented
normal.

## Definitions

The label-run centre at a site is the midpoint of the labelled run through that
site along the oriented normal. The CT ridge is the position of maximum CT
intensity along the same line within a corridor of ±4 voxels of the centre,
found on a linearly interpolated profile sampled every 0.25 voxels. A site is
discarded, and the discard counted, when the corridor leaves the crop, when the
labelled run has no interior, or when the profile is flat to within one grey
level.

## Endpoints

Primary is the median of per-crop median offsets, with a bootstrap interval
over crops. Secondary are the median `|offset|`, which is sign-invariant and
therefore comparable with the historical magnitude, and the fraction of crops
whose median is negative. The corridor is swept over ±3, ±4 and ±8 voxels and
all three are reported, because a magnitude that moves with the corridor is an
estimator property rather than a physical one.

## Declared limits, written before the numbers exist

These are repaired automatic labels rather than hand annotation, so a nonzero
offset says the repair and the CT disagree, not that either is wrong. One
scroll. The corridor bounds the offset by construction, so the estimator cannot
report a displacement larger than the corridor and the sweep is the only guard
against that mattering. Crops come from contact sites, so this samples the
regime where sheets are close rather than the scroll as a whole.

## What would make this wrong

If the discard rate exceeds a third of sites, the estimator is not measuring
the population it claims to. If the median offset tracks the corridor width
proportionally across the sweep, the number is an artefact of the corridor and
this document says so in the results instead of a displacement being claimed.
If the sign distribution is near even, there is no systematic relation and the
correct report is that.
