# Preregistration: how often does the public surface prediction merge two touching sheets, as a function of the measured gap?

Frozen 2026-08-23, before any merge rate was computed.

## Why

The team's own First Letters workflow names "segmentation jumps between
neighbouring layers or drifts away from the papyrus" as one of three reasons a
rendered surface shows no letters. That failure has a precursor which nobody has
measured on real CT in physical units: the surface **prediction** merging two
sheets that are close but distinct, before any tracer touches it.

The tight-contact validation set exists to make that measurable. It ships 300
PHerc1218 crops, each centred on a ray-validated split between two instances,
each carrying both instance ids and a gap measured along the label normal, binned
across five bands from under 34.6 um to above 172.8 um. Nobody has yet used it
for anything, including me. This is its first consumer result, and I am the
consumer, which is worth saying plainly because it means the set and the
measurement share an author and share any systematic error.

The prediction under test is the published `m7` surface prediction for PHerc1218,
`20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr`, read at level 1,
which is the same grid the labels and the crops sit on. Verified: its level-1
shape is 11,624 by 3,797 by 3,797, identical to the CT level-1 and the label grid.

## The measurement

For each contact crop, the two ids `A_id` and `B_id` are the instances a repair
separated. The normal is taken from the A centroid to the B centroid inside the
crop, which is a construction with a definite direction and not an eigenvector.

Twelve points are sampled on instance A, seeded, without replacement. From each
point the walk proceeds along the normal toward B and records three positions:
where the labelled run of A ends, where the labelled run of B begins, and every
step between them. That interval between the two labelled runs is the **gap
interval** and is the only place this measurement looks.

A sampled point is scored **merged** when the m7 prediction is positive at every
0.5-voxel step across its gap interval, and **separated** when m7 is zero at one
or more steps. Points whose gap interval is empty, whose walk leaves the crop, or
which never reach B within 40 voxels are discarded and the discards counted.

## Endpoints

Primary is the merge rate per gap band, computed as the median over crops of the
per-crop fraction of sampled points scored merged. Secondary are the pooled
point-level merge rate per band, and the merge rate as a function of the crop's
own measured gap in microns rather than by band.

The 60 single-sheet control crops are run through the identical code as a
negative control. They have no second instance, so every point there must be
discarded for an empty gap interval. Any control crop that produces a merge
score at all means the estimator is scoring something other than what it claims.

## Units, stated because I got this wrong before

The gap is in level-1 voxels of 17.28 um. The preregistration of the
tight-contact set called them level-0 voxels and that was wrong; the correction
is in results/tightgap1218. Every distance in this document and in its results is
reported in microns as well as voxels, on both sides of any comparison.

## Declared limits, written before the numbers exist

One scroll, one prediction, one repair. The instance labels are repaired
automatic labels rather than hand annotation, so a merge is disagreement between
m7 and a repaired labelling, not proof that m7 is wrong about the papyrus. The
crops sample sites where a repair fired, so this is the contact regime and not
the scroll as a whole; the merge rate here is not the scroll's merge rate. The
normal is a straight line between two centroids and will be a poor normal
wherever the two sheets are not locally parallel.

## What would make this wrong

If the discard rate exceeds a third of sampled points, the estimator is not
measuring the population it claims. If any control crop returns a merge score,
the gap-interval definition is not doing what this document says and the run is
void. If the merge rate is flat across all five bands, the measurement has no
resolving power on this data and that is the result, stated as such rather than
dressed as a null. If the merge rate is at a floor or a ceiling in every band,
the same applies.

## What I expect, so it cannot be constructed afterwards

I expect the merge rate to fall as the gap widens, because a prediction has more
room to resolve two sheets that are further apart. I do not know the shape or
whether the effect clears the band-to-band noise. A flat curve would be the more
interesting outcome, because it would say the prediction's failure to separate
sheets is not driven by how close they are, which would point at something other
than resolution.
