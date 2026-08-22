# Amendment 1 to the ridge-versus-centre preregistration

Amends PREREGISTRATION.md, sha256 prefix `ad186ffbf127216b`. Written 2026-08-22,
after version 5 of the tight-contact set was published and before this
measurement was re-run on it.

## What could not stand as written

The data clause names "the 254 contact crops of the published tight-contact set,
version 2" and its 209-crop both-instances subset. Both numbers are gone. That
set has been through four defect corrections since, recorded in
results/tightgap1218: the crop CT was read from the wrong pyramid level, the
label stitching let overlapping blocks overwrite the split pair, and the block
whose numbering defines A_id and B_id was chosen by an arbitrary tie-break. The
current release is version 5, 300 contact crops, all 300 carrying both split
ids, with the shipped gap reproducing exactly from each crop's own labels.

So the measurement is being re-run on version 5. Sample size goes from 209 to
300, and more importantly the labels are different objects: the earlier run
measured runs through instance arrays that had lost part of the pair to the
overwrite.

## The winding centre, which was never sourced

The frozen document orients the normal outward from "the volume's centre in the
crop's own z slice". The original run supplied half of 7,593 for that, which is
the level-0 width while the sites are level-1 coordinates, so no site could ever
reach it. That constant is retired.

This run reports the primary under an axis fitted from the crops' own normals,
solving for the point that best explains each instance-to-instance direction as
radial. Two sensitivities are reported alongside. One is the level-1 volume centre at
1,898.5 in y and x. The other is the axis-free test, which asks whether the roughly 24
point measurements inside a single crop agree on a sign before any orientation
is applied. The axis-free test is the one that cannot be gamed by an axis
choice, and it is stated as such rather than as a supplement.

## What does not change

The estimator. That means per-sheet sampling of 12 points on each instance, the
labelled run's midpoint along the oriented normal, and the CT maximum on a
linearly interpolated profile at 0.25-voxel steps. Also the ±3, ±4 and ±8
corridor sweep, the discard rules, the per-crop median as the unit, and a
bootstrap over crops. The sign convention stands, negative means the ridge sits
inward of the label-run centre. So do the three failure conditions.

The CT is the crop's own shipped intensity, which version 5 verified byte
identical to CT level 1 at each crop's site on a 30-crop sample, so this run
reads no volume of its own.

## What this predicts, frozen here

The earlier corrected run on 209 crops gave medians of 0.00 at every corridor,
intervals containing zero, 43.1 to 48.3 percent negative, and median absolute
offsets of 0.50, 0.50 and 1.13 voxels. Suppose the null there was real rather than an artefact of damaged labels.
Version 5 should then reproduce it. Medians at or within a quarter voxel of
zero, intervals containing zero, a sign split between 40 and 60 percent, and
median absolute offsets within about 0.2 voxels of those figures.

A systematic displacement appearing now would be the more interesting outcome,
because it would mean the earlier null came from labels that had lost part of the
pair. That reading is written here so it cannot be constructed afterwards.

## What would make this run wrong

The three conditions from the frozen document stand unchanged. If the discard
rate exceeds a third of sites, the estimator is not measuring the population it
claims. If the median offset tracks the corridor width proportionally, the number
is an artefact of the corridor. If the sign distribution is near even, there is
no systematic relation and that is the correct report.

One addition. If the answer differs between the fitted axis and the level-1
volume centre by more than the bootstrap interval, no axis-dependent claim is
made at all and only the axis-free test is reported.
