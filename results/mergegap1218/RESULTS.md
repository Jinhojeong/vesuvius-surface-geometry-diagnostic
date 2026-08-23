# Does the published surface prediction merge touching sheets, and does the gap predict it?

Executed against PREREGISTRATION.md (`a8ae5f7740a5079a`) and its AMENDMENT_1
(`1c8c9419c9ee69ba`). Measured on all 300 contact crops of the tight-contact
validation set version 5, with the 60 single-sheet crops as a negative control.
The prediction is the published PHerc1218 `m7` surface, read at level 1, the same
grid the labels and the crops sit on.

## Answer

**Inside the contact regime the merge rate does not depend on the gap.**

| gap band | microns | crops | merge rate, pooled | 95 percent interval | median interval length |
| --- | --- | --- | --- | --- | --- |
| 0 to 2 vox | under 34.6 | 58 | 0.164 | 0.117 to 0.212 | 19.0 steps |
| 2 to 4 | 34.6 to 69.1 | 58 | 0.236 | 0.174 to 0.304 | 19.0 |
| 4 to 6 | 69.1 to 103.7 | 60 | 0.235 | 0.170 to 0.300 | 19.0 |
| 6 to 10 | 103.7 to 172.8 | 60 | 0.229 | 0.172 to 0.292 | 17.5 |
| 10 and above | above 172.8 | 60 | 0.116 | 0.075 to 0.159 | 27.0 |

Every pair among the four contact bands overlaps every other. About one sampled
contact point in five has m7 predicting continuously across the whole labelled
gap between two sheets a repair had separated, and that is as true at 30 microns
as at 150.

**The preregistration named this as the more interesting outcome and it is the
one that happened.** The frozen text says a flat curve "would say the
prediction's failure to separate sheets is not driven by how close they are,
which would point at something other than resolution". My own recorded
expectation was a merge rate falling with gap, and it is wrong.

## The loosest band is not evidence, because of a length bias I had to check

The widest band looks like a drop, 0.116 against about 0.23. It is not
comparable, and the reason is in the last column above. "Merged" requires m7
positive at every 0.5-voxel step across the gap interval, and that interval is 27
steps wide in the widest band against about 19 in the other four. A longer
interval has more chances to contain a zero, so part of that drop is geometry.

The length-free version is the per-step positive rate. It asks what share of
steps inside the gap interval have m7 positive at all, which no interval length
can inflate or deflate.

| band | all-steps rate | per-step rate |
| --- | --- | --- |
| under 34.6 um | 0.164 | 0.390 |
| 34.6 to 69.1 | 0.236 | 0.395 |
| 69.1 to 103.7 | 0.235 | 0.393 |
| 103.7 to 172.8 | 0.229 | 0.385 |
| above 172.8 | 0.116 | 0.302 |

The per-step rate is flat to within 0.01 across the four contact bands, which is
the same conclusion by a measure that cannot be confounded by interval length.
The widest band falls on both measures, so some of its drop is real, but its
all-steps figure overstates it and this file does not use it to claim a trend.

The per-crop correlation between merge fraction and gap in microns is weak and
negative, r = −0.136 with a bootstrap interval of −0.216 to −0.059. It excludes
zero and it is carried by the widest band.

## The other number, which is the one I did not go looking for

The per-step rate is about 0.39. Across the four contact bands, m7 is positive at
roughly two of every five steps inside a gap the labels call empty. That is a
large amount of predicted surface standing in the space between two sheets, and
it is a different statement from the merge rate, which asks only whether the
prediction bridges the gap end to end.

## The failure conditions, with their numbers

Four were frozen. Three did not fire and one is uncomfortably close.

The negative control passed cleanly. Zero of the 60 single-sheet crops produced a
merge score, as the frozen document requires, so the gap-interval definition is
scoring what it claims.

The discard rate is 33.1 percent against a gate of one third. That is a pass by
two tenths of a point, which is close enough that it should be read as a warning
rather than a clearance. The discards are 966 points whose walk never reached B,
203 with an empty gap interval, 54 crops with too few eligible points, and 9
leaving the crop.

No band fell below 40 crops; the smallest is 58.

The flat-curve condition is the one that fired, and the frozen document says a
flat curve is the result rather than a failure.

## The first run fired the discard gate, and that was my geometry

The original run discarded 55.1 percent. The frozen document did not say where on
instance A to sample, and I sampled uniformly over the whole instance. Each crop
is centred on the split site, so the sheets face each other near the centre and
diverge away from it, and one centroid-to-centroid normal cannot point at B from
a distant point. Reach rate fell monotonically from 91.5 percent within 10 voxels
of the centre to 6.1 percent beyond 70. AMENDMENT_1 restricted sampling to 20
voxels of the centre, 346 microns, chosen where reach was still above 80 percent
and before any merge rate was computed under it.

That is the second time a preregistered discard gate has caught an anchoring
error of mine on this data. The first was the ridge measurement.

## Limits

One scroll, one prediction, one repair. The labels are repaired automatic labels,
so a merge is m7 disagreeing with a repaired labelling and not proof that m7 is
wrong about the papyrus. The crops sample sites where a repair fired, so this is
the contact regime and the rate here is not the scroll's merge rate. The normal
is a straight line between two centroids and is a poor normal wherever the sheets
are not locally parallel, which is part of what the reach failure was measuring.
The set and this measurement share an author, so any systematic error in the set
is inherited here rather than tested.

`p15_merge.py` produces the table, `p15_stats.py` the intervals, `p15_len.py` the
length check, `p15_diag.py` the reach diagnosis.
