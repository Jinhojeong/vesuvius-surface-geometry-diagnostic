# Amendment 1 to the merge-versus-gap preregistration

Amends PREREGISTRATION.md, sha256 prefix `a8ae5f7740a5079a`. Written 2026-08-23,
after the first run fired the document's own discard gate and before any
corrected merge rate was computed.

## What fired

The frozen document says the run is wrong "if the discard rate exceeds a third of
sampled points". The first run discarded 2,358 of 4,276, a rate of 55.1 percent,
and 1,979 of those discards were points whose walk never reached instance B.

The negative control passed. Zero of the 60 single-sheet control crops produced a
merge score, which is what the document requires, so the gap-interval definition
itself is doing what it claims.

## Why, and it is my sampling rather than the data

The document says twelve points are sampled on instance A, and does not say
where on A. I sampled uniformly over the whole instance. Each crop is centred on
the split site, so the two sheets face each other near the crop centre and
diverge away from it, and a single centroid-to-centroid normal cannot point at B
from a point far out on A.

Reach rate against distance from the crop centre, over 4,800 points on 120 crops:

| distance from crop centre | points | reached B |
| --- | --- | --- |
| 0 to 10 vox | 470 | 91.5 percent |
| 10 to 20 vox | 1,301 | 83.1 percent |
| 20 to 30 vox | 1,056 | 67.5 percent |
| 30 to 40 vox | 800 | 49.0 percent |
| 40 to 50 vox | 558 | 27.4 percent |
| 50 to 60 vox | 315 | 17.8 percent |
| 60 to 70 vox | 188 | 14.4 percent |
| 70 to 80 vox | 82 | 6.1 percent |

Monotone, and steep. This is the same class of error as the ridge measurement's
first pass, where the estimator looked like it was failing on the data and was in
fact anchored wrongly, and it is the second time a preregistered discard gate has
caught it.

## What changes

Sampling is restricted to points on instance A within 20 voxels of the crop
centre, which is 346 um. That is the region the crop was built to contain, the
neighbourhood of the ray-validated split site, and it is where the two sheets
actually face each other. Twelve points as before, seeded as before, and a crop
with fewer than twelve eligible points inside that radius is discarded and
counted rather than backfilled from further out.

The radius is fixed at 20 voxels because that is where the reach rate is still
above 80 percent in the table above, not because it produced any particular merge
rate. No merge rate has been computed under it at the time of writing.

## What does not change

Everything else. That means the centroid-to-centroid normal, the gap interval
between the end of A's labelled run and the start of B's, and the merged
criterion of m7 positive at every 0.5-voxel step across that interval. Also the
twelve points, the seed, the per-crop fraction as the unit, and the band medians
and pooled rates as endpoints. The 60 controls stay as the negative control and
all four failure conditions stand.

## What this predicts, frozen here

The discard rate should fall below the one-third gate, since the dominant discard
was the reach failure and the new radius sits where reach is above 80 percent. If
it does not, the reach failure was not the whole story and this amendment is
incomplete. The band medians may move in either direction, and I have no
expectation about which, because the restricted sample is drawn from the contact
neighbourhood rather than the whole sheet and those are different populations.

The expectation recorded in the frozen document stands unchanged. I expect the
merge rate to fall as the gap widens, and a flat curve would be the more
interesting outcome.

## What would make this amendment wrong

If the discard rate stays above a third, the diagnosis is incomplete. If the
control crops start producing merge scores, the restriction has broken something
the first run had right. If fewer than 40 crops survive in any band, the
restricted sample is too thin to carry a band median and the band is reported as
underpowered rather than as a number.
