# Amendment 1 to the tight-contact validation set preregistration

Amends PREREGISTRATION.md, sha256 prefix `9d5d71dbaf45ab85`, written
2026-08-19 after the contact arm was extracted and before any control crop
existed.

## What could not be executed as written

Rule 6 said the control arm would be "drawn the same way from sites the
census examined and did not flag". That population is not recoverable from
the artifacts the rule points at. The repair record files log only flagged
candidates, so every site in them already carries a thickness ratio above
the flag threshold: across 119,468 sampled record sites the minimum ratio is
1.61 and the median is 2.33, and no site sits below 1.5. A first control run
written to the rule returned zero candidates, which is how this was found.

## Replacement rule, fixed here before any control crop is cut

Control sites are drawn from the repaired label blocks directly rather than
from the record files. A control site is a voxel on a labelled sheet whose
locally measured run thickness along the label normal is at or below the
tile median thickness, so it is single-sheet by the same measurement the
census uses to flag. Sites are drawn in block order with a fixed seed of 0
for the within-block pick, and every other rule from the frozen document is
unchanged, including the 128 cubed crop, the eight-octant CT requirement,
the quarter-volume overlap exclusion and the target of 60.

## What this changes about the claim

The control arm is now a re-derivation from labels rather than a leftover
population from the census, so it answers "what does a single-sheet crop
look like here" and not "what did the census decline to flag". Any
comparison that depends on the second reading is not supported by this set,
and the results file repeats this where the control is used.

## Also recorded here, from the contact arm

The frozen rule 4 requires CT in all eight octants of a crop. That rule
removes 15,840 candidate sites overall, and it falls hardest on the tightest
band: of 120 sampled sites at gaps under 2 voxels, 117 fail it and 3 pass.
The 0 to 2 voxel band therefore ships 14 crops rather than 60. This is
reported rather than repaired, per the frozen document's instruction that
realised counts are reported rather than padded. The reading is that the
tightest contacts sit disproportionately in crops that reach outside the
masked volume, which is itself a measurement about where tight contact
occurs and not only a sampling nuisance.
