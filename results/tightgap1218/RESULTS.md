# Tight-contact validation set on PHerc1218: what was built and what it costs

Executed against PREREGISTRATION.md (`9d5d71dbaf45ab85`) and AMENDMENT_1.md
(`a07c2f86f39b3bd8`). Every count below is realised, not targeted.

> **Defect notice, 2026-08-20.** The CT array in every published crop was read
> from the wrong pyramid level and is not the CT at the crop's own labels. The
> labels, the gaps and the gap population are unaffected. The crop intensities,
> the eight-octant acceptance filter and everything downstream of it are.
> Kaggle versions 1 and 2 should not be used for anything that reads the
> intensity array. See the third correction at the end of this file.


## Why it exists

Surface-model validation currently runs on data whose sheets are not in
contact. A contributor training on the public patch set measured his held-out
geometry as median inter-sheet gap 15.6 voxels with 0.02 percent of near-band
positives under 4 voxels, and reported that those labels merge below about 5
voxels, so nothing about touching sheets is measurable there. This set
supplies that regime from real CT.

## The gap population

The normal-direction gap was measured at every ray-validated split site in the
repair records, 49,295 sites in total:

| gap band, level-0 voxels | sites |
| --- | --- |
| 0 to 2 | 429 |
| 2 to 4 | 5,929 |
| 4 to 6 | 12,578 |
| 6 to 10 | 17,545 |
| 10 and above | 12,814 |

12.9 percent of sites sit under 4 voxels, against 0.02 percent in the set that
prompted this. The two failure conditions the preregistration named both pass:
sites spread over 52 slabs with the largest holding 4.4 percent, and the most
repeated instance pair appears 3 times.

## The shipped crops

128 cubed at level 0, centred on the site, with intensity, the instance label
map, the binary surface label, the measured gap, the band, the coordinates and
the CT emptiness of the crop.

| band | crops |
| --- | --- |
| 0 to 2 | 14 |
| 2 to 4 | 60 |
| 4 to 6 | 60 |
| 6 to 10 | 60 |
| 10 and above | 60 |
| control, single sheet | 60 |

## What it costs, measured

**The tightest band is short and the reason is structural.** The frozen rule
requiring CT in all eight octants removes 15,840 candidate sites overall, and
it falls hardest at the tight end. The 0 to 2 band offers 429 candidate sites,
the octant rule rejects 409 of them and the overlap rule another 6, so the band
ships 14 crops rather than 60. A pilot figure of 117 rejections in 120 sampled
sites appeared here and in AMENDMENT_1.md and is superseded, see the second
correction below.
The reading is that the tightest contacts sit disproportionately in crops that
reach outside the masked volume, which is a fact about where tight contact
occurs rather than only a sampling nuisance.

**Four in five contact crops verifiably contain both sheets.** Instance ids in
the source blocks are block-local, so a crop spanning blocks needs a
disambiguated id space; each contributing block is offset by a unique base and
the split pair is remapped through the base of the site's own block. With that
done, 209 of 254 contact crops contain both ids the split separated, 3 contain
one, and 42 contain neither. The 42 are a boundary effect, they sit on a median
of 4 contributing blocks against 2 for the crops that pass. Two earlier counts
are superseded. A pass without the id disambiguation reported 231, wrong because
an id from a neighbouring block can collide with the pair. A later count of 213
of 260 was also wrong, because six crops left over from a pilot extraction were
still in the output directory; see the correction below.

**The octant rule is weak and the emptiness field is how to filter.** It only
asks for 1 percent CT per octant, so a crop with a large empty wedge can pass.
Measured across the contact arm, median emptiness is 0.0002, the ninetieth
percentile is 0.295, 20.5 percent of crops are more than a tenth empty and 11.4
percent are more than a quarter empty. The control arm behaves the same way,
20.0 and 13.3 percent. Consumers wanting full cubes should filter on
`ct_empty_frac`; the selection rule was not changed after the fact.

## Limits

The instance labels are the repaired v2.0 labels, not hand annotation, so this
set measures agreement with a repaired automatic labelling and not with truth.
The gap definition is normal-direction and local and will disagree with pairwise
or connected-component definitions by construction. The control arm is a
re-derivation from labels rather than a leftover census population, per
amendment 1, so it answers what a single-sheet crop looks like here and not
what the census declined to flag. One scroll.

## Correction, 2026-08-19

aviad12g audited the published set against this document and found a count
mismatch: the Kaggle v1 manifest carried 260 contact crops with 17 in the 0 to 2
band and 63 in the 10 and above band, against the 254 and 14 and 60 recorded
here. He also named the cause correctly, `p11_crops.py` writes into an output
directory it does not require to be empty, and the later label and figure passes
glob everything in it, so six crops from an earlier pilot extraction with a
different per-band target were carried along. Those six are not members under
the frozen rule, because the pilot's smaller accepted list changed which later
sites the overlap exclusion rejected.

Membership was re-derived by replaying the frozen rule in one clean pass
(`replay_rule.py`, membership_replay.json). The authoritative set is 254 contact
crops, 14 / 60 / 60 / 60 / 60 by band, plus the 60 controls, and the six orphans
are listed in that file by name. Kaggle version 2 ships that set with a
recomputed manifest. Every derived statistic above has been recomputed on the
254, which moved the both-instances count from 213 of 260 to 209 of 254 and the
emptiness percentiles slightly. The band table and the gap population are
unchanged, since the gap measurement never depended on crop membership.

## Correction, 2026-08-20

The rejection figure for the tightest band was wrong. This file and
AMENDMENT_1.md both said that of 120 sampled sites under 2 voxels, 117 failed
the eight-octant rule and 3 passed, and then concluded that the band therefore
ships 14 crops. Three passing sites cannot yield fourteen crops, so the
sentence never followed from its own numbers. The figure traces to no committed
script and to no line of the extraction log, which is the same unsourced-number
problem the ridge measurement turned up in its own orientation rule.

The real breakdown is a replay rather than a sample. `p11_crops.py` tests in a
fixed order, the band target gate then out of bounds then overlap then read then
octant, and the recorded read failure count is zero. So once the accepted list
is known, the rejection reason for every census row is fixed by arithmetic and
needs no CT. The accepted list is the 254 published crops in census order.
Replaying
that reproduces the frozen `crops_summary.json` exactly, 254 accepted with
15,840 octant, 679 overlap and 24 out of bounds, which is what makes the
per-band split below trustworthy rather than merely plausible.

| band | census sites | accepted | octant | overlap | out of bounds |
| --- | --- | --- | --- | --- | --- |
| 0 to 2 | 429 | 14 | 409 | 6 | 0 |
| 2 to 4 | 5,929 | 60 | 2,154 | 116 | 4 |
| 4 to 6 | 12,578 | 60 | 4,145 | 181 | 11 |
| 6 to 10 | 17,545 | 60 | 5,046 | 170 | 9 |
| 10 and above | 12,814 | 60 | 4,086 | 206 | 0 |

The 0 to 2 band is the only one that exhausted its census population, since the
other four hit the 60 target and stopped early. So its 95.3 percent octant
rejection rate is a rate over everything available rather than over a sample.

The table is an exact account of what the code did. It is not an account of
where tight contacts sit, because the octant test it replays ran against the
wrong CT, which the next section explains. The reading that used to sit here,
that the tightest contacts reach outside the masked volume, is withdrawn along
with the filter that suggested it.

AMENDMENT_1.md is left as written. It is hash chained at `a07c2f86f39b3bd8` and
referenced by `control_summary.json` and `p11_control2.py`, so it is corrected
here rather than edited. `p11_bandreplay.py` reproduces the table.

## Correction, 2026-08-20, second entry: the CT is from the wrong grid

The repaired instance labels live on the CT's level-1 grid. The repair blocks
run to z0 11,368 with 256-voxel blocks, an extent of 11,624, which is the
level-1 array's z extent exactly, and 3,797 matches in y and x the same way.
`p11_crops.py` and `p11_control2.py` take site coordinates from those level-1
block names and then open `ct["0"]`, shape 23,247 by 7,593 by 7,593. Every one
of the 254 contact crops and 60 control crops therefore ships CT from about
half the true offset. It is a real region of PHerc1218. It is not the region
its own labels describe.

Two things fix it, and a third that I first reached for does not.

The block arithmetic. The repair blocks run to z0 11,368 with 256-voxel blocks,
an extent of 11,624, which is the level-1 z extent exactly, and 3,797 matches
the same way in y and x. No site can reach the midpoint of a 7,593-wide axis.

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
shift test.

The shift test is the one that settles alignment. On the corrected crops, reading the CT at the same site but
offset by 20 or 40 voxels scores worse than reading it at the site itself on 30
of 30 sampled crops. The mean CT under the label exceeds the mean under the
background on 97.0 percent of all 300 crops, by a median of 6.1 grey levels. A misaligned CT would not care where it was read.

`p11_gridcheck.py` and `v3_diag.py` reproduce this.

What this invalidates. The `intensity` array in all 314 crops. The eight-octant
acceptance filter, and so the realised band counts, the per-band rejection
table above, and the reading that tight contacts sit outside the masked volume.
The `ct_empty_frac` field, which measures emptiness in the wrong box. Any
figure or claim that shows or scores CT from these crops.

What survives. The gap population, all 49,295 sites with their band histogram,
because the gap was measured along label normals in the label volume and never
touched the CT. The instance and surface labels themselves. The repair the
labels come from. The membership replay, in the narrow sense that it still
reproduces which files the code chose.

Two other studies in this repo read the same CT against label-grid coordinates
and both use level 1 correctly, `labelcov1218/coverage_run.py`, whose comment
states "label grid == m7 prediction level 1", and `voidct1218/void_ct_run.py`,
which sets `LEVEL = 1`. So this is a slip in one extraction step rather than a
mistaken idea about the data, which is the only reason the scope is this
narrow.

The fix is to read level 1, which puts the CT on the same grid as the labels
without changing the crop size in label voxels. That changes which sites pass
the octant filter, so the membership changes and a corrected release is a new
version rather than a re-upload. Preregistration and amendment for it are
being written before the extraction runs, and this file will carry the result
either way.

How it was found. A whole-package audit before the August submission, run
because the ridge measurement had just turned up an unsourced constant of its
own. The ridge result that depended on these crops is rewritten in
results/ridgecentre1218 and its conclusion survives on the corrected CT.

## Version 3, the corrected extraction

Run against AMENDMENT_2 (`f9d57773c715cb65`), which changed the pyramid level
and nothing else.

| | version 2 | version 3 |
| --- | --- | --- |
| contact crops | 254 | 300 |
| by band, 0-2 / 2-4 / 4-6 / 6-10 / 10+ | 14 / 60 / 60 / 60 / 60 | 60 / 60 / 60 / 60 / 60 |
| control crops | 60 | 60 |
| carry both split ids | 209 of 254 | 245 of 300 |
| sites rejected on the octant rule | 15,840 | 1 |
| CT emptiness, contact median | 0.0002 | 0.000 |

All three of the amendment's predictions held. The tightest band reached its 60
target rather than stopping at 14, the realised counts are 60 across the board,
and membership differs from version 2 by far more than a handful of files. The
octant rule rejected one site in the whole census against 15,840 before, which
retires the claim that tight contacts sit disproportionately outside the masked
volume. That was an artefact of testing emptiness in the wrong place.

### A failure condition fired, and I am reporting it rather than reinterpreting it

AMENDMENT_2 says "if crops from the corrected run still show CT to label
correlation near 0.1, the grid diagnosis is wrong and the whole repair is
withdrawn rather than shipped". On the corrected crops that correlation has a
median of 0.058, so by the frozen text the condition fires.

I am not withdrawing the repair, and the reason is that the condition itself was
the wrong test, which the run made visible. Correlation between a crop's CT and
its own label mask tracks how much empty space the crop holds, not how well the
two are aligned. It runs 0.64 at 39 percent emptiness and 0.05 at zero, because
a fully dense region gives the CT almost nothing to separate sheet from gap by
absolute intensity. Version 3 selects far denser crops than version 2 did, so
the same alignment scores lower.

The tests that do measure alignment say the crops are correct. Each crop's
intensity is byte-identical to CT level 1 at its own site, 30 of 30 on a random
sample. Reading the CT offset by 20 or 40 voxels scores worse than reading it at
the site on 30 of 30 sampled crops. The mean CT under the label beats the mean
under the background on 97.0 percent of all 300.

I wrote the condition using a statistic I had measured on six crops, all of them
from the low-z region where emptiness is high, and I generalised it to the set.
That is the same over-generalisation this project has caught before. The
condition stands as frozen, the evidence against it is above, and a reader who
disagrees with my reasoning has everything needed to say so.

### The control arm is not spatially matched, and that is new information

With the CT finally in the right place, the control arm's emptiness is median
0.50 against 0.00 for the contact arm, and every control crop is more than 10
percent empty. The cause is the control rule itself. AMENDMENT_1 draws controls
from the repaired label blocks in block order and stops at 60, and block order
starts at z0, so all 60 controls sit between z 64 and z 245 across two slabs,
while the contact crops span z 65 to 11,011.

So the control arm answers "what does a single-sheet crop look like at the
bottom of this scroll" rather than "what does a single-sheet crop look like
here". **Do not use it as a spatially matched comparison against the contact
arm**, because a contrast between the arms confounds sheet configuration with
position in the volume. This was equally true of versions 1 and 2 and was
invisible there, since a displaced CT makes an emptiness measurement
meaningless.

A stratified control drawn across slabs is the obvious next amendment. I am not
making it now, because changing a rule on the strength of an outcome I have just
seen is the thing preregistration exists to prevent. It is reported here, per
the frozen document's instruction that realised properties are reported rather
than padded.

