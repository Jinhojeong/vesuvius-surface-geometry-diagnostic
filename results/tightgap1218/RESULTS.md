# Tight-contact validation set on PHerc1218: what was built and what it costs

Executed against PREREGISTRATION.md (`9d5d71dbaf45ab85`) and AMENDMENT_1.md
(`a07c2f86f39b3bd8`). Every count below is realised, not targeted.

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
rejection rate is a rate over everything available rather than over a sample,
and it is the honest form of the claim the wrong number was reaching for. The
reading does not change. The tightest contacts sit disproportionately in crops
that reach outside the masked volume.

AMENDMENT_1.md is left as written. It is hash chained at `a07c2f86f39b3bd8` and
referenced by `control_summary.json` and `p11_control2.py`, so it is corrected
here rather than edited. `p11_bandreplay.py` reproduces the table.

