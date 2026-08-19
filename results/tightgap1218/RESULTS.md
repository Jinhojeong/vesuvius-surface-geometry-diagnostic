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
it falls hardest at the tight end: of 120 sampled sites under 2 voxels, 117
fail it and 3 pass. The 0 to 2 band therefore ships 14 crops rather than 60.
The reading is that the tightest contacts sit disproportionately in crops that
reach outside the masked volume, which is a fact about where tight contact
occurs rather than only a sampling nuisance.

**Four in five contact crops verifiably contain both sheets.** Instance ids in
the source blocks are block-local, so a crop spanning blocks needs a
disambiguated id space; each contributing block is offset by a unique base and
the split pair is remapped through the base of the site's own block. With that
done, 213 of 260 contact crops contain both ids the split separated, 3 contain
one, and 44 contain neither. The 44 are a boundary effect, they sit on a median
of 4 contributing blocks against 2 for the crops that pass. A first pass
without the id disambiguation reported 231, and that number was wrong because
an id from a neighbouring block can collide with the pair; it is superseded
here.

**The octant rule is weak and the emptiness field is how to filter.** It only
asks for 1 percent CT per octant, so a crop with a large empty wedge can pass.
Measured across the contact arm, median emptiness is 0.0002, the ninetieth
percentile is 0.29, 20.4 percent of crops are more than a tenth empty and 11.2
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
