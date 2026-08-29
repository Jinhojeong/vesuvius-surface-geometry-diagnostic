# Cross-check of the recovery claims against the pre-repair tree

Requested on PR #1. Two input files from
`Diego-dcv/vesuvius-topological-grid/archives/results/recovery`,
`unlabeled_material_1218.csv` (9,034 claims of papyrus-bright material with no
label in the pre-repair tree) and `wedge_crossings_1218.csv` (3,625 crossings
beyond column truncation, 281 at conf 2). Every point was read at its exact
rounded voxel in the full pre-repair tree, and again in the v2 repaired
blocks.

## Results, in xcheck_result.json

The v2 repair changes label support nowhere. Over the 952 blocks the query
points touch, the repaired and pre-repair support masks differ at zero voxels,
so the two requested checks have one answer.

2,418 of the 9,034 unlabeled-material claims carry a pre-repair label at the
exact voxel. The labelled fraction of the surrounding 5-cube has a median of
0.42 over those hits, and only 27 sit in fully labelled neighbourhoods.
Almost all are therefore edge grazes within the stated 32 to 60 micron
prediction tolerance, and 27 are claims of no sheet where a labelled sheet
exists. For
the wedge crossings the rates are 93 of 3,344 at conf 1 and 7 of 281 at
conf 2.

`xcheck_violations.csv` lists every exact-voxel hit, sorted so the fully
labelled neighbourhoods come first, with the pre-repair global id and the
5-cube labelled fraction.

## The raw-CT glance, in ctglance_result.json

Bright-run thickness at the two theta-0 direct-search loci against pooled
controls at 15 to 30 degrees, over 17 to 21 z planes each. L 1.569 m reads a
median 5.0 voxels against 4.0 pooled control, L 2.838 m reads 4.0 against
4.0. No doubling at either locus. The estimator is a 40th-percentile
threshold on the radial profile with the contiguous bright run nearest the
locus radius, which is crude, and the sampling is 21 planes, so this is a
glance rather than a search.
