# Preregistration: seed-placement test of the GrowPatch hazard-weight effect

Date frozen: 2026-08-12 (runs begin only after this commit is on origin/main).
Author: Jinho Jeong.

## Background and why this design

The hazard-weight A/B on PHerc1218 (villa #191, comments 5200050557 and
5200680531) established, on eight large flagged clusters: area +234k vx2 per
seed against the no-weight control (p 0.003), +236k against a matched
random-field placebo (p 0.003), double-thickness minus 3.12 points against
placebo (p 0.0174), a dead-null placebo, and a kill control showing blanket
release buys area at catastrophic quality cost. The stated mechanism reading
is that placement is the effect. That inference so far rests on randomizing
the FIELD (placebo) and saturating it (kill). It has not been tested by
varying the SEED position, which is the direct form of the claim.

A sixteen-site extension (results/growpatch16) found the per-site effect
scales with flagged-cluster size (Spearman minus 0.67, p 0.0044) and named a
size-stratified design as the next step. That design is not runnable at
honest power: the fresh eligible pool (same eligibility as the original
selection, ratio >= 2.0, tile-interior margin 48, excluding the sixteen used
tiles) contains exactly ONE tile whose best cluster reaches 300 raw sites,
measured 2026-08-12. The original selection consumed the scroll's large
clusters. The placement test below is the confirmatory design the data can
support, and it tests the mechanism the thread actually asserted.

## Hypothesis

The hazard weight's effects require seeding on the flagged cluster. On the
same eight tiles, with the same whole-scroll weight volume active, seeds
placed away from the flagged clusters should show no hazard-weight effect.

## Design

Sites: the eight original demo sites (demo_sites.json, verbatim), each with
its archived on-cluster A and B cells.

New runs: for each site, one OFF-cluster seed in the same tile, two arms
(A = no weight, B = hazard weight, the shipped p1218_conf_v2_amp4.zarr),
five fresh replicates per arm per site. 80 new runs. Same GrowPatch
invocation, direction fields, voxelsize 8.64, MAXJOBS 3 discipline as the
powered block. No new direction fields are needed because off-seeds sit in
the same regions.

Off-seed rule (deterministic, no RNG): within the tile interior (tile-local
48 <= z <= 208, 48 <= y <= 464, 48 <= x <= 464), on a stride-16 lattice,
keep lattice points whose Chebyshev distance in L1 voxels is at least 96
from every census cluster centroid listed for that tile in census8k. Order
surviving points by the md5 hash of "slab/tile:off:z:y:x" and take the first
whose m7 prediction L1 value is >= 128 (same feasibility rule and reader as
select_sites16.py). If no point survives at 96, relax the distance to 64,
then 48, recording the relaxation. Selection script:
scripts/placement_ab/select_offseeds.py in this repo; its emitted seed list
is part of the frozen record.

Archived cells reused: the on-cluster A and B per-run scores are the powered
block's, pinned by sha256
fd9deb2c3b7a3e3f0e5e55b258a5c3230b3ce845cb66e7d24425596472a65876
(powered_scores.json, as shipped in the Kaggle hazard dataset). They predate
this hypothesis's off-arm and are declared as reused, not rerun. The novel
falsifiable content of this study is the off-cluster data and the contrast.

## Endpoints and analysis

Unit of analysis: the site (n = 8). Per site, each cell (arm x placement) is
summarized as the mean over its five replicate runs of the run-level metric,
computed by the unmodified score_ab.py scorer. The analysis script must
first reproduce the archived powered contrasts from the archived JSON to
full float precision before touching new data, as powered16_analyze.py did.

Per site define the double difference DD = (B minus A at the on seed) minus
(B minus A at the off seed).

- PRIMARY: area_vx2. Prediction DD > 0. Paired one-sided t across the eight
  sites, alpha 0.05.
- SECONDARY (declared underpowered): frac_double. Prediction DD < 0, same
  test. Power at the observed on-cell effect (minus 2.65 points, per-site sd
  3.47) is roughly 55 to 60 percent, so this endpoint is reported with its
  interval and carries no confirmatory weight on its own.
- GUARDRAIL: on_sheet_rate DD, two-sided at 0.05. A significant guardrail
  difference flags all endpoints as possibly denominator-driven.
- Sensitivity: leave-one-site-out p range for the primary; the s2/s4
  same-tile pair drops reported explicitly.
- Also reported, not tested: the off-seed (B minus A) itself with its CI,
  which the mechanism predicts is near zero.

Power for the primary: archived per-site on-cell (B minus A) area has mean
+234k and sd 143k. Allowing the off cell comparable sampling spread, DD sd
is at most about 200k, SE at n 8 about 71k, one-sided MDE80 about 198k,
below the +234k point prediction, so power exceeds 80 percent if the
mechanism claim is true.

## Outcome buckets (frozen)

1. Primary confirms and secondary agrees in sign with p < 0.05: placement
   confirmed on both endpoints.
2. Primary confirms, secondary inconclusive: placement confirmed for growth;
   the quality claim stays at its current placebo-level evidence.
3. Primary fails with the off-seed (B minus A) area comparable to on-seed
   (point estimate above half the on-seed effect): placement REFUTED; the
   published mechanism reading gets a public correction in the same thread.
4. Primary fails with wide intervals and off-seed estimate near zero:
   underpowered no-verdict, disclosed as such.
Any guardrail flag is reported alongside whichever bucket obtains.

## Exclusions and amendments

A run that crashes (rc nonzero) is retried once with the same configuration;
a site that cannot complete its cells is dropped with disclosure. No other
exclusions. Any post-freeze change to scripts or parameters is recorded in
an AMENDMENTS section appended to this file with file hashes before and
after, following the retrain A/B practice.

## AMENDMENTS

### Amendment 1 (2026-08-12, before any run was launched)

The design section stated "No new direction fields are needed because
off-seeds sit in the same regions." That premise was false and the mandatory
pre-launch coverage check caught it: the direction fields were built with a
256-voxel write window around the ORIGINAL seeds, not around the tiles, and
five of the eight frozen off-seeds fall outside that window (Chebyshev
distances to their original seeds 292 to 411; measurement in
demo_out/placement/field_extent_check.json). Zero of the 80 runs had
started.

Change: build direction-field regions with the same builder and parameters
(make_direction_field.py, radius 256, band 6.0, sigma 1.0) centered on the
seven unique off-seed positions before running. The frozen off-seed rule,
the emitted seed list (offseeds_placement.json, sha256
b07e905b6dee57c7cc2fec5a755cb099a2357d887933d5c8bbad56ec481290dc), the
arms, endpoints, analysis and buckets are unchanged. This choice preserves
the off-seed semantics (distance from clusters up to 411 voxels) and the
on/off symmetry, since every on-seed likewise grows in a field centered on
itself. The alternative, reselecting off-seeds inside the existing windows,
would have confined them to a 96-to-256 ring around the clusters where a
growing patch reaches the cluster mid-run and blurs the off condition.

Two disclosures recorded with this amendment. First, sites s2 and s4 share a
tile and the deterministic rule therefore gave them the identical off-seed
point; their off cells are independent replicate sets at the same
configuration, and the leave-one-site-out table covers their coupling.
Second, PREREGISTRATION.md sha256 before this amendment:
7dd5c387410f9af35980af5d3518ccb21369064f4ae0dd76bc36780bbbb6c967.

### Amendment 2 (2026-08-12, before any run was launched)

After the Amendment 1 field builds, the mandatory voxel check failed on
four of eight sites (o0, o2, o3, o4): the exact off-seed voxel reads
unwritten. Measured cause: the field writes a band of six voxels around the
prediction surface level set, and these points sit deep inside thick
predicted mass, so the m7 value at the voxel is 255 while the band does not
reach it. Every original on-seed reads a written vector, so an unwritten
off-seed would break the on/off comparability condition, not just a
convenience.

Change, the nudge rule: an off-seed whose voxel is unwritten moves to the
field-written voxel at minimal Chebyshev distance from the frozen point,
preferring voxels with m7 prediction value at or above 128 within that
minimal shell, ties broken by the md5 hash of "nudge:z:y:x". The 96-voxel
cluster-distance rule is re-verified at the nudged point and the run stops
if it fails. Measured nudges: o0 moves 1 voxel, o3 moves 10, o2 and o4
(shared point) move 26; the smallest re-verified cluster distance after
nudging cannot fall below 109 minus 1, 134 minus 10, and 145 minus 26, all
above 96. Arms, endpoints, analysis and buckets remain unchanged. The
frozen and nudged seed lists both ship in the results record.

PREREGISTRATION.md sha256 before this amendment:
e48811c04c4796948ce154555063fa75d651bee97d1b312bad06eafba05137f7.

### Amendment 3 (2026-08-12, after the runs, before the amended analysis ran)

The scored off cells carry essentially zero on-sheet rays, 5 of about
16,000 sampled, and the measured cause is not the scorer. Every off-seed
sits in label-void space. Occupancy of the repaired label volume in a
49-cube around each of the eight off-seeds is 0.000, the label at every
off-seed voxel is 0, and 100 percent of sampled off-run quads have no label
within 12 voxels, against 40 percent within 4 voxels for the matched on
run. The m7 prediction carries surface at all eight points, so these are
regions where the prediction exists and the labels do not, the
label-incompleteness population my July QC measured. The census clusters
that anchor the on-seeds are detected in the labels, so the frozen 96-voxel
exclusion systematically pushed off-seeds into unlabeled space. One further
measured fact bounds the off arms, the hazard volume reads full confidence
at every off-seed and everywhere within 64 voxels of one, so A and B
present identical local inputs at the off cells.

Change: label-referenced endpoints, frac_double, the on-sheet guardrail,
per-thousand and thickness ratios, are unmeasurable at cells with zero
on-sheet rays and are reported with raw counters and no test. The primary
endpoint, area, is unaffected and its analysis runs as frozen. The
secondary and guardrail entries in the outcome buckets therefore read
through the primary alone, and this limitation ships in the results record.

PREREGISTRATION.md sha256 before this amendment:
f4be9d78ddd126375397673efecb0e9c1ae48c07bccb8da51e810843d5390c2c.
