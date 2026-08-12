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

(none)
