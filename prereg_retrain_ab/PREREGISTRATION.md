# Preregistration: does training on the v2.0-repaired labels change what a model learns at weld sites?

Matched-budget A/B fine-tuning on PHerc1218, v1 labels against v2.0 labels, seed
as the inferential unit. Committed in villa #193 comment 5226769475. This
document freezes the design before the first inferential seed runs. The
pre-freeze pilot permitted below is the only training that precedes it.

## 1. Question and scope

Whether fine-tuning on the v2.0-repaired instance labels measurably changes what
a model learns at fused-sheet weld sites, against fine-tuning on the v1 labels
with everything else held fixed. The claim is scoped to learned generalization
of the repair behaviour at held-out sites. It is not a claim about whole-scroll
segmentation correctness, because PHerc1218 has no human ground truth. A null
with a tight confidence interval is a committed deliverable equal to a positive.

## 2. Arms

Arm A trains on the v1 instance labels, the tree the v2.0 release declares as
its diff base, pinned by a manifest-of-manifests sha256 in frozen/hashes.json.
Arm B trains on the v2.0 labels (Kaggle release, MANIFEST.sha256 in the
release). The two trees share a
byte-identical binary mask and instance-id set and differ only in
voxel-to-instance assignment at repaired sites. Any mask-reducible training
target is therefore vacuous across the arms, which is why the target below is
the instance-contact boundary channel.

## 3. Training target and recipe

Single-channel boundary map, derived deterministically from each arm's labels.
A voxel is positive iff it is foreground and has a 6-connected foreground
neighbour carrying a different nonzero instance id. The recipe is the July
FT191 fine-tune unchanged. CE plus soft dice on this channel, 160³ crops,
batch 1 with gradient accumulation 2, AMP fp16, AdamW 1e-4 with cosine decay to
1e-6 over 6,000 steps, both arms initialized from the same frozen checkpoint
(ckpt_ft_full, sha256 in frozen/hashes.json). Seed s in arm A and seed s in
arm B share one seeded shuffle of one frozen crop list, so paired runs differ
only in which label tree the targets were derived from.

## 4. Patch sampling, with the intervention measured in it

One frozen list of 4,000 crops shared by both arms. 2,000 centred on SPLIT
sites of training tiles with a deterministic jitter of up to 16 voxels, 1,000
on ONE_SIDED sites, 1,000 at fixed grid positions. All jitter comes from md5 of
the site key, so there is no RNG anywhere in the data path. Training and
evaluation tiles are disjoint under frozen/tile_split.json, and the 200
eval-list tiles and every site inside them are excluded from all training
streams.

The manipulation check is precomputed over the whole list (frozen/flips.json).
Every one of the 2,000 split-anchored crops differs between the arms, median
28,526.5 differing target voxels per crop. The onesided-anchored crops differ at
a median of 22,301 voxels as well, because sites cluster spatially and those
crops contain neighbouring repaired material. The ONE_SIDED sites themselves
are untouched by the repair, so the control below is a site-level control, not
a crop-level one. 390 of the 1,000 uniform crops carry no differing voxel,
which is the measured sparsity of the repair. This measurement exists because a
July fine-tune with uniform sampling diluted its intervention to nothing and
moved the endpoint at most 0.005. That run is published as a null and is a
prior here, not a surprise waiting to happen.

## 5. Primary endpoint

Held-out-site separation rate at matched budget. Over the 5,000 frozen primary
sites, SPLIT sites on evaluation tiles, a site counts as separated when the
model's thresholded boundary prediction intersects the single-id thick run of
its flagging ray. The ray is re-fired with the shipped recast geometry on the
v1 mask, which is arm-independent because the mask is byte-identical. The
implementation compares the maximum predicted probability along the run points
to the threshold, which is equivalent and runs in one pass. One number per
seed per arm. The test statistic is the mean of the six per-seed paired deltas,
arm B minus arm A.

## 6. Budget rule and tie handling

For each of the twelve models independently, the threshold is the smallest
value on a 2,000-bin grid whose pooled predicted-positive fraction over the
scored foreground of the 191 site-carrying evaluation tiles is at or below b*.
b* is the v2.0 boundary-target fraction over those same voxels, recomputed
inside every scoring run by scripts/retrain_ab/instrument.py and equal to
0.0843838748625497 in the frozen oracle scores. Sites whose maximum run
probability lands exactly on the threshold count as not separated. No other
per-arm or per-seed tuning exists.

## 7. Instrument validation, measured before freezing

Oracle passes with each arm's own boundary target as the prediction:

| site list | n | v1 oracle | v2 oracle |
|---|---|---|---|
| primary | 5,000 | 0.1690 | 0.6142 |
| onesided control | 2,000 | 0.1285 | 0.1310 |
| background control | 2,265 | 0.1532 | 0.1567 |

The read-out has a floor of 0.169 on primary sites, from pre-existing
interfaces that lateral sheet contact places inside mask runs, and a ceiling of
0.614, so the repair signal spans 44.5 points where the planning MDE is 0.75.
The labels move the two negative controls by 0.25 and 0.35 points against 44.5
on primary, a site specificity of better than one hundred to one. An offset
control displaced 30 voxels in y was built, measured, and retired. 69 percent
of displaced sites leave the mask and the survivors respond to the repair at
+0.062 because displaced points land near other repaired material. Its numbers
stay in scores/ as a record of why it is not a gate.

## 8. Controls and strata

1. ONE_SIDED control, 2,000 held-out sites the solver examined and declined to
   cut. Covariates match the primary sites closely already, median thickness
   9.0 against 8.5 and ratio 2.5 against 2.5, and 1,833 one-to-one pairs exist
   within 0.5 pooled standard deviations (frozen/day5_strata.json). The gate
   uses the full control; the matched subset is reported alongside.
2. Background control, 2,265 on-mask points at least 40 voxels from every
   census site of their tile. Gates site specificity the way the retired
   offset control was meant to.
3. A-to-A empirical null, the 15 pairwise deltas among the six arm-A seeds,
   measuring the end-to-end harness noise at no extra GPU cost.
4. was_mega stratum. Only 15 of the 5,000 primary sites join to the v1.4
   agreement table, so by the predeclared rule (under 500) it is reported
   descriptively and carries no inference. A Dataset059 concordance slice is
   dropped because Dataset059 is a different scroll; there is nothing to
   intersect.

## 9. Secondary endpoint

Instance VOI with VOI_merge reported separately (unmerge-cli, villa PR #1301),
computed on connected components of foreground minus thresholded boundary,
against both v1 and v2 references, on the first 20 evaluation tiles in hash
order, seeds 40 to 42 only. This is the endpoint that does not descend from the
ray family that produced the repair.

## 10. Seeds, MDE, and the two predeclared adaptations

Six seeds per arm, seeds 40 to 45. Planning MDE 0.0075 absolute separation at
alpha 0.05 and 80 percent power, from the matched-budget across-seed sd of
0.0038 measured by TAUIL-Abd-Elilah on the m7 A/B, adopted with the seed as
the unit. That sd is borrowed from a different endpoint, so it is checked
twice. After the first two completed pairs, a paired sd above 0.008 triggers two
more seeds per arm, seeds 46 and 47, and nothing else. Before unblinding, the final MDE is recomputed from the 15 A-to-A
deltas and both values are reported. Per-site binomial noise at n = 5,000 is
about 0.007 per seed and is absorbed into seed-level variance.

## 11. Trainability gate, run pre-freeze

One pilot, arm B seed 40 at the full recipe, scored on 2,000 SPLIT sites
drawn from 60 hash-ordered training tiles, which land on 52 distinct tiles
(frozen/sites_trainsplit.json), never on evaluation sites. The v1 oracle floor on those sites is 0.1505
(scores/oracle-v1_ts.json), so the gate passes if pilot separation at matched
budget reaches 0.2005. The pilot ran before this freeze and scored 0.268
(scores/pilot_gate.json), so the gate is passed and the twelve seed runs
launch at this commit. One failure
permits one doubling to 12,000 steps. A second failure is published as
TARGET-UNLEARNABLE and no seeds run. The pilot seed is re-run inside the seed
grid afterward, so nothing from the pilot is reused as an inferential result.

## 12. Outcome buckets, all declared now

- POSITIVE, site-specific: |mean paired delta| at or above MDE, at least 5 of
  6 paired deltas sharing its sign, paired two-sided t under 0.05, and both
  negative-control deltas below MDE.
- SITE-UNSPECIFIC PRIOR: primary delta at or above MDE with sign stability,
  but a negative-control delta also at or above MDE. Reported as a
  placement-unspecific cut prior, with the site-specific claim withheld.
- NULL: |mean paired delta| under MDE with the A-to-A null width also under
  MDE. Published with the 95 percent CI.
- INCONCLUSIVE-UNDERPOWERED: delta at or above MDE without sign stability, at
  or above MDE with sign stability but t at or above 0.05, or under MDE with
  an A-to-A width above MDE. Published as the full seed table.
- TARGET-UNLEARNABLE: the gate in section 11 fails twice.
- UNINTERPRETABLE: frozen/flips.json fails its own recomputation from the
  target caches on a 50-crop spot re-check run before unblinding, or reads
  zero differing voxels on any split-anchored crop there, which would mean the
  arms did not differ in the data. No verdict, fault reported.

## 13. Priors disclosed

The July ignore-mask null (uniform sampling, at most 0.005 movement, published)
motivates the site-concentrated sampling and the manipulation check. The V20
oracle measurement (ground truth substituted at true merge sites moved the
official blend +0.0001 to +0.0026, control level) is why no mask-level or
blend-level endpoint appears anywhere in this design.

## 14. Stop rule

This document freezes at the commit that adds it, before seed 41 of either arm
exists. The only permitted adaptations afterward are the two in sections 10
and 11. All twelve runs complete unless a gate trips, the verdict is computed
once by section 12, and it is published whichever bucket it lands in.

## 15. Frozen artifacts

frozen/tile_split.json, frozen/crops.json, frozen/flips.json,
frozen/sites_primary.json, frozen/sites_onesided.json,
frozen/sites_background.json, frozen/sites_trainsplit.json,
frozen/sites_offset.json (retired, kept as record), frozen/day5_strata.json,
frozen/counts.json, frozen/ab_check.json, scores/oracle-v1.json,
scores/oracle-v2.json, scores/oracle-v1_ts.json, and frozen/hashes.json carrying sha256 of the init checkpoint, both label-tree
manifests, and every script in the path. The scripts ship under
scripts/retrain_ab/ in this repository.
