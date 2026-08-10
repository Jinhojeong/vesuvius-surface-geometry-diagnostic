# Results: the preregistered v1-vs-v2.0 weld-boundary retraining A/B

Verdict first. **NULL**, by the bucket definitions frozen in
[PREREGISTRATION.md](PREREGISTRATION.md) before the first inferential seed ran.
Fine-tuning on the v2.0-repaired labels does not change held-out weld
separation by 0.0075 or more at matched budget, against fine-tuning on the v1
labels with everything else held fixed. The 95 percent confidence interval on
the paired delta is [-0.0075, +0.0036]. The registered criterion is the
measured mean against the MDE with the A-to-A width under it, not the
confidence interval. The preregistration committed to
publishing this bucket exactly as readily as a positive, and this document is
that commitment kept.

## Primary endpoint

Separation rate at matched budget over the 5,000 frozen primary sites, one
number per seed per arm, seed as the inferential unit, delta is arm B minus
arm A.

| seed | arm A (v1) | arm B (v2.0) | delta |
|---|---|---|---|
| 40 | 0.2852 | 0.2772 | -0.0080 |
| 41 | 0.2844 | 0.2814 | -0.0030 |
| 42 | 0.2852 | 0.2836 | -0.0016 |
| 43 | 0.2736 | 0.2762 | +0.0026 |
| 44 | 0.2794 | 0.2848 | +0.0054 |
| 45 | 0.2806 | 0.2734 | -0.0072 |

Mean paired delta -0.0020, paired sd 0.0053, t -0.91, two-sided p 0.40, 4 of 6
deltas negative. The A-to-A empirical null, the 15 pairwise gaps among the six
arm-A seeds, has median width 0.0050, under the MDE, so the harness is quiet
enough to have seen the effect it was powered for. The MDE recomputed from the
measured A-arm sd is 0.00735 against the planning value 0.0075, so the
across-seed sd TAUIL-Abd-Elilah measured on a different endpoint came out
close here, 0.0045 against his 0.0038, and the predeclared two-pair sd check
passed at 0.0035 with no seed extension.

Both arms separate about 28 percent of held-out sites against the v1 oracle
floor of 16.9 percent, so both arms learned the cutting behaviour. This is a
no-difference null, not a no-learning null. The label trees differ at a median
of 28,526.5 boundary voxels per split-anchored training crop and the oracle
window between the two label states spans 44.5 points on this endpoint, so
there was signal to transfer and room to see it. It did not transfer.

## Controls

The two negative controls also drifted slightly negative, onesided -0.0080 and
background -0.0077 mean paired delta. The frozen buckets consult controls only
when the primary reaches MDE, which it did not, but the drift is worth stating
plainly. Arm B models separate marginally less everywhere, at flagged sites,
at solver-declined sites and at background sheet alike, with no site-specific
structure. The drift is consistent rather than noisy, negative in all six
seeds on the onesided control and five of six on background, and its size
sits at the MDE, so it reads as a real placement-unspecific offset. That is
the opposite shape of a site-specific repair effect.

## Secondary endpoint

Instance VOI on the first 20 evaluation tiles in hash order, seeds 40 to 42,
components of foreground minus thresholded boundary, VOI_merge reported
separately (the convention merged upstream in villa PR #1301), against both
label references.

Median voi_merge sits near 5.65 to 5.81 nats for every model, meaning the
thresholded boundary at the matched budget cuts few closed surfaces and the
components stay heavily merged relative to either reference. Within that, against
the v2 reference arm B models run slightly higher voi_split, 0.0205 against
0.0170 median, and lower voi_merge by 0.016, with the merge delta at 0.009
against the v1 reference, consistent in direction across all three seeds and
both references. At n = 3 seeds this is descriptive, and at 0.2 to 0.3
percent of the scale of the merge term it echoes the primary verdict rather
than contradicting it. The endpoint that does not
descend from the ray family agrees there is nothing here.

## Post-freeze amendments, all disclosed

Three engineering events happened after the freeze. The first is recorded in
frozen/hashes.json under amendments_post_freeze with original and amended
hashes. The other two changed no frozen script, so there is no hash to amend,
and they are disclosed here instead.

1. instrument.py crashed on edge evaluation tiles in model mode, because the
   CT cache is zero-padded to the standard block while labels keep the true
   tile shape. The fix crops the cached CT back to the label shape before
   prediction. Oracle numbers are unaffected, no threshold or endpoint or
   site definition changed.
2. The CT fetch was re-run twice to pull tiles the first pass had missed, for
   the background control and the trainsplit gate. The frozen fetch script
   already listed all five site files and was not changed.
3. The scoring crash plus a GPU memory collision between a scoring backfill
   and the training loop meant five trainings were re-run and the predeclared
   two-pair sd check executed while seed 43 was already training rather than
   strictly before it. The only permitted adaptation was adding seeds, which
   the check did not trigger, so no inference was affected. The concurrency
   mistake and its recovery are in the private log.

## What this does and does not say

The v2.0 release's geometry claims are direct measurements on the labels and
are untouched by this result. What this experiment tested is the further
hypothesis that training on the repaired labels changes what a model learns at
weld sites, and at this recipe, a 6,000-step fine-tune of ckpt_ft_full with
site-concentrated sampling, the answer is no at an MDE of 0.0075. A July
fine-tune with QC-derived ignore masks moved nothing either, so this is the
second preregistered-or-published null on label-side interventions reaching
model behaviour in this project. Scope limits are one recipe, one
architecture, one scroll and a site-local endpoint. From-scratch training,
other architectures or longer schedules are untested, and the result says
nothing about downstream tracing consumers of the labels, which never see
training at all.

## Artifacts

results/VERDICT.json, results/VOI_SECONDARY.json and the twelve per-model
score files sit beside this document. The frozen inputs, gates, oracle
validation and hashes are under frozen/ and scores/ as the preregistration
lists them, and every training, scoring, verdict and
secondary-endpoint script is under scripts/retrain_ab/ in this repository.
