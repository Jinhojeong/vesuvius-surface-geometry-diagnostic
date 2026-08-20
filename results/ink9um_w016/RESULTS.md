# What inference-time choices are worth on the released 9 um ink checkpoints

Executed against PREREGISTRATION.md (`8ceda1e7a0536613`). Scoring happens only
inside the `validation_mask` of `pherc0139-w016`, 178,146 voxels at one
annotated depth, 23.16 percent ink positive. w016 is one of the three
validation cases the label README says the released checkpoints report on.

## The bar, set before any ensemble was scored

The two seeds at matched step differ by a median of 0.033 F1 and 0.030 AUC.
That spread is the bar the frozen document fixed: an arm beating the best
single checkpoint by less than it is reported as indistinguishable.

## Arm 1, the fourteen released checkpoints alone

F1 ranges 0.611 to 0.733 and AUC 0.845 to 0.922, so the mask does separate
checkpoints and the preregistered "too easy to separate anything" failure
condition does not fire. Best single is seed 42 at step 20,000, F1 0.7325,
AUC 0.9179, at threshold 0.46.

The best step is not the last one and it differs by seed, 20,000 for seed 42
against 75,000 for seed 43. The tutorial's advice to try a few checkpoints is
measurably right rather than merely cautious.

## Arms 2 to 4, ensembling

| arm | F1 | AUC | vs best single |
| --- | --- | --- | --- |
| all fourteen averaged | 0.7641 | 0.9328 | +0.0316 F1, +0.0150 AUC |
| step ensemble, seed 43 | 0.7487 | 0.9307 | +0.0163 |
| six late checkpoints | 0.7442 | 0.9260 | +0.0117 |
| seed pair at step 75,000 | 0.7436 | 0.9252 | +0.0112 |
| seed pair at step 10,000 | 0.6456 | 0.8854 | −0.0869 |

**By the frozen bar this is indistinguishable, not a gain.** The full ensemble's
+0.0316 F1 sits just under the 0.033 seed spread, and its +0.0150 AUC is half
the 0.030 AUC spread. On w016 alone a rank statement also held, the fourteen-model
average beating every one of its members on both metrics. That statement did not
survive the other two validation segments, see below.

Ensembling only helps when the members are decent. Averaging the two 10,000
step checkpoints lands below either seed's best single run.

## Arm 5, depth window, and two corrections to how it was run

Two things the logs showed after the first pass, both of which changed the
question. The checkpoint fixes its own window at 17 of the 21 slices, indices 2
to 18, so a first setting written as "17 centred" was byte-identical to the
default and measured nothing. And `--layer-start` with `--layer-end` cannot
widen past that window, a request for all 21 slices still selects 2 to 18. The
flag narrows or shifts, it does not widen.

| window | F1 | AUC |
| --- | --- | --- |
| default, 17 slices at 2 to 18 | 0.7325 | 0.9179 |
| 17 shifted up, 0 to 16 | 0.7085 | 0.9106 |
| 17 shifted down, 4 to 20 | 0.6938 | 0.9064 |
| 13 centred | 0.4963 | 0.7594 |
| 9 centred | 0.4424 | 0.6875 |
| 5 centred | 0.4163 | 0.6254 |

The default window is the best of the six. Shifting it costs a little, 0.024
and 0.039 F1, and narrowing it costs a lot and monotonically. So depth window
is a constraint on these checkpoints rather than a knob offering gains.

## Arm 6, direction, where the useful finding is

`--direction both` writes two files rather than fusing them, forward to the
named output and reverse alongside it. Scored separately:

| direction | F1 | AUC |
| --- | --- | --- |
| forward | 0.7325 | 0.9179 |
| reverse | 0.3761 | 0.5383 |
| the two averaged | 0.6813 | 0.8845 |
| elementwise max | 0.6819 | 0.8877 |

Reverse alone is close to chance here, so surface orientation dominates
everything else measured in this document. And fusing the two directions costs
0.051 F1 against just taking forward, because averaging a near-chance map into
a good one dilutes it. Run both directions to find out which way the surface
faces, then use that one. Fusing them blind is a five point penalty.

## What this adds up to

No inference-time setting tested here beats the best single released checkpoint
by the preregistered bar. Ensembling is the only arm that improves anything at
all and it lands just under the bar, while the two settings that look like
tunables can only lose. The practical guidance that follows is cheap and
concrete: try several checkpoints because the best one differs by seed, leave
the depth window alone, and determine direction rather than fusing it.

## Limits

w016 is in the training corpus, so this is held-out region behaviour inside a
trained segment, not transfer to an unseen scroll. One region, one segment, one
annotated depth. Threshold selection on the scoring mask inflates F1 for every
arm equally, which is why the bar is the seed spread rather than zero. The input
is my re-implementation of the documented preprocessing, the public 2.399 um
volume at XY level 2, centred 84 Z planes, 4x Z mean pooling, so a preprocessing
difference would move all arms together.

## Replication on the other two validation segments

The label README names three segments carrying a `validation_mask`, the cases
the released checkpoints report on. Running the same fourteen checkpoints and
the same fourteen-model average on all three, each scored against its own seed
spread:

| segment | ink rate | best single | ensemble | difference | that segment's bar | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| pherc0139-w016 | 0.232 | 0.7325, seed 42 step 20k | 0.7641 | +0.0316 | 0.0330 | indistinguishable |
| pherc0814-46527 | 0.373 | 0.7627, seed 43 step 10k | 0.7651 | +0.0025 | 0.0084 | indistinguishable |
| pherc1667-w029 | 0.230 | 0.7472, seed 43 step 30k | 0.7391 | −0.0080 | 0.0365 | indistinguishable |

**Three of three say indistinguishable, so the null is the finding.** The
ensemble never clears the bar, and on PHerc1667 it is slightly worse than the
best single checkpoint.

**One claim from the single-segment run does not replicate.** "The average beats
every one of its fourteen members" holds on F1 for two of three segments and
fails on PHerc1667. On AUC it holds three of three. So the threshold-free rank
statement is the robust one and the F1 version was a one-segment observation
that should not have been stated as though it were general.

**The guidance that gets stronger is checkpoint choice.** The best checkpoint is
different on all three segments, seed 42 at step 20k, seed 43 at step 10k, and
seed 43 at step 30k, spanning both seeds and three widely separated steps. The
bar itself is segment-dependent too, 0.008 on PHerc0814 against 0.037 on
PHerc1667, so seed-to-seed stability is not a fixed property of these weights.

Depth window and direction were measured on w016 only and are not re-run here.
