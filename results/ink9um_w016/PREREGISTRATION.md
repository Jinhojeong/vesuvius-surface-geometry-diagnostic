# Preregistration: how much do the released 9 um ink checkpoints gain from
# inference-time choices alone?

Frozen 2026-08-19, before any scoring run. Written because the answer is a
comparison between settings, and the temptation once numbers exist is to pick
the setting that wins and call it the finding.

## Why

The team released cross-scroll 9 um ink checkpoints and named the open
directions itself, better augmentation, longer training, other architectures,
ensembling, and more data. Everything on that list except ensembling needs
training runs I cannot afford on one 11 GB GPU. Ensembling and the other
inference-time knobs are affordable, and nobody has published what they are
worth on the released weights. That is the question here, and it is a question
about a core pipeline stage rather than about a diagnostic.

## Data

`pherc0139-w016` from the published `ink_9um` label set, one of exactly three
segments carrying a `validation_mask`, which the label README states are the
online-validation cases the released checkpoints report metrics on. Scoring
happens only inside that mask, 178,146 voxels at a single annotated depth,
23.16 percent ink positive.

Input is built by the documented recipe, the public 2.399 um volume of segment
`20250108000004-w029_2025010827` read at XY pyramid level 2, centred 84 Z
planes, 4x Z mean pooling, giving 21 slices on the label grid. A 128 voxel
margin is added around the mask bounding box so sliding-window inference has
context. No normalisation is applied by me; each checkpoint carries its own.

## Primary endpoint, fixed here

F1 on validation-mask voxels at the threshold that maximises F1 for that run,
reported together with the threshold. Secondary endpoints are AUC, which is
threshold free, and F1 at a fixed 0.5 after the tutorial's display rescaling.
The primary comparison is the best single released checkpoint against the best
ensemble, both selected on the same mask.

## Arms

1. Every released checkpoint alone, two seeds by seven steps, defaults
   otherwise, `--direction both`.
2. Seed ensemble, mean probability of seed 42 and seed 43 at matched step.
3. Step ensemble, mean probability across steps within a seed.
4. Full ensemble, mean across all fourteen.
5. Depth window, `--layer-start` and `--layer-end` swept over the 21 slices in
   five windows, on the single best checkpoint from arm 1.
6. Direction, forward only against reverse only against both, on the same
   checkpoint.

## What counts as a real gain

The mask is one region, so run-to-run variation cannot be estimated by
resampling it. The bar is therefore a margin larger than the spread between
the two seeds at the same step, measured in arm 1 and reported before any
ensemble number is quoted. An ensemble that beats the best single checkpoint by
less than that spread is reported as indistinguishable, not as an improvement.

## Declared limits, written before the numbers exist

w016 is in the training corpus, so this measures held-out region behaviour
inside a trained segment and not transfer to an unseen scroll. One region on
one segment. Threshold selection on the same mask that scores the run inflates
F1 for every arm equally, which is why the seed spread and not zero is the bar.
The recipe is my re-implementation of the documented preprocessing rather than
the team's script run end to end, so a preprocessing difference would shift all
arms together.

## What would make this wrong

If arm 1 shows the fourteen checkpoints scoring within noise of each other,
then the mask is too easy or too small to separate anything and no ensemble
claim follows. If the best F1 lands near the 23.16 percent positive rate's
trivial floor, the run is measuring calling rate rather than ink and this
document says so in the results instead of the set being reissued.
