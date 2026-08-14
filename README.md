# Geometry diagnostics and label repair for Vesuvius surface models

Measurements, preregistered experiments and repaired labels for Herculaneum
scroll surfaces, all on real scroll data. The work runs in the open on
[ScrollPrize/villa#191](https://github.com/ScrollPrize/villa/issues/191),
[#193](https://github.com/ScrollPrize/villa/issues/193) and
[#1114](https://github.com/ScrollPrize/villa/issues/1114), and every number
below has a commit and a script behind it.

## August 2026 in one page

Pictures first, since several of these are easier to see than to read. The
full gallery with sources and selection rules is in
[docs/evidence_gallery](docs/evidence_gallery/INDEX.md).

| what | where |
|---|---|
| The census as a GrowPatch weight volume: same seed, with and without, traced and flattened | [figure](docs/evidence_gallery/render_1218_hazard_ab.png) · [data](results/growpatch16) |
| PHerc1218 fused-sheet repair, before and after on three validated sites | [figure](docs/evidence_gallery/repair_1218_before_after.png) |
| Surfaces traced on PHerc1218 before and after CT-support masking, flattened through the standard renderer | [figure](docs/evidence_gallery/render_1218_support_mask.png) · [data](results/p8_support_mask) |
| Labels cover 45.2 percent of the m7 predicted surface, and the void is structured | [figure](docs/evidence_gallery/labelcov_1218_map.png) · [data](results/labelcov1218) |
| Nearly half of m7 positives sit beyond masked-CT support (villa #1114) | [figure](docs/evidence_gallery/m7_beyond_support.png) · [data](results/voidct1218) |
| PHerc0332 sheet separation on real CT, second scroll | [figure](docs/evidence_gallery/separation_0332_slices.png) · [data](results/repair0332) |
| PHerc1203 thickness zone map, regenerated after a self-caught read bug | [figure](docs/evidence_gallery/thickness1203_v2_map.png) · [data](results/thickness1203) |
| Preregistered placement A/B, off-seeds sit in label-void space | [figure](docs/evidence_gallery/placement_offseed_context.png) · [data](prereg_placement_ab) |

Two experiments were frozen in public before they ran, and both are published
as they came out. [prereg_retrain_ab](prereg_retrain_ab) asked whether
fine-tuning on repaired labels moves held-out weld separation and returned a
null with a tight interval. [prereg_placement_ab](prereg_placement_ab) asked
whether the hazard field's growth effect requires seeding on the flagged
cluster and returned a clean positive, plus three amendments recorded with the
hash of the document each amended.

## The original diagnostic (July 2026)

Quantifies **where and how much** surface-detection models degrade in curved
and compressed regions, using the public
`Dataset059_s1_s4_s5_patches_frangiedt` ground truth (200 random patches,
seed-fixed) and two public baselines. Posted to
[villa#191](https://github.com/ScrollPrize/villa/issues/191#issuecomment-4930879524).

## Headline results

Per-geometry-bin **AUC** (P(surface) at GT surface points vs near-sheet
background points 2–10 vox away — firing-rate-free, ~6M+6M points per model):

**Curvature** (orientation-structure-tensor dispersion, decile bins) — both
models degrade monotonically through all 10 deciles:

| model | flattest decile | most curved decile |
|---|---|---|
| `surface_recto_059_redo` | 0.953 | **0.655** |
| `surface_recto` (nnUNetv2) | 0.970 | **0.846** |

**Compression** (distance to next sheet along the local normal):

| model | <4 vox | >22 vox |
|---|---|---|
| `surface_recto_059_redo` | **0.413** (below chance) | 0.886 |
| `surface_recto` (nnUNetv2) | 0.771 | 0.943 |

Two additional observations. Recall-only metrics flatter high-firing-rate
models in these hard regions (they made `surface_recto` look immune to
compression, which it isn't). And in the lowest-dice patches the predictions
visibly track fine CT sheet structure that the sparse GT bands miss, so part of
the low tail is plausibly **label incompleteness** (see the CT-brightness
check in `results/oracle.csv`).

One limit to keep in mind: patches likely overlap both models' training distribution, so the
curves are optimistic bounds. This is not the official Kaggle Surface
Detection metric; it characterizes the failure mode rather than closing #191.
Geometry proxies derive from GT, so label quality confounds the strata. The
<6 vox spacing bins are thin (~2k background points).

## Follow-up: label-completeness QC (`scripts/qc192_labels.py`)

Turns the label-incompleteness hypothesis above into a measurement. For each
patch it detects **candidate unlabeled sheets** — high-confidence predictions
>3 vox outside the GT that are sheet-shaped (planar PCA) and CT-bright — and
scores the patch's label incompleteness.

- Mean incompleteness across the 200 patches: **6.4%** of sheet voxels;
  the worst decile of patches exceeds **14%**.
- **Cross-model validation:** the score computed from model A's predictions
  correlates with the *independent* model B's per-patch dice at
  **Spearman −0.45 / Pearson −0.49**, so label incompleteness is a real,
  model-independent contributor to low evaluation scores (~24% of variance).
  The same-model correlation is −0.91 but partially circular (the score counts
  a subset of that model's own false positives), so the cross-model number is
  the one to use.
- Visual check: in the worst patches the detected candidates are long coherent
  sheets clearly present in the CT and absent from the GT
  (`results/label_qc/worst_incomplete.png`, red = candidate, green = GT).

## Follow-up 2: why compression fails, and a repair that costs nothing

I traced the compressed-region collapse of `surface_recto_059_redo` to its
cause by elimination, then tested a fix against the official metric:

- Not information limits: a model-free phase-demodulation score (local period
  by autocorrelation along the sheet normal, phase-locked matched filter) sits
  at chance (0.50-0.53) in the tight-spacing bins where a fine-tuned network
  reaches 0.65 (`scripts/phase_prior.py`, `results/finetune/phase_auc.json`).
- Not label incompleteness: fine-tuning with QC-derived ignore masks changes
  nothing vs an identical run without them (max +0.005 in compression bins;
  `scripts/qc_make_ignore.py`, `scripts/ft059.py`).
- Not the loss: fine-tuning with the original combined loss (DC+SkelREC+CE)
  reproduces the repair within ±0.003 of a plain CE+dice arm
  (`scripts/arm_msr.py`, `results/finetune/msr_auc.json`).

What fixes it is continued training itself: 3000 steps on 500 fresh Dataset059
patches (~50 min, one 11GB 2080 Ti) lifts <4 vox spacing AUC from **0.41 to 0.65**,
with the whole degradation curve shifting up. Cost on the official
`topometrics` leaderboard blend: none detectable — paired over n=123 patches,
delta −0.003 ± 0.004 (t p=0.43, Wilcoxon p=0.80), fine-tuned model ahead on
68/123. An n=40 sample had suggested −0.01; it did not survive the larger
sample. The operating threshold shifts from 0.4 to 0.6.

Scaling exposure further (full 1554-patch pool, 6000 steps) gives the best
numbers in 7 of 8 spacing bins and surface dice 0.850, but the compressed
bin saturates around **0.65** — that plateau looks like the real remaining
ceiling for this recipe (`results/finetune/full_auc.json`,
`full_official.json`).

That ceiling looks like a resolution limit, not a thresholding gap a
post-process could clear. Sampling the predicted probability along the sheet
normal at 200 compressed GT points (neighbour sheet <4 vox away), the released
model shows a single broad peak covering both sheets in **78%** of them and two
resolved peaks a splitter could separate in only **0.5%**. Fine-tuning does not
change that shape (still **83%** one peak, 0.5% two peaks); it helps this bin by
detecting sheets it used to miss, not by separating merged ones. A direct
normal-direction split of thick predictions drops the topology score without
helping, which fits a field with one mode to split. So the AUC can reach 0.65
while the topology stays wrong (one connected component where GT has two), and
getting past it needs resolution rather than post-processing
(`scripts/bimodal191.py`, `scripts/peel191.py`,
`results/finetune/bimodal_taxonomy.json`).

**A second diagnostic axis fell out of this work:** discrimination is also
systematically worse where the sheet normal is aligned with a volume axis.
Baseline AUC is 0.90 for oblique normals vs 0.80 for axis-aligned ones, and
the gap survives controlling for curvature (low-curvature half: 0.955 vs
0.853) and for compression (persists within every spacing band, −0.05 to
−0.09). A plausible mechanism is partial-volume stair-stepping for
grid-parallel sheets. Fine-tuning narrows this gap too (0.80 → 0.87).
Tables: `results/finetune/orient_strata.json`, `orient_x_spacing.json`.

Repaired checkpoints (`ckpt_ft_full.pth` recommended; `ckpt_ft_ctrl.pth` and
`ckpt_ft_msr.pth` for the ablation arms) are in the GitHub Release. Known
issue: the official metric library segfaults on one of the s4 evaluation
patches (C++ Betti matching); reproducer available.

## Follow-up 3: does the repair transfer to an unseen scroll?

The repair above was fine-tuned and evaluated on the same dataset, so it shows
more exposure repairs this distribution, not that it generalizes. I checked that
directly with a leave-one-scroll-out test: fine-tune the released 059 on Scroll 1
patches only (same recipe, 500-patch pool, 3000 steps, no s4 or s5 seen), then
evaluate on Scroll 4, a scroll the fine-tune never touched.

Two metrics move in opposite directions.

Per-voxel discrimination transfers strongly. On Scroll 4 the tight-spacing bin
(<4 vox) rises from about 0.52 to about 0.89 AUC, and every wider bin improves
too. A model that only saw Scroll 1 reads Scroll 4's compressed geometry much
better by this measure.

The official topometrics blend regresses. Sweeping the threshold for both arms
and comparing each at its own best, base scores 0.693 (at 0.4) and the fine-tune
0.624 (at 0.5), a drop of about 0.067, driven by topology (0.414 to 0.301). It is
not a threshold artifact: the fine-tune is below base's best at every threshold.

The two disagree because of the resolution limit from Follow-up 2. At <4 vox the
two sheets already read as one probability peak. The bimodality diagnostic on
Scroll 4 shows the fine-tune converts missed sheets into detected ones (no-peak
share 15.8% to 7.6%) but detects them as merged blobs, not resolved pairs
(bimodal stays at 1.1%, single-peak rises 82.1% to 88.6%). Both ground-truth
sheets then get high probability, so per-voxel AUC rises, while the merged blob is
one connected component where there should be two, so topology falls.

This replicates across three fine-tune seeds (1234, 2025, 42): <4 vox AUC lands
at 0.888 / 0.896 / 0.892 and the fair official delta at -0.069 / -0.067 / -0.066,
so it is not a single-run effect. Scope is one architecture (059) and one held-out
scroll (Scroll 4); model- and multi-scroll generality is not tested here.

Takeaway: per-voxel AUC is misleading in compressed geometry, and a cross-scroll
fine-tune that improves it can quietly degrade the topology-aware quality the
official metric measures. Judge cross-scroll transfer with geometry-stratified,
topology-aware metrics, not AUC alone.

Note on baselines: Scroll 4 is intrinsically easier than the Scroll 1 dominated
mix used earlier (base <4 vox AUC 0.52 vs 0.41), so the 0.65 ceiling is a Scroll 1
property, not a universal constant, and the 0.52 to 0.89 jump should be read
against this easier baseline.

Scripts: `scripts/loso_ft.py` (Scroll-1-only fine-tune, `SEED` env),
`scripts/loso_eval.py` (per-bin AUC), `scripts/loso_sweep.py` (fair threshold
sweep), `scripts/loso_seed_val.py` (compact per-seed validation),
`scripts/bimodal_loso.py` (the merged-blob check on Scroll 4). Results in
`results/finetune/loso_*.json`.

## Contents

- `scripts/diag2_191.py` — recall stratification: sliding-window inference +
  GT geometry (structure-tensor curvature, along-normal sheet spacing)
- `scripts/diag4_auc.py` — the two-model per-bin AUC comparison (main result)
- `scripts/diag3_d058.py` — nnUNetv2 baseline runner (trainer-fallback recipe)
- `scripts/loader059.py` — loading recipe for `Model_epoch499.pth`
  (`_orig_mod.` prefix + duplicate per-task encoder keys + `separate_decoders`)
- `scripts/oracle191.py` — CT-brightness label-incompleteness probe
- `scripts/bimodal191.py` — along-normal peak analysis of the compressed-bin
  ceiling (released vs fine-tuned), the resolution-limit vs post-process test
- `scripts/peel191.py` — normal-direction split attempt and its official-metric
  cost (the post-process that does not help)
- `scripts/loso_ft.py` — leave-one-scroll-out fine-tune (Scroll 1 only, `SEED` env)
- `scripts/loso_eval.py`, `scripts/loso_sweep.py`, `scripts/loso_seed_val.py` —
  cross-scroll transfer AUC, fair threshold sweep, and per-seed replication
- `scripts/bimodal_loso.py` — merged-blob mechanism check on the unseen scroll
- `results/` — per-patch CSVs, strata tables (JSON), figures

## Reproduce

```bash
# data: 200-patch random sample (seed-fixed) of Dataset059 from dl.ash2txt.org
# models: huggingface.co/scrollprize/surface_recto_059_redo and /surface_recto
pip install torch tifffile scipy nnunetv2 vesuvius
python scripts/diag2_191.py   # recall strata (per model)
python scripts/diag4_auc.py   # two-model AUC strata
```

Paths at the top of each script point at local data/model directories — adjust
to your layout. GPU: everything runs on a single 11 GB RTX 2080 Ti.

## Reproduce the geometric-channel check

The measurement most worth re-running is the one that kills the obvious next
idea. Feeding a model geometry instead of intensity looks promising until the
control shows the signal is band thickness:

```bash
# eight Dataset059 crops of your own; see the note below on patches-v1
python scripts/reproduce_geometry_check.py --data path/to/patches
```

Expected, on eight patches (n around 674):

| classifier input | AUC |
|---|---|
| full feature set | 0.79 |
| thickness alone | 0.98 |
| structure tensor and curvature only | 0.51 |
| thickness-matched subset | 0.48 |

Local CT geometry carries no weld information; the headline 0.79 is measuring
how thick the band is, which a model already has. Runs on CPU in about half an
hour.

The table needs eight crops, four from Scroll 1 and four from Scroll 4. The
`patches-v1` release holds two, so it is a smoke test for the code rather than a
reproduction of the table. On two crops the per-volume cap of 60 sites per class
leaves n around 170, and the thickness-matched control falls under the 60-pair
minimum the script requires, so that row prints `nan` while the other three
still print. The matched row is the one the conclusion rests on, so treat the
table as reproduced only on eight crops.

## Patch-mode eval for external splitters

`scripts/eval_patch.py` scores any splitter/segmenter output against a GT
surface patch, reporting the official blend (TopoScore / SurfaceDice@2 / VOI)
plus each term, connected-component stats and, when instance labels are
given, instance count and instance-level VOI. A two-input mode prints an
A -> B delta table (e.g. splitter OFF vs ON), which is the number that
matters for split evaluation — see the metric-response notes in the #191
thread (the blend's VOI term rewards predicted mass, so read topology
changes on the component/Betti side).

```bash
# single prediction (mask tif/npy, or npz with an int 'labels' array)
python scripts/eval_patch.py --gt labelsTr/s1_z10240_y2560_x2560.tif        --pred my_instances.npz
# OFF -> ON delta
python scripts/eval_patch.py --gt gt.tif --pred off.npz --pred-b on.npz
```

Patch layout (two sample patches are attached to the `patches-v1` release:
one Scroll 1 crop at 300^3, one Scroll 4 / PHerc1667 crop at 236^3):

```
imagesTr/{scroll}_z{Z}_y{Y}_x{X}_0000.tif  uint8 raw CT, zyx order
labelsTr/{scroll}_z{Z}_y{Y}_x{X}.tif       uint8 binary GT surface, same grid
```

Filename coordinates are the crop origin in scroll voxel space. The raw CT
channel is a direct re-extraction from the released scroll volumes; no
Frangi/EDT preprocessing anywhere in the eval loop.

## License

MIT. Data and models are Vesuvius Challenge resources (CC BY-NC 4.0 — credit
to the Vesuvius Challenge team and the `Dataset059` / model authors).

## Twin validation (2026-08)

The ray flag and the continuity splitter were validated against exact ground
truth on a finite-thickness, instance-labelled build of Diego-dcv's synthetic
scroll twin (harness: `scripts/twin_harness.py`; the build and its two bugs in
the shipped painter are documented in
https://github.com/Diego-dcv/vesuvius-topological-grid/issues/1).

Flag: precision 98.2 percent against true turn ids in both sampling seeds,
recall 73.1 and 72.7 percent, control contamination 0 of roughly 11,400 sites
in both. Splitter (production repair path on merged pairs): 427 of 450
true-interface sites split, identity accuracy 99.2 and 98.9 percent over 2.2
million reassigned voxels on disjoint site subsamples; the 23 refusals are
explicit outcomes concentrated in the three-sheet weld sector. Verdict files
and the build parameters are under `results/twin/`; measurement scripts are
`scripts/twin_measure.py` and `scripts/twin_split.py` (paths inside point at
the box layout and take the run directory as argument).
