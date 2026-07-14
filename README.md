# Geometry-stratified diagnostic for Vesuvius surface models

Quantifies **where and how much** Herculaneum-scroll surface-detection models
degrade in curved and compressed regions, using the public
`Dataset059_s1_s4_s5_patches_frangiedt` ground truth (200 random patches,
seed-fixed) and two public baselines. Posted to
[ScrollPrize/villa#191](https://github.com/ScrollPrize/villa/issues/191#issuecomment-4930879524).

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

Caveats: patches likely overlap both models' training distribution, so the
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
patches (~5h, one 11GB 2080 Ti) lifts <4 vox spacing AUC from **0.41 to 0.65**,
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

## License

MIT. Data and models are Vesuvius Challenge resources (CC BY-NC 4.0 — credit
to the Vesuvius Challenge team and the `Dataset059` / model authors).
