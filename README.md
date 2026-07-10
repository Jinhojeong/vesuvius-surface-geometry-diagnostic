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

## Contents

- `scripts/diag2_191.py` — recall stratification: sliding-window inference +
  GT geometry (structure-tensor curvature, along-normal sheet spacing)
- `scripts/diag4_auc.py` — the two-model per-bin AUC comparison (main result)
- `scripts/diag3_d058.py` — nnUNetv2 baseline runner (trainer-fallback recipe)
- `scripts/loader059.py` — loading recipe for `Model_epoch499.pth`
  (`_orig_mod.` prefix + duplicate per-task encoder keys + `separate_decoders`)
- `scripts/oracle191.py` — CT-brightness label-incompleteness probe
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
