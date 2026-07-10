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

Additional observations: recall-only metrics flatter high-firing-rate models in
exactly these hard regions (they made `surface_recto` look immune to
compression — it isn't); and in the lowest-dice patches the predictions
visibly track fine CT sheet structure that the sparse GT bands miss, so part of
the low tail is plausibly **label incompleteness** (moderate evidence — see the
CT-brightness oracle in `results/oracle.csv`).

Caveats: patches likely overlap both models' training distribution (curves are
optimistic bounds); this is not the official Kaggle Surface Detection metric —
it characterizes the failure mode, it does not claim to close #191; geometry
proxies derive from GT, so label quality confounds the strata; the <6 vox
spacing bins are thin (~2k background points).

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
