# Placement A/B: seed-placement test of the hazard-weight effect

Preregistered design: prereg_placement_ab/PREREGISTRATION.md at commit 8531276; Amendment 1 (off-seed direction-field regions, same builder and parameters) at commit 083f424. Generated 2026-08-12 20:30:36 on the GPU box (CPU only).

## Gate

Archived powered_scores.json sha256 verified (fd9deb2c3b7a3e3f...). All archived 8-site contrasts (5 endpoints: cell means, mean_diff, sd, t, p, CI95, MDE80) reproduced to full float precision before any new data was read.

## Selection (frozen record)

offseeds_placement.json sha256 b07e905b6dee57c7. All 8 off-seeds resolved at the primary distance rule 96 (no relaxation). Sites s2/s4 share a tile; the deterministic rule gives them the identical off-seed point (disclosed in Amendment 1).

| site | slab/tile | off L1 (z,y,x) | m7 | minCheb to cluster | rule |
|---|---|---|---|---|---|
| o0 | z4032/tile_y1344_x448 | 4240,1584,752 | 255 | 109 | 96 |
| o1 | z224/tile_y1792_x1344 | 432,1840,1392 | 255 | 107 | 96 |
| o2 | z0/tile_y448_x3136 | 160,512,3488 | 255 | 145 | 96 |
| o3 | z11200/tile_y448_x2240 | 11264,544,2672 | 255 | 134 | 96 |
| o4 | z0/tile_y448_x3136 | 160,512,3488 | 255 | 145 | 96 |
| o5 | z9184/tile_y1792_x3136 | 9376,2016,3280 | 255 | 103 | 96 |
| o6 | z224/tile_y1792_x2240 | 304,2016,2704 | 255 | 146 | 96 |
| o7 | z672/tile_y2240_x2240 | 768,2288,2672 | 255 | 257 | 96 |

## Off-seed field build verification (Amendment 1)

| site | field u8 (x,y,z) | |v| decoded | written |
|---|---|---|---|
| o0 | (255, 128, 128) | 1.000 | True |
| o1 | (128, 128, 1) | 1.000 | True |
| o2 | (128, 128, 255) | 1.000 | True |
| o3 | (255, 128, 128) | 1.000 | True |
| o4 | (128, 128, 255) | 1.000 | True |
| o5 | (128, 128, 1) | 1.000 | True |
| o6 | (128, 1, 128) | 1.000 | True |
| o7 | (128, 128, 1) | 1.000 | True |

## Completion

80/80 runs scored; 0 missing.

## Headline: double difference DD = (B-A)_on - (B-A)_off

| site | area (B-A)_on | area (B-A)_off | area DD | dbl DD | on_sheet DD |
|---|---|---|---|---|---|
| s0 | 353k | -1k | 354k | n/a | -6.88pp |
| s1 | 271k | 18k | 253k | n/a | +5.00pp |
| s2 | 332k | 25k | 307k | n/a | +1.29pp |
| s3 | 423k | -111k | 534k | n/a | +1.64pp |
| s4 | 177k | 34k | 143k | n/a | -2.81pp |
| s5 | 241k | -71k | 311k | n/a | -4.34pp |
| s6 | -52k | 4k | -56k | n/a | -3.37pp |
| s7 | 127k | -3k | 130k | n/a | -0.35pp |

PRIMARY area_vx2 DD (H1 > 0, one-sided, df 7): mean 247k, sd 176k, t=3.96, one-sided p=0.0027, 95% CI [99k, 394k] -> CONFIRMS

SECONDARY frac_double DD: UNMEASURABLE per Amendment 3 (off cells carry zero on-sheet rays; 8 of 8 sites missing).

GUARDRAIL on_sheet_rate DD (two-sided): mean -1.23pp, t=-0.91, p=0.3952, 95% CI [-4.43pp, +1.98pp] -> not flagged

Reported, not tested: crossing frac_quads_near DD mean -0.58pp (95% CI [-1.15pp, -0.01pp]); per-1000-quads DD mean -7.02 (95% CI [-15.28, 1.24]).

## Off-seed (B-A) itself (mechanism predicts near zero)

| endpoint | mean | 95% CI |
|---|---|---|
| area_vx2 | -13k | [-56k, 30k] |
| frac_double | unmeasurable (Amendment 3) | n/a |
| on_sheet | -0.00pp | [-0.01pp, +0.01pp] |
| crossing near | +0.00pp | [-0.00pp, +0.00pp] |

Archived on-seed (B-A) area mean 234k; off-seed point estimate is -5.5% of it.

## Leave-one-site-out (primary, one-sided p, df 6)

| dropped | p |
|---|---|
| s0 | 0.0080 |
| s1 | 0.0071 |
| s2 (same-tile pair member) | 0.0078 |
| s3 | 0.0045 |
| s4 (same-tile pair member) | 0.0048 |
| s5 | 0.0078 |
| s6 | 0.0007 |
| s7 | 0.0045 |

Range [0.0007, 0.0080]. s2 and s4 share a tile and one off-seed point; their drops are listed above.

## Verdict (frozen buckets)

BUCKET 2: placement confirmed for growth; the quality claim stays at its current placebo-level evidence
