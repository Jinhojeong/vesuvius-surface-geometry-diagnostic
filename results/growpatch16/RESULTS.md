# Powered16: GrowPatch hazard A/B/C pooled over 16 sites

Generated 2026-08-12 11:53:10. 8 archived sites (powered + randctl blocks, 2026-08-06) + 8 new sites (this block, 2026-08-12), every site on its own census tile. Unit = run; paired t across the 16 site pairs of replicate means (df=15); scorer score_ab.score_run UNMODIFIED (2000 rays, near_r 6). Gate: this analyzer reproduced the archived 8-site contrasts exactly before pooling (see gate block in RESULTS.json).

## Paired contrasts (16 sites)

### B-A

| endpoint | mean diff | sd | t(15) | p | 95% CI | MDE80 | sites>0 |
|---|---|---|---|---|---|---|---|
| area_vx2 | +120.3k vx2 | 215.9k | 2.23 | 0.0416 | [+5.2k, +235.3k] | 161.8k | 11/16 |
| crossing frac_quads_near (pp) | -0.42 pp | 0.64 | -2.63 | 0.0188 | [-0.76, -0.08] | 0.48 | 4/16 |
| ray frac_double (pp of on-sheet rays) | -1.01 pp | 3.40 | -1.19 | 0.2528 | [-2.83, +0.80] | 2.55 | 8/16 |
| on_sheet_rate (guardrail, pp) | -1.52 pp | 3.52 | -1.73 | 0.1047 | [-3.39, +0.36] | 2.64 | 6/16 |
| n_double per 1000 sampled quads | -4.55 per-1000 | 11.01 | -1.65 | 0.1191 | [-10.42, +1.32] | 8.25 | 3/16 |

### B-C

| endpoint | mean diff | sd | t(15) | p | 95% CI | MDE80 | sites>0 |
|---|---|---|---|---|---|---|---|
| area_vx2 | +123.0k vx2 | 197.9k | 2.49 | 0.0252 | [+17.5k, +228.4k] | 148.3k | 12/16 |
| crossing frac_quads_near (pp) | -0.42 pp | 0.55 | -3.09 | 0.0075 | [-0.72, -0.13] | 0.41 | 4/16 |
| ray frac_double (pp of on-sheet rays) | -1.59 pp | 2.75 | -2.32 | 0.0352 | [-3.05, -0.13] | 2.06 | 5/16 |
| on_sheet_rate (guardrail, pp) | -1.43 pp | 3.01 | -1.90 | 0.0767 | [-3.04, +0.17] | 2.26 | 4/16 |
| n_double per 1000 sampled quads | -5.68 per-1000 | 9.22 | -2.46 | 0.0263 | [-10.59, -0.77] | 6.91 | 4/16 |

### C-A

| endpoint | mean diff | sd | t(15) | p | 95% CI | MDE80 | sites>0 |
|---|---|---|---|---|---|---|---|
| area_vx2 | -2.7k vx2 | 72.0k | -0.15 | 0.8835 | [-41.0k, +35.7k] | 53.9k | 7/16 |
| crossing frac_quads_near (pp) | +0.00 pp | 0.45 | 0.04 | 0.9703 | [-0.24, +0.24] | 0.34 | 8/16 |
| ray frac_double (pp of on-sheet rays) | +0.58 pp | 1.94 | 1.19 | 0.2506 | [-0.45, +1.61] | 1.45 | 9/16 |
| on_sheet_rate (guardrail, pp) | -0.09 pp | 2.24 | -0.15 | 0.8800 | [-1.28, +1.11] | 1.68 | 8/16 |
| n_double per 1000 sampled quads | +1.12 per-1000 | 6.03 | 0.75 | 0.4672 | [-2.09, +4.33] | 4.52 | 9/16 |

## Per-site replicate means (dbl pp / area k-vx2)

| site | A dbl | B dbl | C dbl | A area | B area | C area |
|---|---|---|---|---|---|---|
| s0 | 19.9 | 14.8 | 20.5 | 636k | 989k | 663k |
| s1 | 18.9 | 13.1 | 18.2 | 842k | 1113k | 800k |
| s2 | 13.2 | 5.4 | 12.4 | 628k | 961k | 620k |
| s3 | 16.7 | 16.8 | 15.1 | 653k | 1076k | 641k |
| s4 | 7.7 | 4.4 | 6.8 | 902k | 1079k | 994k |
| s5 | 14.0 | 16.6 | 17.1 | 910k | 1151k | 862k |
| s6 | 18.3 | 17.2 | 19.9 | 963k | 911k | 901k |
| s7 | 23.3 | 22.6 | 25.8 | 904k | 1031k | 942k |
| s8 | 14.4 | 15.1 | 13.1 | 894k | 835k | 924k |
| s9 | 24.1 | 23.1 | 23.5 | 942k | 1115k | 1066k |
| s10 | 22.9 | 22.9 | 25.8 | 1259k | 1037k | 1066k |
| s11 | 21.5 | 17.5 | 19.0 | 1107k | 821k | 1057k |
| s12 | 15.4 | 20.6 | 19.4 | 687k | 867k | 730k |
| s13 | 25.4 | 26.7 | 25.7 | 914k | 784k | 906k |
| s14 | 17.7 | 19.1 | 18.0 | 639k | 984k | 638k |
| s15 | 16.9 | 18.2 | 19.2 | 887k | 938k | 914k |

## Leave-one-site-out sensitivity (dbl)

### B-A: p range [0.0944, 0.4805]

| dropped site | mean diff (pp) | p |
|---|---|---|
| s0 | -0.74 | 0.4060 |
| s1 | -0.69 | 0.4251 |
| s2 | -0.56 | 0.4805 |
| s3 | -1.09 | 0.2504 |
| s4 | -0.86 | 0.3529 |
| s5 | -1.25 | 0.1736 |
| s6 | -1.01 | 0.2865 |
| s7 | -1.03 | 0.2759 |
| s8 | -1.13 | 0.2315 |
| s9 | -1.02 | 0.2828 |
| s10 | -1.08 | 0.2516 |
| s11 | -0.81 | 0.3739 |
| s12 | -1.43 | 0.0944 |
| s13 | -1.16 | 0.2146 |
| s14 | -1.17 | 0.2108 |
| s15 | -1.17 | 0.2133 |

### B-C: p range [0.0191, 0.0709]

| dropped site | mean diff (pp) | p |
|---|---|---|
| s0 | -1.31 | 0.0709 |
| s1 | -1.36 | 0.0695 |
| s2 | -1.23 | 0.0691 |
| s3 | -1.81 | 0.0209 |
| s4 | -1.54 | 0.0543 |
| s5 | -1.66 | 0.0390 |
| s6 | -1.52 | 0.0567 |
| s7 | -1.48 | 0.0604 |
| s8 | -1.83 | 0.0191 |
| s9 | -1.67 | 0.0379 |
| s10 | -1.51 | 0.0577 |
| s11 | -1.59 | 0.0476 |
| s12 | -1.77 | 0.0251 |
| s13 | -1.76 | 0.0264 |
| s14 | -1.77 | 0.0254 |
| s15 | -1.63 | 0.0431 |

## Block heterogeneity (archived 8 vs new 8, B-A)

The new block is all mid-scroll by construction (the selection mid-band pass filled all 8 slots); the archived block mostly sat near the scroll ends. Same paired method within each block (df=7).

| endpoint | old8 diff | old8 p | new8 diff | new8 p |
|---|---|---|---|---|
| area_vx2 | +234.0k vx2 | 0.0031 | +6.6k vx2 | 0.9348 |
| crossing frac_quads_near (pp) | -0.58 pp | 0.0469 | -0.26 pp | 0.2537 |
| ray frac_double (pp of on-sheet rays) | -2.65 pp | 0.0682 | +0.62 pp | 0.5178 |
| on_sheet_rate (guardrail, pp) | -1.23 pp | 0.3948 | -1.81 pp | 0.1772 |
| n_double per 1000 sampled quads | -7.02 per-1000 | 0.0842 | -2.08 per-1000 | 0.6434 |

## Guardrail

on_sheet_rate does not differ between B and A (p=0.105); the dbl denominator is arm-comparable.
