# Scoring the facing-pairs CC guard against exact truth

flummoxjr's harness (`flummoxjr/facing-pairs-harness`), run on the eight cells
of aviad12g's physical-fusion pilot, which carry a probability map, a surface
label and a per-voxel `turn_id`, exact sheet truth by construction. He asked,
before anyone ran anything, for two numbers, the fraction of step-3 CC-guard
rejections whose two points lie on different sheets in truth, and the
false-accept rate on the sites the guard keeps.

## Chain of custody, all checked on this machine before scoring

The harness self-test passes, three stages. His port is byte-identical to the
archived original by his own verifier, re-run here, PASS. The eight cells pass
aviad's shipped validator, 156 of 156 checks, and the manifest sha256 equals
the value aviad published on villa #191, `506f3353...`. Batch settings were
`--gt-instance turn_id`, `--log-rejected` and `--seed-mode caseid`.

## The two numbers he asked for

| d12 bin | rejected pairs | truly different sheets | accepted sites | false accepts |
| --- | --- | --- | --- | --- |
| 2 to 3 vox | 93,205 | 97.6 percent | 790 | 0 |
| 3 to 4 | 913 | 99.2 percent | 954 | 0 |
| 4 to 5 | 16 | 0.0 percent | 308 | 0 |
| 5 to 6.01 | 0 | none | 348 | 0 |
| all | 94,134 | 97.6 percent | 2,400 | **0.00 percent** |

Ground truth resolved for every pair; no unresolved ids in either file.

**The guard's rejections are real.** Overall 97.6 percent of rejected pairs
are two distinct sheets in truth, pairs the label volume merges locally, which
is exactly the reading his Dataset059 survival rates (0.8, 5.2, 11.4 percent
at 2-3, 3-4, 4-5 vox) assume. And the guard never falsely accepts: all 2,400
kept sites are true two-sheet pairs, in every bin, on every cell.

## The structure worth knowing, which is geometry rather than guard error

Per cell, the truly-different fraction is 95.9 to 99.7 percent at pitches 160
to 250, and **0.0 percent at pitch 300 on both seeds**, 470 and 471 rejects,
every one carrying the same turn_id on both ends. At the loosest pitch the
sheets sit too far apart to produce tight facing pairs, so the only tight
candidates are a single sheet curving back on itself, and the guard rejects
them, correctly, since they are not two sheets. So the rejected population
changes meaning with pitch: at tight pitch it measures label merging, at loose
pitch it measures self-facing curvature. The 4-5 vox row is the same effect at
n=16, all from the pitch-300 cells.

Seeds 0 and 1 agree to within a few tenths of a percent per cell.

## Limits

Synthetic cells, one contrast (papyrus 90), one noise level, four pitches, two
seeds. Turn ids as sheet truth, so a turn that touches itself counts as one
sheet, which is precisely what the pitch-300 rows show. Tight sites only, per
his stated ask. The full rejected and site CSVs are regenerated exactly by his
harness on aviad's cells with the settings above; `fp_score.py` here computes
the table from them, and `facing_pairs_case_diag.csv` is included as shipped.
