# Route A: per-point winding assignment on the full PHerc1218 label tree

Outputs of a full-resolution run requested on
[PR #1](https://github.com/Jinhojeong/vesuvius-surface-geometry-diagnostic/pull/1),
placed here next to the crossing table they were derived from so they can be
pulled by raw URL the same way `positions_merged.csv.gz` is.

## Who did what

The script is `route_a_winding_fullres.py` by Diego-dcv, in
[vesuvius-topological-grid](https://github.com/Diego-dcv/vesuvius-topological-grid).
The labels and the crossing table are mine. The run is mine, on my hardware.
Anything derived from these files should credit the method to him.

## What was run

Input tree is the pre-repair stitch, 1,369 npz blocks, the same source the
step-8 point list was exported from and the same source the crossing table was
cast from. Block visit order `paths`, which is plain lexicographic on the block
paths and matches the exporter. Single process, 920 seconds.

## What came out

216,000,329 labelled voxels lie on the crossing table's 323 planes once the
first-writer rule has resolved block overlaps. Of those, 185,331,060 assign to
a winding without ambiguity, 85.8 percent against a preregistered floor of 60,
and the median radial distance to the assigned crossing is 2.625 voxels against
a preregistered ceiling of 3. The step-8 run gave 85.8 percent and 2.61 on
3,376,037 points, so a sixty-four-fold increase in sampling moved neither
endpoint.

Seam agreement is measured here on all 106,629,002 overlap voxels of those
planes and comes out at 88.8 percent. That supersedes the 82.9 percent I
published from a narrow sample of my own; the scope of that older figure is
described in `results/pointlist_export/NOTE_seam_agreement.md`.

Instance QA flags 139,083 of 475,845 instances whose assigned windings show a
gap of 2 or more, 29.2 percent, against 27.7 percent at step-8. An instance
like that cannot be a single spiral sheet, so the list is a candidate set for
either a label fault on my side or an assignment error on his.

## Files

`winding_maps_1218.npz` carries n, r_p10, r_p90 and median dr per winding, ray
and plane, so `r_p90 - r_p10` is a per-cell sheet-thickness proxy.
`qa_instances_fullres.csv` is the gap-2 instance list. `summary.json` carries
the counts, the exam verdicts, the seam statistics and the parameters.

One limit carried over from the step-8 run and stated by its author. A switch
to the adjacent winding, a gap of exactly 1, is invisible to this test.
