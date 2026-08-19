# Preregistration: a tight-contact validation set on real CT (PHerc1218)

Frozen 2026-08-19, before any crop was extracted. Written because the set is
meant to referee other people's models, so its selection rule has to be
fixed in advance rather than tuned once the contents are visible.

## Why

Training-set validation for surface models currently runs on data whose
sheets are not in contact. A contributor training on the public patch set
measured its held-out geometry as median inter-sheet gap 15.6 voxels with
0.02 percent of near-band positives under 4 voxels, and reported that its
labels merge below about 5 voxels. Any effect attributed to compressed or
touching sheets is therefore unmeasurable on that data. This set exists to
supply the missing regime from real CT rather than from a phantom.

## Source

PHerc1218, masked CT volume `20250521120456-8.640um-1.2m-116keV-masked.zarr`,
level 0, paired with the repaired instance labels of the published
continuity-repair v2.0 release. Contact sites come from the whole-scroll
fused-site census, 644,148 candidate sites of which 109,309 were split, with
the 9,717-row fused-suspect table as the ranked entry point.

## Selection rule, fixed here

1. Candidate sites are census rows with `kind == fused` and a ray-validated
   split in the v2.0 repair records. Ray validation is the same test used in
   the release, unchanged.
2. For each candidate, the local gap is measured along the label normal
   between the two instance ids the split separated, at the site voxel, in
   level-0 voxels. This is the only gap definition used anywhere in the set.
3. Sites are binned by that gap into bands 0 to 2, 2 to 4, 4 to 6, 6 to 10
   and above 10 voxels.
4. Crops are 128 cubed, centred on the site, taken at level 0, and are
   accepted only if the cube lies wholly inside the volume and carries CT in
   every one of its eight octants, so no crop is half empty.
5. Crops are drawn in census row order within each band, not by score, and a
   site is skipped if its cube overlaps an already accepted cube by more than
   one quarter of its volume. Target is 60 crops per band where the band has
   that many, and every band's realised count is reported rather than padded.
6. A control band of 60 crops is drawn the same way from sites the census
   examined and did not flag, to give a matched no-contact arm.

## What ships per crop

Intensity as uint8, the binary surface label, the instance label map, the
measured gap, the band, the census row id, and the coordinates. Nothing is
smoothed, inpainted or relabelled.

## Declared limits, written before the numbers exist

The instance labels are the repaired v2.0 labels, not hand annotation, so
this set measures agreement with a repaired automatic labelling and not with
truth. The repair's own limits carry over, including the measured per-site
rates published with the release. The gap definition is normal-direction and
local, so it will disagree with any pairwise or connected-component
definition by construction, and the disagreement is the point of comparison
rather than an error in either. Sites come from one scroll.

## What would make this set wrong

If the accepted crops turn out to concentrate in a few slabs, or if the
tight bands are dominated by one instance pair repeated across neighbouring
sites, the set is not a sample of the contact regime and this document will
say so in the results rather than the set being reissued quietly.
