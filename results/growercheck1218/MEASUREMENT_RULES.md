# Measurement rules, frozen before the run

Frozen 2026-09-04, after a six-site probe and before any full-set number exists.
This is a measurement, not an intervention test, so it fixes what is measured
and what is discarded rather than predicting an effect.

## The question

When the official surface tracer is seeded on the published PHerc1218 m7
surface prediction, how often does the surface it grows leave the prediction it
was grown from, and does that depend on whether the seed sits at a fused
contact?

## Fixed before the run

**Binary.** `vc_grow_seg_from_seed` from the published AppImage
`VC3D-23adee0-2026-09-03-linux-x86_64.AppImage`, extracted, run through its own
`AppRun`. No local build, no container.

**Volume.** The published prediction
`PHerc1218/representations/predictions/surfaces/20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr`,
read at level 0, shape 23247 by 7593 by 7593, streamed over https. Nothing is
mirrored or preprocessed.

**Seeds.** All 360 rows of the tight-contact validation set version 5, which is
300 contact crops carrying both split instances and 60 single-sheet controls.
Site coordinates are level-1 and are doubled to reach level 0. No site is
selected or dropped on anything measured here.

**Parameters.** mode seed, generations 60, step size 20, voxelsize 8.64,
min_area_cm 0, thread limit 4, no direction field. The probe established that
adding the one direction field available for this scroll, `dir: normal` against
`m7_normals_L1.zarr` at scale 1, moves neither endpoint, so the main run omits
it and a 30-site subset repeats with it as a check.

**Endpoints, in this order.**
1. On-prediction fraction. 300 mesh vertices per surface, sampled without
   replacement under a fixed seed, each rounded to the nearest voxel and read in
   the same prediction the surface was grown from. The fraction reading nonzero.
2. Median absolute z-component of the unit surface normal, taken over quads with
   a non-degenerate cross product.
Both are reported per arm and per gap band, as medians over sites with a
bootstrap interval over sites.

**Reference.** The same two numbers computed by the same code on the six
verified patches published in `IyanDopico/vesuvius-sheet-tools`, which are
level-0 tifxyz on this scroll. The probe put them at 65 to 94 percent
on-prediction and 0.14 to 0.35 normal, and those figures are the comparison the
run is read against.

## Discards, declared

A site is discarded and counted if the grower exits non-zero, if it writes no
tifxyz, if fewer than 30 grid points are valid, or if fewer than 20 quads have a
non-degenerate normal. A vertex is discarded and counted if its rounded
coordinate falls outside the volume. Discards are reported as counts, not
silently dropped.

## What would make this run wrong

If the discard rate exceeds a third of sites, the run is not measuring the
population it claims. If the six reference patches do not reproduce the probe's
65 to 94 percent under the same code in this run, the reader is not stable and
no comparison may be drawn. If the control arm and the contact arm differ by
less than their bootstrap intervals, the answer is that this failure is not
specific to fused contacts, and that is the result rather than a null to work
around.

## What I expect, written so it cannot be constructed afterwards

From the probe I expect on-prediction in the teens to low thirties on both arms
and normals near 0.9 on both, that is, a large gap against the reference and no
separation between contact and control. The probe was six sites, so that prior is wide. A contact-versus-control separation appearing at 360 sites would
be the more interesting outcome and would reopen the question the probe closed.
