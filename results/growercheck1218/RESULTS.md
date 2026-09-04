# Does the grown surface follow a sheet?

Executed against MEASUREMENT_RULES.md, sha256 prefix `8a6d9369bc7167ce`, frozen
and published before the run. 360 seeds on PHerc1218, 68 minutes, no site
excluded. A second scroll was then added to test whether the answer travels.

## Answer

The surfaces lie flat where a sheet stands up, and that holds on both scrolls.

| scroll | surface | n | median normal | 95 percent interval |
| --- | --- | --- | --- | --- |
| PHerc1218 | grown here | 360 | **0.942** | 0.938 to 0.946 |
| PHerc1218 | verified patches, outside project | 6 | **0.211** | not computed, n=6 |
| PHerc1447 | grown here | 15 | **0.927** | 0.857 to 0.953 |
| PHerc1447 | published segments | 15 | **0.220** | 0.211 to 0.232 |

The number is the median absolute z-component of the unit surface normal. A
sheet in a standing scroll is close to vertical, so its normal is close to
horizontal and this reads low. A surface lying flat across the windings reads
near one. Independently produced surfaces sit at 0.21 to 0.22 on both scrolls.
Surfaces grown here on the published m7 predictions sit at 0.93 to 0.94.

## The second metric does not travel, and that is worth saying

The frozen rules put on-prediction fraction first, the share of mesh vertices
that land on a nonzero voxel of the very prediction the surface was grown from.
On PHerc1218 it separated cleanly, 27.0 percent for the grown surfaces against
82.8 percent for the verified patches. It does not survive the second scroll.

| scroll | surface | median on-prediction |
| --- | --- | --- |
| PHerc1218 | grown | 27.0 percent |
| PHerc1218 | verified patches | 82.8 percent |
| PHerc1447 | grown | 21.0 percent |
| PHerc1447 | published segments | **23.3 percent** |

On PHerc1447 the published segments read 23.3 percent, which is where the grown
surfaces read, while their normals say they are correct sheets. So this quantity
depends on how densely the prediction covers the material, not only on whether
the surface is right. It is reported here because the rules named it first, and
it is retired as a cross-scroll discriminator on this evidence. The normal is
the number that carries.

## It is not about contact

On PHerc1218 every gap band comes out the same.

| band | microns | seeds | on-prediction, median | 95 percent interval |
| --- | --- | --- | --- | --- |
| 0 to 2 vox | under 34.6 | 60 | 27.0 percent | 25.7 to 28.0 |
| 2 to 4 | 34.6 to 69.1 | 60 | 27.5 percent | 26.7 to 28.8 |
| 4 to 6 | 69.1 to 103.7 | 60 | 27.7 percent | 26.2 to 28.7 |
| 6 to 10 | 103.7 to 172.8 | 60 | 26.3 percent | 25.2 to 27.2 |
| 10 and above | above 172.8 | 60 | 27.3 percent | 27.0 to 28.3 |

Contact seeds read 27.3 percent and single-sheet controls 24.7 percent, and the
intervals overlap. The frozen rules name that outcome in advance as the result
rather than a null to work around. I expected tight contacts to be worse. They
are not.

## What the failure conditions did

None fired against the run. The discard rate is **0 of 360** against a one third
gate. The reader reproduced the six reference patches at 65.2 to 94.1 percent,
the same range and median as the six-site probe that preceded it. The third
condition did fire in the sense that contact and control did not separate, and
that is reported above.

## What was ruled out

The only term in this tracer that constrains a patch's in-plane basis lives
behind the `direction_fields` array and is a no-op when that array is empty. So
four seeds were rerun with the one direction field published for this scroll,
`dir: normal` against the m7 normals at scale 1, with the log confirming the
field loaded. Neither endpoint moved, on-prediction 21.3 to 23.3, 28.3 to 28.0,
28.3 to 29.0 and 13.2 to 16.7 percent, normals unchanged to two decimals.

One error was caught by the reference arm rather than by me. The first reading
of the verified patches returned 0.0 percent on every one, because those patches
carry level-0 coordinates and I read them against the prediction's level 1.
Without a positive control that would have been read as the tracer failing
completely, and the whole result would have been wrong in the same direction it
now points.

## Limits

Two scrolls, one prediction family, one parameter set, one binary. The PHerc1218
seeds are the tight-contact validation set, so they sample sites where a repair
fired plus single-sheet controls, not the scroll uniformly. The PHerc1447 arm is
fifteen published segments with their own centre vertex as the seed, which makes
it a paired comparison but a small one. Neither reference arm records the
parameters it was produced with, so the comparison is against surfaces of
unknown provenance. Nothing here establishes a cause. In particular no
horizontal fibre direction field was tested, which is what the project's own
evaluation config supplies for seeding and which is not published for either
scroll.

## Reproducing

`MEASUREMENT_RULES.md` carries the binary, the volume URLs, the parameters and
the discard rules. PHerc1218 seeds are all 360 rows of
`jhjeong0815/pherc1218-tight-contact-val` version 5, doubled from level 1 to
level 0. PHerc1447 seeds are the centre valid vertex of each published segment's
`mesh/intermediate/tifxyz_original`. `summary.json` in each run directory carries
the per-arm figures.
