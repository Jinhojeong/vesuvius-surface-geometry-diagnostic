# Amendment 4 to the tight-contact validation set preregistration

Amends PREREGISTRATION.md, sha256 prefix `9d5d71dbaf45ab85`. AMENDMENT_1
(`a07c2f86f39b3bd8`), AMENDMENT_2 (`f9d57773c715cb65`) and AMENDMENT_3
(`db61e452e4cce9ca`) stand unchanged. Written 2026-08-22, after version 4 was
built and before it was published, and before any version 5 crop exists.

## What is wrong

AMENDMENT_3 made the block containing the site authoritative and wrote it first.
The site falls inside more than one block for 76 of the 300 crops, so that rule
had to pick one, and it picked the first in coordinate order. That is the wrong
one.

`A` and `B` in the census are block-local ids **in the numbering of the block the
census read**, which is the block named in the crop's own filename. The
authoritative block has to be that one, because `A_id = base + A` is only
meaningful when `base` belongs to the numbering `A` came from. For the 37 crops
where coordinate order picks a different block, `A_id` and `B_id` name objects
that do not exist.

| | |
| --- | --- |
| crops where the coordinate-order block differs from the census block | 37 of 300 |
| of those, carrying both ids under version 4's manifest | 0 |
| both-instances count, version 4 as built | 263 of 300 |
| both-instances count, ids resolved through the census block's base | 276 of 300 |

Version 4's own verification reported "A_id/B_id are built on the authoritative
site block, 300 of 300" and passed. That check was tautological. It confirmed
that the ids were built on whichever block the code had called authoritative,
which is true by construction and says nothing about whether that block was the
right one.

Versions 1 to 3 chose the site block by taking whichever candidate came last in
filename order, which is a third arbitrary rule and also not the census block.

## A second assumption removed

The block index was built from filenames with the shape assumed to be 256 by 512
by 512. Seventeen of the 1,321 blocks are trimmed at the volume edge, ten at
(256, 512, 213) and seven at (256, 213, 512). Sixteen crops touch one of those.
Overlap detection happened to agree under the assumed and the real shapes, 28
pairs either way, and a mismatched slice would have raised rather than passed
silently, so version 4 is not corrupted by this. Version 5 reads the real shape
from each block's header anyway, so the assumption is gone rather than checked.

## What changes

The authoritative block is the census block, identified by the slab and tile the
census row carries and encoded in the crop filename. It is written first, and
every other contributing block fills only voxels nothing has written yet. Block
shapes come from the stored array headers.

## What does not change

Every selection rule, every band edge, the per-band target, the octant rule, the
overlap exclusion, census order, the control-arm rule and seed, and the level-1
CT read. The 300 crops and 60 controls are the same crops. Gap values are carried
over and not recomputed.

## What this predicts, frozen here

Both-instances should come out at 276 of 300. That number is not a guess, it is
what resolving version 4's own crops through the census base already gives, so
the corrected run should reproduce it exactly rather than approximately. The 37
crops whose authoritative block changes should be exactly the crops whose ids
change; no other crop's `A_id` or `B_id` should move. Pair-voxel loss should stay
at zero. No gap, band or membership should change.

AMENDMENT_3 predicted a rise of about 3 and got 18, because it described only the
overwrite defect and the site-block ambiguity had not been found when it was
frozen. This document is written after both are understood, so a miss here has no
such excuse.

## What would make this amendment wrong

If both-instances comes out at anything other than 276 of 300, the account above
is incomplete. If any crop outside the 37 changes its ids, the fix reaches
further than it should. If any gap value, band or crop identity changes, the run
is void. If pair-voxel loss is non-zero anywhere, AMENDMENT_3's fix has been
undone.

## Disclosure

The defect was found by an adversarial pre-publication audit of the version-4
draft comment and card, run because five errors had already been found in this
material over two days. Version 4 was never published; it exists only on the
working box. Versions 1 to 3 remain on Kaggle, marked defective, so that anyone
who pulled them can tell what they have.
