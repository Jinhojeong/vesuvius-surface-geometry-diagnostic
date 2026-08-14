#!/bin/bash
# Render the archived hazard A/B patches (no weight vs hazard weight) so the
# measured +236k / -3.12pp result exists as pictures.
set -u
ROOT=/mnt/vesuvius/vcbuild; P8=/mnt/vesuvius/p8_sprint
CT="https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr"
mkdir -p $P8/render_hazard
for site in s0_z4032_tile_y1344_x448 s1_z224_tile_y1792_x1344 s3_z11200_tile_y448_x2240; do
  for arm in A B; do
    src=$ROOT/demo_out/powered/$site/${arm}_r1/out
    [ -d "$src" ] || { echo "missing $src"; continue; }
    o=$P8/render_hazard/${site%%_*}_${arm}
    [ -f "$o/00.tif" ] && { echo "skip $site $arm"; continue; }
    mkdir -p "$o"
    (cd $ROOT && timeout 2400 env OMP_NUM_THREADS=4 nice -n 12 ./vcrun.sh vc_render_tifxyz \
       -v "$CT" --remote-url "$CT" --scale 0.5 -g 1 -s "$src" \
       --tif-output "$o" --flatten --auto-crop --cache-gb 8 --pyramid 0) \
       > $P8/logs/rh_${site%%_*}_${arm}.log 2>&1
    echo "$site $arm rc=$? $(ls $o 2>/dev/null | head -1)"
  done
done
echo HAZARD_RENDER_DONE
