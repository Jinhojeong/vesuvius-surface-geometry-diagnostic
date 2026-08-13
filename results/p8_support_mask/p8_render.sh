#!/bin/bash
set -u
ROOT=/mnt/vesuvius/vcbuild; P8=/mnt/vesuvius/p8_sprint
CT="https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr"
mkdir -p $P8/render
for i in 0 1 2; do for arm in before after; do
  o=$P8/render/s${i}_${arm}
  [ -f "$o/00000.tif" ] && { echo "skip s$i $arm"; continue; }
  mkdir -p "$o"
  (cd $ROOT && timeout 2400 env OMP_NUM_THREADS=4 nice -n 12 ./vcrun.sh vc_render_tifxyz \
     -v "$CT" --remote-url "$CT" --scale 0.5 -g 1 \
     -s $P8/trace/r1_s${i}_${arm} --tif-output "$o" --flatten --auto-crop \
     --cache-gb 8 --pyramid 0) > $P8/logs/render_s${i}_${arm}.log 2>&1
  echo "s$i $arm rc=$? $(ls $o 2>/dev/null | head -1)"
done; done
echo RENDER_DONE
