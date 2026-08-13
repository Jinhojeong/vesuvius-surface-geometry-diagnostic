#!/bin/bash
set -u
ROOT=/mnt/vesuvius/vcbuild; P8=/mnt/vesuvius/p8_sprint; TAG=$1
R=512; L0Z=23247; L0YX=7593
clamp(){ v=$1; lo=$2; hi=$3; [ "$v" -lt "$lo" ] && v=$lo; [ "$v" -gt "$hi" ] && v=$hi; echo "$v"; }
n=0
while read -r sz sy sx; do
  for arm in before after; do
    out=$P8/trace/${TAG}_s${n}_${arm}; mkdir -p "$out"
    [ -f "$out/meta.json" ] && { echo "skip $n $arm"; continue; }
    zmin=$(clamp $((sz-R)) 0 $L0Z); zmax=$(clamp $((sz+R)) 0 $L0Z)
    ymin=$(clamp $((sy-R)) 0 $L0YX); ymax=$(clamp $((sy+R)) 0 $L0YX)
    xmin=$(clamp $((sx-R)) 0 $L0YX); xmax=$(clamp $((sx+R)) 0 $L0YX)
    sed -e "s/@GENERATIONS@/60/" -e "s/@ZMIN@/$zmin/" -e "s/@ZMAX@/$zmax/" \
        -e "s/@YMIN@/$ymin/" -e "s/@YMAX@/$ymax/" -e "s/@XMIN@/$xmin/" -e "s/@XMAX@/$xmax/" \
        $ROOT/params_ab_A.template.json > "$out/params.json"
    vol=$P8/${TAG}_before.zarr; [ "$arm" = after ] && vol=$P8/${TAG}_after.zarr
    echo "trace s$n $arm seed=($sx,$sy,$sz)"
    (cd "$ROOT" && timeout 2400 env OMP_NUM_THREADS=4 nice -n 12 ./vcrun.sh vc_grow_seg_from_seed \
       -v "$vol" -t "$out" --segment-name "${TAG}_s${n}_${arm}" -p "$out/params.json" \
       -s "$sx" "$sy" "$sz" --skip-overlap-check) > "$P8/logs/trace_${TAG}_s${n}_${arm}.log" 2>&1
    echo "  rc=$? $([ -f "$out/meta.json" ] && echo metaOK || echo NOMETA)"
  done
  n=$((n+1))
done < "$P8/${TAG}_seeds.txt"
echo TRACE_DONE
