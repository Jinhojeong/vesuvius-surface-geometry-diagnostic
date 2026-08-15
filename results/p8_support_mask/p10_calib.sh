#!/bin/bash
# Wait for the v2 build, then run the dilation calibration axiosdevs asked for:
# 3 seeds x dilation {0,1,2,4}, support_volume = the masked-CT store we built.
set -u
P8=/mnt/vesuvius/p8_sprint; ROOT=/mnt/vesuvius/vcbuild
B=$ROOT/vc3d-rootfs/src/build/bin/vc_grow_seg_from_seed
R=512; L0Z=23247; L0YX=7593
# wait for ninja to finish and the binary to carry the v2 symbol
while pgrep -f "ninja -C /src/build" >/dev/null; do sleep 60; done
strings $B | grep -q "volume_support_dilation\|support volume" || { echo "V2 SYMBOL MISSING, abort"; exit 1; }
echo "[$(date +%H:%M)] v2 build done, starting calibration"
clamp(){ v=$1; lo=$2; hi=$3; [ "$v" -lt "$lo" ] && v=$lo; [ "$v" -gt "$hi" ] && v=$hi; echo "$v"; }
for d in 0 1 2 4; do
  n=0
  while read -r sz sy sx; do
    out=$P8/calib1451/d${d}_s${n}; mkdir -p "$out"
    if [ -f "$out/meta.json" ]; then echo "skip d$d s$n"; n=$((n+1)); continue; fi
    zmin=$(clamp $((sz-R)) 0 $L0Z); zmax=$(clamp $((sz+R)) 0 $L0Z)
    ymin=$(clamp $((sy-R)) 0 $L0YX); ymax=$(clamp $((sy+R)) 0 $L0YX)
    xmin=$(clamp $((sx-R)) 0 $L0YX); xmax=$(clamp $((sx+R)) 0 $L0YX)
    sed -e "s/@GENERATIONS@/60/" -e "s/@ZMIN@/$zmin/" -e "s/@ZMAX@/$zmax/" \
        -e "s/@YMIN@/$ymin/" -e "s/@YMAX@/$ymax/" -e "s/@XMIN@/$xmin/" -e "s/@XMAX@/$xmax/" \
        $ROOT/params_ab_A.template.json > "$out/params.json"
    ~/Vesuvius/.venv/bin/python3 - "$out/params.json" "$d" <<PY
import json,sys
p=sys.argv[1]; d=json.load(open(p))
d["require_volume_support"]=True
d["volume_support_threshold"]=0
d["support_volume"]="/mnt/vesuvius/p8_sprint/r1_support.zarr"
d["volume_support_dilation"]=int(sys.argv[2])
json.dump(d,open(p,"w"),indent=1)
PY
    echo "[$(date +%H:%M)] trace d=$d s$n seed=($sx,$sy,$sz)"
    (cd $ROOT && timeout 2400 env OMP_NUM_THREADS=4 nice -n 12 ./vcrun1451.sh \
       -v $P8/r1_after.zarr -t "$out" --segment-name "d${d}_s${n}" \
       -p "$out/params.json" -s "$sx" "$sy" "$sz" --skip-overlap-check) \
       > $P8/logs/calib_d${d}_s${n}.log 2>&1
    echo "  rc=$? $([ -f "$out/meta.json" ] && echo metaOK || echo NOMETA)"
    n=$((n+1))
  done < "$P8/r1_seeds.txt"
done
echo CALIB_DONE
