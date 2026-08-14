#!/bin/bash
# PR #1451 at-source support test: same 3 seeds, support rule ON, masked m7
# as the volume so the test is "keep a quad if CT-supported prediction exists".
set -u
P8=/mnt/vesuvius/p8_sprint; ROOT=/mnt/vesuvius/vcbuild
R=512; L0Z=23247; L0YX=7593
clamp(){ v=$1; lo=$2; hi=$3; [ "$v" -lt "$lo" ] && v=$lo; [ "$v" -gt "$hi" ] && v=$hi; echo "$v"; }
n=0
while read -r sz sy sx; do
  out=$P8/trace1451/s${n}; mkdir -p "$out"
  if [ -f "$out/meta.json" ]; then echo "skip s$n"; n=$((n+1)); continue; fi
  zmin=$(clamp $((sz-R)) 0 $L0Z); zmax=$(clamp $((sz+R)) 0 $L0Z)
  ymin=$(clamp $((sy-R)) 0 $L0YX); ymax=$(clamp $((sy+R)) 0 $L0YX)
  xmin=$(clamp $((sx-R)) 0 $L0YX); xmax=$(clamp $((sx+R)) 0 $L0YX)
  sed -e "s/@GENERATIONS@/60/" -e "s/@ZMIN@/$zmin/" -e "s/@ZMAX@/$zmax/" \
      -e "s/@YMIN@/$ymin/" -e "s/@YMAX@/$ymax/" -e "s/@XMIN@/$xmin/" -e "s/@XMAX@/$xmax/" \
      $ROOT/params_ab_A.template.json > "$out/params.json"
  ~/Vesuvius/.venv/bin/python3 - "$out/params.json" <<PY
import json,sys
p=sys.argv[1]; d=json.load(open(p))
d["require_volume_support"]=True; d["volume_support_threshold"]=0
json.dump(d,open(p,"w"),indent=1)
PY
  echo "trace s$n support-rule ON seed=($sx,$sy,$sz)"
  (cd $ROOT && timeout 2400 env OMP_NUM_THREADS=4 nice -n 12 ./vcrun1451.sh \
     -v $P8/r1_after.zarr -t "$out" --segment-name "s${n}_rule" \
     -p "$out/params.json" -s "$sx" "$sy" "$sz" --skip-overlap-check) \
     > $P8/logs/t1451_s${n}.log 2>&1
  echo "  rc=$? $([ -f "$out/meta.json" ] && echo metaOK || echo NOMETA) rejects=$(grep -c "support test rejected" $P8/logs/t1451_s${n}.log)"
  n=$((n+1))
done < "$P8/r1_seeds.txt"
echo RULE_TRACE_DONE
