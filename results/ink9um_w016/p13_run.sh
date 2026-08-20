#!/bin/bash
# Arm 1 replicated on the two additional validation segments.
set -u
R=/mnt/vesuvius/ink9um
W=/mnt/vesuvius/ink9um_w016
export PYTHONPATH=$R/villa-main/vesuvius/src
for seg in pherc0814-46527 pherc1667-w029; do
  IN=$W/${seg}_val_input.zarr
  OUT=$W/preds_$seg
  mkdir -p "$OUT"
  for ck in "$R"/checkpoints/*.pth; do
    n=$(basename "$ck" .pth)
    o="$OUT/${n}.tif"
    [ -f "$o" ] && { echo "skip $seg $n"; continue; }
    nice -n 12 ~/Vesuvius/.venv/bin/python3 -m vesuvius.ink_detection.inference.infer \
      "$IN" "$ck" "$o" --overlap 0.5 --blend-mode hann --batch-size 2 \
      --direction forward --no-compile > "$R/logs/${seg}_${n}.log" 2>&1
    echo "$seg $n rc=$?"
  done
  echo "$seg done: $(ls $OUT/*.tif 2>/dev/null | wc -l) maps"
done
echo ARM1x2_DONE
