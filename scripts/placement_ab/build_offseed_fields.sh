#!/bin/bash
# Amendment 1 (commit 083f424): build direction-field regions, same builder
# and parameters (make_direction_field.py, radius 256, band 6.0, sigma 1.0),
# centered on the SEVEN unique off-seed positions from
# offseeds_placement.json (o2/o4 share one point). Sequential, one region at
# a time (memory ~3.5 GB per region), nice 12, capped threads.
# Run: nohup bash /mnt/vesuvius/vcbuild/build_offseed_fields.sh \
#        > /dev/null 2>&1 < /dev/null &
set -u
PLACE=/mnt/vesuvius/vcbuild/demo_out/placement
LOG=$PLACE/dirfield_offseeds.log
PY=$HOME/Vesuvius/.venv/bin/python3
MK=/mnt/vesuvius/hazard_zarr_smoke/make_direction_field.py
echo "[$(date '+%F %T')] offseed field builds start (7 regions, r256 band6 sigma1)" >> "$LOG"
while read -r cz cy cx; do
  echo "[$(date '+%F %T')] build --center $cz $cy $cx" >> "$LOG"
  nice -n 12 env OMP_NUM_THREADS=4 "$PY" "$MK" \
      --center "$cz" "$cy" "$cx" --radius 256 --verify >> "$LOG" 2>&1 \
      || { echo "[$(date '+%F %T')] BUILD_FAILED $cz $cy $cx" >> "$LOG"; exit 1; }
done <<'CENTERS'
4240 1584 752
432 1840 1392
160 512 3488
11264 544 2672
9376 2016 3280
304 2016 2704
768 2288 2672
CENTERS
echo "[$(date '+%F %T')] ALL_BUILDS_DONE" >> "$LOG"
