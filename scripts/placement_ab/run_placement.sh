#!/bin/bash
# run_placement.sh - prereg_placement_ab off-seed runs (frozen at 8531276,
# Amendment 1 at 083f424). 8 off-seeds x 2 arms x 5 fresh reps = 80 runs:
#   A = direction field only              (params_ab_A.template.json)
#   B = field + p1218_conf_v2_amp4.zarr   (template B, v1 -> v2_amp4 swap,
#       exactly the powered/powered16 construction)
# Same knobs as the powered block: GENERATIONS=60, RADIUS_L0=512,
# voxelsize 8.64 + min_area_cm 0.0 in the templates, MAXJOBS=3,
# OMP_NUM_THREADS=4, nice 12, timeout 2400s.
# Per the frozen exclusion rule, a run that crashes (no out/meta.json after
# the attempt) is retried ONCE with the same configuration; the retry is
# logged. Resumable: a run with out/meta.json is skipped.
# Run: nohup bash /mnt/vesuvius/vcbuild/run_placement.sh >> \
#        /mnt/vesuvius/vcbuild/demo_out/placement/nohup.out 2>&1 < /dev/null &
set -u
ROOT=/mnt/vesuvius/vcbuild
PLACE=$ROOT/demo_out/placement
SITES=$PLACE/off_sites.txt
FIELD=/mnt/vesuvius/hazard_zarr_smoke/m7_normals_L1.zarr
HAZ_B=/mnt/vesuvius/hazard_zarr/p1218_conf_v2_amp4.zarr
PY=$HOME/Vesuvius/.venv/bin/python
M7URL="https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/representations/predictions/surfaces/20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr"
MAXJOBS=3
GENERATIONS=60
RADIUS_L0=512
RUN_TIMEOUT=2400
L0_Z=23247; L0_YX=7593
mkdir -p "$PLACE"
LOG=$PLACE/placement.log
log() { echo "[$(date '+%F %T')] $*" >> "$LOG"; }

log "placement start pid=$$ (8 off-seeds x 2 arms x 5 reps, MAXJOBS=$MAXJOBS)"

# ------------------------------------------------------- gates (single shot)
[ -f "$HAZ_B/1/.zarray" ] || { log "FATAL: B zarr missing ($HAZ_B)"; exit 1; }
[ -f "$SITES" ] || { log "FATAL: off_sites.txt missing"; exit 1; }
[ -f "$PLACE/field_verify.json" ] || {
    log "FATAL: field_verify.json missing (Amendment 1 post-build voxel check must pass first)"; exit 1; }
grep -q '"all_written": true' "$PLACE/field_verify.json" || {
    log "FATAL: field_verify.json does not certify all off-seed voxels written"; exit 1; }
log "gates OPEN (B zarr, off_sites, field_verify all_written)"

clamp() {
    local v=$1
    [ "$v" -lt "$2" ] && v=$2
    [ "$v" -gt "$3" ] && v=$3
    echo "$v"
}

run_one() {
    local site=$1 arm=$2 rep=$3 sx=$4 sy=$5 sz=$6
    local sdir=$PLACE/$site/${arm}_r$rep
    [ -f "$sdir/out/meta.json" ] && { log "skip $site/$arm r$rep"; return 0; }
    local attempt
    for attempt in 1 2; do
        rm -rf "$sdir/out"
        mkdir -p "$sdir/out"
        cp "$PLACE/$site/params_$arm.json" "$sdir/params.json"
        local t0=$(date +%s)
        (cd "$ROOT" && timeout $RUN_TIMEOUT env OMP_NUM_THREADS=4 nice -n 12 \
            ./vcrun.sh vc_grow_seg_from_seed \
            -v "$M7URL" -t "$sdir/out" --segment-name "${site}_${arm}_r${rep}" \
            -p "$sdir/params.json" -s "$sx" "$sy" "$sz" --skip-overlap-check) \
            > "$sdir/run.log" 2>&1
        local rc=$?
        log "done $site/$arm r$rep attempt=$attempt rc=$rc $(( $(date +%s) - t0 ))s $([ -f "$sdir/out/meta.json" ] && echo meta_OK || echo NO_META)"
        [ -f "$sdir/out/meta.json" ] && return 0
        [ "$attempt" = 1 ] && log "RETRY $site/$arm r$rep (crashed, retrying once per frozen exclusion rule)"
    done
    log "EXCLUDE-CANDIDATE $site/$arm r$rep failed twice"
}

# ------------------------------------------------------- params + run loop
while read -r i slab tile gz gy gx; do
    site="o${i}_${slab}_${tile}"
    sx=$((gx * 2)); sy=$((gy * 2)); sz=$((gz * 2))
    sdir=$PLACE/$site
    mkdir -p "$sdir"
    zmin=$(clamp $((sz - RADIUS_L0)) 0 $L0_Z);  zmax=$(clamp $((sz + RADIUS_L0)) 0 $L0_Z)
    ymin=$(clamp $((sy - RADIUS_L0)) 0 $L0_YX); ymax=$(clamp $((sy + RADIUS_L0)) 0 $L0_YX)
    xmin=$(clamp $((sx - RADIUS_L0)) 0 $L0_YX); xmax=$(clamp $((sx + RADIUS_L0)) 0 $L0_YX)
    for cfg in A B; do
        sed -e "s/@GENERATIONS@/$GENERATIONS/" \
            -e "s/@ZMIN@/$zmin/" -e "s/@ZMAX@/$zmax/" \
            -e "s/@YMIN@/$ymin/" -e "s/@YMAX@/$ymax/" \
            -e "s/@XMIN@/$xmin/" -e "s/@XMAX@/$xmax/" \
            "$ROOT/params_ab_${cfg}.template.json" > "$sdir/params_$cfg.tmp"
    done
    mv "$sdir/params_A.tmp" "$sdir/params_A.json"
    grep -q "weight_zarr" "$sdir/params_A.json" && { log "FATAL A params has weight_zarr $site"; exit 1; }
    sed "s#/mnt/vesuvius/hazard_zarr/p1218_conf_v1\.zarr#$HAZ_B#" \
        "$sdir/params_B.tmp" > "$sdir/params_B.json"
    grep -q "p1218_conf_v2_amp4\.zarr" "$sdir/params_B.json" || { log "FATAL bad B params $site"; exit 1; }
    rm -f "$sdir/params_B.tmp"
    log "params $site seed=($sx,$sy,$sz) bbox z[$zmin,$zmax] y[$ymin,$ymax] x[$xmin,$xmax]"
    for rep in 1 2 3 4 5; do
        for arm in A B; do
            while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n; done
            run_one "$site" "$arm" "$rep" "$sx" "$sy" "$sz" &
        done
    done
done < "$SITES"
wait
log "placement runs complete"

# ------------------------------------------------------- scoring
log "scoring"
nice -n 12 "$PY" "$ROOT/placement_score.py" >> "$LOG" 2>&1
log "scorer rc=$?"
log "CHAIN_COMPLETE"
