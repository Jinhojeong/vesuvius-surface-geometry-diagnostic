#!/bin/bash
# run_powered16.sh - extension of the powered hazard A/B to 8 NEW sites.
# 8 new sites (indices 8-15 of demo_sites16.json, one per tile, no tile
# shared with the original 8) x 3 arms x 5 fresh reps = 120 runs:
#   A = direction field only              (params_ab_A.template.json)
#   B = field + p1218_conf_v2_amp4.zarr   (template B, v1 -> v2_amp4 swap,
#       exactly the ab2/powered construction)
#   C = field + p1218_conf_rand_amp4_ext8.zarr (random placebo, same
#       census clusters redrawn uniformly in the NEW region cubes,
#       radius 288 L1, rng 20260812; analogous to randctl's C)
# Same knobs as rounds 1/2/powered/randctl: GENERATIONS=60, RADIUS_L0=512,
# voxelsize 8.64 + min_area_cm 0.0 in the templates, MAXJOBS=3,
# OMP_NUM_THREADS=4, nice 12, timeout 2400s.
# Resumable: a run with out/meta.json is skipped.
# Run: nohup bash /mnt/vesuvius/vcbuild/run_powered16.sh >> \
#        /mnt/vesuvius/vcbuild/demo_out/powered16/nohup.out 2>&1 < /dev/null &
set -u
ROOT=/mnt/vesuvius/vcbuild
PW16=$ROOT/demo_out/powered16
SITES=/mnt/vesuvius/hazard_zarr_smoke/demo_sites16.json
FIELD=/mnt/vesuvius/hazard_zarr_smoke/m7_normals_L1.zarr
HAZ_B=/mnt/vesuvius/hazard_zarr/p1218_conf_v2_amp4.zarr
HAZ_C=/mnt/vesuvius/hazard_zarr/p1218_conf_rand_amp4_ext8.zarr
RANDLOG=/mnt/vesuvius/hazard_zarr/rand_ext8_build.log
PY=$HOME/Vesuvius/.venv/bin/python
M7URL="https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/representations/predictions/surfaces/20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr"
MAXJOBS=3
GENERATIONS=60
RADIUS_L0=512
RUN_TIMEOUT=2400
L0_Z=23247; L0_YX=7593
mkdir -p "$PW16"
LOG=$PW16/powered16.log
log() { echo "[$(date '+%F %T')] $*" >> "$LOG"; }

log "powered16 start pid=$$ (8 new sites x 3 arms x 5 reps, MAXJOBS=$MAXJOBS)"

# ------------------------------------------------------- gates (single shot)
[ -f "$HAZ_B/1/.zarray" ] || { log "FATAL: B zarr missing ($HAZ_B)"; exit 1; }
grep -q "\[done\].*p1218_conf_rand_amp4_ext8\.zarr" "$RANDLOG" 2>/dev/null || {
    log "FATAL: rand_ext8 build has no [done] marker"; exit 1; }
[ -f "$HAZ_C/1/.zarray" ] || { log "FATAL: C zarr missing ($HAZ_C)"; exit 1; }

"$PY" - "$SITES" "$FIELD" <<'PYEOF' >> "$LOG" 2>&1
import json, os, sys
sites = json.load(open(sys.argv[1])); base = sys.argv[2]
R, C = 256, 128
missing = []
for i in range(8, 16):
    s = sites[i]
    ok = True
    for comp in "xyz":
        found = False
        for cz in range(max(0, s["gz"] - R) // C, (s["gz"] + R - 1) // C + 1):
            for cy in range(max(0, s["gy"] - R) // C, (s["gy"] + R - 1) // C + 1):
                for cx in range(max(0, s["gx"] - R) // C, (s["gx"] + R - 1) // C + 1):
                    if os.path.isfile(os.path.join(base, comp, "1",
                                                   str(cz), str(cy), str(cx))):
                        found = True
                        break
                if found: break
            if found: break
        ok = ok and found
    if not ok:
        missing.append(i)
if missing:
    print("FATAL: new sites missing field chunks: %s" % missing)
    sys.exit(1)
print("field check OK: sites 8-15 have x/y/z chunks")
PYEOF
[ $? -eq 0 ] || { log "FATAL: direction-field check failed"; exit 1; }
log "gates OPEN (B zarr, C ext8 zarr, field chunks 8-15)"

clamp() {
    local v=$1
    [ "$v" -lt "$2" ] && v=$2
    [ "$v" -gt "$3" ] && v=$3
    echo "$v"
}

run_one() {
    local site=$1 arm=$2 rep=$3 sx=$4 sy=$5 sz=$6
    local sdir=$PW16/$site/${arm}_r$rep
    [ -f "$sdir/out/meta.json" ] && { log "skip $site/$arm r$rep"; return 0; }
    mkdir -p "$sdir/out"
    cp "$PW16/$site/params_$arm.json" "$sdir/params.json"
    local t0=$(date +%s)
    (cd "$ROOT" && timeout $RUN_TIMEOUT env OMP_NUM_THREADS=4 nice -n 12 \
        ./vcrun.sh vc_grow_seg_from_seed \
        -v "$M7URL" -t "$sdir/out" --segment-name "${site}_${arm}_r${rep}" \
        -p "$sdir/params.json" -s "$sx" "$sy" "$sz" --skip-overlap-check) \
        > "$sdir/run.log" 2>&1
    local rc=$?
    log "done $site/$arm r$rep rc=$rc $(( $(date +%s) - t0 ))s $([ -f "$sdir/out/meta.json" ] && echo meta_OK || echo NO_META)"
}

# ------------------------------------------------------- params + run loop
while read -r i slab tile gz gy gx; do
    site="s${i}_${slab}_${tile}"
    sx=$((gx * 2)); sy=$((gy * 2)); sz=$((gz * 2))
    sdir=$PW16/$site
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
    sed "s#/mnt/vesuvius/hazard_zarr/p1218_conf_v1\.zarr#$HAZ_C#" \
        "$sdir/params_B.tmp" > "$sdir/params_C.json"
    grep -q "p1218_conf_rand_amp4_ext8\.zarr" "$sdir/params_C.json" || { log "FATAL bad C params $site"; exit 1; }
    rm -f "$sdir/params_B.tmp"
    log "params $site seed=($sx,$sy,$sz) bbox z[$zmin,$zmax] y[$ymin,$ymax] x[$xmin,$xmax]"
    for rep in 1 2 3 4 5; do
        for arm in A B C; do
            while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n; done
            run_one "$site" "$arm" "$rep" "$sx" "$sy" "$sz" &
        done
    done
done < "$PW16/new_sites.txt"
wait
log "powered16 runs complete"

# ------------------------------------------------------- scoring
log "scoring"
nice -n 12 "$PY" "$ROOT/powered16_score.py" >> "$LOG" 2>&1
log "scorer rc=$?"
log "CHAIN_COMPLETE"
