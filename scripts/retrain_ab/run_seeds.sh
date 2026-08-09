#!/bin/bash
# Twelve inferential runs, alternating arms so a mid-sequence loss hits both
# arms evenly. Each run scores itself on the frozen eval lists right after
# training, so the seed table accumulates in scores/ as it goes.
cd /mnt/vesuvius/experiments/retrain_ab
for SEED in 40 41 42 43 44 45; do
  for ARM in v1 v2; do
    T="ckpts/ckpt_${ARM}_s${SEED}.pth"
    if [ ! -f "$T" ]; then
      env ARM=$ARM SEED=$SEED PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
        ~/Vesuvius/.venv/bin/python train_arm.py > "train_${ARM}_s${SEED}.log" 2>&1
    fi
    S="scores/ckpts_ckpt_${ARM}_s${SEED}.json"
    if [ -f "$T" ] && [ ! -f "$S" ]; then
      ~/Vesuvius/.venv/bin/python instrument.py "model:ckpts/ckpt_${ARM}_s${SEED}.pth" \
        > "score_${ARM}_s${SEED}.log" 2>&1
    fi
  done
done
echo "ALL_SEEDS_DONE" > seeds_done.marker
