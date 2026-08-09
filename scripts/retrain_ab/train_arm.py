"""Prereg build, day 3b: one arm, one seed, boundary-channel fine-tune.

The ft059/FT191 recipe with three swaps and nothing else changed:
  data    frozen crop list + local CT cache instead of Dataset059 tifs
  labels  the arm's instance tree -> 6-connected instance-contact boundary
  init    ckpt_ft_full.pth, identical for both arms

env: ARM=v1|v2  SEED=40..45  SMOKE=1 (200 steps, no ckpt overwrite protection)
The per-seed crop order is a seeded shuffle of the SAME frozen list, so seed s
in arm A and seed s in arm B see byte-identical crop sequences. The only
difference between paired runs is which label tree the reader opens.
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

sys.path.insert(0, "/mnt/vesuvius")
from loader059 import load_059

ARM = os.environ.get("ARM", "v1")
SEED = int(os.environ.get("SEED", "40"))
SMOKE = os.environ.get("SMOKE") == "1"

TREES = {"v1": Path("/mnt/vesuvius/p1218_full/blocks"),
         "v2": Path("/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired")}
AB = Path("/mnt/vesuvius/experiments/retrain_ab")
CACHE = AB / "ct_cache"
CKPT_INIT = "/mnt/vesuvius/experiments/FT191/ckpt_ft_full.pth"
OUT = AB / "ckpts"
OUT.mkdir(exist_ok=True)

PS = 160
ACCUM = 2
STEPS = 200 if SMOKE else 6000
LR = 1e-4


def boundary_target(lab: np.ndarray) -> np.ndarray:
    out = np.zeros(lab.shape, bool)
    fg = lab > 0
    for ax in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[ax] = slice(None, -1)
        b[ax] = slice(1, None)
        a, b = tuple(a), tuple(b)
        la, lb = lab[a], lab[b]
        diff = (la > 0) & (lb > 0) & (la != lb)
        out[a] |= diff
        out[b] |= diff
    return out & fg


def main() -> None:
    torch.manual_seed(SEED)
    rng = np.random.default_rng(SEED)
    crops = json.loads((AB / "frozen/crops.json").read_text())["crops"]
    order = rng.permutation(len(crops))

    net, norm, props = load_059()
    ck = torch.load(CKPT_INIT, map_location="cpu", weights_only=False)
    net.load_state_dict(ck["model"], strict=True)
    net.cuda().train()
    mean, std = props["mean"], props["std"]
    tree = TREES[ARM]  # recorded in the ckpt; targets come precomputed from it

    tcache = AB / f"target_cache_{ARM}"

    opt = torch.optim.AdamW(net.parameters(), lr=LR, weight_decay=1e-5)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(
        opt, T_max=STEPS, eta_min=1e-6)
    scaler = torch.amp.GradScaler()
    print(f"ARM={ARM} SEED={SEED} steps={STEPS} tree={tree}", flush=True)

    t0, run_loss, pos_seen = time.time(), 0.0, 0
    first_losses = []
    for step in range(1, STEPS + 1):
        for a in range(ACCUM):
            c = crops[order[(step * ACCUM + a) % len(crops)]]
            ct = np.load(CACHE / f"{c['key']}.npy")
            y = np.unpackbits(np.load(tcache / f"{c['key']}.npy"))[
                :PS ** 3].reshape(PS, PS, PS).astype(bool)
            pos_seen += int(y.sum())
            x = torch.from_numpy(
                ((ct.astype(np.float32) - mean) / std)[None, None]).cuda()
            t = torch.from_numpy(y.astype(np.int64))[None].cuda()
            with torch.amp.autocast("cuda", dtype=torch.float16):
                out = net(x)
                if isinstance(out, dict):
                    out = out["surface"]
                if isinstance(out, (list, tuple)):
                    out = out[0]
                out = out.float()
                ce = F.cross_entropy(out, t)
                p1 = torch.softmax(out, 1)[:, 1]
                t1 = (t == 1).float()
                dice = 1 - (2 * (p1 * t1).sum() + 1) / (p1.sum() + t1.sum() + 1)
                loss = (ce + dice) / ACCUM
            scaler.scale(loss).backward()
            run_loss += float(loss) * ACCUM
        scaler.step(opt)
        scaler.update()
        opt.zero_grad()
        sched.step()
        if step <= 5:
            first_losses.append(run_loss if step == 1 else None)
        if step % 50 == 0:
            print(f"[{ARM}/s{SEED}] {step}/{STEPS} loss={run_loss/50:.4f} "
                  f"pos/crop={pos_seen/(50*ACCUM):.0f} "
                  f"{(time.time()-t0)/step:.2f}s/step", flush=True)
            run_loss, pos_seen = 0.0, 0

    tag = f"{ARM}_s{SEED}" + ("_smoke" if SMOKE else "")
    torch.save({"model": net.state_dict(), "arm": ARM, "seed": SEED,
                "steps": STEPS, "init": CKPT_INIT,
                "normalization_scheme": norm, "intensity_properties": props},
               OUT / f"ckpt_{tag}.pth")
    print(f"SAVED ckpt_{tag}.pth  {(time.time()-t0)/60:.1f}min", flush=True)


if __name__ == "__main__":
    main()
