"""LOSO external-validity test: fine-tune on SCROLL 1 patches ONLY, so that
Scroll 4 is a genuinely unseen scroll at eval time.

Identical recipe to the published repair (500-patch pool, 3000 steps, 160^3
crops, batch 1, accum 2, AMP fp16, AdamW 1e-4 cosine, seed 1234) -- the ONLY
change is that the training pool is restricted to s1.
Question: does the compressed-region repair transfer to a scroll never trained on?
"""
import os, glob, sys, time, random, json
import numpy as np, tifffile, torch
import torch.nn.functional as F

FT = "/mnt/vesuvius/surf191_ft"
LBL_DIR = os.path.join(FT, "labelsTr")
OUTDIR = "/mnt/vesuvius/experiments/FT191"; os.makedirs(OUTDIR, exist_ok=True)
PS = 160; ACCUM = 2; EPOCHS = 6; ITERS = 500; LR = 1e-4
SEED = int(os.environ.get('SEED','1234'))
POOL = 500                     # match the published recipe's pool size
TRAIN_SCROLL = "s1_"           # held out: everything else (s4, s5)
ARM = "loso_s1"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loader059 import load_059

def crop_pair(img, lbl, rng, ps=PS):
    Z, Y, X = img.shape
    z = rng.integers(0, Z - ps + 1); y = rng.integers(0, Y - ps + 1)
    x = rng.integers(0, X - ps + 1)
    return img[z:z+ps, y:y+ps, x:x+ps], lbl[z:z+ps, y:y+ps, x:x+ps]

def main():
    random.seed(SEED); np.random.seed(SEED); torch.manual_seed(SEED)
    rng = np.random.default_rng(SEED)
    net, norm, props = load_059()
    net.train()
    mean, std = props["mean"], props["std"]
    imgs = sorted(glob.glob(os.path.join(FT, "imagesTr", "*_0000.tif")))
    imgs = [f for f in imgs if os.path.basename(f).startswith(TRAIN_SCROLL)]
    pairs = []
    for f in imgs:
        lf = os.path.join(LBL_DIR, os.path.basename(f).replace("_0000.tif", ".tif"))
        if os.path.exists(lf): pairs.append((f, lf))
    # seeded subsample to the published pool size
    sel = np.random.default_rng(SEED).permutation(len(pairs))[:POOL]
    pairs = [pairs[i] for i in sel]
    scrolls = sorted({os.path.basename(f).split("_")[0] for f, _ in pairs})
    print(f"ARM={ARM}  train scrolls={scrolls}  n_pairs={len(pairs)} "
          f"(HELD OUT: s4, s5)", flush=True)
    assert scrolls == ["s1"], f"training pool leaked non-s1 scrolls: {scrolls}"
    assert len(pairs) >= 400, "pool too small"

    opt = torch.optim.AdamW(net.parameters(), lr=LR, weight_decay=1e-5)
    total_steps = EPOCHS * ITERS
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=total_steps, eta_min=1e-6)
    scaler = torch.amp.GradScaler()
    step = 0; t0 = time.time(); run_loss = 0.0
    for ep in range(EPOCHS):
        for it in range(ITERS):
            for a in range(ACCUM):
                f, lf = pairs[rng.integers(0, len(pairs))]
                img = tifffile.imread(f); lbl = tifffile.imread(lf)
                ci, cl = crop_pair(img, lbl, rng)
                x = torch.from_numpy(((ci.astype(np.float32) - mean) / std)[None, None]).cuda()
                y = torch.from_numpy(cl.astype(np.int64))[None].cuda()
                with torch.amp.autocast("cuda", dtype=torch.float16):
                    out = net(x)
                    if isinstance(out, dict): out = out["surface"]
                    if isinstance(out, (list, tuple)): out = out[0]
                    out = out.float()
                    ce = F.cross_entropy(out, y.clamp(max=1), reduction="none")
                    valid = (y != 2).float()
                    ce = (ce * valid).sum() / valid.sum().clamp(min=1)
                    p1 = torch.softmax(out, 1)[:, 1]
                    t1 = (y == 1).float()
                    inter = (p1 * t1 * valid).sum()
                    dice = 1 - (2*inter + 1) / ((p1*valid).sum() + t1.sum() + 1)
                    loss = (ce + dice) / ACCUM
                scaler.scale(loss).backward()
                run_loss += float(loss)
            scaler.step(opt); scaler.update(); opt.zero_grad(); sched.step()
            step += 1
            if step % 100 == 0:
                print(f"[{ARM}] ep{ep} step {step}/{total_steps} loss={run_loss/100:.4f} "
                      f"lr={sched.get_last_lr()[0]:.2e} ({(time.time()-t0)/step:.1f}s/step)",
                      flush=True)
                run_loss = 0.0
    ck = {"model": net.state_dict(), "arm": ARM, "steps": total_steps,
          "train_scrolls": scrolls, "pool": len(pairs), "seed": SEED,
          "normalization_scheme": norm, "intensity_properties": props}
    path = os.path.join(OUTDIR, f"ckpt_ft_{ARM}_s{SEED}.pth")
    torch.save(ck, path)
    print(f"SAVED {path}", flush=True)
    print("LOSO FT COMPLETE", flush=True)

if __name__ == "__main__":
    main()
