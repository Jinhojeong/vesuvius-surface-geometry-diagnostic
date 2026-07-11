"""A-arm: fine-tune with the ORIGINAL combined loss (DC+SkelREC+CE) — disambiguates
"repair came from extra training" vs "repair came from dropping the srec term".

Stage 1 (SKEL env unset): precompute tubed medial-surface skeletons per patch.
Stage 2 (SKEL=trained): fine-tune, same seed/steps/crops as ft059.py ctrl arm.
"""
import os, glob, sys, time
import numpy as np, tifffile, torch
from skimage.morphology import skeletonize, dilation, closing
from scipy import ndimage as ndi

FT = "/mnt/vesuvius/surf191_ft"
SKELDIR = os.path.join(FT, "skelsTr")
OUTDIR = "/mnt/vesuvius/experiments/FT191"
PS = 160; ACCUM = 2; EPOCHS = 6; ITERS = 500; LR = 1e-4; SEED = 1234

sys.path.insert(0, "/mnt/vesuvius")

def make_skels():
    os.makedirs(SKELDIR, exist_ok=True)
    lbls = sorted(glob.glob(os.path.join(FT, "labelsTr", "*.tif")))
    ball = np.ones((3, 3, 3), bool)
    for i, lf in enumerate(lbls):
        of = os.path.join(SKELDIR, os.path.basename(lf))
        if os.path.exists(of): continue
        t0 = time.time()
        gt = tifffile.imread(lf) > 0
        sk = skeletonize(closing(gt, ball))
        tube = dilation(dilation(sk, ball), ball)        # ~2px tube
        tube = (tube & ~(~gt)).astype(np.uint8) if False else tube.astype(np.uint8)
        tifffile.imwrite(of, tube, compression="zlib")
        if (i + 1) % 25 == 0:
            print(f"[skel {i+1}/{len(lbls)}] ({time.time()-t0:.0f}s/patch)", flush=True)
    print("SKELS DONE", flush=True)

def train():
    from loader059 import load_059
    from vesuvius.models.training.loss.skeleton_recall import DC_SkelREC_and_CE_loss
    rng = np.random.default_rng(SEED); torch.manual_seed(SEED)
    net, norm, props = load_059(); net.train()
    mean, std = props["mean"], props["std"]
    loss_fn = DC_SkelREC_and_CE_loss(
        soft_dice_kwargs={"batch_dice": False, "smooth": 1e-5, "do_bg": False, "ddp": False},
        soft_skelrec_kwargs={"batch_dice": False, "smooth": 1e-5, "do_bg": False, "ddp": False},
        ce_kwargs={}, weight_ce=1, weight_dice=1, weight_srec=1)
    imgs = sorted(glob.glob(os.path.join(FT, "imagesTr", "*_0000.tif")))
    pairs = []
    for f in imgs:
        b = os.path.basename(f).replace("_0000.tif", ".tif")
        lf = os.path.join(FT, "labelsTr", b); sf = os.path.join(SKELDIR, b)
        if os.path.exists(lf) and os.path.exists(sf): pairs.append((f, lf, sf))
    print(f"n_pairs={len(pairs)}", flush=True)
    assert len(pairs) >= 400
    opt = torch.optim.AdamW(net.parameters(), lr=LR, weight_decay=1e-5)
    total = EPOCHS * ITERS
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=total, eta_min=1e-6)
    scaler = torch.amp.GradScaler()
    step = 0; t0 = time.time(); run = 0.0
    for ep in range(EPOCHS):
        for it in range(ITERS):
            for a in range(ACCUM):
                f, lf, sf = pairs[rng.integers(0, len(pairs))]
                img = tifffile.imread(f); lbl = tifffile.imread(lf); sk = tifffile.imread(sf)
                Z, Y, X = img.shape
                z = rng.integers(0, Z-PS+1); y = rng.integers(0, Y-PS+1); x = rng.integers(0, X-PS+1)
                ci = img[z:z+PS, y:y+PS, x:x+PS]
                cl = lbl[z:z+PS, y:y+PS, x:x+PS]
                cs = sk[z:z+PS, y:y+PS, x:x+PS]
                xb = torch.from_numpy(((ci.astype(np.float32)-mean)/std)[None, None]).cuda()
                yb = torch.from_numpy(cl.astype(np.int64))[None, None].cuda()
                sb = torch.from_numpy(cs.astype(np.float32))[None, None].cuda()
                with torch.amp.autocast("cuda", dtype=torch.float16):
                    out = net(xb)
                    if isinstance(out, dict): out = out["surface"]
                    if isinstance(out, (list, tuple)): out = out[0]
                    l = loss_fn(out.float(), yb, sb) / ACCUM
                scaler.scale(l).backward()
                run += float(l) * ACCUM
            scaler.step(opt); scaler.update(); opt.zero_grad(); sched.step()
            step += 1
            if step % 50 == 0:
                print(f"[msr] ep{ep} step {step}/{total} loss={run/50/ACCUM:.4f} "
                      f"({(time.time()-t0)/step:.1f}s/step)", flush=True)
                run = 0.0
    ck = {"model": net.state_dict(), "arm": "msr", "steps": total,
          "normalization_scheme": norm, "intensity_properties": props}
    torch.save(ck, os.path.join(OUTDIR, "ckpt_ft_msr.pth"))
    print("SAVED ckpt_ft_msr.pth", flush=True)

if __name__ == "__main__":
    if os.environ.get("STAGE") == "train":
        train()
    else:
        make_skels()
