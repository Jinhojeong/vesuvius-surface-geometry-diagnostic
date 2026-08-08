"""Two-sided harness control for the leg-2 range check.

Before I tell TAUIL that leg 2 cannot cast an independent vote, this script has to
show that the conclusion is the checkpoint's and not my code's. It runs BOTH models
through the SAME predict-and-score path used on the 892, over the 80 Dataset059
patches that produced the July numbers:

  surface_recto_3dunet     July: pooled dice 0.23, precision ceiling ~0.22
  surface_recto_059_redo   July: pooled dice 0.81 @th0.4

If the negative reproduces and the positive reproduces, the path is sound in both
directions and the 892 result stands. If the positive fails, the 892 result is about
my code and has to be thrown away.

Estimators and normalisation are identical to leg2_892.py. Both checkpoints declare
z-score with no intensity properties, so both get per-volume z-scoring.
"""
from __future__ import annotations
import glob, json, os, sys
from pathlib import Path

import numpy as np
import tifffile
import torch

PATCHES = "/mnt/vesuvius/surf191"
NBIN = 2000
MODELS = {
    "surface_recto_3dunet": ("/mnt/vesuvius/models/surface_recto_3dunet/"
                             "checkpoint_inference_ready.pth"),
    "surface_recto_059_redo": ("/mnt/vesuvius/models/surface_recto_059_redo/"
                               "Model_epoch499.pth"),
}


def load(ckpt_path: str):
    """Both checkpoints instantiate through NetworkFromConfig. The 059 one stores a
    duplicate per-task encoder that the current builder does not build, which is why
    it needs separate_decoders and the task-encoder keys dropped. That recipe is
    probe059.py's and is not being changed here."""
    from vesuvius.models.build.build_network_from_config import NetworkFromConfig
    ck = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    mc = dict(ck.get("model_config") or {})
    if not mc:
        mc = json.load(open(os.path.join(os.path.dirname(ckpt_path), "config.json")))
    has_task_enc = any(k.replace("_orig_mod.", "").startswith(
        "task_decoders.") for k in ck.get("model", {}))
    if has_task_enc:
        mc["separate_decoders"] = True

    class CM:
        def __init__(s, mc):
            s.model_config = mc
            s.targets = mc.get("targets", {})
            s.train_patch_size = mc.get("train_patch_size", mc.get("patch_size"))
            s.train_batch_size = mc.get("train_batch_size", mc.get("batch_size", 2))
            s.in_channels = mc.get("in_channels", 1)
            s.autoconfigure = mc.get("autoconfigure", False)
            s.model_name = mc.get("model_name", "Model")
            s.spacing = [1] * len(s.train_patch_size)

    net = NetworkFromConfig(CM(mc))

    def strip(k):
        for p in ("module.", "_orig_mod."):
            while k.startswith(p):
                k = k[len(p):]
        return k

    sd = {strip(k): v for k, v in (ck.get("model", ck)).items()}
    dropped = [k for k in sd if k.startswith("task_decoders.surface.encoder.")]
    for k in dropped:
        del sd[k]
    net.load_state_dict(sd, strict=True)
    ps = mc.get("patch_size", [128, 128, 128])
    return net.cuda().eval(), int(ps[0]), ck.get("normalization_scheme"), len(dropped)


@torch.no_grad()
def predict(net, vol: np.ndarray, ps: int) -> np.ndarray:
    """Same normalisation, same softmax channel, same uniform blending as the 892 run.
    Sliding window with 50% step because these patches are 300^3 and the declared
    patch can be smaller."""
    a = vol.astype(np.float32)
    a = (a - a.mean()) / (a.std() + 1e-6)
    Z, Y, X = a.shape
    if ps >= max(Z, Y, X):
        pad = [(0, max(ps - Z, 0)), (0, max(ps - Y, 0)), (0, max(ps - X, 0))]
        b = np.pad(a, pad, mode="reflect") if any(p[1] for p in pad) else a
        t = torch.from_numpy(b[None, None]).cuda()
        with torch.autocast("cuda", torch.float16):
            out = net(t)
        if isinstance(out, dict):
            out = out.get("surface", list(out.values())[0])
        if isinstance(out, (list, tuple)):
            out = out[0]
        return torch.softmax(out.float(), 1)[0, 1].cpu().numpy()[:Z, :Y, :X]
    prob = np.zeros((Z, Y, X), np.float32)
    wsum = np.zeros((Z, Y, X), np.float32)
    step = ps // 2

    def starts(n):
        return sorted(set(list(range(0, n - ps + 1, step)) + [n - ps]))

    for z0 in starts(Z):
        for y0 in starts(Y):
            for x0 in starts(X):
                cube = a[z0:z0 + ps, y0:y0 + ps, x0:x0 + ps]
                t = torch.from_numpy(cube[None, None]).cuda()
                with torch.autocast("cuda", torch.float16):
                    out = net(t)
                if isinstance(out, dict):
                    out = out.get("surface", list(out.values())[0])
                if isinstance(out, (list, tuple)):
                    out = out[0]
                p = torch.softmax(out.float(), 1)[0, 1].cpu().numpy()
                prob[z0:z0 + ps, y0:y0 + ps, x0:x0 + ps] += p
                wsum[z0:z0 + ps, y0:y0 + ps, x0:x0 + ps] += 1.0
    return prob / np.maximum(wsum, 1e-6)


sys.path.insert(0, "/mnt/vesuvius/experiments/leg2")
from leg2_finalize import auc_from_hist, curve  # noqa: E402  the same estimators


def main() -> None:
    imgs = sorted(glob.glob(os.path.join(PATCHES, "imagesTr", "*_0000.tif")))
    print(f"MARKER patches {len(imgs)}", flush=True)
    out = {}
    for name, path in MODELS.items():
        net, ps, norm, dropped = load(path)
        print(f"MARKER {name} loaded ps={ps} norm={norm} dropped_task_enc={dropped}",
              flush=True)
        agg_s = np.zeros(NBIN)
        agg_t = np.zeros(NBIN)
        rows = []
        for f in imgs:
            lf = os.path.join(PATCHES, "labelsTr",
                              os.path.basename(f).replace("_0000.tif", ".tif"))
            vol = np.asarray(tifffile.imread(f))
            lab = np.asarray(tifffile.imread(lf))
            p = predict(net, vol, ps)
            scored = lab != 2
            sheet = lab == 1
            if not sheet.sum():
                continue
            hs = np.histogram(p[sheet], bins=NBIN, range=(0.0, 1.0))[0]
            ht = np.histogram(p[scored], bins=NBIN, range=(0.0, 1.0))[0]
            agg_s += hs
            agg_t += ht
            c = curve(hs, ht)
            rec, prec, dice, tp, tt, P = c
            i = int(np.nanargmax(dice))
            rows.append({"patch": os.path.basename(f),
                         "auc": auc_from_hist(hs, ht),
                         "dice_best": float(dice[i]),
                         "base_rate": float(sheet.sum() / max(scored.sum(), 1))})
        c = curve(agg_s, agg_t)
        rec, prec, dice, tp, tt, P = c
        i = int(np.nanargmax(dice))
        ok = rec >= 0.05
        out[name] = {
            "n_patches": len(rows),
            "pooled_auc": auc_from_hist(agg_s, agg_t),
            "pooled_dice_best": float(dice[i]),
            "pooled_dice_best_threshold": round(i / NBIN, 4),
            "pooled_dice_at_0.40": float(dice[int(0.40 * NBIN)]),
            "pooled_dice_at_0.50": float(dice[int(0.50 * NBIN)]),
            "pooled_precision_ceiling_at_recall_05": float(np.nanmax(prec[ok])),
            "pooled_base_rate": float(P / agg_t.sum()),
            "median_auc": round(float(np.median([r["auc"] for r in rows])), 5),
            "median_dice_best": round(float(np.median([r["dice_best"] for r in rows])), 5),
        }
        print(f"MARKER {name} {json.dumps(out[name])}", flush=True)
        del net
        torch.cuda.empty_cache()
    Path("/mnt/vesuvius/experiments/leg2/control_059.json").write_text(
        json.dumps({"what": "two-sided harness control on the 80 Dataset059 patches, "
                            "same predict-and-score path as the 892 run",
                    "july_reference": {"surface_recto_3dunet": "pooled dice 0.23, "
                                       "precision ceiling ~0.22",
                                       "surface_recto_059_redo": "pooled dice 0.81 "
                                       "@th0.4"},
                    "measured": out}, indent=1))
    print("MARKER STATUS=DONE", flush=True)


if __name__ == "__main__":
    main()
