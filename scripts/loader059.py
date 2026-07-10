"""Validated loader + sliding-window inference for surface_recto_059_redo.
(Extracted from probe059.py as a proper module — complete functions.)"""
import os, json
import numpy as np
import torch

CKPT = "/mnt/vesuvius/models/surface_recto_059_redo/Model_epoch499.pth"

def load_059():
    from vesuvius.models.build.build_network_from_config import NetworkFromConfig
    ck = torch.load(CKPT, map_location="cpu", weights_only=False)
    mc = dict(ck["model_config"])
    mc["separate_decoders"] = True  # ckpt keys are task_decoders.surface.*
    print("normalization_scheme:", ck.get("normalization_scheme"))
    class M:
        def __init__(s, mc):
            s.model_config = mc; s.targets = mc.get("targets", {})
            s.train_patch_size = mc.get("train_patch_size", mc.get("patch_size"))
            s.train_batch_size = mc.get("train_batch_size", mc.get("batch_size", 2))
            s.in_channels = mc.get("in_channels", 1)
            s.autoconfigure = mc.get("autoconfigure", False)
            s.model_name = mc.get("model_name", "Model")
            s.spacing = [1]*len(s.train_patch_size)
    net = NetworkFromConfig(M(mc))
    def strip(k):
        for p in ("module.", "_orig_mod."):
            while k.startswith(p): k = k[len(p):]
        return k
    sd = {strip(k): v for k, v in ck["model"].items()}
    sd = {k: v for k, v in sd.items()
          if not k.startswith("task_decoders.surface.encoder.")}
    net.load_state_dict(sd, strict=True)
    print("loaded OK (duplicate per-task encoder keys dropped)")
    return net.cuda().eval(), ck.get("normalization_scheme"), ck.get("intensity_properties")

@torch.no_grad()
def predict(net, img, ps, norm, props):
    """Sliding-window softmax P(surface); ct-norm (global mean/std) or instance."""
    a = img.astype(np.float32)
    if norm and "instance" not in str(norm) and props and "mean" in props:
        a = (a - props["mean"]) / (props["std"] + 1e-6)
        inst = False
    else:
        inst = True
    Z, Y, X = a.shape
    prob = np.zeros((Z, Y, X), np.float32); wsum = np.zeros((Z, Y, X), np.float32)
    step = ps - ps//4
    ax = lambda L: sorted(set(list(range(0, max(L-ps, 0)+1, step)) + [max(L-ps, 0)]))
    for z0 in ax(Z):
        for y0 in ax(Y):
            for x0 in ax(X):
                cube = a[z0:z0+ps, y0:y0+ps, x0:x0+ps]
                pz, py, px = cube.shape
                if inst:
                    cube = (cube - cube.mean()) / (cube.std() + 1e-6)
                pad = [(0, ps-pz), (0, ps-py), (0, ps-px)]
                if any(p[1] for p in pad): cube = np.pad(cube, pad, mode="reflect")
                t = torch.from_numpy(cube[None, None]).cuda()
                out = net(t)
                if isinstance(out, dict): out = out.get("surface", list(out.values())[0])
                if isinstance(out, (list, tuple)): out = out[0]
                p = torch.softmax(out.float(), 1)[0, 1].cpu().numpy()[:pz, :py, :px]
                prob[z0:z0+pz, y0:y0+py, x0:x0+px] += p
                wsum[z0:z0+pz, y0:y0+py, x0:x0+px] += 1
    return prob / np.maximum(wsum, 1e-6)
