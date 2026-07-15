"""#191 frontier diagnostic: is the compressed-bin (<4vox spacing) ceiling a
POST-PROCESS problem or a RESOLUTION problem?
At each GT surface voxel whose nearest neighbour sheet is <4vox away, sample the
PREDICTED prob field along the GT normal across both sheets:
  - 2 prob peaks (valley between) => sheets ARE resolved in prob, merged only
    after threshold => SPLITTABLE (a post-process could win Betti/topo).
  - 1 broad peak => genuinely UNRESOLVED => needs higher input resolution.
  - no peak (<0.5) => MISSED sheet.
Decides the next #191 research direction. Reuses the released 059 model.
"""
import os, glob, sys, json
import numpy as np, tifffile, torch
from scipy import ndimage as ndi
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loader059 import load_059, predict

OUT = "/mnt/vesuvius/surf191_rand"
FULL = "/mnt/vesuvius/experiments/FT191/ckpt_ft_loso_s1.pth"
NP = 36
SPAN = 6.0                    # +/- vox along normal
K = int(2*SPAN/0.5) + 1       # 0.5-vox sampling
COMP = 4.0                    # compressed if neighbour sheet < 4 vox

def gt_normals(gt):
    surf = gt > 0
    d_out = ndi.distance_transform_edt(~surf).astype(np.float32)
    d_in = ndi.distance_transform_edt(surf).astype(np.float32)
    sdf = ndi.gaussian_filter(d_out - d_in, 1.0)
    n = np.stack(np.gradient(sdf), 0).astype(np.float32)
    n /= (np.sqrt((n**2).sum(0)) + 1e-6)
    return surf, n

def runs(b):
    out = []; i = 0
    while i < len(b):
        if b[i]:
            j = i
            while j < len(b) and b[j]: j += 1
            out.append((i+j-1)/2.0); i = j
        else: i += 1
    return out

def taxonomy(net, norm, props, sel, imgs, rng):
    cats = {"bimodal_splittable": 0, "unimodal_unresolved": 0, "missed": 0, "other": 0}
    valley = []; step = SPAN*2/(K-1)
    offs = np.linspace(-SPAN, SPAN, K)
    for c, si in enumerate(sel):
        ip = imgs[si]; lp = ip.replace("imagesTr", "labelsTr").replace("_0000.tif", ".tif")
        img = tifffile.imread(ip); gt = tifffile.imread(lp).astype(np.uint8)
        prob = predict(net, img, 192, norm, props)
        surf, n = gt_normals(gt)
        pts = np.argwhere(surf)
        pts = pts[rng.permutation(len(pts))[:2000]]
        nv = n[:, pts[:, 0], pts[:, 1], pts[:, 2]].T
        coords = pts[:, None, :] + offs[None, :, None]*nv[:, None, :]
        cc = coords.reshape(-1, 3).T
        gl = ndi.map_coordinates(gt.astype(np.float32), cc, order=0, mode="constant").reshape(-1, K)
        pr = ndi.map_coordinates(prob, cc, order=1, mode="constant").reshape(-1, K)
        for m in range(len(pts)):
            rc = runs(gl[m] > 0.5)
            if len(rc) < 2: continue
            spac = min(abs(rc[k+1]-rc[k]) for k in range(len(rc)-1))*step
            if spac >= COMP: continue
            prs = ndi.uniform_filter1d(pr[m], 3)
            pk = [k for k in range(1, K-1) if prs[k] > 0.5 and prs[k] >= prs[k-1] and prs[k] >= prs[k+1]]
            mp = []
            for k in pk:
                if mp and k-mp[-1] <= 2: continue
                mp.append(k)
            if len(mp) >= 2:
                valley.append(min(prs[mp[0]], prs[mp[1]]) - prs[mp[0]:mp[1]+1].min())
                cats["bimodal_splittable"] += 1
            elif len(mp) == 1: cats["unimodal_unresolved"] += 1
            elif prs.max() < 0.5: cats["missed"] += 1
            else: cats["other"] += 1
        print(f"  [{c+1}/{len(sel)}] compressed voxels so far {sum(cats.values())}", flush=True)
    return cats, valley

def report(tag, cats, valley):
    tot = sum(cats.values())
    print(f"\n=== [{tag}] #191 COMPRESSED-BIN (<4vox) FAILURE TAXONOMY  (n={tot}) ===")
    for k, v in cats.items():
        print(f"  {k:22s} {v:6d}  {100*v/max(tot,1):5.1f}%")
    if valley:
        print(f"  bimodal valley depth mean {np.mean(valley):.3f}")
    return {"cats": cats, "n": tot, "valley_mean": float(np.mean(valley)) if valley else None,
            "bimodal_frac": cats["bimodal_splittable"]/max(tot, 1)}

def main():
    imgs = sorted(glob.glob(f"{OUT}/imagesTr/s4_*_0000.tif"))  # UNSEEN scroll
    sel = np.random.default_rng(0).choice(len(imgs), min(NP, len(imgs)), replace=False)
    out = {}
    # released baseline
    net, norm, props = load_059()
    print("### RELEASED 059 ###", flush=True)
    cB, vB = taxonomy(net, norm, props, sel, imgs, np.random.default_rng(0))
    out["released"] = report("RELEASED", cB, vB)
    # fine-tuned (ceiling model)
    netF, _, _ = load_059()
    ck = torch.load(FULL, map_location="cpu", weights_only=False)
    netF.load_state_dict(ck["model"], strict=True); netF.cuda().eval()
    print("\n### FINE-TUNED (ceiling) ###", flush=True)
    cF, vF = taxonomy(netF, norm, props, sel, imgs, np.random.default_rng(0))
    out["loso_s1"] = report("LOSO_S1 (unseen s4)", cF, vF)
    print("\nVERDICT: ceiling is",
          "RESOLUTION-limited (few bimodal even after FT)" if out["loso_s1"]["bimodal_frac"] < 0.2
          else "partly post-processable (FT resolves sheets into 2 peaks)")
    json.dump(out, open(f"{OUT}/bimodal_loso_s4.json", "w"), indent=1)
    print("BIMODAL DIAG COMPLETE")

if __name__ == "__main__":
    main()
