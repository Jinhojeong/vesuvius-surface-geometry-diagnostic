"""Build the 9um aligned input window for pherc0139-w016's validation region.

Documented recipe (vesuvius/docs/ink_detection.md, prepare_9um_isotropic_input):
read the 2.399um volume at XY pyramid level 2, select the centred 84 Z planes,
then 4x mean-pool in Z, giving 21 slices on the label grid.
"""
import json, os
import numpy as np, zarr

SRC=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/PHerc0139/segments/"
     "20250108000004-w029_2025010827/surface-volumes/2.399um-0.22m-78keV-volume-20260102150214.zarr")
LBL="/mnt/vesuvius/ink9um_labels/w016/pherc0139-w016_%s.zarr"
OUT="/mnt/vesuvius/ink9um_w016"
os.makedirs(OUT, exist_ok=True)
MARGIN=128

g=zarr.open(SRC, mode="r"); a=g["2"]
Z,Y,X=a.shape
z0=(Z-84)//2
print("level2", a.shape, "| centred 84 planes start at z", z0, flush=True)

val=zarr.open(LBL%"validation_mask", mode="r"); val=val["0"] if "0" in val else val
V=np.asarray(val)
nz=np.argwhere(V>0); lo=nz.min(0); hi=nz.max(0)+1
y0,y1=max(0,int(lo[1])-MARGIN), min(V.shape[1],int(hi[1])+MARGIN)
x0,x1=max(0,int(lo[2])-MARGIN), min(V.shape[2],int(hi[2])+MARGIN)
print("validation bbox y %d-%d x %d-%d | window with margin y %d-%d x %d-%d"
      %(lo[1],hi[1],lo[2],hi[2],y0,y1,x0,x1), flush=True)

raw=np.asarray(a[z0:z0+84, y0:y1, x0:x1]).astype(np.float32)
pooled=raw.reshape(21,4,raw.shape[1],raw.shape[2]).mean(1)
vol=np.clip(np.rint(pooled),0,255).astype(np.uint8)
print("pooled input", vol.shape, "nonzero %.3f"%float((vol>0).mean()), flush=True)

z=zarr.open(OUT+"/w016_val_input.zarr", mode="w", shape=vol.shape,
            chunks=(vol.shape[0],256,256), dtype=vol.dtype)
z[:]=vol

ink=zarr.open(LBL%"inklabels", mode="r"); ink=ink["0"] if "0" in ink else ink
sup=zarr.open(LBL%"supervision_mask", mode="r"); sup=sup["0"] if "0" in sup else sup
I=np.asarray(ink[:, y0:y1, x0:x1]); Vw=np.asarray(val[:, y0:y1, x0:x1]); S=np.asarray(sup[:, y0:y1, x0:x1])
np.savez_compressed(OUT+"/w016_val_labels.npz", ink=I, validation=Vw, supervision=S,
                    window=np.array([0,21,y0,y1,x0,x1]))

# alignment check: material presence under the annotated slice vs elsewhere
zc=int(lo[0]); m=Vw[zc]>0
print("annotated z:", zc, "| masked voxels in window:", int(m.sum()))
print("input material under mask: %.3f  vs whole slice: %.3f"
      %(float((vol[zc][m]>0).mean()), float((vol[zc]>0).mean())))
print("ink frac under mask: %.4f"%float((I[zc][m]>0).mean()))
json.dump(dict(level2_shape=[Z,Y,X], z_start=z0, window=[y0,y1,x0,x1],
               input_shape=list(vol.shape), annotated_z=zc,
               masked_voxels=int(m.sum()),
               ink_frac=round(float((I[zc][m]>0).mean()),5)),
          open(OUT+"/input_probe.json","w"), indent=1)
print("saved", OUT)
