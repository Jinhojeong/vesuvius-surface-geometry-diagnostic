"""Build 9um inputs for the two additional validation segments, same documented recipe."""
import json, os, sys
import numpy as np, zarr

VOL={
 "pherc0814-46527": ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/PHerc0814/segments/"
   "20260226000000-46527_2um_try2/surface-volumes/2.399um-0.22m-78keV-volume-20260309142202.zarr"),
 "pherc1667-w029": ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/PHerc1667/segments/"
   "20251212185248-w029_20251212185248662_flatboi/surface-volumes/2.399um-0.22m-78keV-volume-20251217075048.zarr"),
}
MARGIN=128; OUT="/mnt/vesuvius/ink9um_w016"
res={}
for seg,url in VOL.items():
    d="/mnt/vesuvius/ink9um_labels/"+seg
    L=lambda k: (lambda z: z["0"] if "0" in z else z)(zarr.open("%s/%s_%s.zarr"%(d,seg,k),mode="r"))
    val=L("validation_mask"); V=np.asarray(val)
    nz=np.argwhere(V>0); lo,hi=nz.min(0),nz.max(0)+1
    y0,y1=max(0,int(lo[1])-MARGIN),min(V.shape[1],int(hi[1])+MARGIN)
    x0,x1=max(0,int(lo[2])-MARGIN),min(V.shape[2],int(hi[2])+MARGIN)
    g=zarr.open(url,mode="r"); a=g["2"]
    Z=a.shape[0]; z0=(Z-84)//2
    if z0<0: print(seg,"SKIP: level2 z",Z,"< 84"); continue
    print(seg,"level2",a.shape,"z start",z0,"| window y%d-%d x%d-%d"%(y0,y1,x0,x1),flush=True)
    if a.shape[1]!=V.shape[1] or a.shape[2]!=V.shape[2]:
        print("   GRID MISMATCH: volume %s vs labels %s"%(a.shape[1:],V.shape[1:]))
        res[seg]=dict(error="grid mismatch", vol=list(a.shape), lab=list(V.shape)); continue
    raw=np.asarray(a[z0:z0+84, y0:y1, x0:x1]).astype(np.float32)
    vol=np.clip(np.rint(raw.reshape(21,4,raw.shape[1],raw.shape[2]).mean(1)),0,255).astype(np.uint8)
    p=OUT+"/%s_val_input.zarr"%seg
    z=zarr.open(p,mode="w",shape=vol.shape,chunks=(21,256,256),dtype=vol.dtype); z[:]=vol
    I=np.asarray(L("inklabels")[:,y0:y1,x0:x1]); Vw=np.asarray(val[:,y0:y1,x0:x1])
    np.savez_compressed(OUT+"/%s_val_labels.npz"%seg, ink=I, validation=Vw,
                        window=np.array([0,21,y0,y1,x0,x1]))
    zc=int(lo[0]); m=Vw[zc]>0
    res[seg]=dict(input_shape=list(vol.shape), window=[y0,y1,x0,x1], annotated_z=zc,
                  masked=int(m.sum()), ink_frac=round(float((I[zc][m]>0).mean()),5),
                  material=round(float((vol[zc]>0).mean()),4))
    print("   input",vol.shape,"| masked",int(m.sum()),"| ink %.3f"%res[seg]["ink_frac"],flush=True)
json.dump(res,open(OUT+"/input_probe2.json","w"),indent=1)
