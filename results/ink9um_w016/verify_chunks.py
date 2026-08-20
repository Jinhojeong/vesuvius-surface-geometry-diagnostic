"""Find truncated or unreadable chunks, since a matching file count proved nothing."""
import glob, os
import numpy as np, zarr
for seg in ("pherc0814-46527","pherc1667-w029"):
    d="/mnt/vesuvius/ink9um_labels/"+seg
    for kind in ("validation_mask","inklabels","supervision_mask"):
        p="%s/%s_%s.zarr"%(d,seg,kind)
        files=sorted(glob.glob(p+"/0/*"))
        zero=[f for f in files if os.path.getsize(f)==0]
        z=zarr.open(p,mode="r"); a=z["0"] if "0" in z else z
        bad=[]
        for f in files:
            name=os.path.basename(f)
            if name.startswith("."): continue
            try:
                idx=tuple(int(x) for x in name.split("."))
                sl=tuple(slice(i*c,(i+1)*c) for i,c in zip(idx,a.chunks))
                _=np.asarray(a[sl])
            except Exception:
                bad.append(name)
                if len(bad)>6: break
        print("%-18s %-18s files %4d | zero-byte %2d | unreadable %s"%(
            seg,kind,len(files),len(zero),bad[:5] if bad else "none"))
