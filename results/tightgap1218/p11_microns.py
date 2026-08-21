"""Restate the gap population in microns and give the like-for-like figure.

Our gaps were measured inside repaired label blocks, which sit on PHerc1218 CT
level 1. Level 0 is 8.640 um, so our voxel is 17.28 um.
Dataset059 patches are level-0 coordinates on 7.91 um volumes, so flummoxjr's
"under 4 vox" is 4 x 7.91 = 31.64 um, not 4 x 17.28 = 69.12 um.
"""
import json
import numpy as np
OURS_UM=17.28; HIS_UM=7.91
rows=json.load(open("/mnt/vesuvius/tightgap1218/sites_gaps.json"))
g=np.array([r["gap"] for r in rows],dtype=float)
print("sites: %d"%len(g))
print("\nour gap distribution")
print("  in our voxels (level 1): median %.2f  p10 %.2f  p90 %.2f"%(np.median(g),np.percentile(g,10),np.percentile(g,90)))
um=g*OURS_UM
print("  in microns            : median %.1f  p10 %.1f  p90 %.1f"%(np.median(um),np.percentile(um,10),np.percentile(um,90)))
print("\nTHE PUBLISHED COMPARISON, and the like-for-like one")
pub=100*(g<4).mean()
print("  published : 'under 4 voxels' = under %.2f um -> %.2f%% of sites"%(4*OURS_UM,pub))
thr=4*HIS_UM
print("  his 4 vox : under %.2f um = %.2f of OUR voxels"%(thr,thr/OURS_UM))
like=100*(um<thr).mean()
print("  like-for-like at his physical threshold          -> %.2f%% of sites"%like)
print("  his figure                                        -> 0.02%%")
print("  ratio published %.0fx   ratio like-for-like %.0fx"%(pub/0.02,like/0.02))
print("\nhis median gap 15.6 vox = %.1f um = %.2f of our voxels"%(15.6*HIS_UM,15.6*HIS_UM/OURS_UM))
print("our median %.2f vox = %.1f um = %.2f of HIS voxels"%(np.median(g),np.median(um),np.median(um)/HIS_UM))
print("\nband table restated")
for lo,hi,lab in ((0,2,"0-2"),(2,4,"2-4"),(4,6,"4-6"),(6,10,"6-10"),(10,1e9,"10+")):
    n=int(((g>=lo)&(g<hi)).sum())
    print("  %-5s  %6d sites   %6.1f to %6.1f um"%(lab,n,lo*OURS_UM,min(hi,50)*OURS_UM))
