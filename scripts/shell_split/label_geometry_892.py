"""The shell geometry of all 892 volumes, from labels alone.

Two of the descriptive numbers in his comment are properties of the labels, not of m7:
the share of scored non-sheet volume that sits beyond distance k, and the volume share of
each one-voxel shell. Those cost no GPU and do not have to be estimated off 60 volumes,
so compute them on the whole 892 at the same crop his analysis uses.
"""
import json, sys, time
from pathlib import Path
import numpy as np, tifffile
from scipy.ndimage import distance_transform_edt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from geom import TRIM, SIZE

LB = Path("/mnt/vesuvius/kaggle892/labels")
names = sorted(p.stem for p in LB.glob("sample_*.tif"))
out = Path("/mnt/vesuvius/experiments/shell_split/label_geometry_892.jsonl")
done = set()
if out.exists():
    done = {json.loads(l)["sample"] for l in out.read_text().splitlines() if l.strip()}
fh = out.open("a")
t0 = time.time()
for k, nm in enumerate(names):
    if nm in done:
        continue
    lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
    off = (lab.shape[0] - SIZE) // 2
    lo, hi = off + TRIM, off + SIZE - TRIM
    c = lab[lo:hi, lo:hi, lo:hi]
    sheet, scored = c == 1, c != 2
    ns = scored & ~sheet
    n_ns = int(ns.sum())
    row = {"sample": nm, "n_sheet": int(sheet.sum()), "n_scored": int(scored.sum()),
           "n_nonsheet_scored": n_ns, "shape": list(lab.shape)}
    if n_ns > 0 and sheet.sum() > 0:
        d = distance_transform_edt(~sheet)
        row["shell_vox"] = [int((((d > j - 1) & (d <= j)) & ns).sum()) for j in range(1, 6)]
        for j in (1, 2, 3, 4):
            row[f"beyond{j}_vox"] = int(((d > j) & ns).sum())
    fh.write(json.dumps(row) + "\n"); fh.flush()
    if k % 50 == 0:
        print(f"MARKER [{k}/{len(names)}] {nm} {time.time()-t0:.0f}s", flush=True)
fh.close()
print("MARKER done", flush=True)
