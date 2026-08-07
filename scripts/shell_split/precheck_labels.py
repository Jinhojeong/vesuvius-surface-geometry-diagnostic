"""Pre-gate: label-only. Recompute the margin from CT+labels and check it against the
counts TAUIL-Abd-Elilah's shipped rows carry (n_margin_scored, margin_share_of_nonsheet).
Those two numbers do not involve m7 at all, so a match proves (a) our kaggle892 CT and
labels are the same bytes his were and (b) our margin_relabel re-implementation is his.
"""
import json, sys, time
from pathlib import Path
import numpy as np, tifffile
sys.path.insert(0, str(Path(__file__).parent))
from geom import margin_mask, TRIM, SIZE

IM = Path("/mnt/vesuvius/kaggle892/images")
LB = Path("/mnt/vesuvius/kaggle892/labels")
his = json.loads(Path(sys.argv[1]).read_text())["rows"]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 6
out = []
for r in his[:n]:
    if r.get("status") != "ok":
        continue
    nm = r["sample"]
    t0 = time.time()
    ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
    lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
    off = (ct.shape[0] - SIZE) // 2
    sl = (slice(off + TRIM, off + SIZE - TRIM),) * 3
    mg = margin_mask(ct, lab)[sl]
    labc = lab[sl]
    scored = labc != 2
    ns = scored & (labc != 1)
    n_margin = int((mg & scored).sum())
    share = n_margin / float(ns.sum())
    row = {"sample": nm, "shape": list(ct.shape),
           "ours_n_margin_scored": n_margin, "his_n_margin_scored": r["n_margin_scored"],
           "ours_margin_share": round(share, 6), "his_margin_share": r["margin_share_of_nonsheet"],
           "match": n_margin == r["n_margin_scored"], "sec": round(time.time() - t0, 1)}
    out.append(row)
    print(json.dumps(row), flush=True)
print("EXACT MATCHES: %d/%d" % (sum(r["match"] for r in out), len(out)))
