"""P9: quantify what a support-aware stopping rule would buy.

The tracer has no support constraint, and we measured that unsupported quads
concentrate at the mesh fringe. This applies the rule the tracer lacks, as a
post-hoc trim on the emitted tifxyz, and measures the cost and the benefit:

  keep a quad if CT support exists at its position, then erode the surviving
  mask so that isolated survivors past the frontier do not remain,
  and report supported share and retained area before and after.

Run on the P8 traces (both arms) so the effect is measured on real patches.
"""
import json
import sys

import numpy as np
import zarr
from scipy import ndimage as ndi

sys.path.insert(0, "/mnt/vesuvius/vcbuild")
import score_ab  # noqa: E402

ROOT = "/mnt/vesuvius/p8_sprint"
rep = json.load(open("%s/r1_repair.json" % ROOT))
lo, hi = rep["bounds_l0"]
sup = zarr.open("%s/r1_support.zarr/0" % ROOT, mode="r")

out = []
for i in range(3):
    for arm in ("before", "after"):
        run = score_ab.load_tifxyz("%s/trace/r1_s%d_%s" % (ROOT, i, arm))
        v = run["valid"]
        pts = np.stack([run["z"], run["y"], run["x"]], -1)
        inb = ((pts >= np.array(lo)) & (pts < np.array(hi))).all(-1) & v
        idx = np.argwhere(inb)
        if len(idx) < 50:
            continue
        s = np.array([int(sup[int(round(pts[a, b, 0])),
                              int(round(pts[a, b, 1])),
                              int(round(pts[a, b, 2]))]) for a, b in idx])
        supported = np.zeros(v.shape, bool)
        for (a, b), val in zip(idx, s):
            supported[a, b] = val > 0

        # rule: keep supported quads, then drop survivors that lost their
        # neighbourhood (opening removes fringe spurs the trim leaves behind)
        keep = ndi.binary_opening(supported, structure=np.ones((3, 3)))
        n_in = int(inb.sum())
        n_sup = int(supported.sum())
        n_keep = int((keep & inb).sum())
        # fringe profile of what the rule removes
        d = ndi.distance_transform_edt(v)
        removed = inb & ~keep
        out.append({
            "seed": i, "arm": arm,
            "quads_in_region": n_in,
            "supported_frac_before": round(n_sup / n_in, 4),
            "retained_frac": round(n_keep / n_in, 4),
            "supported_frac_after_rule": 1.0 if n_keep == 0 else round(
                float((supported & keep & inb).sum()) / n_keep, 4),
            "removed_mean_edge_dist": round(float(d[removed].mean()), 2)
            if removed.any() else None,
            "kept_mean_edge_dist": round(float(d[keep & inb].mean()), 2)
            if (keep & inb).any() else None,
        })
        print(json.dumps(out[-1]))

json.dump({"rule": "keep quads with CT>0 at their position, then 3x3 binary "
                   "opening to drop isolated survivors; measured on P8 r1 "
                   "traces, both arms, quads inside the built support region",
           "rows": out}, open("%s/r1_trim.json" % ROOT, "w"), indent=1)
print("saved r1_trim.json")
