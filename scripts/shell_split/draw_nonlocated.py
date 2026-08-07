"""The extension cohort, drawn once and written down so it cannot be re-drawn.

pool  = the 703 volumes groups892.json calls `nonlocated`, sorted by sample id
draw  = numpy.random.default_rng(0).choice(len(pool), size=60, replace=False), sorted
No screening on label content: whatever the draw gives is what gets run, and volumes that
come back degenerate are reported as degenerate rather than replaced.
"""
import json
from pathlib import Path
import numpy as np

g = json.loads(Path("/mnt/vesuvius/kaggle892/groups892.json").read_text())
pool = sorted(k for k, v in g.items() if v == "nonlocated")
assert len(pool) == 703, len(pool)
idx = np.random.default_rng(0).choice(len(pool), size=60, replace=False)
draw = sorted(pool[i] for i in idx)
Path("/mnt/vesuvius/experiments/shell_split/nonlocated60.json").write_text(json.dumps(draw))
Path("/mnt/vesuvius/experiments/shell_split/located189.json").write_text(
    json.dumps(sorted(k for k, v in g.items() if v in ("located", "intersecting", "iou1"))))
print(json.dumps({"pool": len(pool), "seed": 0, "n": 60, "draw": draw}, indent=1))
