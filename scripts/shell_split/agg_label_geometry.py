"""Label-only shell geometry, by population, both conventions."""
import json
from pathlib import Path
import numpy as np

rows = [json.loads(l) for l in Path("/mnt/vesuvius/experiments/shell_split/label_geometry_892.jsonl").read_text().splitlines() if l.strip()]
g = json.loads(Path("/mnt/vesuvius/kaggle892/groups892.json").read_text())
his60 = set(json.loads(Path("/tmp/his60.json").read_text()))
nl60 = set(json.loads(Path("/mnt/vesuvius/experiments/shell_split/nonlocated60.json").read_text()))

def pop(r):
    return "located" if g.get(r["sample"]) in ("located", "intersecting", "iou1") else "nonlocated"

def summ(sel, tag):
    sel = [r for r in sel if r.get("shell_vox") and r["n_nonsheet_scored"] > 0]
    if not sel: return {"tag": tag, "n": 0}
    ns = np.array([r["n_nonsheet_scored"] for r in sel], float)
    o = {"tag": tag, "n": len(sel)}
    for k in (1, 2, 3, 4):
        v = np.array([r[f"beyond{k}_vox"] for r in sel], float)
        o[f"beyond{k}_vol_share_pct"] = {
            "per_volume_median": round(float(np.median(v / ns)) * 100, 2),
            "pooled": round(float(v.sum() / ns.sum()) * 100, 2)}
    sv = np.array([r["shell_vox"] for r in sel], float)
    o["shell_vol_share_pct"] = {
        "per_volume_median": [round(float(np.median(sv[:, i] / ns)) * 100, 2) for i in range(5)],
        "pooled": [round(float(sv[:, i].sum() / ns.sum()) * 100, 2) for i in range(5)]}
    o["median_base_rate"] = round(float(np.median([r["n_sheet"] / r["n_scored"] for r in sel])), 4)
    o["median_sheet_voxels"] = int(np.median([r["n_sheet"] for r in sel]))
    o["n_with_no_sheet"] = sum(1 for r in rows if r in sel and r["n_sheet"] == 0)
    return o

out = {
    "all892": summ(rows, "all 892"),
    "located189": summ([r for r in rows if pop(r) == "located"], "located 189"),
    "nonlocated703": summ([r for r in rows if pop(r) == "nonlocated"], "nonlocated 703"),
    "his60": summ([r for r in rows if r["sample"] in his60], "his cohort 60"),
    "nonlocated60": summ([r for r in rows if r["sample"] in nl60], "extension draw 60"),
    "degenerate": {"n_no_sheet_all892": sum(1 for r in rows if r["n_sheet"] == 0),
                   "n_no_sheet_nonlocated": sum(1 for r in rows if pop(r) == "nonlocated" and r["n_sheet"] == 0),
                   "n_no_sheet_in_draw60": sum(1 for r in rows if r["sample"] in nl60 and r["n_sheet"] == 0),
                   "n_sheet_under_400_in_draw60": sum(1 for r in rows if r["sample"] in nl60 and r["n_sheet"] < 400)},
}
Path("/mnt/vesuvius/experiments/shell_split/label_geometry_summary.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
