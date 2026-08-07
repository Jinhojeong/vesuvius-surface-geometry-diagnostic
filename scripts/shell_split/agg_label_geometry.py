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

def summ(pop, tag):
    sel = [r for r in pop if r.get("shell_vox") and r["n_nonsheet_scored"] > 0]
    head = {"tag": tag, "n_population": len(pop), "n_scored": len(sel),
            "n_excluded_no_sheet": len(pop) - len(sel)}
    if not sel: return head
    ns = np.array([r["n_nonsheet_scored"] for r in sel], float)
    o = dict(head)
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
    return o

out = {
    "legend": ("n_population is how many volumes the label carries, n_scored is how many of "
               "them have labelled sheet inside the centred crop, and every figure in the "
               "block is over n_scored. The excluded volumes are exactly the ones with no "
               "sheet in the crop, where the shell geometry is undefined."),
    "all": summ(rows, "all public volumes"),
    "located": summ([r for r in rows if pop(r) == "located"], "located population"),
    "nonlocated": summ([r for r in rows if pop(r) == "nonlocated"], "nonlocated population"),
    "his60": summ([r for r in rows if r["sample"] in his60], "his cohort"),
    "nonlocated_draw60": summ([r for r in rows if r["sample"] in nl60], "extension draw"),
    "degenerate": {"n_no_sheet_all892": sum(1 for r in rows if r["n_sheet"] == 0),
                   "n_no_sheet_nonlocated": sum(1 for r in rows if pop(r) == "nonlocated" and r["n_sheet"] == 0),
                   "n_no_sheet_in_draw60": sum(1 for r in rows if r["sample"] in nl60 and r["n_sheet"] == 0),
                   "n_sheet_under_400_in_draw60": sum(1 for r in rows if r["sample"] in nl60 and r["n_sheet"] < 400)},
}
Path("/mnt/vesuvius/experiments/shell_split/label_geometry_summary.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
