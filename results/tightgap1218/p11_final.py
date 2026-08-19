"""Add a measured CT-emptiness field to every crop, and rebuild the figure
choosing, per band, the crops at the median emptiness rather than by eye.

The frozen rule 4 only requires 1 percent CT in each octant, which a crop with
a large empty wedge can still pass. That rule is not changed here; the crops it
selected are the crops that ship. What is added is the measurement a consumer
needs to filter on, plus the distribution reported honestly.
"""
import glob, json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/mnt/vesuvius/tightgap1218"
bands = ["0-2", "2-4", "4-6", "6-10", "10+"]

stats = {"contact": [], "control": []}
for arm, key in (("crops", "contact"), ("control", "control")):
    for f in sorted(glob.glob("%s/%s/*.npz" % (OUT, arm))):
        d = dict(np.load(f, allow_pickle=True))
        v = d["intensity"]
        empty = float((v == 0).mean())
        d["ct_empty_frac"] = np.float32(empty)
        np.savez_compressed(f, **d)
        stats[key].append(empty)

rep = {}
for k, v in stats.items():
    a = np.array(v)
    rep[k] = dict(n=len(a), median=round(float(np.median(a)), 4),
                  p90=round(float(np.percentile(a, 90)), 4),
                  frac_over_10pct=round(float((a > 0.10).mean()), 4),
                  frac_over_25pct=round(float((a > 0.25).mean()), 4))
json.dump(rep, open(OUT + "/emptiness_summary.json", "w"), indent=1)
print(json.dumps(rep, indent=1))


def best_view(v):
    cands = [v[64], v[:, 64], v[:, :, 64]]
    e = [float(np.abs(np.gradient(c.astype(np.float32))).mean()) for c in cands]
    return cands[int(np.argmax(e))]


byband = {}
for f in sorted(glob.glob(OUT + "/crops/*.npz")):
    d = np.load(f, allow_pickle=True)
    if int(d.get("A_id", -1)) < 0:
        continue
    lab = d["instance"]
    ids = set(np.unique(lab).tolist())
    if not (int(d["A_id"]) in ids and int(d["B_id"]) in ids):
        continue
    byband.setdefault(str(d["band"]), []).append((float(d["ct_empty_frac"]), float(d["gap"]), f))

fig, axes = plt.subplots(2, 6, figsize=(17, 6.2))
for j, b in enumerate(bands):
    items = sorted(byband.get(b, []))
    for i in range(2):
        ax = axes[i, j]
        if items:
            k = len(items) // 2 + (i - 1) * max(1, len(items) // 6)
            k = max(0, min(len(items) - 1, k))
            _, g, f = items[k]
            v = np.load(f, allow_pickle=True)["intensity"]
            ax.imshow(best_view(v), cmap="gray")
            ax.set_title("%s vox, gap %.1f" % (b, g), fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
ctl = sorted(glob.glob(OUT + "/control/*.npz"))
for i in range(2):
    ax = axes[i, 5]
    v = np.load(ctl[20 + i * 10], allow_pickle=True)["intensity"]
    ax.imshow(best_view(v), cmap="gray")
    ax.set_title("control, single sheet", fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])

c = json.load(open(OUT + "/crops_summary.json"))["realised"]
s = json.load(open(OUT + "/labels_summary.json"))["contact"]
fig.suptitle("PHerc1218 tight-contact validation set, real CT at level 0\n"
             "crops per band  0-2:%d  2-4:%d  4-6:%d  6-10:%d  10+:%d  control:60   "
             "(%d of %d contact crops carry both split instances; examples drawn at the band median CT emptiness, view chosen to show sheets edge-on)"
             % (c.get("0-2", 0), c.get("2-4", 0), c.get("4-6", 0), c.get("6-10", 0),
                c.get("10+", 0), s["both"], s["n"]), fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(OUT + "/tightgap_bands.png", dpi=130)
print("figure saved")
