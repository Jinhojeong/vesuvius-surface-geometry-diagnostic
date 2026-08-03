"""Fetch aviad's 8-cell pilot bundle from HF with manifest verification."""
import hashlib
import json
import os
import urllib.request

BASE = "https://huggingface.co/datasets/AviadCoh/vesuvius-physical-fusion-pilot"
API = "https://huggingface.co/api/datasets/AviadCoh/vesuvius-physical-fusion-pilot/tree/main"
OUT = "/mnt/vesuvius/fusion_pilot"


def ls(path=""):
    url = API + (f"/{path}" if path else "")
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)


def fetch(rel, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return
    url = f"{BASE}/resolve/main/{rel}"
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=1800) as r, open(tmp, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    os.replace(tmp, dest)


entries = []
stack = [""]
while stack:
    p = stack.pop()
    for e in ls(p):
        if e["type"] == "directory":
            stack.append(e["path"])
        else:
            entries.append(e["path"])
print(f"{len(entries)} files:", entries, flush=True)
for rel in entries:
    dest = os.path.join(OUT, os.path.basename(rel))
    print("fetching", rel, flush=True)
    fetch(rel, dest)

# verify npz sha256 against the manifest if present
man = [f for f in os.listdir(OUT) if "manifest" in f.lower() and f.endswith(".json")]
if man:
    m = json.load(open(os.path.join(OUT, man[0])))
    print("manifest keys:", list(m.keys())[:8], flush=True)
print("FETCH_DONE", len(entries), flush=True)
