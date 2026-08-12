"""PHerc1203 thickness zone map, the deliverable TAUIL accepted in villa #191
comment 5241945617, on his terms: exact L1 origin, shape and hashes, usable as
a predeclared secondary stratum only.

Input: 7jycwjmbfn-eng's physical-truth labels (uint8 bit flags on the lo
volume 20250820131727 L1 grid, origin_l1 [3936,0,0], shape [2016,3456,3456]).
No instance ids exist there, so this is the thickness-and-centerline half of
the census instrument. Per sampled ray run on the material bit:

  run thickness   material run length along the local normal, in L1 voxels
  n_centerlines   centerline-bit points crossed within the run
  fused signature thickness ratio >= 1.6 x the global median AND exactly one
                  centerline crossed. A resolved double sheet crosses two.

Zone map: 96^3-voxel cells (1.8 mm at 18.724 um). Per cell: n_runs, median
thickness ratio, fused-signature share, boundary_poor voxel share, valid
share. Rays are sampled on a deterministic per-cell md5 grid, no RNG.

Chunk-stream memory shape: cells are processed in z-slabs of 96 with a 16-vox
halo, so peak memory stays near one slab of the five bit planes.
"""
from __future__ import annotations
import hashlib, json, sys, time
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

SRC = Path("/mnt/vesuvius/audit0139/labels1203_L1.zarr")
OUT = Path("/mnt/vesuvius/audit0139/thickness_map_1203")
OUT.mkdir(exist_ok=True)

CELL = 96
HALO = 16
K, STEP = 41, 0.75
SPAN = (K // 2) * STEP
OFFS = np.arange(K, dtype=np.float32) * STEP - SPAN
RAYS_PER_CELL = 24
THICK_F = 1.6


_ZARR = None


def read_slab(z0, z1, meta):
    # v2: decode through zarr. v1 copied raw chunk-file bytes into the label
    # planes, but the store is blosc-zstd, so v1 measured compression noise.
    global _ZARR
    if _ZARR is None:
        import zarr
        _ZARR = zarr.open(str(SRC), mode="r")
    return np.asarray(_ZARR[z0:z1])


def cell_normals(mask, zlo, zhi, ys, xs, w=12):
    """One smoothed SDF per cell window, gradients returned for lookup. The
    per-ray version recomputed this EDT 24 times per cell."""
    sl = (slice(max(0, zlo - w), min(mask.shape[0], zhi + w)),
          slice(max(0, ys - w), min(mask.shape[1], ys + CELL + w)),
          slice(max(0, xs - w), min(mask.shape[2], xs + CELL + w)))
    win = mask[sl]
    sm = (ndi.distance_transform_edt(win)
          - ndi.distance_transform_edt(~win)).astype(np.float32)
    sm = ndi.gaussian_filter(sm, 1.5)
    g = np.gradient(sm)
    orig = (sl[0].start, sl[1].start, sl[2].start)
    return g, orig


def normal_from(g, orig, p0):
    c = (p0[0] - orig[0], p0[1] - orig[1], p0[2] - orig[2])
    n = np.array([g[a][c] for a in range(3)], np.float32)
    nn = np.linalg.norm(n)
    return (n / nn) if nn > 1e-4 else None


def run_stats(matf, cen, p0, nv):
    cc = (np.asarray(p0, float)[None, :] + OFFS[:, None] * nv[None, :]).T
    on = ndi.map_coordinates(matf, cc, order=0, mode="constant") > 0.5
    c = K // 2
    if not on[c]:
        return None
    a = c
    while a > 0 and on[a - 1]:
        a -= 1
    b = c
    while b < K - 1 and on[b + 1]:
        b += 1
    if a == 0 or b == K - 1:
        return None                     # run leaves the probe, thickness unknown
    kk = np.arange(a, b + 1)
    q = np.round(np.asarray(p0, float)[None, :]
                 + (kk * STEP - SPAN)[:, None] * nv[None, :]).astype(int)
    ok = ((q >= 0) & (q < np.array(matf.shape))).all(1)
    q = q[ok]
    ncen = int((cen[q[:, 0], q[:, 1], q[:, 2]]).sum() > 0) if False else \
        int(np.count_nonzero(np.diff(np.r_[0, cen[q[:, 0], q[:, 1], q[:, 2]]
                                           .astype(np.int8)]) == 1))
    return (b - a + 1) * STEP, ncen


def main() -> None:
    meta = json.loads((SRC / ".zarray").read_text())
    attrs = json.loads((SRC / ".zattrs").read_text())
    sh = meta["shape"]
    ncz, ncy, ncx = [int(np.ceil(s / CELL)) for s in sh]
    grids = {k: np.zeros((ncz, ncy, ncx), np.float32)
             for k in ("n_runs", "median_th", "fused_share",
                       "boundary_poor_share", "valid_share")}
    all_th = []
    cellruns: dict[tuple, list] = {}
    t0 = time.time()
    for zc in range(ncz):
        z0 = zc * CELL
        z1 = min(z0 + CELL, sh[0])
        lo, hi = max(0, z0 - HALO), min(sh[0], z1 + HALO)
        slab = read_slab(lo, hi, meta)
        mat = (slab & 2) > 0
        matf = mat.astype(np.float32)      # once per slab; per-ray cast was 5 GB
        cen = (slab & 4) > 0
        bp = (slab & 16) > 0
        va = (slab & 1) > 0
        nv_, nm_ = int(va.sum()), int(mat.sum())
        if nv_ > 1_000_000:
            r = nm_ / nv_
            assert 0.3 < r <= 1.0, (
                f"slab {zc}: material/valid {r:.3f} outside label band, "
                f"read path broken (noise reads ~1.0 with mat~=valid~=bpoor)")
        for yc in range(ncy):
            for xc in range(ncx):
                ys, xs = yc * CELL, xc * CELL
                cell_mat = mat[z0 - lo:z1 - lo, ys:ys + CELL, xs:xs + CELL]
                nmat = int(cell_mat.sum())
                grids["valid_share"][zc, yc, xc] = \
                    float(va[z0 - lo:z1 - lo, ys:ys + CELL, xs:xs + CELL].mean())
                grids["boundary_poor_share"][zc, yc, xc] = \
                    float(bp[z0 - lo:z1 - lo, ys:ys + CELL, xs:xs + CELL].mean())
                if nmat < 2000:
                    continue
                idx = np.argwhere(cell_mat)
                key0 = f"{zc}:{yc}:{xc}"
                g, orig = cell_normals(mat, z0 - lo, z1 - lo, ys, xs)
                runs = []
                for r in range(RAYS_PER_CELL):
                    h = int(hashlib.md5(f"{key0}:{r}".encode()).hexdigest()[:8],
                            16)
                    p_local = idx[h % len(idx)]
                    p0 = (int(p_local[0]) + z0 - lo, int(p_local[1]) + ys,
                          int(p_local[2]) + xs)
                    nv = normal_from(g, orig, p0)
                    if nv is None:
                        continue
                    st = run_stats(matf, cen, p0, nv)
                    if st is not None:
                        runs.append(st)
                if runs:
                    cellruns[(zc, yc, xc)] = runs
                    all_th += [t for t, _ in runs]
        print(f"slab {zc + 1}/{ncz} {(time.time() - t0) / 60:.1f}min "
              f"runs={len(all_th)}", flush=True)

    med = float(np.median(all_th))
    for (zc, yc, xc), runs in cellruns.items():
        th = np.array([t for t, _ in runs])
        nc = np.array([c for _, c in runs])
        grids["n_runs"][zc, yc, xc] = len(runs)
        grids["median_th"][zc, yc, xc] = float(np.median(th) / med)
        grids["fused_share"][zc, yc, xc] = \
            float(((th >= THICK_F * med) & (nc == 1)).mean())

    np.savez_compressed(OUT / "zones.npz", **grids)

    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for b in iter(lambda: f.read(1 << 20), b""):
                h.update(b)
        return h.hexdigest()

    prov = {
        "what": "PHerc1203 thickness zone map for use as a predeclared "
                "secondary stratum, per villa #191 comment 5241945617",
        "input_labels": "labels1203_L1.tar, 7jycwjmbfn-eng/pherc0139-"
                        "physical-audit release v1.0",
        "input_tar_sha256": sha("/mnt/vesuvius/audit0139/labels1203_L1.tar"),
        "label_grid": attrs["grid"],
        "origin_l1": attrs["origin_l1"],
        "label_shape": sh,
        "cell_vox": CELL,
        "cell_um": CELL * 18.724,
        "zone_grid_shape": [ncz, ncy, ncx],
        "rays_per_cell": RAYS_PER_CELL,
        "probe": {"K": K, "step_vox": STEP,
                  "runs_touching_probe_end_discarded": True},
        "read_path": "zarr blosc decode (v2); supersedes retracted "
                     "v1 raw chunk-byte read",
        "global_median_run_vox": med,
        "thick_factor": THICK_F,
        "fused_signature": "run thickness >= 1.6 x global median AND exactly "
                           "one centerline segment crossed; a resolved double "
                           "sheet crosses two",
        "boundary_poor_carried": "per-cell voxel share of bit 16, reported "
                                 "alongside rather than filtered",
        "sampling": "24 rays per 96^3 cell at md5-hash-selected material "
                    "voxels, no RNG",
        "n_runs_total": int(sum(len(r) for r in cellruns.values())),
        "output_sha256": sha(OUT / "zones.npz"),
    }
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=1))
    print(json.dumps({k: v for k, v in prov.items()
                      if k not in ("input_tar_sha256", "output_sha256")},
                     indent=1)[:800], flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
