#!/usr/bin/env python3
"""Does the saturation plateau survive first-order multipath?

The sensing matrix is a free-space line-of-sight model, and the fair objection
is that indoor WiFi is defined by multipath. Two effects compete: specular
reflections act as image transmitters, adding sensitivity along reflected paths
(more structure, potentially more rank), while every reflected path overlaps
the direct paths near the endpoints (more correlation, potentially less).
Which wins at the plateau decides whether the saturation result is a property
of the propagation model or of the geometry.

First-order image-source check: each link's row gains four reflected-path
Fresnel kernels, one per wall, built by mirroring the transmitter across that
wall, with amplitude Gamma relative to the direct path. Gamma = 0 recovers the
paper's model; Gamma = 0.6 is a strongly reflective room.

Memory (Rule 7.5): rows are built voxel-chunk-wise and accumulated into the
L x L Gram, so no dense L x V matrix ever exists. Peak resident is the Gram
(128 MB at N = 90) plus one chunk.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import resource
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = Path(os.environ.get(
    "NDF_RESULTS_DIR", _REPO_ROOT / "docs/4-research/results"))
RESULTS_JSON = "ndf_multipath_check.json"

TAU = 0.01
ROOM_W_M, ROOM_H_M = 5.0, 4.0
WAVELENGTH_M = 0.125
WALL_INSET_M = 0.10
GRID_RES_M = 0.035           # resolves sigma_F of the shortest perimeter link
N_VALUES = (20, 40, 60, 90)
GAMMAS = (0.0, 0.3, 0.6)
VOXEL_CHUNK = 4000
FLAT_TOL = 0.02


def perimeter(n, w=ROOM_W_M, h=ROOM_H_M):
    iw, ih = w - 2 * WALL_INSET_M, h - 2 * WALL_INSET_M
    per = 2 * (iw + ih)
    pts = np.empty((n, 2))
    for i in range(n):
        s = (i + 0.5) * per / n
        if s < iw:
            pts[i] = (WALL_INSET_M + s, WALL_INSET_M)
        elif s < iw + ih:
            pts[i] = (WALL_INSET_M + iw, WALL_INSET_M + (s - iw))
        elif s < 2 * iw + ih:
            pts[i] = (WALL_INSET_M + iw - (s - iw - ih), WALL_INSET_M + ih)
        else:
            pts[i] = (WALL_INSET_M, WALL_INSET_M + ih - (s - 2 * iw - ih))
    return pts


def mirror(p, wall):
    x, y = p
    if wall == "left":
        return (-x, y)
    if wall == "right":
        return (2 * ROOM_W_M - x, y)
    if wall == "bottom":
        return (x, -y)
    return (x, 2 * ROOM_H_M - y)


def kernel(tx, rx, gx, gy):
    """Fresnel weight of the path tx->rx over voxel coordinates (gx, gy)."""
    tx = np.asarray(tx, float)
    rx = np.asarray(rx, float)
    d = rx - tx
    length = float(np.hypot(*d))
    if length < 1e-9:
        return np.zeros_like(gx)
    ux, uy = d / length
    relx, rely = gx - tx[0], gy - tx[1]
    t = (relx * ux + rely * uy) / length
    d_perp = np.abs(relx * uy - rely * ux)
    sigma = np.sqrt(WAVELENGTH_M * length / 4.0) / 2.0
    w = np.exp(-d_perp ** 2 / (2.0 * sigma ** 2)) * 4.0 * t * (1.0 - t)
    w[(t < 0) | (t > 1)] = 0.0
    return w


def gram_for(pos, gamma):
    n = len(pos)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    L = len(pairs)
    xs = np.arange(GRID_RES_M / 2, ROOM_W_M, GRID_RES_M)
    ys = np.arange(GRID_RES_M / 2, ROOM_H_M, GRID_RES_M)
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    gx, gy = gx.ravel(), gy.ravel()
    G = np.zeros((L, L))
    for lo in range(0, len(gx), VOXEL_CHUNK):
        cx, cy = gx[lo:lo + VOXEL_CHUNK], gy[lo:lo + VOXEL_CHUNK]
        blk = np.empty((L, len(cx)))
        for r, (i, j) in enumerate(pairs):
            w = kernel(pos[i], pos[j], cx, cy)
            if gamma > 0:
                for wall in ("left", "right", "bottom", "top"):
                    w = w + gamma * kernel(mirror(pos[i], wall), pos[j], cx, cy)
            blk[r] = w
        G += blk @ blk.T
        del blk
        gc.collect()
    return G, L


def ndf_of(G):
    ev = np.clip(np.linalg.eigvalsh(G)[::-1], 0.0, None)
    s = np.sqrt(ev)
    return int((s > TAU * s[0]).sum()) if s[0] > 0 else 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / RESULTS_JSON)
    args = ap.parse_args()
    out = {"gammas": list(GAMMAS), "n_values": list(N_VALUES),
           "grid_res_m": GRID_RES_M, "ndf": {}}
    print(f"{'gamma':>6}" + "".join(f"{'N='+str(n):>8}" for n in N_VALUES)
          + f"{'plateau growth':>16}")
    for g in GAMMAS:
        row = []
        for n in N_VALUES:
            G, L = gram_for(perimeter(n), g)
            k = ndf_of(G)
            row.append(k)
            del G
            gc.collect()
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
            print(f"    .. gamma={g} N={n}: NDF={k} (RSS {rss:.1f} GB)",
                  flush=True)
        growth = (row[-1] - row[-2]) / max(row[-2], 1)
        out["ndf"][str(g)] = row
        out.setdefault("plateau_growth", {})[str(g)] = growth
        print(f"{g:>6}" + "".join(f"{k:>8}" for k in row)
              + f"{growth:>15.1%}", flush=True)
    flat = all(abs(v) < FLAT_TOL for v in out["plateau_growth"].values())
    print(f"\nVERDICT: plateau {'SURVIVES' if flat else 'DOES NOT SURVIVE'} "
          f"first-order multipath up to Gamma={max(GAMMAS)}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
