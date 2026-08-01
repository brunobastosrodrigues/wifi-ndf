#!/usr/bin/env python3
"""Sensing efficiency of deliberately mis-placed meshes, with replicates.

The placement claim in the TMC paper rests on designed layouts rather than on
uniformly sampled ones, because uniform sampling without a minimum separation
produces node pairs closer than MIN_LINK_LENGTH_M whose links are degenerate
rather than merely redundant (36% of configurations in our earlier sweep).

Every layout here is physically realisable: all node pairs are held at least
MIN_SEPARATION_M apart. Replicates jitter each node within a small radius so
the reported efficiency carries an uncertainty rather than being a single
deterministic draw.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "opencsi" / "src"))

from opencsi.geometry.ndf import (  # noqa: E402
    build_sensing_matrix,
    DEFAULT_NDF_THRESHOLD,
    MIN_LINK_LENGTH_M,
)

RESULTS_DIR = Path(os.environ.get(
    "NDF_RESULTS_DIR", _REPO_ROOT / "docs/4-research/results"))
RESULTS_JSON = "ndf_failure_modes.json"

ROOM_W_M, ROOM_H_M = 5.0, 4.0
WAVELENGTH_M = 0.125
N_NODES = 14
WALL_INSET_M = 0.10
MIN_SEPARATION_M = 0.25
JITTER_M = 0.05
N_REPLICATES = 12
BASE_SEED = 20260731
VOXEL_LINK_RATIO = 4.0
MIN_GRID_RES_M = 0.015
SIGMA_F_CELLS = 2.0            # cells per Fresnel std dev
MAX_GRID_RES_M = 0.25


def perimeter(n: int) -> np.ndarray:
    w, h = ROOM_W_M - 2 * WALL_INSET_M, ROOM_H_M - 2 * WALL_INSET_M
    per = 2 * (w + h)
    pts = np.empty((n, 2))
    for i in range(n):
        s = (i + 0.5) * per / n
        if s < w:
            pts[i] = (WALL_INSET_M + s, WALL_INSET_M)
        elif s < w + h:
            pts[i] = (WALL_INSET_M + w, WALL_INSET_M + (s - w))
        elif s < 2 * w + h:
            pts[i] = (WALL_INSET_M + w - (s - w - h), WALL_INSET_M + h)
        else:
            pts[i] = (WALL_INSET_M, WALL_INSET_M + h - (s - 2 * w - h))
    return pts


def corner_cluster(n: int) -> np.ndarray:
    side = int(np.ceil(np.sqrt(n)))
    xs = np.linspace(WALL_INSET_M, WALL_INSET_M + 1.8, side)
    ys = np.linspace(WALL_INSET_M, WALL_INSET_M + 1.8, side)
    return np.array([(x, y) for y in ys for x in xs])[:n]


def two_groups(n: int) -> np.ndarray:
    half = n // 2
    out = []
    for k, x0 in ((half, WALL_INSET_M), (n - half, ROOM_W_M - WALL_INSET_M - 0.8)):
        side = int(np.ceil(np.sqrt(k)))
        xs = np.linspace(x0, x0 + 0.8, side)
        ys = np.linspace(WALL_INSET_M, ROOM_H_M - WALL_INSET_M, side)
        out += [(x, y) for y in ys for x in xs][:k]
    return np.array(out)


def single_wall(n: int) -> np.ndarray:
    xs = np.linspace(WALL_INSET_M, ROOM_W_M - WALL_INSET_M, n)
    return np.array([(x, WALL_INSET_M) for x in xs])


LAYOUTS = {
    "Well-spread": perimeter,
    "One wall": single_wall,
    "Corner cluster": corner_cluster,
    "Two groups": two_groups,
}


def efficiency(pos: np.ndarray) -> tuple[float, int, int]:
    n_links = len(pos) * (len(pos) - 1) // 2
    # Physical grid criterion: resolve the narrowest Fresnel zone in the layout.
    # The voxel-per-link rule under-resolves exactly the near-degenerate layouts
    # measured here, and since coarsening merges near-degenerate directions it
    # depresses their rank, biasing the failure modes to look worse than they are.
    d_all = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    d_min = d_all[np.triu_indices(len(pos), k=1)].min()
    res = float(max(np.sqrt(WAVELENGTH_M * d_min / 4.0) / 2.0 / SIGMA_F_CELLS,
                    MIN_GRID_RES_M))
    w, _, _ = build_sensing_matrix(
        pos, grid_resolution_m=res, room_bounds=((0.0, 0.0), (ROOM_W_M, ROOM_H_M)),
        wavelength_m=WAVELENGTH_M)
    dense = w.toarray()
    del w
    sigma = np.linalg.svd(dense, compute_uv=False)
    del dense
    ndf = 0 if sigma[0] <= 0 else int((sigma > DEFAULT_NDF_THRESHOLD * sigma[0]).sum())
    return ndf / n_links, ndf, n_links


def min_separation(pos: np.ndarray) -> float:
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    return float(d[np.triu_indices(len(pos), k=1)].min())


def main() -> None:
    out = {}
    for name, fn in LAYOUTS.items():
        base = fn(N_NODES)
        sep = min_separation(base)
        if sep < MIN_SEPARATION_M:
            raise RuntimeError(
                f"{name}: min separation {sep:.3f} m below "
                f"{MIN_SEPARATION_M} m; layout is not physically realisable")
        vals = []
        for r in range(N_REPLICATES):
            rng = np.random.default_rng(BASE_SEED + r)
            pos = base + rng.uniform(-JITTER_M, JITTER_M, base.shape)
            pos[:, 0] = np.clip(pos[:, 0], WALL_INSET_M, ROOM_W_M - WALL_INSET_M)
            pos[:, 1] = np.clip(pos[:, 1], WALL_INSET_M, ROOM_H_M - WALL_INSET_M)
            if min_separation(pos) < MIN_LINK_LENGTH_M:
                continue
            eff, ndf, n_links = efficiency(pos)
            vals.append(eff)
        out[name] = {
            "mean": float(np.mean(vals)), "std": float(np.std(vals)),
            "n_replicates": len(vals), "min_separation_m": sep,
            "n_links": N_NODES * (N_NODES - 1) // 2,
        }
        print(f"{name:<16} NDF/L = {out[name]['mean']:.3f} "
              f"+/- {out[name]['std']:.3f}  (n={len(vals)}, "
              f"min sep {sep:.2f} m)")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / RESULTS_JSON).write_text(json.dumps({
        "config": {
            "room": [ROOM_W_M, ROOM_H_M], "n_nodes": N_NODES,
            "wavelength_m": WAVELENGTH_M, "jitter_m": JITTER_M,
            "min_separation_m": MIN_SEPARATION_M,
            "ndf_threshold": DEFAULT_NDF_THRESHOLD,
        },
        "layouts": out,
    }, indent=2))
    print(f"\nwrote {RESULTS_DIR / RESULTS_JSON}")


if __name__ == "__main__":
    main()
