#!/usr/bin/env python3
"""Two questions the paper's practical claims rest on but never computes.

TMC R3 (practitioner), M1: the paper tells vendors of 2-5 node kits that
placement does not matter because they sit below N*. But N* was measured for
WELL-SPREAD layouts only, and Table II shows that at N=14, below N*, layouts
differ by 3.5x. So N* is itself placement-dependent and the claim is untested
exactly where the installed base lives. A real 3-node kit is frequently
collinear: router at the ONT, node two on the TV wall, node three at the desk,
all along the long axis. Q1 computes NDF for those layouts.

TMC R3, M3: NDF counts independent measurements available at one INSTANT. The
paper never models what it costs to acquire them. A mesh sounds L = C(N,2)
links out of one shared airtime budget, so per-link update rate falls as
1/L = O(1/N^2) while NDF grows at best linearly. Q2 computes the spatial
information RATE, which is the quantity a deployer is actually buying.

Q2's algebra collapses to something worth stating plainly. Under a fixed
sounding budget of R packets/s shared over L links, each link is sampled at
R/L Hz, so independent measurements per second is

    I(N) = NDF(N) x R/L(N) = R x (NDF(N)/L(N)),

i.e. the total sounding rate times sensing efficiency. Below N* efficiency is
exactly 1 and I(N) = R regardless of node count. Above N* efficiency falls, so
I(N) strictly DECREASES. Under a fixed airtime budget, nodes beyond N* do not
buy diminishing returns, they buy negative returns.

This treats temporal samples of one link as exchangeable with simultaneous
samples of different links, which holds only when the scene decorrelates
between samples. For a moving target that is optimistic, so I(N) is an upper
bound. The Nyquist floor is separate and binds first: a link sampled below
f_min aliases whatever it sees, capping L at R/f_min outright.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "opencsi" / "src"))

from opencsi.geometry.ndf import build_sensing_matrix  # noqa: E402

RESULTS_DIR = Path(os.environ.get(
    "NDF_RESULTS_DIR", _REPO_ROOT / "docs/4-research/results"))
RESULTS_JSON = "ndf_small_and_airtime.json"

TAU = 0.01
WAVELENGTH_M = 0.125
WALL_INSET_M = 0.10
VOXEL_LINK_RATIO = 4.0
MIN_GRID_RES_M = 0.02
MAX_GRID_RES_M = 0.25

# ---- Q1: small-N layouts, the actual installed base -------------------------
# Reference room plus a long apartment living space, where the collinear failure
# mode is most natural: everything lands along the one long wall with the outlets.
Q1_ROOMS = {"reference 5x4 m": (5.0, 4.0), "long room 8x4 m": (8.0, 4.0)}
Q1_N_VALUES = (2, 3, 4, 5, 6, 8)
# The voxel-per-link heuristic hits its coarse clip when L is small, which is
# exactly the small-N regime Q1 asks about, and at 0.25 m the grid barely
# resolves a Fresnel half-width (sigma_F ~ 0.2 m at d = 5 m). That under-counted
# rank at N=8 (17 vs 18). These matrices are tiny, so Q1 uses a fixed fine grid.
# Verified grid-converged: NDF is identical at 0.25, 0.15, 0.10, 0.05, 0.03 m
# for N <= 6, and stable from 0.10 m down at N = 8.
Q1_GRID_RES_M = 0.05
COLLINEAR_JITTER_M = 0.0     # exactly collinear: the worst realistic case
NEAR_COLLINEAR_M = 0.5       # what "roughly along one wall" actually looks like

# ---- Q2: airtime ------------------------------------------------------------
# The paper's own testbed: ~100 Hz aggregate CSI over 91 links.
TESTBED_AGGREGATE_HZ = 100.0
# Per-link rate a detector needs. Presence/slow motion at the low end, gait at
# the high end. Both are order-of-magnitude figures, quoted as such.
F_MIN_PRESENCE_HZ = 1.0
F_MIN_GAIT_HZ = 10.0
Q2_N_VALUES = tuple(range(4, 61, 2))


def _spectrum(pos, w, h, lam=WAVELENGTH_M, res=None):
    pos = np.asarray(pos, dtype=float)
    n_links = len(pos) * (len(pos) - 1) // 2
    if res is None:
        res = float(np.clip(np.sqrt(w * h / (VOXEL_LINK_RATIO * n_links)),
                            MIN_GRID_RES_M, MAX_GRID_RES_M))
    m, _, _ = build_sensing_matrix(pos, grid_resolution_m=res,
                                   room_bounds=((0.0, 0.0), (w, h)),
                                   wavelength_m=lam)
    actual_links = m.shape[0]
    dense = m.toarray()
    del m
    s = np.linalg.svd(dense, compute_uv=False)
    del dense
    gc.collect()
    return s, actual_links


def _ndf(pos, w, h, res=None):
    s, n_links = _spectrum(pos, w, h, res=res)
    n = int((s > TAU * s[0]).sum()) if s[0] > 0 else 0
    return n, n_links


# ---- small-N layout generators ---------------------------------------------
def collinear(n, w, h, offset_m=COLLINEAR_JITTER_M):
    """All nodes along the long wall. offset_m alternates them off the line."""
    xs = np.linspace(WALL_INSET_M, w - WALL_INSET_M, n)
    ys = np.full(n, WALL_INSET_M) + offset_m * (np.arange(n) % 2)
    return np.column_stack([xs, ys])


def l_shaped(n, w, h):
    """Two adjacent walls: the natural layout when the room has one door."""
    half = (n + 1) // 2
    a = [(WALL_INSET_M + i * (w - 2 * WALL_INSET_M) / max(half - 1, 1),
          WALL_INSET_M) for i in range(half)]
    b = [(WALL_INSET_M, WALL_INSET_M + (j + 1) * (h - 2 * WALL_INSET_M)
          / max(n - half, 1)) for j in range(n - half)]
    return np.array(a + b)


def spread(n, w, h):
    """Best case: nodes distributed around the perimeter."""
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


Q1_LAYOUTS = {
    "Collinear": lambda n, w, h: collinear(n, w, h, COLLINEAR_JITTER_M),
    "Near-collinear": lambda n, w, h: collinear(n, w, h, NEAR_COLLINEAR_M),
    "L-shaped": l_shaped,
    "Spread": spread,
}


def run_q1():
    out = {}
    for room, (w, h) in Q1_ROOMS.items():
        out[room] = {}
        for name, fn in Q1_LAYOUTS.items():
            rec = {}
            for n in Q1_N_VALUES:
                ndf, n_links = _ndf(fn(n, w, h), w, h, res=Q1_GRID_RES_M)
                rec[str(n)] = {"ndf": ndf, "links": n_links,
                               "efficiency": ndf / n_links if n_links else 0.0}
            out[room][name] = rec
    return out


def run_q2():
    """Spatial information rate under a fixed sounding budget."""
    w, h = Q1_ROOMS["reference 5x4 m"]
    rows = {}
    for n in Q2_N_VALUES:
        ndf, n_links = _ndf(spread(n, w, h), w, h)
        eff = ndf / n_links if n_links else 0.0
        f_link = TESTBED_AGGREGATE_HZ / n_links
        rows[str(n)] = {
            "ndf": ndf, "links": n_links, "efficiency": eff,
            "per_link_hz": f_link,
            "info_rate": TESTBED_AGGREGATE_HZ * eff,
            "feasible_presence": f_link >= F_MIN_PRESENCE_HZ,
            "feasible_gait": f_link >= F_MIN_GAIT_HZ,
        }
    caps = {}
    for label, f_min in (("presence", F_MIN_PRESENCE_HZ), ("gait", F_MIN_GAIT_HZ)):
        l_max = TESTBED_AGGREGATE_HZ / f_min
        n_max = int(np.floor((1 + np.sqrt(1 + 8 * l_max)) / 2))
        caps[label] = {"f_min_hz": f_min, "max_links": l_max, "max_nodes": n_max}
    return {"aggregate_hz": TESTBED_AGGREGATE_HZ, "per_n": rows, "caps": caps}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / RESULTS_JSON)
    ap.add_argument("--dry-run", action="store_true",
                    help="report the largest matrix that will be built")
    args = ap.parse_args()

    if args.dry_run:
        worst = max(Q2_N_VALUES)
        n_links = worst * (worst - 1) // 2
        w, h = Q1_ROOMS["reference 5x4 m"]
        res = float(np.clip(np.sqrt(w * h / (VOXEL_LINK_RATIO * n_links)),
                            MIN_GRID_RES_M, MAX_GRID_RES_M))
        v = int(w / res) * int(h / res)
        print(f"largest dense matrix {n_links} x {v} = "
              f"{n_links * v * 8 / 1e6:.1f} MB; SVD workspace ~3x that")
        return

    q1, q2 = run_q1(), run_q2()

    print("=" * 78)
    print("Q1  Does placement matter BELOW N*, where the installed base lives?")
    print("=" * 78)
    for room, layouts in q1.items():
        print(f"\n  {room}    NDF/L by layout")
        print(f"    {'N':>3}{'L':>5}" + "".join(f"{k:>16}" for k in Q1_LAYOUTS))
        for n in Q1_N_VALUES:
            links = layouts["Spread"][str(n)]["links"]
            row = "".join(f"{layouts[k][str(n)]['efficiency']:>16.3f}"
                          for k in Q1_LAYOUTS)
            print(f"    {n:>3}{links:>5}" + row)
        worst = min(layouts[k][str(3)]["efficiency"] for k in Q1_LAYOUTS)
        best = max(layouts[k][str(3)]["efficiency"] for k in Q1_LAYOUTS)
        print(f"    3-node kit: worst {worst:.3f} vs best {best:.3f}"
              f"  ({'no penalty' if best - worst < 1e-9 else f'{best/max(worst,1e-9):.2f}x'})")

    print("\n" + "=" * 78)
    print("Q2  Spatial information RATE under a fixed sounding budget")
    print("=" * 78)
    print(f"\n  aggregate budget {q2['aggregate_hz']:.0f} Hz (the paper's testbed)")
    for label, c in q2["caps"].items():
        print(f"  {label:>9}: needs {c['f_min_hz']:>4.0f} Hz/link -> at most "
              f"{c['max_links']:.0f} links -> N <= {c['max_nodes']}")
    print(f"\n    {'N':>3}{'L':>6}{'NDF':>6}{'NDF/L':>8}{'Hz/link':>9}"
          f"{'indep/s':>10}  feasible")
    for n in Q2_N_VALUES:
        r = q2["per_n"][str(n)]
        flags = ("presence " if r["feasible_presence"] else "") + \
                ("gait" if r["feasible_gait"] else "")
        print(f"    {n:>3}{r['links']:>6}{r['ndf']:>6}{r['efficiency']:>8.3f}"
              f"{r['per_link_hz']:>9.2f}{r['info_rate']:>10.1f}  {flags or '-'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"q1_small_layouts": q1, "q2_airtime": q2,
                                    "tau": TAU}, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
