#!/usr/bin/env python3
"""Threshold ablation: do the paper's CONCLUSIONS survive the choice of tau?

GLOBECOM 2026 reviewer R5: "The threshold used to define NDF, namely counting
singular values above 1% of the largest singular value, is not fully justified.
Although a limited sensitivity comparison is provided, a more systematic
threshold analysis would improve confidence in the metric."

Reporting that NDF itself varies with tau answers a weaker question than the one
that matters. tau is the single free parameter in the definition, so the test
that earns confidence is whether every claim the paper makes is invariant to it.
We therefore re-derive three conclusions across two decades of tau:

  C1  Layout ordering. Is one-wall always worst and well-spread always best?
      (Sec. 5.2: collinearity destroys efficiency, clustering does not.)
  C2  Saturation asymmetry. Do perimeter meshes always plateau while interior
      meshes always keep growing? (Sec. 5.3.)
  C3  Deployment classification. Do the published deployments stay above the
      first boundary? (Sec. 6.2.)

A conclusion that flips inside the sweep is reported as flipped.
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
RESULTS_JSON = "ndf_threshold_ablation.json"

TAUS = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10)
NOMINAL_TAU = 0.01

ROOM_W_M, ROOM_H_M = 5.0, 4.0
WAVELENGTH_M = 0.125
WALL_INSET_M = 0.10
MIN_SEPARATION_M = 0.25          # physically realisable layouts only
VOXEL_LINK_RATIO = 4.0
MIN_GRID_RES_M = 0.015
SIGMA_F_CELLS = 2.0            # cells per Fresnel std dev (physical grid rule)
MAX_GRID_RES_M = 0.25

C1_N = 14
C2_N_VALUES = (20, 40, 60, 90)
C2_FLAT_TOL = 0.02               # relative growth below this counts as flat
REGIME1_RATIO = 0.95


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


def one_wall(n, w=ROOM_W_M, h=ROOM_H_M):
    return np.array([(x, WALL_INSET_M)
                     for x in np.linspace(WALL_INSET_M, w - WALL_INSET_M, n)])


def corner_cluster(n, w=ROOM_W_M, h=ROOM_H_M):
    side = int(np.ceil(np.sqrt(n)))
    xs = np.linspace(WALL_INSET_M, WALL_INSET_M + 1.8, side)
    ys = np.linspace(WALL_INSET_M, WALL_INSET_M + 1.8, side)
    return np.array([(x, y) for y in ys for x in xs])[:n]


def two_groups(n, w=ROOM_W_M, h=ROOM_H_M):
    half = n // 2
    out = []
    for k, x0 in ((half, WALL_INSET_M), (n - half, w - WALL_INSET_M - 0.8)):
        side = int(np.ceil(np.sqrt(k)))
        xs = np.linspace(x0, x0 + 0.8, side)
        ys = np.linspace(WALL_INSET_M, h - WALL_INSET_M, side)
        out += [(x, y) for y in ys for x in xs][:k]
    return np.array(out)


def poisson_disk(n, rng, w=ROOM_W_M, h=ROOM_H_M):
    pts = []
    for _ in range(20000):
        if len(pts) >= n:
            break
        c = np.array([rng.uniform(WALL_INSET_M, w - WALL_INSET_M),
                      rng.uniform(WALL_INSET_M, h - WALL_INSET_M)])
        if not pts or np.min(np.linalg.norm(np.array(pts) - c, axis=1)) >= MIN_SEPARATION_M:
            pts.append(c)
    return np.array(pts)


C1_LAYOUTS = {"Well-spread": perimeter, "Two groups": two_groups,
              "Corner cluster": corner_cluster, "One wall": one_wall}


def spectrum(pos, w=ROOM_W_M, h=ROOM_H_M, lam=WAVELENGTH_M):
    """Singular values once; every tau is then just a different count.

    Grid set by the physical criterion res <= sigma_F(d_min)/2, not by the
    voxel-per-link ratio. The ratio rule leaves near-degenerate layouts
    under-resolved, and since coarsening merges near-degenerate directions it
    depresses their rank: it is exactly the failure-mode layouts in C1 that it
    would distort.
    """
    pos = np.asarray(pos, dtype=float)
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    d_min = d[np.triu_indices(len(pos), k=1)].min()
    res = float(max(np.sqrt(lam * d_min / 4.0) / 2.0 / SIGMA_F_CELLS,
                    MIN_GRID_RES_M))
    m, _, _ = build_sensing_matrix(pos, grid_resolution_m=res,
                                   room_bounds=((0.0, 0.0), (w, h)),
                                   wavelength_m=lam)
    actual = m.shape[0]
    # Never densify: L < V here, so eigenvalues of the L x L sparse Gram give
    # sigma^2 at a fraction of the memory. At N=90 that is 128 MB against
    # ~1.3 GB dense plus SVD workspace.
    gram = (m @ m.T).toarray(); del m; gc.collect()
    ev = np.clip(np.linalg.eigvalsh(gram)[::-1], 0.0, None); del gram
    gc.collect()
    return np.sqrt(ev), actual


def count(sig, tau):
    return int((sig > tau * sig[0]).sum()) if sig[0] > 0 else 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / RESULTS_JSON)
    args = ap.parse_args()
    out = {"taus": list(TAUS), "nominal_tau": NOMINAL_TAU}

    # ---- C1: does the layout ordering hold at every tau? -------------------
    spectra = {k: spectrum(fn(C1_N)) for k, fn in C1_LAYOUTS.items()}
    c1 = {}
    for tau in TAUS:
        eff = {k: count(s, tau) / L for k, (s, L) in spectra.items()}
        order = sorted(eff, key=lambda k: -eff[k])
        c1[str(tau)] = {"efficiency": eff, "order": order}
    ref_order = c1[str(NOMINAL_TAU)]["order"]
    c1_holds = all(v["order"] == ref_order for v in c1.values())
    out["C1"] = {"per_tau": c1, "reference_order": ref_order, "invariant": c1_holds}

    # ---- C2: does the saturation asymmetry hold at every tau? --------------
    rng = np.random.default_rng(20260731)
    per_spec = {n: spectrum(perimeter(n)) for n in C2_N_VALUES}
    int_spec = {n: spectrum(poisson_disk(n, rng)) for n in C2_N_VALUES}
    c2 = {}
    for tau in TAUS:
        p = [count(per_spec[n][0], tau) for n in C2_N_VALUES]
        i = [count(int_spec[n][0], tau) for n in C2_N_VALUES]
        p_growth = (p[-1] - p[-2]) / max(p[-2], 1)
        i_growth = (i[-1] - i[-2]) / max(i[-2], 1)
        c2[str(tau)] = {"perimeter": p, "interior": i,
                        "perimeter_final_growth": p_growth,
                        "interior_final_growth": i_growth,
                        "asymmetry_holds": p_growth < C2_FLAT_TOL <= i_growth}
    c2_holds = all(v["asymmetry_holds"] for v in c2.values())
    out["C2"] = {"n_values": list(C2_N_VALUES), "per_tau": c2, "invariant": c2_holds}

    # ---- C3: do published deployments stay above the boundary? ------------
    DEPLOY = {"Wilson2010": (28, 4.27, 4.27), "Wilson2011": (34, 8.51, 8.51),
              "Utah-open": (30, 8.37, 8.37), "Utah-apartment": (33, 7.0, 8.25),
              "Utah-office": (32, 8.19, 8.19), "McGill-Trottier": (24, 8.0, 8.0),
              "McGill-CNLab": (24, 9.0, 9.0), "BUPT": (28, 5.2, 6.7),
              "Savazzi(14)": (14, 4.0, 3.0)}
    dep_spec = {k: spectrum(perimeter(n, w, h), w, h) for k, (n, w, h) in DEPLOY.items()}
    c3 = {}
    for tau in TAUS:
        cls = {k: ("1" if count(s, tau) / L >= REGIME1_RATIO else "2")
               for k, (s, L) in dep_spec.items()}
        c3[str(tau)] = cls
    ref_cls = c3[str(NOMINAL_TAU)]
    c3_holds = all(v == ref_cls for v in c3.values())
    out["C3"] = {"per_tau": c3, "reference": ref_cls, "invariant": c3_holds}

    # ---- report -----------------------------------------------------------
    print("=" * 74)
    print("Threshold ablation: are the paper's conclusions invariant to tau?")
    print("=" * 74)
    print(f"\nC1  Layout ordering ({C1_N} nodes)      "
          f"{'INVARIANT' if c1_holds else 'FLIPS'}")
    print(f"    {'tau':>7}" + "".join(f"{k[:13]:>15}" for k in C1_LAYOUTS))
    for tau in TAUS:
        e = c1[str(tau)]["efficiency"]
        mark = "  <- nominal" if tau == NOMINAL_TAU else ""
        print(f"    {tau:>7.3f}" + "".join(f"{e[k]:>15.3f}" for k in C1_LAYOUTS) + mark)

    print(f"\nC2  Saturation asymmetry              "
          f"{'INVARIANT' if c2_holds else 'FLIPS'}")
    print(f"    {'tau':>7}{'perim growth':>15}{'interior growth':>17}  holds?")
    for tau in TAUS:
        v = c2[str(tau)]
        print(f"    {tau:>7.3f}{v['perimeter_final_growth']:>15.3f}"
              f"{v['interior_final_growth']:>17.3f}  "
              f"{'yes' if v['asymmetry_holds'] else 'NO'}")

    print(f"\nC3  Deployment classification         "
          f"{'INVARIANT' if c3_holds else 'FLIPS'}")
    print(f"    at every tau in [{TAUS[0]}, {TAUS[-1]}]: "
          f"{sum(1 for v in ref_cls.values() if v=='2')}/{len(ref_cls)} in regime 2")
    if not c3_holds:
        for tau in TAUS:
            diff = {k: v for k, v in c3[str(tau)].items() if v != ref_cls[k]}
            if diff:
                print(f"    tau={tau}: reclassified {diff}")

    print("\n" + "=" * 74)
    n_ok = sum([c1_holds, c2_holds, c3_holds])
    print(f"VERDICT: {n_ok}/3 conclusions invariant across a 100x range in tau")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
