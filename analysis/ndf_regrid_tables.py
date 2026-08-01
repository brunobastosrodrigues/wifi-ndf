#!/usr/bin/env python3
"""Recompute the paper's small-N tables under a PHYSICAL grid criterion.

TMC R2: the paper's grid rule, V/L >= 4, is not a convergence criterion. It
ties resolution to link count, so it hands the largest meshes the finest grids
and leaves Tables I, II and IV (N = 14-34) the least resolved results in the
paper. At N = 14 it yields 0.234 m while the narrowest Fresnel zone in the
layout has sigma_F = 0.177 m, i.e. 0.75 cells across the beam. The grid does
not resolve the quantity the model is built on.

The criterion should be absolute and set by the physics: the narrowest zone in
the layout is the one that must be resolved, so

    res <= sigma_F(d_min)/2,   sigma_F(d) = sqrt(lambda*d/4)/2,

with d_min the shortest link. Note the direction of the error this fixes.
Coarsening MERGES near-degenerate directions, so it is rank-LOWERING: it can
manufacture saturation but it cannot manufacture full rank. The layouts it
damages are therefore the near-degenerate ones, which is precisely the failure
modes of Table II, while the well-spread rows are near-orthogonal and barely
move. The paper's convergence check was run on the configurations least able
to expose the problem.

We also report two diagnostics the paper currently omits, both of which turn a
threshold crossing into a falsifiable statement:

  sigma_min/sigma_1   margin by which the last retained mode clears tau.
                      NDF = L is only as strong as this ratio.
  near-tau density    how many singular values sit within a factor of 2 of tau.
                      A layout with many is one whose NDF is a tau artifact;
                      a layout with none has a genuine spectral gap.
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
RESULTS_JSON = "ndf_regrid_tables.json"

TAU = 0.01
NEAR_TAU_FACTOR = 2.0
WAVELENGTH_M = 0.125
ROOM_W_M, ROOM_H_M = 5.0, 4.0
WALL_INSET_M = 0.10
TABLE2_N = 14

# the paper's rule, kept only to quantify what it cost
LEGACY_VOXEL_LINK_RATIO = 4.0
LEGACY_MIN_RES_M, LEGACY_MAX_RES_M = 0.02, 0.25

# the physical rule
SIGMA_F_CELLS = 2.0              # cells per Fresnel std dev
ABS_MIN_RES_M = 0.015            # floor, to bound cost
TAU_INTERVAL = (0.005, 0.02)     # the interval the ablation already concedes


def sigma_f(d, lam=WAVELENGTH_M):
    return np.sqrt(lam * d / 4.0) / 2.0


def shortest_link(pos):
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    return d[np.triu_indices(len(pos), k=1)].min()


def legacy_res(n_links, w=ROOM_W_M, h=ROOM_H_M):
    return float(np.clip(np.sqrt(w * h / (LEGACY_VOXEL_LINK_RATIO * n_links)),
                         LEGACY_MIN_RES_M, LEGACY_MAX_RES_M))


def physical_res(pos):
    return float(max(sigma_f(shortest_link(pos)) / SIGMA_F_CELLS, ABS_MIN_RES_M))


# ---- layouts ---------------------------------------------------------------
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


TABLE2_LAYOUTS = {"Well-spread": perimeter, "Two groups": two_groups,
                  "Corner cluster": corner_cluster, "One wall": one_wall}


def spectrum(pos, res, w=ROOM_W_M, h=ROOM_H_M, lam=WAVELENGTH_M):
    m, _, _ = build_sensing_matrix(pos, grid_resolution_m=res,
                                   room_bounds=((0.0, 0.0), (w, h)),
                                   wavelength_m=lam)
    n_links, n_vox = m.shape
    gram = (m @ m.T).toarray()
    del m
    gc.collect()
    ev = np.clip(np.linalg.eigvalsh(gram)[::-1], 0.0, None)
    del gram
    gc.collect()
    return np.sqrt(ev), n_links, n_vox


def diagnose(sig, n_links, tau=TAU):
    rel = sig / sig[0]
    k = int((rel > tau).sum())
    sig_min_ratio = float(rel[k - 1]) if k else 0.0
    near = int(((rel >= tau / NEAR_TAU_FACTOR) &
                (rel <= tau * NEAR_TAU_FACTOR)).sum())
    swing = [float((rel > t).sum()) / n_links for t in TAU_INTERVAL]
    return {"ndf": k, "links": n_links, "efficiency": k / n_links,
            "sigma_min_over_sigma1": sig_min_ratio,
            "margin_over_tau": sig_min_ratio / tau,
            "near_tau_count": near,
            "efficiency_at_tau_interval": swing,
            "efficiency_swing": abs(swing[0] - swing[1])}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / RESULTS_JSON)
    args = ap.parse_args()
    out = {"tau": TAU, "criterion": "res <= sigma_F(d_min)/2"}

    # ---- Table II under both grid rules ------------------------------------
    print("=" * 82)
    print(f"Table II (N={TABLE2_N}): what the grid rule cost")
    print("=" * 82)
    t2 = {}
    print(f"\n  {'layout':>15}{'d_min':>7}{'sigma_F':>9}{'legacy':>8}{'phys':>7}"
          f"{'NDF/L legacy':>14}{'NDF/L phys':>12}")
    for name, fn in TABLE2_LAYOUTS.items():
        pos = fn(TABLE2_N)
        n_links = TABLE2_N * (TABLE2_N - 1) // 2
        r_leg, r_phy = legacy_res(n_links), physical_res(pos)
        dmin = shortest_link(pos)
        s_l, L, _ = spectrum(pos, r_leg)
        s_p, _, _ = spectrum(pos, r_phy)
        d_l, d_p = diagnose(s_l, L), diagnose(s_p, L)
        t2[name] = {"d_min_m": float(dmin), "sigma_f_m": float(sigma_f(dmin)),
                    "res_legacy_m": r_leg, "res_physical_m": r_phy,
                    "legacy": d_l, "physical": d_p}
        print(f"  {name:>15}{dmin:>7.2f}{sigma_f(dmin):>9.3f}{r_leg:>8.3f}"
              f"{r_phy:>7.3f}{d_l['efficiency']:>14.3f}{d_p['efficiency']:>12.3f}")
    out["table2"] = t2

    ws = t2["Well-spread"]["physical"]["efficiency"]
    ow_l = t2["One wall"]["legacy"]["efficiency"]
    ow_p = t2["One wall"]["physical"]["efficiency"]
    print(f"\n  abstract claim 'one wall wastes {100*(1-ow_l):.0f}% of links' "
          f"-> converged {100*(1-ow_p):.0f}%")
    print(f"  ordering preserved: best {ws:.3f} vs worst {ow_p:.3f} "
          f"({ws/ow_p:.2f}x)")

    print(f"\n  threshold diagnostics (the paper reports neither)")
    print(f"  {'layout':>15}{'sig_min/sig1':>14}{'x over tau':>12}"
          f"{'near-tau':>10}{'eff swing':>11}")
    for name in TABLE2_LAYOUTS:
        d = t2[name]["physical"]
        print(f"  {name:>15}{d['sigma_min_over_sigma1']:>14.4f}"
              f"{d['margin_over_tau']:>12.1f}{d['near_tau_count']:>10}"
              f"{d['efficiency_swing']:>11.3f}")

    # ---- N*: the exactness boundary, per layout family ----------------------
    print("\n" + "=" * 82)
    print("N*: largest N with NDF == L exactly, under the physical grid")
    print("=" * 82)
    families = {"Well-spread": perimeter, "One wall": one_wall,
                "Corner cluster": corner_cluster}
    nstar = {}
    for name, fn in families.items():
        rec, last_exact = {}, None
        for n in range(3, 21):
            pos = fn(n)
            sig, L, _ = spectrum(pos, physical_res(pos))
            d = diagnose(sig, L)
            rec[str(n)] = d
            if d["ndf"] == L:
                last_exact = n
        nstar[name] = {"per_n": rec, "n_star_exact": last_exact}
        print(f"\n  {name}: NDF == L exactly up to N = {last_exact}")
        print(f"    {'N':>4}{'L':>6}{'NDF':>6}{'NDF/L':>9}{'sig_min/sig1':>14}")
        for n in range(3, 21):
            d = rec[str(n)]
            flag = " *" if d["ndf"] == d["links"] else ""
            print(f"    {n:>4}{d['links']:>6}{d['ndf']:>6}{d['efficiency']:>9.4f}"
                  f"{d['sigma_min_over_sigma1']:>14.4f}{flag}")
    out["n_star"] = nstar

    print("\n" + "=" * 82)
    print("N* is placement-dependent: " + ", ".join(
        f"{k} = {v['n_star_exact']}" for k, v in nstar.items()))
    print("=" * 82)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
