#!/usr/bin/env python3
"""Is the interior-mesh growth real, or is it under-resolution?

The paper's one qualitative claim is an ASYMMETRY: perimeter (1-D curve) meshes
saturate, interior (2-D area) meshes keep growing. Only the perimeter arm was
grid-refined. That is the wrong arm to check. Under-resolution inflates rank on
the arm with the shortest links, because a voxel grid coarser than the Fresnel
half-width cannot represent two overlapping zones as overlapping, so it reports
them as independent. Interior layouts have the shorter links. The untested arm
is therefore the one biased toward the claimed result.

sigma_F = r_F/2 = sqrt(lambda*d/4)/2. At lambda = 0.125 m a 0.5 m link has
sigma_F = 6.3 cm, so a grid at the paper's coarse end genuinely cannot resolve
it. This script refines both arms together until NDF stops moving.

Memory (Rule 7.5): the sensing matrix is sparse and L < V here, so we never
densify it. Singular values come from eigenvalues of the L x L Gram matrix
G = M M^T, computed sparsely. At N = 90, G is 4005^2 x 8 B = 128 MB, versus
707 MB for the dense M and ~2 GB of SVD workspace. Eigenvalues of G are
sigma^2, so a relative threshold tau on sigma is tau^2 on the eigenvalues;
with tau = 0.01 that is 1e-4, far above double-precision noise.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import resource
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "opencsi" / "src"))

from opencsi.geometry.ndf import build_sensing_matrix  # noqa: E402

RESULTS_DIR = Path(os.environ.get(
    "NDF_RESULTS_DIR", _REPO_ROOT / "docs/4-research/results"))
RESULTS_JSON = "ndf_interior_convergence.json"

TAU = 0.01
WAVELENGTH_M = 0.125
ROOM_W_M, ROOM_H_M = 5.0, 4.0
WALL_INSET_M = 0.10
MIN_SEPARATION_M = 0.25
RNG_SEED = 20260731

# Refinement ladder. The coarse end is what the paper used.
GRID_RESOLUTIONS_M = (0.10, 0.07, 0.05, 0.035, 0.025)
N_VALUES = (20, 40, 60, 90)
CONVERGENCE_TOL = 0.02           # relative NDF change counted as converged
FLAT_TOL = 0.02                  # relative growth below this counts as flat

RSS_WARN_FRACTION = 0.50
_TOTAL_RAM_B = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")


def _rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def _check_rss(where):
    rss = _rss_gb()
    if rss > RSS_WARN_FRACTION * _TOTAL_RAM_B / 1e9:
        print(f"  !! peak RSS {rss:.1f} GB at {where} "
              f"({100*rss*1e9/_TOTAL_RAM_B:.0f}% of RAM)", flush=True)
    return rss


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


def interior(n, rng, w=ROOM_W_M, h=ROOM_H_M):
    """Poisson-disk interior placement with a physical minimum separation."""
    pts = []
    for _ in range(200000):
        if len(pts) >= n:
            break
        c = np.array([rng.uniform(WALL_INSET_M, w - WALL_INSET_M),
                      rng.uniform(WALL_INSET_M, h - WALL_INSET_M)])
        if not pts or np.min(np.linalg.norm(np.array(pts) - c, axis=1)) >= MIN_SEPARATION_M:
            pts.append(c)
    if len(pts) < n:
        raise RuntimeError(f"could only place {len(pts)}/{n} nodes")
    return np.array(pts)


def ndf_sparse(pos, res, w=ROOM_W_M, h=ROOM_H_M, lam=WAVELENGTH_M):
    """NDF without ever densifying: eigenvalues of the sparse Gram matrix."""
    m, _, _ = build_sensing_matrix(pos, grid_resolution_m=res,
                                   room_bounds=((0.0, 0.0), (w, h)),
                                   wavelength_m=lam)
    n_links, n_vox = m.shape
    gram = (m @ m.T).toarray()
    del m
    gc.collect()
    ev = np.linalg.eigvalsh(gram)
    del gram
    gc.collect()
    ev = np.clip(ev[::-1], 0.0, None)
    sig = np.sqrt(ev)
    k = int((sig > TAU * sig[0]).sum()) if sig[0] > 0 else 0
    return k, n_links, n_vox


def shortest_link_sigma_f(pos, lam=WAVELENGTH_M):
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    dmin = d[np.triu_indices(len(pos), k=1)].min()
    return np.sqrt(lam * dmin / 4.0) / 2.0, dmin


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / RESULTS_JSON)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rng = np.random.default_rng(RNG_SEED)
    layouts = {"perimeter": {n: perimeter(n) for n in N_VALUES},
               "interior": {n: interior(n, rng) for n in N_VALUES}}

    if args.dry_run:
        print(f"machine RAM {_TOTAL_RAM_B/1e9:.1f} GB")
        n = max(N_VALUES)
        L = n * (n - 1) // 2
        res = min(GRID_RESOLUTIONS_M)
        v = int(ROOM_W_M / res) * int(ROOM_H_M / res)
        print(f"worst case N={n}: L={L}, V={v}")
        print(f"  dense M would be {L*v*8/1e9:.2f} GB  (avoided)")
        print(f"  sparse Gram L x L is {L*L*8/1e9:.3f} GB  (used)")
        for arm, per_n in layouts.items():
            s, d = shortest_link_sigma_f(per_n[n])
            print(f"  {arm:>9}: shortest link {d:.2f} m -> sigma_F {s*100:.1f} cm")
        return

    out = {"taus": TAU, "resolutions_m": list(GRID_RESOLUTIONS_M),
           "n_values": list(N_VALUES), "arms": {}}

    for arm, per_n in layouts.items():
        out["arms"][arm] = {}
        sf, dmin = shortest_link_sigma_f(per_n[max(N_VALUES)])
        out["arms"][arm]["shortest_link_m"] = float(dmin)
        out["arms"][arm]["sigma_f_m"] = float(sf)
        for n in N_VALUES:
            rec = {}
            for res in GRID_RESOLUTIONS_M:
                k, L, V = ndf_sparse(per_n[n], res)
                rec[str(res)] = {"ndf": k, "links": L, "voxels": V,
                                 "efficiency": k / L}
                print(f"  {arm:>9} N={n:>3} res={res:.3f}  NDF={k:>5}/{L:<5}"
                      f" V={V:>6}  RSS={_check_rss(f'{arm} N={n}'):.1f} GB",
                      flush=True)
            out["arms"][arm][str(n)] = rec

    # ---- verdicts ----------------------------------------------------------
    print("\n" + "=" * 78)
    print("Does NDF converge as the grid refines?")
    print("=" * 78)
    finest = str(GRID_RESOLUTIONS_M[-1])
    second = str(GRID_RESOLUTIONS_M[-2])
    for arm in layouts:
        print(f"\n  {arm}")
        print(f"    {'N':>4}" + "".join(f"{r:>10}" for r in GRID_RESOLUTIONS_M)
              + "   converged?")
        for n in N_VALUES:
            rec = out["arms"][arm][str(n)]
            a, b = rec[second]["ndf"], rec[finest]["ndf"]
            drift = abs(b - a) / max(a, 1)
            print(f"    {n:>4}" + "".join(f"{rec[str(r)]['ndf']:>10}"
                                          for r in GRID_RESOLUTIONS_M)
                  + f"   {'yes' if drift <= CONVERGENCE_TOL else f'NO ({drift:.1%})'}")

    print("\n" + "=" * 78)
    print("Does the ASYMMETRY survive refinement?")
    print("=" * 78)
    print(f"\n    {'res':>7}{'perim growth':>15}{'interior growth':>17}  asymmetry?")
    verdict = {}
    for res in GRID_RESOLUTIONS_M:
        g = {}
        for arm in layouts:
            a = out["arms"][arm][str(N_VALUES[-2])][str(res)]["ndf"]
            b = out["arms"][arm][str(N_VALUES[-1])][str(res)]["ndf"]
            g[arm] = (b - a) / max(a, 1)
        holds = g["perimeter"] < FLAT_TOL <= g["interior"]
        verdict[str(res)] = {"perimeter_growth": g["perimeter"],
                             "interior_growth": g["interior"], "holds": holds}
        print(f"    {res:>7.3f}{g['perimeter']:>15.3f}{g['interior']:>17.3f}"
              f"  {'yes' if holds else 'NO'}")
    out["asymmetry_by_resolution"] = verdict
    survives = verdict[finest]["holds"]
    out["asymmetry_survives_refinement"] = survives
    print(f"\nVERDICT at the finest grid ({finest} m): "
          f"asymmetry {'SURVIVES' if survives else 'DOES NOT SURVIVE'}")
    print(f"peak RSS {_rss_gb():.2f} GB")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
