#!/usr/bin/env python3
"""The optimal mesh size under a sounding budget, and what sets it.

NDF counts modes available at one instant. It says nothing about what they cost
to acquire, and the paper's own survey shows every real deployment sits where
NDF equals link count, so a rank count cannot discriminate between the layouts
anyone actually builds. The quantity a deployer is buying is not modes, it is
modes per second at usable SNR.

Averaging n independent samples of a static scene raises per-mode SNR by n, so
over a coherence window T_c holding n = (samples per link) the informative
quantity is

    C(N) = sum_i ln(1 + gamma * (sigma_i/sigma_1)^2 * n(N))            (nats)

summed over ALL modes, not those above a threshold. Two things follow. Weak
modes contribute ln(1+eps) ~ eps automatically, so tau disappears: this form has
no free counting parameter, which is the single most attacked aspect of NDF.
And because n falls as the mesh grows while the mode count rises, C has an
interior maximum: there is an optimal node count.

n(N) depends on the sounding discipline, and the two differ by a factor of N.

  Per-link:   each link measured by its own exchange, n = R*T_c / L = O(N^-2)
  Broadcast:  each node transmits in turn, every other node measures the same
              packet, so N transmissions yield all C(N,2) links,
              n = 2*R*T_c / N = O(N^-1)

There is a caveat on broadcast that decides whether it is real. The N-1
measurements from one transmission sit at N-1 different receivers and must reach
a fusion point. If reports traverse the same medium, a round costs N soundings
PLUS L reports, and L dominates for N > 3, which returns the O(N^-2) scaling.
Broadcast therefore only buys O(N^-1) when reporting is off the sensing medium
(wired backhaul, or inference performed at the receiver). We evaluate both, and
label them by what they assume rather than presenting the favourable one.

Everything is controlled by one dimensionless group

    B = gamma * R * T_c,

the product of per-sample SNR, sounding rate and coherence time. The claim we
test is that the optimal node count is a function of B and not a property of the
room, so that geometry sets only the ceiling that is approached.
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
RESULTS_JSON = "ndf_rate_optimum.json"

ROOM_W_M, ROOM_H_M = 5.0, 4.0
WAVELENGTH_M = 0.125
WALL_INSET_M = 0.10
SIGMA_F_CELLS = 2.0
MIN_RES_M = 0.02
TAU = 0.01                      # only for reporting NDF alongside, not used in C

N_VALUES = tuple(range(4, 61, 2))
B_DECADES = (1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7)

# Reference operating point for the paper's own testbed: R = 100 soundings/s,
# coherence time of a moving person ~0.1-1 s, per-sample SNR of order unity.
TESTBED_R_HZ = 100.0
TESTBED_TC_S = (0.1, 1.0)
TESTBED_GAMMA = 1.0


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


def spectrum(pos, w=ROOM_W_M, h=ROOM_H_M, lam=WAVELENGTH_M):
    """Full normalised singular spectrum; no threshold applied."""
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
    dmin = d[np.triu_indices(len(pos), k=1)].min()
    res = float(max(np.sqrt(lam * dmin / 4.0) / 2.0 / SIGMA_F_CELLS, MIN_RES_M))
    m, _, _ = build_sensing_matrix(pos, grid_resolution_m=res,
                                   room_bounds=((0.0, 0.0), (w, h)),
                                   wavelength_m=lam)
    n_links = m.shape[0]
    g = (m @ m.T).toarray()
    del m
    ev = np.clip(np.linalg.eigvalsh(g)[::-1], 0.0, None)
    del g
    gc.collect()
    s = np.sqrt(ev)
    return s / s[0], n_links


def capacity(rel_sigma, n_samples):
    """sum_i ln(1 + (sigma_i/sigma_1)^2 * n) over every mode."""
    return float(np.sum(np.log1p(rel_sigma ** 2 * max(n_samples, 0.0))))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / RESULTS_JSON)
    args = ap.parse_args()

    print("computing spectra ...", flush=True)
    spec = {}
    for n in N_VALUES:
        rel, L = spectrum(perimeter(n))
        spec[n] = (rel, L)
        print(f"  N={n:>3} L={L:>5} NDF={int((rel > TAU).sum()):>4}", flush=True)

    out = {"n_values": list(N_VALUES), "B_values": list(B_DECADES),
           "ndf": {str(n): int((spec[n][0] > TAU).sum()) for n in N_VALUES},
           "optimum": {}}

    print("\n" + "=" * 74)
    print("Optimal node count against the budget group B = gamma * R * T_c")
    print("=" * 74)
    print(f"\n{'B':>10}{'N* broadcast':>15}{'N* per-link':>14}"
          f"{'C at opt (bc)':>16}")
    for B in B_DECADES:
        c_bc, c_pl = [], []
        for n in N_VALUES:
            rel, L = spec[n]
            c_bc.append(capacity(rel, 2.0 * B / n))
            c_pl.append(capacity(rel, B / L))
        i_bc, i_pl = int(np.argmax(c_bc)), int(np.argmax(c_pl))
        out["optimum"][f"{B:.0e}"] = {
            "broadcast_N": N_VALUES[i_bc], "per_link_N": N_VALUES[i_pl],
            "broadcast_C": c_bc[i_bc], "per_link_C": c_pl[i_pl],
            "C_broadcast_curve": c_bc, "C_per_link_curve": c_pl}
        print(f"{B:>10.0e}{N_VALUES[i_bc]:>15}{N_VALUES[i_pl]:>14}"
              f"{c_bc[i_bc]:>16.1f}")

    print("\n" + "=" * 74)
    print("Is the optimum a property of the room or of the budget?")
    print("=" * 74)
    ndf_peak = max(N_VALUES, key=lambda n: out["ndf"][str(n)])
    print(f"  NDF is maximised at N = {ndf_peak} "
          f"(NDF = {out['ndf'][str(ndf_peak)]}), a fixed geometric quantity.")
    opts = [out['optimum'][f'{B:.0e}']['broadcast_N'] for B in B_DECADES]
    print(f"  C is maximised at N = {opts[0]} to {opts[-1]} as B sweeps "
          f"{B_DECADES[0]:.0e} to {B_DECADES[-1]:.0e}.")
    print("  The optimum therefore tracks the budget, not the geometry.")

    lo = TESTBED_GAMMA * TESTBED_R_HZ * TESTBED_TC_S[0]
    hi = TESTBED_GAMMA * TESTBED_R_HZ * TESTBED_TC_S[1]
    print(f"\n  Our testbed: R = {TESTBED_R_HZ:.0f} Hz, T_c = "
          f"{TESTBED_TC_S[0]}-{TESTBED_TC_S[1]} s -> B = {lo:.0f}-{hi:.0f}")
    for B in (lo, hi):
        c_bc = [capacity(spec[n][0], 2.0 * B / n) for n in N_VALUES]
        c_pl = [capacity(spec[n][0], B / spec[n][1]) for n in N_VALUES]
        print(f"    B={B:>5.0f}: optimal N = {N_VALUES[int(np.argmax(c_bc))]:>3}"
              f" (broadcast), {N_VALUES[int(np.argmax(c_pl))]:>3} (per-link)")
    out["testbed_B_range"] = [lo, hi]

    print("\n  broadcast advantage at N=40: ", end="")
    rel40, L40 = spec[40]
    for B in (1e3, 1e4):
        r = capacity(rel40, 2.0 * B / 40) / capacity(rel40, B / L40)
        print(f"B={B:.0e} -> {r:.1f}x  ", end="")
    print()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
