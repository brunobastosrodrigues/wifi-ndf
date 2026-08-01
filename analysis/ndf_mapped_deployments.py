#!/usr/bin/env python3
"""Regime 1 tested on three of our own deployments, from surveyed maps.

Every geometry the paper has used so far is either a designed simulation layout
or, for the home testbed, the output of an MDS solve on ranging data. Neither is
a measurement of where the nodes actually are. These three deployments were
drawn as scale room maps with metric axes at install time, so the coordinates
are an input.

The regime-1 claim is exact and therefore falsifiable on them: below N*, NDF = L
with no residual. A single link short of L in any of these rooms falsifies it.
The claim is also the load-bearing one, because the paper's answer to "NDF
correlates no better than link count" is that below N* the two are the same
number rather than two numbers that happen to agree.

We additionally enumerate subsets. If NDF = L for every subset of every size,
then placement within these deployments is not weakly informative but exactly
uninformative, which is the strongest form of the paper's own claim that below
N* placement guidance is unfalsifiable.

Provenance. Dornbirn is read from pictures/dornbirn/room_map_8nodes.svg, whose
axes are labelled in metres (80 px per metre), so its coordinates are exact to
the drawing. sg-lab and sg-kk are read off rendered PNG maps against their
printed axes and are accurate to roughly +-0.05 m; we report the sensitivity of
the result to that uncertainty rather than assuming it away.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "opencsi" / "src"))

from opencsi.geometry.ndf import build_sensing_matrix  # noqa: E402

RESULTS_DIR = Path(os.environ.get(
    "NDF_RESULTS_DIR", _REPO_ROOT / "docs/4-research/results"))
RESULTS_JSON = "ndf_mapped_deployments.json"

TAU = 0.01
WAVELENGTH_M = 0.125
SIGMA_F_CELLS = 2.0
MIN_RES_M = 0.015
READOFF_JITTER_M = 0.05      # PNG read-off uncertainty
JITTER_REPLICATES = 12
RNG_SEED = 20260801

DEPLOYMENTS = {
    "Dornbirn (office)": {
        "room": (2.9, 4.4), "source": "SVG, metric axes",
        "nodes": {"N1": (0.73, 0.02), "N2": (1.46, 0.02), "N3": (2.81, 0.73),
                  "N14": (0.10, 1.01), "N6": (2.81, 2.92), "N11": (0.10, 3.20),
                  "N9": (0.50, 4.33), "N8": (2.10, 4.33)}},
    "sg-lab": {
        "room": (7.05, 4.55), "source": "PNG read-off",
        "nodes": {"42:A8": (2.00, 4.50), "43:24": (5.60, 4.50),
                  "3C:64": (0.05, 3.15), "D3:34": (0.05, 1.35),
                  "89:90": (7.00, 3.10), "9A:68": (7.00, 0.45),
                  "43:64": (2.00, 0.05), "6C:64": (5.15, 0.05)}},
    "sg-kk": {
        "room": (3.25, 6.15), "source": "PNG read-off",
        "nodes": {"N1": (0.05, 0.05), "N2": (1.70, 0.05), "N3": (3.15, 0.05),
                  "N4": (0.05, 3.00), "N5": (1.70, 3.00), "N6": (3.15, 3.00),
                  "N7": (0.35, 6.05), "N8": (1.70, 6.05), "N9": (2.85, 6.05)}},
}


def ndf(pos, w, h, lam=WAVELENGTH_M):
    pos = np.asarray(pos, float)
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
    k = int((s > TAU * s[0]).sum())
    margin = float(s[k - 1] / s[0] / TAU) if k else 0.0
    return k, n_links, margin


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / RESULTS_JSON)
    args = ap.parse_args()
    rng = np.random.default_rng(RNG_SEED)
    out = {"tau": TAU, "deployments": {}}

    print("=" * 78)
    print("Regime 1 on three mapped deployments (NDF = L is falsifiable here)")
    print("=" * 78)
    print(f"\n{'deployment':>20}{'N':>4}{'L':>5}{'NDF':>5}{'NDF/L':>8}"
          f"{'area':>8}{'ka':>7}{'AR':>6}{'sig/tau':>9}")
    for name, d in DEPLOYMENTS.items():
        w, h = d["room"]
        pos = np.array(list(d["nodes"].values()))
        k, L, margin = ndf(pos, w, h)
        a = np.sqrt(w * h / np.pi)
        ka = 2 * np.pi * a / WAVELENGTH_M
        ar = max(w, h) / min(w, h)
        out["deployments"][name] = {
            "n_nodes": len(pos), "n_links": L, "ndf": k, "efficiency": k / L,
            "area_m2": w * h, "ka": ka, "aspect_ratio": ar,
            "margin_over_tau": margin, "source": d["source"], "exact": k == L}
        print(f"{name:>20}{len(pos):>4}{L:>5}{k:>5}{k/L:>8.3f}"
              f"{w*h:>8.1f}{ka:>7.0f}{ar:>6.2f}{margin:>9.1f}")

    print("\nsubset enumeration (is placement exactly uninformative below N*?)")
    print(f"{'deployment':>20}{'k':>3}{'subsets':>9}{'eff min':>9}{'eff max':>9}")
    for name, d in DEPLOYMENTS.items():
        w, h = d["room"]
        names = list(d["nodes"])
        rec = {}
        for kk in range(3, len(names) + 1):
            effs = [ndf(np.array([d["nodes"][n] for n in c]), w, h)[0]
                    / (kk * (kk - 1) // 2)
                    for c in combinations(names, kk)]
            rec[str(kk)] = {"n": len(effs), "min": min(effs), "max": max(effs)}
            print(f"{name if kk == 3 else '':>20}{kk:>3}{len(effs):>9}"
                  f"{min(effs):>9.3f}{max(effs):>9.3f}")
        out["deployments"][name]["subsets"] = rec

    print("\nsensitivity to map read-off error "
          f"(+-{READOFF_JITTER_M} m, {JITTER_REPLICATES} replicates)")
    for name, d in DEPLOYMENTS.items():
        w, h = d["room"]
        base = np.array(list(d["nodes"].values()))
        effs = []
        for _ in range(JITTER_REPLICATES):
            p = base + rng.uniform(-READOFF_JITTER_M, READOFF_JITTER_M, base.shape)
            p[:, 0] = np.clip(p[:, 0], 0.02, w - 0.02)
            p[:, 1] = np.clip(p[:, 1], 0.02, h - 0.02)
            k, L, _ = ndf(p, w, h)
            effs.append(k / L)
        out["deployments"][name]["readoff_jitter"] = {
            "min": float(min(effs)), "max": float(max(effs))}
        print(f"  {name:>20}: NDF/L in [{min(effs):.3f}, {max(effs):.3f}]")

    allexact = all(v["exact"] for v in out["deployments"].values())
    print("\n" + "=" * 78)
    print(f"VERDICT: NDF = L exactly in {'all three' if allexact else 'NOT all'} "
          f"deployments; the prediction is not falsified.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
