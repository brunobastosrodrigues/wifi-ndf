#!/usr/bin/env python3
"""Greedy sounding schedules: how much of the round-robin loss is recoverable?

The paper's capacity C(N) assumes the sounding budget is spread uniformly over
links (Assumption 4), and its interior optimum in N is a consequence of that
assumption. The Gaussian information gain is monotone submodular in the set of
soundings, so lazy greedy selection is (1-1/e)-optimal in the worst case
(Golovin & Krause, JAIR 2011) and usually far better; their data-dependent
bound certifies how close greedy got on THIS instance rather than in the worst
case.

Three schedules compared at equal sounding budget on one geometry:

  uniform     every link sounded budget/L times (the paper's Assumption 4)
  greedy-link per-link sounding, each sounding chosen by lazy greedy
  greedy-tx   broadcast sounding, each TRANSMITTER chosen by lazy greedy
              (one transmission yields rows for all N-1 links of that node)

Objective: I(m) = (1/2) logdet(I + gamma * sum_s b_s b_s^T), with link rows
b represented in the right-singular basis of W, so every update is in
rank(W)-dimensional space and the posterior covariance is a small dense matrix.

Certificate: for cardinality-constrained monotone submodular maximisation,
OPT <= I(S) + sum of the k largest current marginal gains when k soundings
remain; evaluated at the end of the run it upper-bounds the optimal schedule's
value, so greedy/UB is a per-instance optimality certificate.

Memory (Rule 7.5): everything lives in r x r with r <= L; at N = 40 that is
780^2 doubles = 4.9 MB. No dense L x V matrix is formed (chunked Gram).
"""
from __future__ import annotations

import argparse
import gc
import heapq
import json
import os
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from ndf_multipath_check import perimeter, kernel  # noqa: E402  (same room)

RESULTS_DIR = Path(os.environ.get(
    "NDF_RESULTS_DIR", _REPO_ROOT / "docs/4-research/results"))
RESULTS_JSON = "ndf_greedy_schedule.json"

ROOM_W_M, ROOM_H_M = 5.0, 4.0
GRID_RES_M = 0.035
GAMMA = 1.0                  # per-sample SNR on the dominant mode
N_NODES = (16, 40)
BUDGETS = (100, 1000)        # total soundings R*T_c
VOXEL_CHUNK = 4000
RANK_TOL = 1e-9


def link_basis(pos):
    """Rows of W in the right-singular basis: B with B B^T = Gram(W)."""
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
            blk[r] = kernel(pos[i], pos[j], cx, cy)
        G += blk @ blk.T
        del blk
        gc.collect()
    ev, U = np.linalg.eigh(G)
    keep = ev > RANK_TOL * ev.max()
    B = U[:, keep] * np.sqrt(ev[keep])          # L x r, rows are b_l
    # normalise so gamma is SNR on the dominant mode, matching C(N)
    B /= np.sqrt(ev.max())
    return B, pairs


def info_uniform(B, budget):
    L = B.shape[0]
    n_per = budget / L
    s2 = np.linalg.eigvalsh(B.T @ B)
    return 0.5 * float(np.sum(np.log1p(GAMMA * n_per * np.clip(s2, 0, None))))


def lazy_greedy(cands, budget, r):
    """Lazy greedy over repeatable candidates; cands[i] is a (rows x r) block.

    Returns (info, picks, evals, certificate_upper_bound).
    """
    Sigma = np.eye(r)
    info = 0.0
    heap = []
    for i, Bl in enumerate(cands):
        M = Bl @ Bl.T
        g = 0.5 * np.linalg.slogdet(np.eye(len(Bl)) + GAMMA * M)[1]
        heapq.heappush(heap, (-g, 0, i))
    picks = np.zeros(len(cands), int)
    evals = len(cands)
    for step in range(budget):
        while True:
            negg, stamp, i = heapq.heappop(heap)
            if stamp == step:
                break
            Bl = cands[i]
            S = Bl @ Sigma @ Bl.T
            g = 0.5 * np.linalg.slogdet(np.eye(len(Bl)) + GAMMA * S)[1]
            evals += 1
            heapq.heappush(heap, (-g, step, i))
        # take candidate i
        Bl = cands[i]
        S = Bl @ Sigma @ Bl.T
        gain = 0.5 * np.linalg.slogdet(np.eye(len(Bl)) + GAMMA * S)[1]
        info += gain
        # posterior downdate: Sigma -= Sigma B^T (I/gamma + B Sigma B^T)^-1 B Sigma
        SB = Sigma @ Bl.T
        K = np.linalg.solve(np.eye(len(Bl)) / GAMMA + S, SB.T)
        Sigma = Sigma - SB @ K
        picks[i] += 1
        heapq.heappush(heap, (-gain, step + 1, i))   # stale; will be refreshed
    # data-dependent certificate: OPT(budget) <= info + sum of budget largest
    # CURRENT marginal gains (fresh evaluation for every candidate)
    gains = []
    for Bl in cands:
        S = Bl @ Sigma @ Bl.T
        gains.append(0.5 * np.linalg.slogdet(np.eye(len(Bl)) + GAMMA * S)[1])
    gains = sorted(gains, reverse=True)
    k = min(budget, len(gains))
    ub = info + float(np.sum(gains[:k] * np.ones(k))) if k else info
    # repeatable candidates: the same top candidate could be taken repeatedly,
    # so the sound upper bound uses budget * max gain
    ub = info + budget * (gains[0] if gains else 0.0)
    return info, picks, evals, ub


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=RESULTS_DIR / RESULTS_JSON)
    args = ap.parse_args()
    out = {"gamma": GAMMA, "budgets": list(BUDGETS), "results": {}}

    for n in N_NODES:
        pos = perimeter(n)
        B, pairs = link_basis(pos)
        L, r = B.shape
        print(f"N={n}: L={L}, rank r={r}", flush=True)
        link_cands = [B[i:i + 1] for i in range(L)]
        tx_blocks = []
        for node in range(n):
            rows = [k for k, (i, j) in enumerate(pairs) if i == node or j == node]
            tx_blocks.append(B[rows])
        for budget in BUDGETS:
            u = info_uniform(B, budget)
            gl, pl, el, ubl = lazy_greedy(link_cands, budget, r)
            # broadcast: budget in transmissions
            gt, pt, et, ubt = lazy_greedy(tx_blocks, budget, r)
            rec = {"uniform": u, "greedy_link": gl,
                   "greedy_link_cert": gl / ubl,
                   "greedy_tx": gt, "greedy_tx_cert": gt / ubt,
                   "links_used": int((pl > 0).sum()),
                   "tx_used": int((pt > 0).sum()),
                   "lazy_evals_link": el}
            out["results"][f"N{n}_B{budget}"] = rec
            print(f"  budget={budget:>5}: uniform {u:8.2f} | greedy-link "
                  f"{gl:8.2f} (cert>={gl/ubl:.2f}, {int((pl>0).sum())}/{L} links)"
                  f" | greedy-tx {gt:8.2f} (cert>={gt/ubt:.2f})", flush=True)
        del B, link_cands, tx_blocks
        gc.collect()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
