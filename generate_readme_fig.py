"""Generate README demo figure showing NDF for good vs bad placement."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from wifi_ndf import compute_ndf

C_BLUE = '#0072B2'
C_RED = '#D55E00'
C_GREEN = '#009E73'

fig, axes = plt.subplots(1, 3, figsize=(9, 3))

# --- Panel 1: Well-spread ---
nodes_good = [
    (0.3, 0.3), (1.5, 0.3), (2.5, 0.3), (3.7, 0.3), (4.7, 0.3),
    (0.3, 2.0), (4.7, 2.0),
    (0.3, 3.7), (1.2, 3.7), (2.0, 3.7), (2.8, 3.7), (3.5, 3.7),
    (4.7, 3.7), (4.7, 1.0),
]
r_good = compute_ndf(nodes_good, room=(5.0, 4.0))

ax = axes[0]
ax.set_xlim(-0.2, 5.2)
ax.set_ylim(-0.2, 4.2)
ax.set_aspect('equal')
ax.add_patch(plt.Rectangle((0, 0), 5, 4, fill=False, edgecolor='black', lw=1.5))
xs, ys = zip(*nodes_good)
ax.scatter(xs, ys, c=C_GREEN, s=60, zorder=5, edgecolor='black', lw=0.5)
ax.set_title(f'Well-spread\nNDF={r_good.ndf}, NDF/L={r_good.efficiency:.2f}', fontsize=10, fontweight='bold', color=C_GREEN)
ax.set_xlabel('5 x 4 m room', fontsize=9)
ax.tick_params(labelsize=7)

# --- Panel 2: Corner cluster ---
nodes_bad = [(0.5 + i*0.3, 0.5 + j*0.3) for i in range(4) for j in range(4)][:14]
r_bad = compute_ndf(nodes_bad, room=(5.0, 4.0))

ax = axes[1]
ax.set_xlim(-0.2, 5.2)
ax.set_ylim(-0.2, 4.2)
ax.set_aspect('equal')
ax.add_patch(plt.Rectangle((0, 0), 5, 4, fill=False, edgecolor='black', lw=1.5))
xs, ys = zip(*nodes_bad)
ax.scatter(xs, ys, c=C_RED, s=60, zorder=5, edgecolor='black', lw=0.5)
ax.set_title(f'Corner cluster\nNDF={r_bad.ndf}, NDF/L={r_bad.efficiency:.2f}', fontsize=10, fontweight='bold', color=C_RED)
ax.set_xlabel('Same 14 nodes, same L=91', fontsize=9)
ax.tick_params(labelsize=7)

# --- Panel 3: NDF vs N scaling ---
ax = axes[2]
ns = range(4, 30)
ndfs = []
for n in ns:
    np.random.seed(42)
    nodes = [(np.random.uniform(0.3, 4.7), np.random.uniform(0.3, 3.7)) for _ in range(n)]
    r = compute_ndf(nodes, room=(5.0, 4.0))
    ndfs.append(r.ndf)

ax.plot(list(ns), ndfs, '-o', color=C_BLUE, markersize=4, lw=1.5)
ax.axvline(x=15, color='gray', linestyle='--', lw=0.8, alpha=0.7)
ax.text(15.5, max(ndfs)*0.3, 'N*~15\ndiminishing\nreturns', fontsize=7, color='gray')
ax.set_xlabel('Number of nodes N', fontsize=9)
ax.set_ylabel('NDF', fontsize=9)
ax.set_title('Scaling (5x4m, 2.4 GHz)', fontsize=10, fontweight='bold')
ax.tick_params(labelsize=7)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('docs/demo.png', dpi=150, bbox_inches='tight')
plt.savefig('docs/demo.svg', bbox_inches='tight')
print("Saved docs/demo.png and docs/demo.svg")
