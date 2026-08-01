# wifi-ndf

**How large should a WiFi sensing mesh be, and where should its nodes sit?**

`wifi-ndf` computes the Number of spatial Degrees of Freedom (NDF) of a WiFi
mesh: the count of independent measurements its pairwise links provide about a
room, from node coordinates, room dimensions, and wavelength alone. No CSI is
needed. It is the reference implementation for the paper *Optimal Mesh Size for
WiFi Sensing Under a Sounding Budget* (under review).

## Install

```bash
pip install wifi-ndf
```

## Quick start

```python
from wifi_ndf import compute_ndf

result = compute_ndf(
    nodes=[(0.5, 0.5), (4.5, 0.5), (0.5, 3.5), (4.5, 3.5),
           (2.5, 0.5), (2.5, 3.5), (0.5, 2.0), (4.5, 2.0)],
    room=(5.0, 4.0),
    freq_ghz=2.4,
)
print(result.summary())
# NDF=28 from 28 links (8 nodes), efficiency=1.00, room=5.0x4.0m @ 2.4 GHz
```

Eight well-spread nodes achieve `NDF = L` *exactly*: every link carries
independent information. This equality holds up to a placement-dependent
boundary (12 nodes well-spread, 4 along a single wall, for the room above) and
was confirmed on three mapped hardware deployments across all 904 node subsets.

## How to read the number

- `0 <= NDF <= L`, larger is better; `NDF/L` is the fraction of links carrying
  unique information.
- Collinearity, not clustering, destroys efficiency: 14 nodes in a 2x2 m corner
  keep 96% of their links, the same 14 along one wall keep 37%.
- Geometry caps NDF absolutely: a wall-confined mesh in the reference room
  saturates near `NDF = 312` at ~40 nodes; interior meshes do not saturate at
  any tested count. Ceiling law: `NDF_inf ~ 1.19*ka + 14.21*sqrt(ka)`
  (fitted on three electrical sizes, predicting two held-out ones at -5.3%
  and -2.6%).
- NDF prices no measurement by acquisition cost. The optimal *node count* is a
  property of the sounding budget `B = SNR * rate * coherence time`, not of the
  room: it moves from 8 to 46 nodes across one geometry as `B` sweeps six
  decades. See the paper and `analysis/ndf_rate_optimum.py`.

## v0.2.0: corrected grid criterion

Earlier releases used a fixed 0.25 m grid, which under-resolves the Fresnel
zones of short links and **under-counts NDF for near-collinear and clustered
layouts** (a 14-node corner cluster reported 0.64 efficiency; the converged
value is 0.96). The default is now automatic: half the Gaussian Fresnel width
of the shortest link, `sigma_F(d_min)/2`. Pass `grid_resolution_m` explicitly
to override.

## API

```python
compute_ndf(
    nodes,               # list of (x, y) tuples or (N, 2) array [m]
    room=(5.0, 4.0),     # (width, height) in meters, or None for auto
    freq_ghz=2.4,        # carrier frequency
    threshold=0.01,      # relative singular-value threshold
    grid_resolution_m=None,  # None = physical criterion sigma_F(d_min)/2
) -> NDFResult
```

**NDFResult**: `ndf`, `num_nodes`, `num_links`, `efficiency`, `ka`,
`singular_values`, `summary()`, `diagnose()`.

## Analysis archive

`analysis/` contains the scripts and committed results behind the paper's
claims: the capacity optimum against the sounding budget, greedy-vs-round-robin
scheduling with data-dependent certificates, the first-order multipath
robustness check, grid-convergence and threshold ablations, and the mapped
hardware deployments. `docs/theory/` contains the per-node rank theorem with
its 240-configuration numerical verification
(`analysis/results/gaussian_count_lemma_verification.json`). Some scripts
import geometry helpers from the authors' testbed repository and are archived
here for inspection alongside their committed outputs.

## Citation

```bibtex
@misc{khamaisi2026ndf,
  title={Optimal Mesh Size for WiFi Sensing Under a Sounding Budget},
  author={Khamaisi, Karim and Rodrigues, Bruno},
  year={2026},
  note={Under review}
}
```

## License

MIT. See [LICENSE](LICENSE).
