"""Core NDF computation for WiFi mesh sensing.

Computes the Number of spatial Degrees of Freedom (NDF) of a WiFi mesh
by constructing a Fresnel-weighted sensing matrix and analyzing its
singular value spectrum.

Reference:
    K. Khamaisi and B. Rodrigues, "How Many Nodes Do We Need? A Spatial
    Coverage Metric for Wireless Sensing," IEEE GLOBECOM, 2026.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, Union

import numpy as np
import scipy.sparse
from scipy.sparse.linalg import svds

logger = logging.getLogger(__name__)

# =========================================================================
# Constants
# =========================================================================

SPEED_OF_LIGHT_M_S = 299_792_458.0

# Default wavelengths
WAVELENGTH_2_4GHZ_M = SPEED_OF_LIGHT_M_S / 2.4e9  # ~0.125 m
WAVELENGTH_5GHZ_M = SPEED_OF_LIGHT_M_S / 5.0e9     # ~0.060 m

# Grid discretization
DEFAULT_GRID_RESOLUTION_M = 0.25

# Fresnel weight model
FRESNEL_SIGMA_FRACTION = 0.5
FRESNEL_TRUNCATION_SIGMA = 3.0
ENDPOINT_WEIGHT_SCALE = 4.0

# NDF threshold
DEFAULT_NDF_THRESHOLD = 0.01

# Bounding box
ROOM_MARGIN_FRACTION = 0.1
DEGENERATE_DIMENSION_M = 1.0
MIN_LINK_LENGTH_M = 0.1

# SVD strategy
DENSE_SVD_THRESHOLD = 200


# =========================================================================
# Data structures
# =========================================================================

@dataclass
class NDFResult:
    """Result of NDF analysis.

    Attributes:
        ndf: Number of spatial degrees of freedom.
        num_nodes: Number of nodes in the mesh.
        num_links: Number of valid bistatic links.
        singular_values: Singular value spectrum (descending).
        room_size: (width, height) of the sensing region in meters.
        freq_ghz: Carrier frequency in GHz.
        threshold: Singular value threshold used.
    """
    ndf: int
    num_nodes: int
    num_links: int
    singular_values: np.ndarray
    room_size: Tuple[float, float]
    freq_ghz: float
    threshold: float

    @property
    def efficiency(self) -> float:
        """NDF/L: fraction of links carrying independent information.
        Values near 1.0 mean every link is useful; below 0.8 indicates
        wasted links from poor placement."""
        if self.num_links == 0:
            return 0.0
        return self.ndf / self.num_links

    @property
    def ka(self) -> float:
        """Electrical size of the room (2*pi*a/lambda)."""
        a = max(self.room_size) / 2
        wavelength = SPEED_OF_LIGHT_M_S / (self.freq_ghz * 1e9)
        return 2 * np.pi * a / wavelength

    def summary(self) -> str:
        """One-line summary string."""
        return (
            f"NDF={self.ndf} from {self.num_links} links "
            f"({self.num_nodes} nodes), "
            f"efficiency={self.efficiency:.2f}, "
            f"room={self.room_size[0]:.1f}x{self.room_size[1]:.1f}m "
            f"@ {self.freq_ghz} GHz"
        )

    def diagnose(self) -> str:
        """Placement quality diagnosis."""
        eff = self.efficiency
        if eff >= 0.9:
            return "Excellent: nearly every link adds independent information."
        elif eff >= 0.8:
            return "Good: minor redundancy, acceptable for most applications."
        elif eff >= 0.6:
            return "Fair: consider spreading nodes across more walls."
        else:
            return "Poor: most links are redundant. Reposition nodes."


# =========================================================================
# Public API
# =========================================================================

def compute_ndf(
    nodes: Union[np.ndarray, Sequence[Tuple[float, float]]],
    room: Optional[Tuple[float, float]] = None,
    freq_ghz: float = 2.4,
    threshold: float = DEFAULT_NDF_THRESHOLD,
    grid_resolution_m: float = DEFAULT_GRID_RESOLUTION_M,
) -> NDFResult:
    """Compute the NDF of a WiFi mesh from node positions.

    Args:
        nodes: Node positions as (N, 2) array or list of (x, y) tuples.
            Coordinates in meters.
        room: (width, height) of the room in meters. If None, auto-computed
            from node positions with 10% margin.
        freq_ghz: Carrier frequency in GHz (default 2.4).
        threshold: Relative singular value threshold (default 0.01 = 1%).
        grid_resolution_m: Grid cell size in meters (default 0.25).

    Returns:
        NDFResult with NDF, efficiency, and diagnostic info.

    Example:
        >>> from wifi_ndf import compute_ndf
        >>> result = compute_ndf(
        ...     nodes=[(0.5, 0.5), (4.5, 0.5), (0.5, 3.5), (4.5, 3.5)],
        ...     room=(5.0, 4.0),
        ... )
        >>> print(result.summary())
        >>> print(result.diagnose())
    """
    positions = np.asarray(nodes, dtype=np.float64)
    if positions.ndim == 1:
        positions = positions.reshape(-1, 2)
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError(f"nodes must be (N, 2), got shape {positions.shape}")
    if positions.shape[0] < 2:
        raise ValueError(f"Need at least 2 nodes, got {positions.shape[0]}")

    wavelength = SPEED_OF_LIGHT_M_S / (freq_ghz * 1e9)

    # Compute room bounds
    if room is not None:
        w, h = room
        cx = (positions[:, 0].min() + positions[:, 0].max()) / 2
        cy = (positions[:, 1].min() + positions[:, 1].max()) / 2
        bounds = ((cx - w / 2, cy - h / 2), (cx + w / 2, cy + h / 2))
        room_size = (w, h)
    else:
        bounds = _compute_room_bounds(positions)
        (x0, y0), (x1, y1) = bounds
        room_size = (x1 - x0, y1 - y0)

    # Build sensing matrix
    W, _, _ = _build_sensing_matrix(
        positions, grid_resolution_m, bounds, wavelength
    )

    # Compute NDF via SVD
    ndf, sigma = _compute_ndf_svd(W, threshold)

    return NDFResult(
        ndf=ndf,
        num_nodes=positions.shape[0],
        num_links=W.shape[0],
        singular_values=sigma,
        room_size=room_size,
        freq_ghz=freq_ghz,
        threshold=threshold,
    )


# =========================================================================
# Internal functions
# =========================================================================

def _compute_room_bounds(positions):
    x_min, y_min = positions.min(axis=0)
    x_max, y_max = positions.max(axis=0)
    width = max(x_max - x_min, DEGENERATE_DIMENSION_M)
    height = max(y_max - y_min, DEGENERATE_DIMENSION_M)
    if x_max - x_min < 1e-6:
        cx = (x_min + x_max) / 2
        x_min, x_max = cx - width / 2, cx + width / 2
    if y_max - y_min < 1e-6:
        cy = (y_min + y_max) / 2
        y_min, y_max = cy - height / 2, cy + height / 2
    mx = (x_max - x_min) * ROOM_MARGIN_FRACTION
    my = (y_max - y_min) * ROOM_MARGIN_FRACTION
    return ((x_min - mx, y_min - my), (x_max + mx, y_max + my))


def _build_sensing_matrix(positions, grid_res, bounds, wavelength):
    (x_min, y_min), (x_max, y_max) = bounds
    grid_x = np.arange(x_min + grid_res / 2, x_max, grid_res)
    grid_y = np.arange(y_min + grid_res / 2, y_max, grid_res)
    if len(grid_x) == 0 or len(grid_y) == 0:
        raise ValueError("Grid is empty — check room bounds and resolution.")
    gx, gy = np.meshgrid(grid_x, grid_y, indexing="ij")
    cells = np.column_stack([gx.ravel(), gy.ravel()])
    n_cells = len(cells)
    n_nodes = positions.shape[0]

    rows, cols, data = [], [], []
    link_idx = 0
    min_weight = np.exp(-0.5 * FRESNEL_TRUNCATION_SIGMA ** 2)

    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            tx, rx = positions[i], positions[j]
            vec = rx - tx
            length = np.linalg.norm(vec)
            if length < MIN_LINK_LENGTH_M:
                link_idx += 1
                continue
            direction = vec / length
            r_f = np.sqrt(wavelength * length / 4.0)
            sigma_f = r_f * FRESNEL_SIGMA_FRACTION
            margin = FRESNEL_TRUNCATION_SIGMA * sigma_f

            # Bounding box filter
            mask = (
                (cells[:, 0] >= min(tx[0], rx[0]) - margin) &
                (cells[:, 0] <= max(tx[0], rx[0]) + margin) &
                (cells[:, 1] >= min(tx[1], rx[1]) - margin) &
                (cells[:, 1] <= max(tx[1], rx[1]) + margin)
            )
            idx = np.nonzero(mask)[0]
            if len(idx) == 0:
                link_idx += 1
                continue

            cands = cells[idx]
            t = (cands - tx) @ direction / length
            valid = (t >= 0) & (t <= 1)
            if not valid.any():
                link_idx += 1
                continue

            vi, vc, vt = idx[valid], cands[valid], t[valid]
            proj = tx + np.outer(vt * length, direction)
            d_perp = np.linalg.norm(vc - proj, axis=1)
            w = np.exp(-d_perp**2 / (2 * sigma_f**2))
            w *= ENDPOINT_WEIGHT_SCALE * vt * (1 - vt)
            sig = w >= min_weight
            if sig.any():
                rows.extend([link_idx] * sig.sum())
                cols.extend(vi[sig].tolist())
                data.extend(w[sig].tolist())
            link_idx += 1

    n_links = link_idx
    if data:
        W = scipy.sparse.csr_matrix((data, (rows, cols)), shape=(n_links, n_cells))
    else:
        W = scipy.sparse.csr_matrix((n_links, n_cells))
    return W, grid_x, grid_y


def _compute_ndf_svd(W, threshold):
    n_links, n_cells = W.shape
    if n_links == 0 or n_cells == 0:
        return 0, np.array([])
    min_dim = min(n_links, n_cells)
    if min_dim < 2:
        return 0, np.array([])

    if min_dim <= DENSE_SVD_THRESHOLD:
        W_d = W.toarray() if scipy.sparse.issparse(W) else np.asarray(W)
        _, sigma, _ = np.linalg.svd(W_d, full_matrices=False)
    else:
        _, sigma, _ = svds(W.astype(np.float64), k=min_dim - 1)
        sigma = sigma[::-1]

    if sigma[0] <= 0:
        return 0, sigma
    ndf = int(np.sum(sigma > threshold * sigma[0]))
    return ndf, sigma
