"""wifi-ndf: Spatial Degrees of Freedom for WiFi Mesh Sensing.

Compute how many independent spatial measurements a WiFi mesh provides,
from node coordinates and room geometry alone. No CSI data needed.

Usage:
    from wifi_ndf import compute_ndf

    result = compute_ndf(
        nodes=[(0.5, 0.5), (4.5, 0.5), (0.5, 3.5), (4.5, 3.5)],
        room=(5.0, 4.0),
    )
    print(f"NDF = {result.ndf}, efficiency = {result.efficiency:.2f}")
"""

from wifi_ndf.ndf import compute_ndf, NDFResult

__version__ = "0.1.0"
__all__ = ["compute_ndf", "NDFResult"]
