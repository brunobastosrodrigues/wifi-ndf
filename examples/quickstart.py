"""Quick start: compute NDF for a 14-node mesh in a 5x4m room."""

from wifi_ndf import compute_ndf

# 14 nodes spread across three walls of a 5x4m room
nodes = [
    (0.3, 0.3), (1.5, 0.3), (2.5, 0.3), (3.7, 0.3), (4.7, 0.3),  # bottom
    (0.3, 2.0), (4.7, 2.0),                                        # sides
    (0.3, 3.7), (1.2, 3.7), (2.0, 3.7), (2.8, 3.7), (3.5, 3.7),  # top
    (4.7, 3.7), (4.7, 1.0),                                        # right
]

result = compute_ndf(nodes, room=(5.0, 4.0), freq_ghz=2.4)

print(result.summary())
print(result.diagnose())
print(f"\nScaling model prediction: C*sqrt(L*ka) = 0.92*sqrt({result.num_links}*{result.ka:.0f}) = {0.92 * (result.num_links * result.ka)**0.5:.0f}")


# Compare: same N=14, but all nodes in one corner
corner_nodes = [(0.5 + i*0.3, 0.5 + j*0.3) for i in range(4) for j in range(4)][:14]
corner = compute_ndf(corner_nodes, room=(5.0, 4.0), freq_ghz=2.4)

print(f"\n--- Corner cluster ---")
print(corner.summary())
print(corner.diagnose())
print(f"Link count is the same ({corner.num_links}), but NDF dropped from {result.ndf} to {corner.ndf}")
