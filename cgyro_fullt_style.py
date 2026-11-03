"""Shared FULLT palette and display lattice, independent of the GUI/readers."""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap


FULLT_DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "fullt_blue_white_red_suppressed_weak",
    [(0.00, "#053061"), (0.18, "#2166ac"), (0.34, "#92c5de"),
     (0.46, "#f7fbff"), (0.50, "#ffffff"), (0.54, "#fff5f0"),
     (0.66, "#f4a582"), (0.82, "#b2182b"), (1.00, "#67001f")],
    N=512,
)


def diamond_lattice(x_axis, y_axis, x_edges, y_edges, values):
    """Fade saved samples to neutral cell edges for the classic FULLT appearance.

    The inserted zero nodes are display geometry, not extra physical samples.
    The caller must retain the original grid for extrema, statistics and export.
    """
    x_axis = np.asarray(x_axis, dtype=float).reshape(-1)
    y_axis = np.asarray(y_axis, dtype=float).reshape(-1)
    x_edges = np.asarray(x_edges, dtype=float).reshape(-1)
    y_edges = np.asarray(y_edges, dtype=float).reshape(-1)
    values = np.ma.asarray(values, dtype=float).filled(np.nan)
    if values.shape != (y_axis.size, x_axis.size):
        return None, None, None
    if x_edges.size != x_axis.size + 1 or y_edges.size != y_axis.size + 1:
        return None, None, None

    x_nodes = np.empty(2 * x_axis.size + 1, dtype=float)
    y_nodes = np.empty(2 * y_axis.size + 1, dtype=float)
    x_nodes[0::2], x_nodes[1::2] = x_edges, x_axis
    y_nodes[0::2], y_nodes[1::2] = y_edges, y_axis
    z_nodes = np.zeros((y_nodes.size, x_nodes.size), dtype=float)
    z_mask = np.zeros_like(z_nodes, dtype=bool)
    finite = np.isfinite(values)
    z_nodes[1::2, 1::2] = np.where(finite, values, 0.)
    z_mask[1::2, 1::2] = ~finite
    x_grid, y_grid = np.meshgrid(x_nodes, y_nodes)
    return x_grid, y_grid, np.ma.array(z_nodes, mask=z_mask)
