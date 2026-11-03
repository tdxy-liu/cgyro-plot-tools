"""Shared FULLT palette and display lattice, independent of the GUI/readers."""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.tri import Triangulation


FULLT_DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "fullt_blue_white_red_suppressed_weak",
    [(0.00, "#053061"), (0.18, "#2166ac"), (0.34, "#92c5de"),
     (0.46, "#f7fbff"), (0.50, "#ffffff"), (0.54, "#fff5f0"),
     (0.66, "#f4a582"), (0.82, "#b2182b"), (1.00, "#67001f")],
    N=512,
)

# Fraction of each sample-to-cell-edge distance retaining the sample's color.
# This is shared display geometry, not an amplitude adjustment or a data kernel.
SMOOTH_CORE_FRACTION = .4


def smooth_diamond_mesh(x_axis, y_axis, x_edges, y_edges, values):
    """Return separate Gouraud glyphs with a full-color core and a shorter fade.

    Every saved sample has a central diamond of constant color, then a smooth
    transition to neutral at its cell edges. Explicit triangles keep neighboring
    samples separate; all triangles belonging to a missing sample are masked.
    Added vertices are display geometry only, never interpolated physical data.
    """
    x_axis = np.asarray(x_axis, dtype=float).reshape(-1)
    y_axis = np.asarray(y_axis, dtype=float).reshape(-1)
    x_edges = np.asarray(x_edges, dtype=float).reshape(-1)
    y_edges = np.asarray(y_edges, dtype=float).reshape(-1)
    values = np.ma.asarray(values, dtype=float).filled(np.nan)
    if (not x_axis.size or not y_axis.size or
            values.shape != (y_axis.size, x_axis.size) or
            x_edges.size != x_axis.size + 1 or y_edges.size != y_axis.size + 1):
        raise ValueError("Smooth FULLT samples and cell edges must match the coordinate grid.")

    xx, yy = np.meshgrid(x_axis, y_axis)
    centres = np.column_stack((xx.ravel(), yy.ravel()))
    offsets = np.zeros((values.size, 4, 2), dtype=float)
    # Counterclockwise corners, using each actual cell edge on nonuniform grids.
    offsets[:, 0, 0] = np.tile(x_edges[1:] - x_axis, y_axis.size)
    offsets[:, 1, 1] = np.repeat(y_edges[1:] - y_axis, x_axis.size)
    offsets[:, 2, 0] = np.tile(x_edges[:-1] - x_axis, y_axis.size)
    offsets[:, 3, 1] = np.repeat(y_edges[:-1] - y_axis, x_axis.size)
    vertices = np.concatenate((centres[:, None, :],
                               centres[:, None, :] + SMOOTH_CORE_FRACTION * offsets,
                               centres[:, None, :] + offsets), axis=1)
    finite = np.isfinite(values.ravel())
    colors = np.zeros((values.size, 9), dtype=float)
    colors[:, :5] = np.where(finite, values.ravel(), 0.)[:, None]
    cell_triangles = []
    for corner in range(4):
        following = (corner + 1) % 4
        cell_triangles.extend(((0, 1 + corner, 1 + following),
                               (1 + corner, 5 + corner, 5 + following),
                               (1 + corner, 5 + following, 1 + following)))
    triangles = (np.asarray(cell_triangles)[None, :, :] +
                 9 * np.arange(values.size)[:, None, None]).reshape(-1, 3)
    mesh = Triangulation(vertices[:, :, 0].ravel(), vertices[:, :, 1].ravel(),
                        triangles=triangles, mask=np.repeat(~finite, len(cell_triangles)))
    return mesh, colors.ravel()


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
