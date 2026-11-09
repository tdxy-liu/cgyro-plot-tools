"""Shared FULLT palette and display lattice, independent of the GUI/readers."""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm, TwoSlopeNorm
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from matplotlib.tri import Triangulation
from matplotlib.transforms import Affine2D, ScaledTranslation


FULLT_COLOR_SCALES = {"Weak signals (asinh)": "asinh", "Linear": "linear", "Symmetric log": "symlog"}
ASINH_WIDTH_FRACTION = .25

FULLT_DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "fullt_blue_white_red",
    [(0.00, "#053061"), (0.22, "#2166ac"), (0.40, "#92c5de"),
     (0.48, "#edf5fa"), (0.50, "#ffffff"), (0.52, "#fceee8"),
     (0.60, "#f4a582"), (0.78, "#b2182b"), (1.00, "#67001f")],
    N=2048,
)

# Fraction of each sample-to-cell-edge distance retaining the sample's color.
# This is shared display geometry, not an amplitude adjustment or a data kernel.
SMOOTH_CORE_FRACTION = .2


class FulltModeMarker(PathPatch):
    """A black diamond sized to the current mode's cell at render time."""

    def __init__(self, kx, ky, x_edges, y_edges):
        ix = np.searchsorted(x_edges[1:-1], kx)
        iy = np.searchsorted(y_edges[1:-1], ky)
        self._cell_points = [(x_edges[ix], y_edges[iy]), (kx, ky),
                             (x_edges[ix + 1], y_edges[iy + 1])]
        self._diamond_transform = Affine2D()
        super().__init__(Path.unit_regular_polygon(4), transform=self._diamond_transform,
                         facecolor="black", edgecolor="white", zorder=6)

    def draw(self, renderer):
        lower, center, upper = self.axes.transData.transform(self._cell_points)
        # Fit the diamond to the nearest cell edge in screen space. A single
        # radius keeps it square even when the two axes are zoomed differently.
        radius = .75 * np.minimum(np.abs(center - lower), np.abs(upper - center)).min()
        self._diamond_transform.clear().scale(radius).translate(*center)
        self.set_linewidth(.12 * radius / renderer.points_to_pixels(1))
        super().draw(renderer)


def apply_fullt_grid(ax):
    """Apply the light major-tick guides used in stacked trace panels."""
    ax._cgyro_custom_grid = True
    ax.set_axisbelow(False)
    ax.grid(True, which="major", color=".5", linewidth=.30, alpha=.10)
    ax.axhline(0, color=".45", linewidth=.35, alpha=.18, zorder=2.7)
    ax.axvline(0, color=".45", linewidth=.35, alpha=.18, zorder=2.7)


def fullt_color_norm(limit, scale):
    """Map both signs identically, with the colorbar in the displayed T units."""
    if scale == "asinh":
        width = max(limit * ASINH_WIDTH_FRACTION, np.finfo(float).tiny)
        try:
            from matplotlib.colors import AsinhNorm
        except ImportError:  # Matplotlib before 3.6 on some compute nodes.
            from matplotlib.colors import FuncNorm
            return FuncNorm((lambda t: np.arcsinh(t / width), lambda u: width * np.sinh(u)),
                            vmin=-limit, vmax=limit)
        return AsinhNorm(linear_width=width, vmin=-limit, vmax=limit)
    if scale == "linear":
        return TwoSlopeNorm(vmin=-limit, vcenter=0., vmax=limit)
    if scale == "symlog":
        return SymLogNorm(linthresh=max(limit * .02, np.finfo(float).tiny),
                          linscale=1., vmin=-limit, vmax=limit, base=10)
    raise ValueError(f"Unknown FULLT color scale: {scale}")


def fullt_color_scale(values, *, normalize=False, color_max=None, logarithmic=False, scale="linear"):
    """Use one absolute maximum for both signs and a symmetric color range."""
    values = np.ma.asarray(values, dtype=float).filled(np.nan)
    finite = values[np.isfinite(values)]
    maximum = float(np.max(np.abs(finite))) if finite.size else 0.
    divisor = maximum if normalize and maximum > 0 else 1.
    limit = maximum / divisor if maximum > 0 else 1.
    if color_max is not None:
        try:
            limit = float(color_max)
        except (TypeError, ValueError):
            raise ValueError("FULLT color max must be finite and positive.") from None
        if not np.isfinite(limit) or limit <= 0:
            raise ValueError("FULLT color max must be finite and positive.")
    return divisor, fullt_color_norm(limit, "symlog" if logarithmic else scale), maximum


def format_fullt_colorbar(colorbar, font_size=8.):
    """Share ticks and place the scientific multiplier above the right edge."""
    colorbar.locator = MaxNLocator(nbins=6)
    if hasattr(colorbar.norm, "linear_width"):
        from matplotlib.ticker import AsinhLocator
        colorbar.locator = AsinhLocator(colorbar.norm.linear_width, numticks=7, base=0)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 3))
    colorbar.formatter = formatter
    colorbar.update_ticks()
    colorbar.ax.tick_params(labelsize=font_size, length=2.5)
    colorbar.ax.yaxis.set_offset_position("right")
    offset = colorbar.ax.yaxis.get_offset_text()
    offset.set_fontsize(font_size)
    offset.set_ha("left")
    offset.set_transform(offset.get_transform() +
                         ScaledTranslation(3 / 72, 0, colorbar.ax.figure.dpi_scale_trans))


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
