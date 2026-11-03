"""Publication layout and reproducible export for the standalone FULLT tracer."""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional, Tuple
import csv
import re
import textwrap

import matplotlib as mpl
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import numpy as np

try:
    from cgyro_io import write_json_atomic
except ImportError:
    from .cgyro_io import write_json_atomic


PUBLICATION_STYLE = {
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9, "mathtext.fontset": "dejavusans", "pdf.fonttype": 42,
    "ps.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .65,
    "axes.unicode_minus": True, "text.usetex": False,
}
PREVIEW_DPI = 110  # Interactive preview only; exported rasters use settings.dpi (>=300).


@dataclass(frozen=True)
class FigureSettings:
    width_mm: float = 183.
    panel_height_mm: float = 32.
    font_pt: float = 9.
    color_max: Optional[float] = None
    normalize: bool = False
    glyph: str = "diamonds"
    xlim: Optional[Tuple[float, float]] = None
    ylim: Optional[Tuple[float, float]] = None
    title: Optional[str] = None
    colorbar_label: str = ""
    dpi: int = 600

    def validate(self):
        for name, low, high in (("width_mm", 89, 300), ("panel_height_mm", 22, 100),
                                 ("font_pt", 8.5, 16), ("dpi", 300, 1200)):
            value = getattr(self, name)
            if not np.isfinite(value) or not low <= value <= high:
                raise ValueError(f"{name} must be between {low} and {high}.")
        if self.color_max is not None and (not np.isfinite(self.color_max) or self.color_max <= 0):
            raise ValueError("The shared color maximum must be finite and positive.")
        if self.glyph not in ("diamonds", "cells"):
            raise ValueError("Choose diamonds or cells.")
        for limits in (self.xlim, self.ylim):
            if limits is not None and (len(limits) != 2 or not np.all(np.isfinite(limits)) or limits[0] >= limits[1]):
                raise ValueError("Axis limits must be a finite increasing pair.")


def _edges(axis):
    axis = np.asarray(axis, dtype=float)
    if axis.size == 1:
        return np.array([axis[0] - .05, axis[0] + .05])
    if np.any(np.diff(axis) <= 0):
        raise ValueError("Plot coordinates must be strictly increasing.")
    return np.r_[axis[0] - .5 * (axis[1] - axis[0]),
                 .5 * (axis[1:] + axis[:-1]), axis[-1] + .5 * (axis[-1] - axis[-2])]


def shared_scale(result, settings):
    """One linear, zero-centred scale for every panel, without percentile clipping."""
    maxima = [float(np.max(np.abs(p.values[np.isfinite(p.values)])))
              for p in result.panels if np.any(np.isfinite(p.values))]
    maximum = max(maxima, default=0.)
    divisor = maximum if settings.normalize and maximum > 0 else 1.
    limit = settings.color_max if settings.color_max is not None else maximum / divisor
    return divisor, (limit if limit > 0 else 1.), maximum


def build_figure(result, settings=None):
    """Stack one step per row; whitespace headers never obscure measured points."""
    settings = settings or FigureSettings()
    settings.validate()
    if not result.panels:
        raise ValueError("There are no trace panels to draw.")
    divisor, limit, maximum = shared_scale(result, settings)
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0., vmax=limit)
    n = len(result.panels)
    # Reserve a physical header strip for current/next kx and a bottom legend.
    narrow = settings.width_mm / (settings.font_pt / 9.) < 130
    header_mm = (17. if narrow else 10.5) * settings.font_pt / 9.
    title = settings.title
    if title is None:
        title = f"{result.case}  |  t = {result.times[0]:g}–{result.times[-1]:g}  |  {result.settings.direction} T"
    title = "\n".join(textwrap.wrap(title, width=max(24, int((settings.width_mm - 20) * 5.2 / settings.font_pt))))
    top_mm = max(4., 4.5 * settings.font_pt / 9. * len(title.splitlines()) + 4.)
    bottom_mm = 20. * settings.font_pt / 9.
    row_mm, gap_mm = settings.panel_height_mm + header_mm, 3.
    height_mm = n * row_mm + (n - 1) * gap_mm + bottom_mm + top_mm
    with mpl.rc_context(PUBLICATION_STYLE):
        fig = Figure(figsize=(settings.width_mm / 25.4, height_mm / 25.4), dpi=PREVIEW_DPI)
        left_mm = 20. * settings.font_pt / 9.
        right_mm = (23. if settings.colorbar_label else 17.) * settings.font_pt / 9.
        outer = fig.add_gridspec(n, 2, width_ratios=[1, .026], hspace=gap_mm / row_mm, wspace=.07,
                                left=left_mm / settings.width_mm,
                                right=1 - right_mm / settings.width_mm,
                                bottom=bottom_mm / height_mm, top=1 - top_mm / height_mm)
        axes = []
        xedges, yedges = _edges(result.kx), _edges(result.ky)
        colors = mpl.colormaps["RdBu_r"](np.linspace(0, 1, 257))
        colors[128] = (1., 1., 1., 1.)
        cmap = LinearSegmentedColormap.from_list("RdBu_white_zero", colors, N=257)
        for i, panel in enumerate(result.panels):
            row = outer[i, 0].subgridspec(2, 1, height_ratios=[header_mm, settings.panel_height_mm], hspace=0.)
            header = fig.add_subplot(row[0])
            header.set_axis_off()
            ax = fig.add_subplot(row[1], sharex=axes[0] if axes else None, sharey=axes[0] if axes else None)
            axes.append(ax)
            values = np.ma.masked_invalid(panel.values / divisor)
            if settings.glyph == "diamonds":
                # Each diamond is a mark for ONE saved grid sample. No smoothing,
                # synthetic zero nodes, filtering or amplitude thresholding.
                xx, yy = np.meshgrid(result.kx, result.ky)
                wx, wy = np.meshgrid(np.diff(xedges) * .48, np.diff(yedges) * .48)
                vertices = np.stack((np.stack((xx - wx, yy), -1), np.stack((xx, yy + wy), -1),
                                     np.stack((xx + wx, yy), -1), np.stack((xx, yy - wy), -1)), axis=-2)
                artist = PolyCollection(vertices.reshape(-1, 4, 2), array=values.ravel(),
                                        cmap=cmap, norm=norm, edgecolors="none", rasterized=True)
                ax.add_collection(artist)
            else:
                artist = ax.pcolormesh(xedges, yedges, values, cmap=cmap, norm=norm,
                                      shading="flat", rasterized=True)
            ax.axhline(0, color=".65", linewidth=.45, zorder=2.7)
            ax.axvline(0, color=".65", linewidth=.45, zorder=2.7)
            ax.grid(color=".9", linewidth=.35)
            ax.set_axisbelow(False)
            ax.plot(panel.mode.kx, panel.mode.ky, "o", mfc="none", mec=".1", mew=1., ms=4., zorder=5)
            if panel.peak:
                ax.plot(panel.peak.kx, panel.peak.ky, "D", mfc="none", mec="black", mew=.9, ms=5., zorder=6)
            ax.set_xlim(settings.xlim or (xedges[0], xedges[-1]))
            ax.set_ylim(settings.ylim or (yedges[0], yedges[-1]))
            ax.set_ylabel(r"$k_y'\rho_s$", fontsize=settings.font_pt)
            ax.tick_params(labelsize=settings.font_pt - 1, length=2.5, width=.6, pad=2)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
            ax.tick_params(labelbottom=i == n - 1)
            if i == n - 1:
                ax.set_xlabel(r"$k_x'\rho_s$", fontsize=settings.font_pt, labelpad=2)
            quantity = r"$T^\Phi(k,k')/\max|T|$" if settings.normalize else r"$T^\Phi(k,k')$"
            label = f"({chr(97 + i)}) " + quantity
            header.text(0, .84 if narrow else .72, label, fontsize=settings.font_pt, va="center")
            current = rf"$k_x\rho_s={panel.mode.kx:.4g},\ k_y\rho_s={panel.mode.ky:.4g}$"
            if panel.peak:
                value_label = "raw T" if settings.normalize else "T"
                following = f"next kx = {panel.peak.next_mode.kx:.4g}; {value_label} = {panel.peak.value:.3g}"
                if panel.peak.conjugated:
                    following += " (*)"
            else:
                following = "No eligible extremum"
            if narrow:
                header.text(0, .5, current, ha="left", va="center", fontsize=settings.font_pt)
                header.text(0, .16, following, ha="left", va="center", fontsize=settings.font_pt - 1)
            else:
                header.text(1, .72, current, ha="right", va="center", fontsize=settings.font_pt)
                header.text(1, .23, following, ha="right", va="center", fontsize=settings.font_pt - 1)
        cax = fig.add_subplot(outer[:, 1])
        position = cax.get_position()
        cax.set_position([position.x0, axes[-1].get_position().y0, position.width,
                          axes[0].get_position().y1 - axes[-1].get_position().y0])
        cbar = fig.colorbar(artist, cax=cax, extend="both" if maximum / divisor > limit else "neither")
        cbar.set_label(settings.colorbar_label, fontsize=settings.font_pt)
        cbar.locator = MaxNLocator(nbins=6)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-3, 3))
        cbar.formatter = formatter
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=settings.font_pt - 1, length=2.5)
        cbar.ax.yaxis.get_offset_text().set_fontsize(settings.font_pt - 1)
        fig.text(left_mm / settings.width_mm, 1 - 2 / height_mm, title,
                 ha="left", va="top", fontsize=settings.font_pt)
        handles = [Line2D([], [], marker="o", color="none", mec=".1", mfc="none", ms=4, label="Current mode"),
                   Line2D([], [], marker="D", color="none", mec=".1", mfc="none", ms=4, label="Selected extremum")]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, 0), ncol=2,
                   frameon=False, fontsize=settings.font_pt - 1, handlelength=1., columnspacing=1.3)
        fig._fullt_trace_axes = axes
        fig._fullt_trace_scale = {"divisor": divisor, "limit": limit}
        return fig


def export_bundle(result, settings, destination, formats=("pdf", "svg", "png")):
    """Create a new bundle directory; existing figures/data are never overwritten."""
    settings.validate()
    if not formats or any(fmt not in ("pdf", "svg", "png", "tiff") for fmt in formats):
        raise ValueError("Export formats must be pdf, svg, png or tiff.")
    destination = Path(destination).expanduser().absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    target = destination
    serial = 1
    while True:
        try:
            target.mkdir()
            break
        except FileExistsError:
            serial += 1
            target = destination.with_name(destination.name + f"_{serial}")
    fig = build_figure(result, settings)
    with mpl.rc_context(PUBLICATION_STYLE):
        if "pdf" in formats:
            fig.savefig(target / "trace.pdf", dpi=settings.dpi)
        if "svg" in formats:
            fig.savefig(target / "trace.svg", dpi=settings.dpi)
        if "png" in formats:
            fig.savefig(target / "trace.png", dpi=settings.dpi)
        if "tiff" in formats:
            fig.savefig(target / "trace.tiff", dpi=settings.dpi)
    divisor, limit, maximum = shared_scale(result, settings)
    with (target / "trace.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "kx", "ky", "selected_kx_prime", "selected_ky_prime", "T_mean",
                         "next_kx", "next_ky", "conjugated", "ties", "excluded_visited"])
        for i, panel in enumerate(result.panels, 1):
            p = panel.peak
            writer.writerow([i, panel.mode.kx, panel.mode.ky] +
                            ([p.kx, p.ky, p.value, p.next_mode.kx, p.next_mode.ky, p.conjugated, p.ties]
                             if p else [""] * 7) + [panel.excluded_visited])
    with (target / "plot_data.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "kx_prime", "ky_prime", "T_mean", "display_value"])
        for i, panel in enumerate(result.panels, 1):
            for iy, ky in enumerate(result.ky):
                writer.writerows((i, float(kx), float(ky), float(panel.values[iy, ix]), float(panel.values[iy, ix] / divisor))
                                 for ix, kx in enumerate(result.kx))
    np.savez_compressed(target / "maps.npz", kx=result.kx, ky=result.ky, times=result.times,
                        **{f"step_{i:02d}": panel.values for i, panel in enumerate(result.panels, 1)})
    metadata = result.metadata()
    metadata["figure"] = asdict(settings)
    metadata["shared_color_scale"] = {"divisor": divisor, "limit": limit, "raw_max_abs": maximum}
    metadata["formats"] = list(formats)
    metadata["version"] = (Path(__file__).resolve().parent / "VERSION").read_text().strip()
    write_json_atomic(target / "trace.json", metadata)
    (target / "README.txt").write_text(
        "Each panel is a successive extremum of the same time-averaged native Re(T) diagnostic.\n"
        "All panels share one zero-centred linear color scale. No percentile clipping, smoothing or per-panel scaling.\n"
        "Diamonds, when selected, are discrete sample marks, not interpolated lobes.\n"
        "(*) means the selected negative-ky partner is represented by (-kx,-ky) at the next step.\n"
        "trace.csv records the path; plot_data.csv and maps.npz retain every grid sample.\n"
        "Nonfinite samples invalidate that grid point for the entire selected average.\n"
        "trace.json records time samples, options, color scaling, input identity and stop reason.\n"
        "The trace is an extremum-following diagnostic, not a temporal trajectory or proof of a unique energy cascade.\n"
        + result.stop_reason + "\n", encoding="utf-8")
    return target


def suggested_bundle_name(result):
    name = re.sub(r"[^\w.-]+", "_", result.case).strip("._") or "case"
    return f"{name}_{result.settings.direction}_trace"
