"""Publication layout and reproducible export for the standalone FULLT tracer."""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional, Tuple
import csv
import re
import textwrap

import matplotlib as mpl
from matplotlib.collections import PolyCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
import numpy as np

try:
    from cgyro_io import write_json_atomic
    from cgyro_fullt_style import (FULLT_DIVERGING_CMAP, SMOOTH_CORE_FRACTION,
                                  smooth_diamond_mesh, fullt_color_scale, format_fullt_colorbar)
    from cgyro_fonts import plot_font_style, configure_plot_fonts
except ImportError:
    from .cgyro_io import write_json_atomic
    from .cgyro_fullt_style import (FULLT_DIVERGING_CMAP, SMOOTH_CORE_FRACTION,
                                   smooth_diamond_mesh, fullt_color_scale, format_fullt_colorbar)
    from .cgyro_fonts import plot_font_style, configure_plot_fonts


PUBLICATION_STYLE = {
    **plot_font_style(), "font.size": 9, "axes.linewidth": .65,
}
PREVIEW_DPI = 110  # Interactive preview only; exported rasters use settings.dpi (>=300).


@dataclass(frozen=True)
class FigureSettings:
    width_mm: float = 183.
    panel_height_mm: float = 32.
    font_pt: float = 9.
    color_max: Optional[float] = None
    normalize: bool = False
    glyph: str = "smooth"
    xlim: Optional[Tuple[float, float]] = None
    ylim: Optional[Tuple[float, float]] = None
    title: Optional[str] = None
    colorbar_label: str = ""
    dpi: int = 600
    panel_scaling: str = "auto"
    panel_multipliers: Optional[Tuple[float, ...]] = None

    def validate(self):
        for name, low, high in (("width_mm", 89, 300), ("panel_height_mm", 22, 100),
                                 ("font_pt", 8.5, 16), ("dpi", 300, 1200)):
            value = getattr(self, name)
            if not np.isfinite(value) or not low <= value <= high:
                raise ValueError(f"{name} must be between {low} and {high}.")
        if self.color_max is not None and (not np.isfinite(self.color_max) or self.color_max <= 0):
            raise ValueError("The shared color maximum must be finite and positive.")
        if self.glyph not in ("smooth", "diamonds", "cells"):
            raise ValueError("Choose smooth, diamonds or cells.")
        if self.panel_scaling not in ("auto", "none", "custom"):
            raise ValueError("Panel scaling must be auto, none or custom.")
        if self.panel_multipliers is not None:
            factors = np.asarray(self.panel_multipliers, dtype=float)
            if (factors.ndim != 1 or not 1 <= factors.size <= 20 or
                    not np.all(np.isfinite(factors)) or np.any(factors <= 0)):
                raise ValueError("Panel multipliers must contain 1-20 finite positive numbers.")
        if self.panel_scaling == "custom" and self.panel_multipliers is None:
            raise ValueError("Enter custom multipliers in panel order, for example: 1, 1, 10.")
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


def _panel_maxima(result):
    return [float(np.max(np.abs(p.values[np.isfinite(p.values)])))
            if np.any(np.isfinite(p.values)) else 0. for p in result.panels]


def panel_factors(result, settings):
    """Display-only multipliers; automatic 1/2/5 decades never clip a panel."""
    count = len(result.panels)
    if settings.panel_scaling == "none":
        return [1.] * count
    if settings.panel_scaling == "custom":
        if len(settings.panel_multipliers) < count:
            raise ValueError(f"Enter at least {count} multipliers, one for each displayed panel.")
        # A trace can stop before the requested number of steps.
        return [float(x) for x in settings.panel_multipliers[:count]]
    maxima = _panel_maxima(result)
    maximum = max(maxima, default=0.)
    divisor = maximum if settings.normalize and maximum > 0 else 1.
    target = settings.color_max if settings.color_max is not None else maximum / divisor
    factors = []
    for peak in maxima:
        ratio = target / (peak / divisor) if peak > 0 else 1.
        if not np.isfinite(ratio):
            raise ValueError("Automatic panel scaling overflowed; choose custom multipliers.")
        if ratio < 2:
            factors.append(1.)
            continue
        decade = 10. ** np.floor(np.log10(ratio))
        if decade > ratio:  # log10 can round a value just below a decade upwards.
            decade /= 10.
        factors.append(float(max(x * decade for x in (1., 2., 5.) if x <= ratio / decade)))
    return factors


def shared_scale(result, settings):
    """Keep the shared range independent of display multipliers.

    Otherwise enlarging a weak panel also enlarges the colorbar and washes out
    the unchanged panels. Explicit color_max still overrides the raw reference.
    """
    maxima = _panel_maxima(result)
    divisor, norm, maximum = fullt_color_scale(maxima, normalize=settings.normalize,
                                              color_max=settings.color_max)
    displayed = [peak / divisor * factor for peak, factor in zip(maxima, panel_factors(result, settings))]
    if not np.all(np.isfinite(displayed)):
        raise ValueError("Panel multipliers overflow the displayed values; use smaller factors.")
    return divisor, norm.vmax, maximum


def build_figure(result, settings=None):
    """Stack one step per row; show only each panel's current mode in its header."""
    configure_plot_fonts()
    settings = settings or FigureSettings()
    settings.validate()
    if not result.panels:
        raise ValueError("There are no trace panels to draw.")
    divisor, limit, maximum = shared_scale(result, settings)
    factors = panel_factors(result, settings)
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0., vmax=limit)
    n = len(result.panels)
    # Reserve physical space for panel headers and the common bottom x label.
    narrow = settings.width_mm / (settings.font_pt / 9.) < 130
    header_mm = (11. if narrow else 6.) * settings.font_pt / 9.
    title = settings.title
    if title is None:
        title = f"{result.case}  |  t = {result.times[0]:g}–{result.times[-1]:g}  |  {result.settings.direction} T"
    title = "\n".join(textwrap.wrap(title, width=max(24, int((settings.width_mm - 20) * 5.2 / settings.font_pt))))
    top_mm = max(4., 4.5 * settings.font_pt / 9. * len(title.splitlines()) + 4.)
    bottom_mm = 12. * settings.font_pt / 9.
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
        cmap = FULLT_DIVERGING_CMAP
        for i, panel in enumerate(result.panels):
            row = outer[i, 0].subgridspec(2, 1, height_ratios=[header_mm, settings.panel_height_mm], hspace=0.)
            header = fig.add_subplot(row[0])
            header.set_axis_off()
            ax = fig.add_subplot(row[1], sharex=axes[0] if axes else None, sharey=axes[0] if axes else None)
            axes.append(ax)
            values = np.ma.masked_invalid(panel.values / divisor * factors[i])
            if settings.glyph == "smooth":
                mesh, nodes = smooth_diamond_mesh(result.kx, result.ky, xedges, yedges, values)
                artist = ax.tripcolor(mesh, nodes, cmap=cmap, norm=norm, shading="gouraud",
                                      edgecolors="none", linewidth=0., antialiased=True, rasterized=True)
            elif settings.glyph == "diamonds":
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
            ax.axhline(0, color=".45", alpha=.18, linewidth=.35, zorder=2.7)
            ax.axvline(0, color=".45", alpha=.18, linewidth=.35, zorder=2.7)
            ax.grid(color=".5", alpha=.10, linewidth=.30)
            ax.set_axisbelow(False)
            ax.plot(panel.mode.kx, panel.mode.ky, "D", mfc="black", mec="white",
                    mew=.7, ms=6.5, zorder=6)
            ax.set_xlim(settings.xlim or (xedges[0], xedges[-1]))
            ax.set_ylim(settings.ylim or (yedges[0], yedges[-1]))
            ax.set_ylabel(r"$k_y'\rho_s$", fontsize=settings.font_pt)
            ax.tick_params(labelsize=settings.font_pt - 1, length=2.5, width=.6, pad=2)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
            ax.tick_params(labelbottom=i == n - 1)
            if i == n - 1:
                ax.set_xlabel(r"$k_x'\rho_s$", fontsize=settings.font_pt, labelpad=2)
            quantity = r"T(k,k')/\max|T|" if settings.normalize else r"T(k,k')"
            multiplier = rf"{factors[i]:g}\times " if factors[i] != 1. else ""
            label = f"({chr(97 + i)}) $" + multiplier + quantity + "$"
            header.text(0, .77 if narrow else .5, label, fontsize=settings.font_pt, va="center")
            current = rf"$k_x\rho_s={panel.mode.kx:.4g},\ k_y\rho_s={panel.mode.ky:.4g}$"
            if narrow:
                header.text(0, .25, current, ha="left", va="center", fontsize=settings.font_pt)
            else:
                header.text(1, .5, current, ha="right", va="center", fontsize=settings.font_pt)
        cax = fig.add_subplot(outer[:, 1])
        position = cax.get_position()
        cax.set_position([position.x0, axes[-1].get_position().y0, position.width,
                          axes[0].get_position().y1 - axes[-1].get_position().y0])
        cbar = fig.colorbar(artist, cax=cax, extend="both", extendfrac=.025)
        cbar.set_label(settings.colorbar_label, fontsize=settings.font_pt)
        format_fullt_colorbar(cbar, settings.font_pt - 1)
        fig.text(left_mm / settings.width_mm, 1 - 2 / height_mm, title,
                 ha="left", va="top", fontsize=settings.font_pt)
        fig._fullt_trace_axes = axes
        fig._fullt_trace_scale = {"divisor": divisor, "limit": limit, "panel_multipliers": factors}
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
    factors = panel_factors(result, settings)
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
        writer.writerow(["step", "kx_prime", "ky_prime", "T_mean", "display_value", "display_multiplier"])
        for i, panel in enumerate(result.panels, 1):
            for iy, ky in enumerate(result.ky):
                writer.writerows((i, float(kx), float(ky), float(panel.values[iy, ix]),
                                  float(panel.values[iy, ix] / divisor * factors[i - 1]), factors[i - 1])
                                 for ix, kx in enumerate(result.kx))
    np.savez_compressed(target / "maps.npz", kx=result.kx, ky=result.ky, times=result.times,
                        **{f"step_{i:02d}": panel.values for i, panel in enumerate(result.panels, 1)})
    metadata = result.metadata()
    metadata["figure"] = asdict(settings)
    metadata["shared_color_scale"] = {"divisor": divisor, "limit": limit, "raw_max_abs": maximum,
                                      "reference": "explicit_color_max" if settings.color_max is not None else "unscaled_global_max_abs",
                                      "saturated_points": [int(np.count_nonzero(
                                          np.isfinite(p.values) & (np.abs(p.values / divisor * factor) > limit)))
                                          for p, factor in zip(result.panels, factors)]}
    metadata["display"] = {"panel_multipliers": factors, "formula": "display_value = T_mean / divisor * panel_multiplier",
                           "palette": FULLT_DIVERGING_CMAP.name, "rendering": settings.glyph,
                           "smooth_lattice": "sample-centred constant-color diamond cores with Gouraud fade to neutral cell edges; display only",
                           "smooth_core_fraction": SMOOTH_CORE_FRACTION if settings.glyph == "smooth" else None}
    metadata["formats"] = list(formats)
    metadata["version"] = (Path(__file__).resolve().parent / "VERSION").read_text().strip()
    write_json_atomic(target / "trace.json", metadata)
    (target / "README.txt").write_text(
        "Each panel is a successive extremum of the same time-averaged native Re(T) diagnostic.\n"
        "All panels share one zero-centred linear color scale for the DISPLAYED values.\n"
        "The automatic color range is fixed by the UNMULTIPLIED maps, not recomputed after panel amplification.\n"
        "Values beyond that range use the endpoint colors; exact values remain in plot_data.csv and maps.npz.\n"
        "display_value = T_mean / global_divisor * panel_multiplier. Multipliers appear before T in each panel header.\n"
        f"Panel multipliers: {', '.join(format(x, 'g') for x in factors)}; global divisor: {divisor:g}.\n"
        "Smooth rendering uses the same full-color diamond core and short Gouraud edge fade in every panel.\n"
        f"The core extends {SMOOTH_CORE_FRACTION:g} of the sample-to-cell-edge distance; the remaining rim fades to neutral.\n"
        "Those extra nodes are display geometry, not physical data. No weak-signal threshold or percentile clipping is applied.\n"
        "Diamonds and cells show discrete samples; a filled black diamond marks the current diagnosed mode.\n"
        "The initial matched ky stays fixed. Extrema are selected only along ky_prime=ky; only kx changes.\n"
        "Other ky_prime rows remain visible but do not participate in tracing.\n"
        "trace.csv records the path; plot_data.csv and maps.npz retain every grid sample.\n"
        "Nonfinite samples invalidate that grid point for the entire selected average.\n"
        "trace.json records time samples, options, color scaling, input identity and stop reason.\n"
        "The trace is an extremum-following diagnostic, not a temporal trajectory or proof of a unique energy cascade.\n"
        + result.stop_reason + "\n", encoding="utf-8")
    return target


def suggested_bundle_name(result):
    name = re.sub(r"[^\w.-]+", "_", result.case).strip("._") or "case"
    return f"{name}_{result.settings.direction}_trace"
