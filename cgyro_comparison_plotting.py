"""
Computation and plotting mixin for CGYRO comparison GUI.
"""

import os
import re
import traceback
import difflib
from tkinter import messagebox

import matplotlib.animation as animation
import matplotlib as mpl
import numpy as np

DEFAULT_POD_Z_WINDOW_PI = 8.0
LINE_COLOR_PALETTE = [
    "#F14040",
    "#1A6FDF",
    "#37AD6B",
    "#B177DE",
    "#CC9900",
    "#00CBCC",
    "#7D4E4E",
    "#8E8E00",
    "#FB6501",
    "#6699CC",
    "#6FB802",
]
GAMMA_LIN_COLOR = "#515151"
CONTOUR_COLOR_PALETTE = [
    "#08306B",
    "#08519C",
    "#2171B5",
    "#4292C6",
    "#6BAED6",
    "#9ECAE1",
    "#C6DBEF",
    "#DEEBF7",
    "#FEE0D2",
    "#FCBBA1",
    "#FC9272",
    "#FB6A4A",
    "#EF3B2C",
    "#CB181D",
    "#A50F15",
    "#67000D",
]
GRID_COLOR = "#808080"
GRID_LINEWIDTH = 0.5
AXIS_BORDER_LINEWIDTH = 1.5
PLOT_LINEWIDTH = 3.0
TICK_LABEL_FONTSIZE = 18
MAIN_FONT_FONTSIZE = 22
FONT_FAMILY = "Arial"

try:
    from cgyro_comparison_plotting_zf import ZfPlotting
except ImportError:
    from .cgyro_comparison_plotting_zf import ZfPlotting

try:
    from cgyro_comparison_plotting_frequency import FrequencyPlotting
    from cgyro_comparison_plotting_fft import FftPlotting
    from cgyro_comparison_plotting_fluctuation import FluctuationPlotting
    from cgyro_comparison_plotting_flux import FluxPlotting
    from cgyro_comparison_plotting_energy import EnergyPlotting
    from cgyro_comparison_plotting_other import OtherPlotting
except ImportError:
    from .cgyro_comparison_plotting_frequency import FrequencyPlotting
    from .cgyro_comparison_plotting_fft import FftPlotting
    from .cgyro_comparison_plotting_fluctuation import FluctuationPlotting
    from .cgyro_comparison_plotting_flux import FluxPlotting
    from .cgyro_comparison_plotting_energy import EnergyPlotting
    from .cgyro_comparison_plotting_other import OtherPlotting

try:
    import scipy.signal as sp_signal
    from scipy.optimize import curve_fit as sp_curve_fit
except Exception:
    sp_signal = None
    sp_curve_fit = None


class Plotting(FrequencyPlotting, FftPlotting, FluctuationPlotting, FluxPlotting, EnergyPlotting, OtherPlotting, ZfPlotting):
    @staticmethod
    def _get_line_color_palette():
        """Unified line-color palette for 1D/line plots."""
        return list(LINE_COLOR_PALETTE)

    @staticmethod
    def _get_gamma_lin_color():
        """Dedicated color for gamma_lin curves."""
        return str(GAMMA_LIN_COLOR)

    @staticmethod
    def _get_2d_contour_cmap():
        """Unified smooth colormap built from the requested 16-color gradient list."""
        try:
            return mpl.colors.LinearSegmentedColormap.from_list(
                "cgyro_cmp_contour", CONTOUR_COLOR_PALETTE, N=256
            )
        except Exception:
            return "viridis"

    def _apply_global_plot_color_style(self):
        """Apply shared line color cycle for current plotting session."""
        try:
            self.fig.set_prop_cycle(color=self._get_line_color_palette())
        except Exception:
            pass
        try:
            mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=self._get_line_color_palette())
        except Exception:
            pass
        style_updates = {
            'font.family': FONT_FAMILY,
            'font.size': MAIN_FONT_FONTSIZE,
            'axes.titlesize': MAIN_FONT_FONTSIZE,
            'axes.labelsize': MAIN_FONT_FONTSIZE,
            'legend.fontsize': MAIN_FONT_FONTSIZE,
            'xtick.labelsize': TICK_LABEL_FONTSIZE,
            'ytick.labelsize': TICK_LABEL_FONTSIZE,
            'axes.linewidth': AXIS_BORDER_LINEWIDTH,
            'grid.color': GRID_COLOR,
            'grid.linewidth': GRID_LINEWIDTH,
            'grid.linestyle': '-',
            'axes.axisbelow': True,
            'lines.linewidth': PLOT_LINEWIDTH,
        }
        for k, v in style_updates.items():
            try:
                mpl.rcParams[k] = v
            except Exception:
                pass

    def _apply_unified_visual_style_to_figure(self):
        """Enforce final plotting style on every axes/colorbar in the current figure."""
        try:
            axes_all = list(self.fig.axes)
        except Exception:
            axes_all = []
        for ax in axes_all:
            try:
                ax.set_axisbelow(True)
            except Exception:
                pass
            try:
                ax.grid(True, color=GRID_COLOR, linewidth=GRID_LINEWIDTH)
            except Exception:
                pass
            for side in ['left', 'right', 'top', 'bottom']:
                try:
                    ax.spines[side].set_linewidth(AXIS_BORDER_LINEWIDTH)
                except Exception:
                    pass
            try:
                ax.tick_params(axis='both', which='both', labelsize=TICK_LABEL_FONTSIZE)
            except Exception:
                pass
            for tick in list(getattr(ax, 'get_xticklabels', lambda: [])()) + list(getattr(ax, 'get_yticklabels', lambda: [])()):
                try:
                    tick.set_fontname(FONT_FAMILY)
                except Exception:
                    pass
            try:
                ax.xaxis.label.set_fontname(FONT_FAMILY)
                ax.yaxis.label.set_fontname(FONT_FAMILY)
                ax.xaxis.label.set_fontsize(MAIN_FONT_FONTSIZE)
                ax.yaxis.label.set_fontsize(MAIN_FONT_FONTSIZE)
            except Exception:
                pass
            try:
                ax.title.set_fontname(FONT_FAMILY)
                ax.title.set_fontsize(MAIN_FONT_FONTSIZE)
            except Exception:
                pass
            try:
                leg = ax.get_legend()
                if leg is not None:
                    for txt in leg.get_texts():
                        txt.set_fontname(FONT_FAMILY)
                        txt.set_fontsize(MAIN_FONT_FONTSIZE)
                    ttl = leg.get_title()
                    if ttl is not None:
                        ttl.set_fontname(FONT_FAMILY)
                        ttl.set_fontsize(MAIN_FONT_FONTSIZE)
            except Exception:
                pass
            try:
                for ln in ax.get_lines():
                    ln.set_linewidth(PLOT_LINEWIDTH)
            except Exception:
                pass

        # Also style figure-level texts (e.g., suptitle).
        try:
            txts = list(getattr(self.fig, 'texts', []))
        except Exception:
            txts = []
        for txt in txts:
            try:
                txt.set_fontname(FONT_FAMILY)
                txt.set_fontsize(MAIN_FONT_FONTSIZE)
            except Exception:
                pass

    def _average_mode_name(self):
        """Return current amplitude-averaging mode name from UI."""
        try:
            mode = str(self.fluc_average_var.get()).strip().lower()
        except Exception:
            mode = "root mean square"
        if mode == "mean absolute":
            return "Mean Absolute"
        return "Root Mean Square"

    def _use_mean_absolute_average(self):
        """True when using Mean Absolute averaging mode."""
        return self._average_mode_name() == "Mean Absolute"

    def _average_mode_short_tag(self):
        """Compact tag for legends/labels."""
        return "MA" if self._use_mean_absolute_average() else "RMS"

    def _reset_figure_layout_defaults(self):
        """
        Restore figure-level subplot/layout state to Matplotlib defaults.

        Important: `fig.clear()` does NOT reset subplot parameters modified by
        previous `tight_layout()` calls (e.g. POD parity 1x2 panel). Without this
        reset, subsequent single-axis plots may look stretched/deformed.
        """
        try:
            if hasattr(self.fig, 'set_layout_engine'):
                self.fig.set_layout_engine(None)
        except Exception:
            pass

        try:
            self.fig.subplots_adjust(
                left=float(mpl.rcParams.get('figure.subplot.left', 0.125)),
                right=float(mpl.rcParams.get('figure.subplot.right', 0.9)),
                bottom=float(mpl.rcParams.get('figure.subplot.bottom', 0.11)),
                top=float(mpl.rcParams.get('figure.subplot.top', 0.88)),
                wspace=float(mpl.rcParams.get('figure.subplot.wspace', 0.2)),
                hspace=float(mpl.rcParams.get('figure.subplot.hspace', 0.2)),
            )
        except Exception:
            try:
                self.fig.subplots_adjust(
                    left=0.125, right=0.9, bottom=0.11, top=0.88, wspace=0.2, hspace=0.2
                )
            except Exception:
                pass

    @staticmethod
    def _read_input_cgyro_scalars(case_dir):
        """
        Parse scalar assignments from input.cgyro into a dict.

        Supports common forms such as:
          KEY = value
          KEY value
        and ignores comments after `#`.
        """
        scalars = {}
        if not case_dir:
            return scalars
        path = os.path.join(case_dir, 'input.cgyro')
        if not os.path.isfile(path):
            return scalars

        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    s = line.split('#', 1)[0].strip()
                    if not s:
                        continue
                    m = re.match(
                        r'^([A-Za-z][A-Za-z0-9_]*)\s*(?:=\s*|\s+)([-+0-9.eEdD]+)\s*$',
                        s
                    )
                    if not m:
                        continue
                    key = m.group(1).upper()
                    val_txt = m.group(2).replace('D', 'E').replace('d', 'e')
                    try:
                        val = float(val_txt)
                    except Exception:
                        continue
                    scalars[key] = val
        except Exception:
            return {}
        return scalars

    def _use_fluc2d_x_rhoe(self):
        """Return True when Fluctuation-2D spatial axes should use electron-scale normalization."""
        try:
            return bool(self.fluc2d_x_elec_var.get())
        except Exception:
            return False

    @staticmethod
    def _positive_ky_axis(ky_like):
        """Return ky array for plotting with non-negative convention."""
        ky = np.asarray(ky_like, dtype=float).reshape(-1)
        if ky.size <= 0:
            return ky
        return np.abs(ky)

    def _use_flux_real_ion_norm(self):
        """Return True when flux should be renormalized by real-ion GyroBohm units."""
        try:
            return bool(self.flux_norm_real_ion_var.get())
        except Exception:
            return False

    def _get_flux_real_ion_norm_context(self, data, label=''):
        """
        Build cached real-ion normalization context from case species metadata.

        Real ion is selected as the ion species (Z>0) with the largest equilibrium
        density. The conversion follows D-normalized CGYRO units (m_D=1) and uses:
          c_s = sqrt(T_e/m_ref), rho_s ~ sqrt(m_ref)
        so for reference change D -> real-ion:
          cs_ri/cs = sqrt(1/m_i), rho_s_ri/rho_s = sqrt(m_i)
          gc = Gamma_GB,ri / Gamma_GB = (cs_ri/cs)*(rho_s_ri/rho_s)^2
          qc = Q_GB,ri / Q_GB = (cs_ri/cs)*(rho_s_ri/rho_s)^2
        so plotted flux values should be multiplied by 1/gc (particle) or 1/qc (energy).
        """
        cached = getattr(data, '_cmp_flux_real_ion_norm_ctx', None)
        if isinstance(cached, dict):
            return cached

        ctx = {
            'valid': False,
            'gc': 1.0,
            'qc': 1.0,
            'vc': 1.0,
            'rhoc': 1.0,
            'ion_index': None,
            'z': None,
            'mass': None,
            'temp': None,
            'dens': None,
        }

        try:
            case_dir = self._resolve_case_dir(data)
            scalars = self._read_input_cgyro_scalars(case_dir)

            z_attr = np.asarray(getattr(data, 'z', []), dtype=float).reshape(-1)
            mass_attr = np.asarray(getattr(data, 'mass', []), dtype=float).reshape(-1)
            temp_attr = np.asarray(getattr(data, 'temp', []), dtype=float).reshape(-1)
            dens_attr = np.asarray(getattr(data, 'dens', []), dtype=float).reshape(-1)

            n_from_data = int(getattr(data, 'n_species', 0))
            try:
                n_from_input = int(round(float(scalars.get('N_SPECIES', 0))))
            except Exception:
                n_from_input = 0

            n_species = max(
                n_from_data,
                n_from_input,
                int(z_attr.size),
                int(mass_attr.size),
                int(temp_attr.size),
                int(dens_attr.size),
            )
            if n_species <= 0:
                raise ValueError("no species metadata available")

            z = np.full(n_species, np.nan, dtype=float)
            mass = np.full(n_species, np.nan, dtype=float)
            temp = np.full(n_species, np.nan, dtype=float)
            dens = np.full(n_species, np.nan, dtype=float)

            for src, dst in (
                (z_attr, z),
                (mass_attr, mass),
                (temp_attr, temp),
                (dens_attr, dens),
            ):
                if src.size > 0:
                    ncopy = min(int(src.size), n_species)
                    dst[:ncopy] = src[:ncopy]

            for i in range(n_species):
                j = i + 1
                if not np.isfinite(z[i]):
                    z_try = scalars.get(f'Z_{j}', None)
                    if z_try is not None:
                        z[i] = float(z_try)
                if not np.isfinite(mass[i]):
                    m_try = scalars.get(f'MASS_{j}', None)
                    if m_try is not None:
                        mass[i] = float(m_try)
                if not np.isfinite(temp[i]):
                    t_try = scalars.get(f'TEMP_{j}', None)
                    if t_try is not None:
                        temp[i] = float(t_try)
                if not np.isfinite(dens[i]):
                    d_try = scalars.get(f'DENS_{j}', None)
                    if d_try is not None:
                        dens[i] = float(d_try)

            # Real ion: ion species with largest density (Z>0).
            ion_candidates = np.where(np.isfinite(z) & (z > 0.0))[0]
            if ion_candidates.size <= 0:
                raise ValueError("no ion species (Z>0) found")

            if np.any(np.isfinite(dens[ion_candidates])):
                idens = np.where(np.isfinite(dens[ion_candidates]), dens[ion_candidates], -np.inf)
                i_main = int(ion_candidates[int(np.argmax(idens))])
            else:
                i_main = int(ion_candidates[0])

            zi = float(z[i_main]) if np.isfinite(z[i_main]) else np.nan
            mi = float(mass[i_main]) if np.isfinite(mass[i_main]) else np.nan
            ti = float(temp[i_main]) if np.isfinite(temp[i_main]) else np.nan
            ni = float(dens[i_main]) if np.isfinite(dens[i_main]) else np.nan

            if (not np.isfinite(zi)) or zi <= 0.0:
                raise ValueError("invalid main-ion charge")
            if (not np.isfinite(mi)) or mi <= 0.0:
                raise ValueError("invalid main-ion mass")
            # Follow user's normalization formulas:
            #   c_s = sqrt(T_e/m_ref), rho_s ∝ sqrt(m_ref)
            # with m_ref:D=1 -> real-ion:m_i
            #   cs_ri/cs = sqrt(1/m_i), rho_s_ri/rho_s = sqrt(m_i)
            vc = float(np.sqrt(1.0 / mi))
            rhoc = float(np.sqrt(mi))
            gc = float(vc * (rhoc ** 2))
            qc = float(vc * (rhoc ** 2))

            if (not np.isfinite(gc)) or gc <= 0.0:
                raise ValueError("invalid gc normalization factor")
            if (not np.isfinite(qc)) or qc <= 0.0:
                raise ValueError("invalid qc normalization factor")

            ctx.update({
                'valid': True,
                'gc': gc,
                'qc': qc,
                'vc': vc,
                'rhoc': rhoc,
                'ion_index': i_main,
                'z': zi,
                'mass': mi,
                'temp': ti,
                'dens': ni,
            })
        except Exception as e:
            if label and (not getattr(data, '_cmp_flux_real_ion_norm_warned', False)):
                print(
                    f"Warning: failed to compute real-ion flux normalization for {label}; "
                    f"fallback to default D normalization ({e})."
                )
                setattr(data, '_cmp_flux_real_ion_norm_warned', True)

        setattr(data, '_cmp_flux_real_ion_norm_ctx', ctx)
        if (
            label
            and ctx.get('valid', False)
            and (not getattr(data, '_cmp_flux_real_ion_norm_logged', False))
        ):
            i_main = int(ctx.get('ion_index', 0)) + 1
            qgb_ratio = float(ctx.get('qc', 1.0))   # QGB_Ri / QGB_D
            rho_ratio = float(ctx.get('rhoc', 1.0)) # rho_ri / rho_s(D)
            cs_ratio = float(ctx.get('vc', 1.0))    # cs_ri / cs
            gc_ratio = float(ctx.get('gc', 1.0))    # GammaGB_Ri / GammaGB_D
            print(
                "Info: {} uses real-ion normalization from species #{} "
                "(Z={:.4g}, M={:.4g}, dens={:.4g}).".format(
                    label,
                    i_main,
                    float(ctx.get('z', 0.0)),
                    float(ctx.get('mass', 0.0)),
                    float(ctx.get('dens', np.nan)),
                )
            )
            print(
                "Info: normalization ratios -> "
                "QGB_Ri/QGB={:.6g}, rho_s_ri/rho_s={:.6g}, cs_ri/cs={:.6g}, "
                "GammaGB_Ri/GammaGB={:.6g}".format(
                    qgb_ratio, rho_ratio, cs_ratio, gc_ratio
                )
            )
            print(
                "Info: applied plot scales -> "
                "Flux(E): x{}, y*1/{:.6g}; Flux(P): x{}, y*1/{:.6g}".format(
                    " (ky,time)", qgb_ratio, " (ky,time)", gc_ratio
                )
            )
            setattr(data, '_cmp_flux_real_ion_norm_logged', True)
        return ctx

    def _get_flux_real_ion_norm_scale(self, data, moment_idx, label=''):
        """Return multiplicative scale from D-normalized flux to real-ion-normalized flux."""
        if not self._use_flux_real_ion_norm():
            return 1.0

        ctx = self._get_flux_real_ion_norm_context(data, label=label)
        if not ctx.get('valid', False):
            return 1.0

        denom_ratio = float(ctx.get('qc', 1.0)) if int(moment_idx) == 1 else float(ctx.get('gc', 1.0))
        if (not np.isfinite(denom_ratio)) or denom_ratio <= 0.0:
            return 1.0
        return float(1.0 / denom_ratio)

    def _get_rhos_to_rhoe_factor(self, data, label=''):
        """
        Compute rho_s/rho_e conversion factor using CGYRO's D-normalized convention.

        x/rho_e = (x/rho_s) * (rho_s/rho_e)
        y/rho_e = (y/rho_s) * (rho_s/rho_e)

        In this module we follow the user's convention:
        - Deuterium is the normalization reference, m_D = 1.
        - Therefore rho_s/rho_e = sqrt(m_D/m_e) = 1/sqrt(m_e_norm).
        - Electron mass priority:
            1) input.cgyro: MASS_{N_SPECIES}
            2) data.mass last species
            3) physical fallback
        - If electron mass is not found from case metadata, fallback to
          physical ratio m_e/m_D ~= 9.109e-28 / 3.345e-24.
        """
        cached = getattr(data, '_cmp_rhos_to_rhoe_factor', None)
        try:
            if cached is not None and np.isfinite(float(cached)) and float(cached) > 0.0:
                return float(cached)
        except Exception:
            pass

        factor = 1.0
        try:
            me_norm_default = float(9.109e-28 / 3.345e-24)
            me_norm = None

            # Priority 1: input.cgyro with MASS_{N_SPECIES}.
            case_dir = self._resolve_case_dir(data)
            scalars = self._read_input_cgyro_scalars(case_dir)
            n_species_raw = scalars.get('N_SPECIES', None)
            if n_species_raw is not None:
                try:
                    n_species = int(round(float(n_species_raw)))
                except Exception:
                    n_species = None
                if isinstance(n_species, int) and n_species > 0:
                    key = f'MASS_{n_species}'
                    me_try = scalars.get(key, None)
                    if me_try is not None:
                        me_try = float(me_try)
                        if np.isfinite(me_try) and me_try > 0.0:
                            me_norm = me_try

            # Priority 2: fallback to last species in data.mass.
            if me_norm is None:
                m_arr = np.asarray(getattr(data, 'mass', []), dtype=float).reshape(-1)
                if m_arr.size > 0:
                    me_try = float(np.abs(m_arr[-1]))
                    if np.isfinite(me_try) and me_try > 0.0:
                        me_norm = me_try

            # Priority 3: physical constant ratio m_e/m_D.
            if me_norm is None:
                me_norm = me_norm_default

            if (not np.isfinite(me_norm)) or me_norm <= 0.0:
                raise ValueError("invalid electron normalized mass")

            # D-normalized convention: m_D = 1.
            factor = float(np.sqrt(1.0 / me_norm))
            if (not np.isfinite(factor)) or factor <= 0.0:
                raise ValueError("non-finite rho_s/rho_e factor")
        except Exception as e:
            factor = 1.0
            if label:
                print(
                    f"Warning: failed to compute rho_s/rho_e for {label}, "
                    f"fallback to rho_s scale ({e})."
                )

        setattr(data, '_cmp_rhos_to_rhoe_factor', float(factor))
        if label and (factor > 1.0) and (not getattr(data, '_cmp_rhos_to_rhoe_logged', False)):
            print(f"Info: {label} uses rho_s/rho_e = {factor:.4g} for Fluctuation 2D spatial axes.")
            setattr(data, '_cmp_rhos_to_rhoe_logged', True)
        return float(factor)

    def _get_time_indices(self, t):
        """Resolve time-window indices from GUI range entries with safe fallbacks."""
        try:
            t_start_str = self.t_start_var.get().strip()
            t_end_str = self.t_end_var.get().strip()

            t_arr = np.asarray(t).reshape(-1)
            if t_arr.size == 0:
                return [], 0, 0

            t_min = float(np.min(t_arr))
            t_max = float(np.max(t_arr))
            default_start = t_min + 0.5 * (t_max - t_min)

            t_start = float(t_start_str) if t_start_str else default_start
            t_end = float(t_end_str) if t_end_str else t_max

            # Clamp user input to available range.
            t_start = max(t_min, min(t_start, t_max))
            t_end = max(t_min, min(t_end, t_max))

            # If invalid (reversed or collapsed), fallback to the last 50%.
            if t_start >= t_end:
                t_start = default_start
                t_end = t_max

            indices = np.where((t_arr >= t_start) & (t_arr <= t_end))[0]
            if len(indices) == 0:
                # Final fallback for sparse/irregular grids.
                start_idx = t_arr.size // 2
                indices = np.arange(start_idx, t_arr.size)
                t_start = float(t_arr[indices[0]])
                t_end = float(t_arr[indices[-1]])

            return indices, t_start, t_end
        except Exception as e:
            print(f"Error parsing time range: {e}")
            return [], 0, 0

    @staticmethod
    def _format_avg_suffix(t_start, t_end, prefix="Avg"):
        """Build a normalized averaging-range suffix like ` (Avg: t0-t1)`."""
        try:
            t0 = float(t_start)
            t1 = float(t_end)
            if (not np.isfinite(t0)) or (not np.isfinite(t1)):
                return ""
            if t1 < t0:
                t0, t1 = t1, t0
            return f" ({prefix}: {t0:.1f}-{t1:.1f})"
        except Exception:
            return ""

    def _append_avg_suffix(self, label, t_start, t_end, prefix="Avg"):
        """Append a formatted averaging-range suffix to one label."""
        return f"{label}{self._format_avg_suffix(t_start, t_end, prefix=prefix)}"

    @staticmethod
    def _format_avg_range_from_axis(x_axis, valid_idx, prefix="Avg"):
        """Build averaging suffix from axis values and selected index array."""
        try:
            x = np.asarray(x_axis, dtype=float).reshape(-1)
            idx = np.asarray(valid_idx, dtype=int).reshape(-1)
            if x.size <= 0 or idx.size <= 0:
                return ""
            idx = idx[(idx >= 0) & (idx < x.size)]
            if idx.size <= 0:
                return ""
            return f" ({prefix}: {float(x[idx[0]]):.1f}-{float(x[idx[-1]]):.1f})"
        except Exception:
            return ""

    def _get_case_linestyle(self, case_label, map_attr, line_styles):
        """Return stable per-case line style from one persistent style map."""
        style_map = getattr(self, map_attr, None)
        if not isinstance(style_map, dict):
            style_map = {}
            setattr(self, map_attr, style_map)
        if case_label not in style_map:
            style_map[case_label] = line_styles[len(style_map) % len(line_styles)]
        return style_map[case_label]

    @staticmethod
    def _maptoreal_fft(nr, nn, nx, ny, c):
        """Map CGYRO spectral coefficients to real-space field via 2D inverse FFT."""
        # Storage for numpy inverse real transform (irfft2)
        d = np.zeros([nx, nn], dtype=complex)

        for i in range(nr):
            p = i - nr // 2
            # k is the "standard FFT index"
            if -p < 0:
                k = -p + nx
            else:
                k = -p
            # Use identity f(p,-n) = f(-p,n)*
            d[k, 0:nn] = np.conj(c[i, 0:nn])

        # 2D inverse real Hermitian transform
        # NOTE: using inverse FFT with convention exp(ipx+iny), so need n -> -n
        # NOTE: need factor of 0.5 to match half-sum method of slow maptoreal()
        f = np.fft.irfft2(d, s=[nx, ny], norm='forward') * 0.5
        return f

    def _ensure_bigfield_loaded(self, data, label):
        """Try loading heavy field data required by fluctuation/FFT plots."""
        try:
            data.getbigfield()
            return True
        except Exception as e:
            print(f"Could not load big field data for {label}: {e}")
            return False

    @staticmethod
    def _reconstruct_x_from_kx(c_kx_t, nx=None):
        """
        Reconstruct real-space x profiles from complex kx coefficients.
        Input shape: [n_radial, n_time], output shape: [n_time, nx].
        """
        coeff = np.asarray(c_kx_t, dtype=complex)
        if coeff.ndim != 2 or coeff.shape[0] <= 0 or coeff.shape[1] <= 0:
            return None, None

        nr, _nt = coeff.shape
        if nx is None or int(nx) <= 0:
            nx = nr + 1
        nx = int(nx)

        x = np.arange(nx, dtype=float) * (2.0 * np.pi / max(nx, 1))

        # Match pygacode/cgyro/plot_rt.py mapping for ky=0 x-t contours.
        d = np.zeros((coeff.shape[1], nx), dtype=complex)
        for ix in range(-nr // 2 + 1, nr // 2):
            i = ix if ix >= 0 else ix + nx
            src = -ix + nr // 2
            if 0 <= i < nx and 0 <= src < nr:
                d[:, i] = coeff[src, :]

        # Shape: [nt, nx]
        f_tx = np.real(np.fft.fft(d, axis=1)) * 0.5

        # Fallback for unusual grids (e.g., odd nr) where FFT mapping can be sparse.
        if not np.any(np.isfinite(f_tx)) or np.nanmax(np.abs(f_tx)) == 0.0:
            p = np.arange(nr, dtype=float) - (nr // 2)
            phase = np.exp(1j * np.outer(p, x))  # [nr, nx]
            f_tx = np.real(np.matmul(coeff.T, phase))

        return x, f_tx

    def _plot_1d(self, x, y, label, plot_type):
        """Shared finite-checking wrapper for 1D line plotting."""
        x_arr = np.asarray(x).reshape(-1)
        y_arr = np.asarray(y).reshape(-1)
        if x_arr.size == 0 or y_arr.size == 0:
            return

        if x_arr.size != y_arr.size:
            print(f"Dimension mismatch for {label}: x={x_arr.size}, y={y_arr.size}")
            return

        finite_mask = np.isfinite(x_arr) & np.isfinite(y_arr)
        if not np.any(finite_mask):
            print(f"No finite points to plot for {label}")
            return

        x_plot = x_arr[finite_mask]
        y_plot = y_arr[finite_mask]
        if "vs Time" in plot_type:
            self.ax.plot(x_plot, y_plot, label=label) # No marker for time traces
        else:
            self.ax.plot(x_plot, y_plot, marker='o', label=label)

    def _get_flux_species_subscript(self):
        """Build TeX-friendly species subscript token for flux y-axis labels."""
        # For multi-species overlay, use a generic species index subscript.
        try:
            if self.plot_all_species_var.get():
                return "s"
        except Exception:
            pass

        selected = ""
        try:
            selected = self.species_var.get().strip()
        except Exception:
            selected = ""

        if not selected:
            return "s"

        if selected == "Main Ion (D+T)":
            return "i"
        if selected == "All Ions":
            return "i"

        species_name = selected.split(" (Z=")[0].strip()
        name_to_sub = {
            "Electron": "e",
            "Hydrogen": "H",
            "Deuterium": "D",
            "Tritium": "T",
            "Helium-3": "He3",
            "Helium-4": "He4",
            "Lithium": "Li",
            "Beryllium": "Be",
            "Boron": "B",
            "Carbon": "C",
            "Nitrogen": "N",
            "Oxygen": "O",
            "Neon": "Ne",
            "Argon": "Ar",
            "Nickel": "Ni",
            "Tungsten": "W",
        }
        if species_name in name_to_sub:
            return name_to_sub[species_name]

        # Fallback from charge sign in the selection string.
        match = re.search(r"Z=([-\d\.]+)", selected)
        if match:
            try:
                z = float(match.group(1))
                return "e" if z < 0 else "i"
            except Exception:
                pass

        # Last resort: compact plain-text token.
        token = re.sub(r"[^A-Za-z0-9]+", "", species_name)
        return token if token else "s"

    def _apply_axis_labels(self, plot_type):
        """Apply default axis labels for plot types that did not set labels explicitly."""
        # Respect labels already set by specialized plotting functions.
        cur_x = self.ax.get_xlabel().strip()
        cur_y = self.ax.get_ylabel().strip()

        x_label = None
        y_label = None

        if plot_type == "Frequency":
            x_label = r"$k_y \rho_s$"
            y_label = r"$\omega\ (c_s/a)$"
        elif plot_type == "Growth Rate":
            x_label = r"$k_y \rho_s$"
            y_label = r"$\gamma\ (c_s/a)$"
        elif plot_type == "ZF ExB Shearing Rate":
            x_label = r"$t\ (a/c_s)$"
            y_label = r"$\omega_{E\times B}^{ZF}\ (c_s/a)$"
        elif plot_type == "ZF ExB Shearing Spectrum":
            x_label = r"$k_x \rho_s$"
            y_label = r"$\omega_{E\times B}^{ZF}\ (c_s/a)$"
        elif plot_type == "ZF Phi Spectrum (theta0)":
            x_label = r"$k_x \rho_s$"
            y_label = r"$\langle |\phi_{ZF}(k_x,\theta=0)|/\rho_s \rangle_t$"
        elif plot_type in ["ZF ExB vs gamma_lin (kx=ky)", "ZF ExB Fig4 (kx=ky)"]:
            x_label = r"$k\rho_s\ (k_x=k_y)$"
            y_label = r"$\langle\omega_{E\times B}^{ZF}\rangle,\ k_x\langle V_{ZF}\rangle,\ \gamma_{lin}\ (c_s/a)$"
        elif plot_type == "Integration Error":
            x_label = r"$t\ (a/c_s)$"
            y_label = r"$\mathrm{Integration~Error}$"
        elif plot_type == "Radial Correlation (rcorr_phi)":
            x_label = r"$\Delta r/\rho_s$"
            y_label = r"$C(\Delta r)$"
        elif "Flux" in plot_type:
            use_real_ion_norm = self._use_flux_real_ion_norm()
            if "vs ky" in plot_type:
                x_label = r"$k_y \rho_{s,ri}$" if use_real_ion_norm else r"$k_y \rho_s$"
            elif "vs kx" in plot_type:
                x_label = r"$k_x \rho_s$"
            else:
                x_label = r"$(c_{s,ri}/a)\,t$" if use_real_ion_norm else r"$t\ (a/c_s)$"
            sub = self._get_flux_species_subscript()
            if "Energy" in plot_type:
                denom = r"Q_{GB,\mathrm{ri}}" if use_real_ion_norm else r"Q_{GB}"
                y_label = rf"$Q_{{{sub}}}/{denom}$"
            else:
                if use_real_ion_norm:
                    y_label = rf"$\Gamma_{{{sub}}}/\Gamma_{{GB,\mathrm{{ri}}}}$"
                else:
                    y_label = rf"$\Gamma_{{{sub}}}/Q_{{GB}}$"
        elif "vs ky" in plot_type:
            x_label = r"$k_y \rho_s$"
        elif "vs kx" in plot_type:
            x_label = r"$k_x \rho_s$"
        elif "vs Time" in plot_type:
            x_label = r"$t\ (a/c_s)$"

        if not cur_x and x_label:
            self.ax.set_xlabel(x_label)
        if not cur_y and y_label:
            self.ax.set_ylabel(y_label)

    def _reset_plot_area(self):
        """Reset figure/canvas state before starting a new plot action."""
        self._stop_animation()
        self._energy_entropy_axes = None
        self._energy_entropy_axes_active = False
        self._energy_entropy_case_style = {}
        self.fig.clear()
        self._reset_figure_layout_defaults()
        self._apply_global_plot_color_style()
        self.ax = self.fig.add_subplot(111)


    def _plot_selected_cases(self, selected_cases, plot_type):
        """Plot all selected cases using current plot mode."""
        for case_name in selected_cases:
            data = self.cases[case_name]
            try:
                self._plot_single_case(data, case_name, plot_type)
            except Exception as e:
                print(f"Error plotting {case_name}: {e}")
                traceback.print_exc()

    def _finalize_plot(self, plot_type_selection, plot_type, display_plot_type):
        """Finalize legends, labels, scales, title, and redraw canvas."""
        if getattr(self, "_energy_entropy_axes_active", False):
            axes = getattr(self, "_energy_entropy_axes", None)
            if isinstance(axes, (list, tuple)) and len(axes) == 2:
                for axi in axes:
                    h, _l = axi.get_legend_handles_labels()
                    if h:
                        axi.legend(loc='best')
                    axi.grid(True, color=GRID_COLOR, linewidth=GRID_LINEWIDTH)
                    try:
                        axi.set_axisbelow(True)
                    except Exception:
                        pass
            self._apply_unified_visual_style_to_figure()
            self.canvas.draw()
            return

        if not self._is_standard_line_plot(plot_type):
            self._apply_unified_visual_style_to_figure()
            self.canvas.draw()
            return

        handles, _labels = self.ax.get_legend_handles_labels()
        if handles:
            self.ax.legend()
        self.ax.grid(True, color=GRID_COLOR, linewidth=GRID_LINEWIDTH)
        self.ax.set_axisbelow(True)

        title = f"Comparison: {display_plot_type}"
        if self.norm_ky_var.get() and plot_type_selection in ["Frequency", "Growth Rate"]:
            title = f"{title}/ky"
        self.ax.set_title(title)

        self._apply_axis_labels(plot_type)

        if self.log_x_var.get():
            self.ax.set_xscale('log')
        if self.log_y_var.get():
            self.ax.set_yscale('log')

        self._apply_unified_visual_style_to_figure()
        self.canvas.draw()

    def plot_comparison(self):
        """
        Main GUI action: render current comparison plot for selected cases.

        This method builds effective plot mode, applies single-case constraints
        for contour-like plots, dispatches plotting, then finalizes figure state.
        """
        self._reset_plot_area()
        plot_type_selection, plot_type, display_plot_type = self._build_effective_plot_type()
        selected_cases = self._get_selected_case_names()

        if self._is_contour_like_plot(plot_type) and len(selected_cases) > 1:
            messagebox.showwarning(
                "Warning",
                f"{plot_type} only supports one case at a time. Plotting the first selected case."
            )
            selected_cases = list(selected_cases)[:1]

        # Hard stop for parity POD when theta resolution is insufficient.
        # Leave plot area blank and show warning dialog.
        if plot_type == "POD Parity" and len(selected_cases) > 0:
            data0 = self.cases.get(selected_cases[0], None)
            if data0 is not None and self._case_theta_resolution_is_insufficient(data0):
                self._show_pod_theta_resolution_warning()
                self.canvas.draw()
                return

        if (
            plot_type_selection == "Flux"
            and ("vs 2D" in plot_type)
            and len(selected_cases) > 0
        ):
            self._plot_flux_vs_2d_selected_cases(selected_cases, plot_type)
        elif (
            plot_type_selection == "Energy balance"
            and (plot_type == "Energy Balance vs 2D")
            and len(selected_cases) > 0
        ):
            self._plot_energy_balance_vs_2d_selected_cases(selected_cases)
        elif self.plot_all_species_var.get() and "Flux" in plot_type and len(selected_cases) > 0:
            self._plot_all_species_flux_first_case(selected_cases, plot_type)
        else:
            self._plot_selected_cases(selected_cases, plot_type)

        self._finalize_plot(plot_type_selection, plot_type, display_plot_type)

    def plot_case_info(self):
        """
        Show one case `out.cgyro.info` in a continuously scrollable view.

        Selection rule:
        - no selected case -> use the first loaded case
        - one or more selected cases -> use the first selected case
        """
        self._reset_plot_area()

        total_cases = int(self.case_listbox.size())
        if total_cases <= 0:
            messagebox.showwarning("Warning", "No loaded case to display out.cgyro.info.")
            self.canvas.draw()
            return

        selected_indices = list(self.case_listbox.curselection())
        if len(selected_indices) == 0:
            case_name = self.case_listbox.get(0)
        else:
            case_name = self.case_listbox.get(selected_indices[0])
            if len(selected_indices) > 1:
                messagebox.showwarning(
                    "Warning",
                "Multiple cases selected. Displaying out.cgyro.info of the first selected case."
                )

        data = self.cases.get(case_name, None)
        if data is None:
            self.ax.axis('off')
            self.ax.text(
                0.01,
                0.99,
                f"Failed to resolve case object: {case_name}",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            self.canvas.draw()
            return

        self._plot_other_case_info(data, case_name)
        self.canvas.draw()

    def plot_input_diff(self):
        """
        Show `input.cgyro` differences in scrollable text view.

        Selection rule:
        - >=2 selected: compare all selected cases
        - 1 selected: compare selected case against all other loaded cases
        - 0 selected: compare all loaded cases

        Rendering mode:
        - always scalar-parameter difference matrix (only changed parameters)
        """
        self._reset_plot_area()

        total_cases = int(self.case_listbox.size())
        if total_cases < 2:
            messagebox.showwarning("Warning", "Need at least two loaded cases to diff input.cgyro.")
            self.canvas.draw()
            return

        all_cases = [self.case_listbox.get(i) for i in range(total_cases)]
        selected_indices = list(self.case_listbox.curselection())
        selected_cases = [self.case_listbox.get(i) for i in selected_indices]

        if len(selected_cases) >= 2:
            compare_cases = list(selected_cases)
        elif len(selected_cases) == 1:
            base = selected_cases[0]
            compare_cases = [base] + [nm for nm in all_cases if nm != base]
        else:
            compare_cases = list(all_cases)

        if len(compare_cases) < 2:
            messagebox.showwarning("Warning", "Need at least two distinct cases to diff.")
            self.canvas.draw()
            return

        self._plot_input_diff_multi_cases(compare_cases)
        self.canvas.draw()

    def _plot_input_diff_two_cases(self, case_a, case_b):
        """Render colored unified diff for two selected cases."""
        data_a = self.cases.get(case_a, None)
        data_b = self.cases.get(case_b, None)
        if data_a is None or data_b is None:
            self.ax.axis('off')
            self.ax.text(
                0.01,
                0.99,
                f"Failed to resolve case objects:\nA={case_a}\nB={case_b}",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            return

        dir_a = self._resolve_case_dir(data_a)
        dir_b = self._resolve_case_dir(data_b)
        path_a = os.path.join(dir_a, "input.cgyro") if dir_a else ""
        path_b = os.path.join(dir_b, "input.cgyro") if dir_b else ""

        if not path_a or not os.path.isfile(path_a) or not path_b or not os.path.isfile(path_b):
            self.ax.axis('off')
            self.ax.text(
                0.01,
                0.99,
                f"input.cgyro not found for one or both cases:\nA: {path_a}\nB: {path_b}",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            return

        try:
            with open(path_a, "r", encoding="utf-8", errors="replace") as fa:
                txt_a = fa.read().replace('\r\n', '\n').replace('\r', '\n')
            with open(path_b, "r", encoding="utf-8", errors="replace") as fb:
                txt_b = fb.read().replace('\r\n', '\n').replace('\r', '\n')
        except Exception as e:
            self.ax.axis('off')
            self.ax.text(
                0.01,
                0.99,
                f"Failed to read input.cgyro for diff:\n{e}",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            return

        lines_a = txt_a.split('\n')
        lines_b = txt_b.split('\n')
        diff_lines = list(
            difflib.unified_diff(
                lines_a,
                lines_b,
                fromfile=f"{case_a}/input.cgyro",
                tofile=f"{case_b}/input.cgyro",
                lineterm="",
            )
        )
        if len(diff_lines) == 0:
            diff_lines = ["No differences found."]

        header_lines = [
            "Diff input.cgyro (2-case unified diff)",
            f"A: {case_a}",
            f"B: {case_b}",
            f"Path A: {path_a}",
            f"Path B: {path_b}",
            "-" * 110,
            "",
        ]

        def line_color(line):
            if line.startswith('+++') or line.startswith('---'):
                return "#6a1b9a"
            if line.startswith('@@'):
                return "#1565c0"
            if line.startswith('+') and not line.startswith('+++'):
                return "#2e7d32"
            if line.startswith('-') and not line.startswith('---'):
                return "#c62828"
            return "#222222"

        colors = ["#0d47a1"] * len(header_lines) + [line_color(ln) for ln in diff_lines]

        self._case_info_case_name = f"Diff {case_a} vs {case_b}"
        self._case_info_file_path = f"{path_a}  <->  {path_b}"
        self._case_info_lines = header_lines + diff_lines
        self._case_info_line_colors = colors
        self._enable_case_info_scroll()

    def _plot_input_diff_multi_cases(self, compare_cases):
        """Render multi-case scalar-parameter diff matrix with highlights."""
        case_infos = []
        missing_cases = []
        for case_name in compare_cases:
            data = self.cases.get(case_name, None)
            if data is None:
                missing_cases.append((case_name, "", "case object missing"))
                continue
            case_dir = self._resolve_case_dir(data)
            path = os.path.join(case_dir, "input.cgyro") if case_dir else ""
            if not path or not os.path.isfile(path):
                missing_cases.append((case_name, path, "input.cgyro not found"))
                continue
            scalars = self._read_input_cgyro_scalars(case_dir)
            case_infos.append((case_name, path, scalars))

        if len(case_infos) < 2:
            self.ax.axis('off')
            self.ax.text(
                0.01,
                0.99,
                "Need at least two valid cases with input.cgyro for multi-case diff.",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            return

        all_keys = set()
        for _case_name, _path, scalars in case_infos:
            all_keys.update(str(k) for k in scalars.keys())

        def norm_val(v):
            if v is None:
                return "<MISSING>"
            return str(v).strip()

        diff_keys = []
        for key in sorted(all_keys):
            vals = [norm_val(sc.get(key, None)) for _nm, _p, sc in case_infos]
            if len(set(vals)) > 1:
                diff_keys.append(key)

        lines = []
        colors = []

        def add_line(text, color="#222222"):
            lines.append(str(text))
            colors.append(str(color))

        add_line("Diff input.cgyro (multi-case scalar matrix)", "#0d47a1")
        add_line(f"Compared cases: {len(case_infos)}", "#0d47a1")
        for idx, (nm, path, _sc) in enumerate(case_infos, start=1):
            add_line(f"  [{idx}] {nm} -> {path}", "#37474f")

        if len(missing_cases) > 0:
            add_line("", "#222222")
            add_line("Skipped cases:", "#c62828")
            for nm, pth, reason in missing_cases:
                add_line(f"  - {nm}: {reason} ({pth})", "#c62828")

        add_line("-" * 110, "#0d47a1")
        add_line("", "#222222")

        if len(diff_keys) == 0:
            add_line("No scalar differences found among valid cases.", "#2e7d32")
        else:
            add_line(f"Differing scalar parameters: {len(diff_keys)}", "#ef6c00")
            add_line("", "#222222")
            for key in diff_keys:
                add_line(f"[{key}]", "#1565c0")
                case_vals = [norm_val(sc.get(key, None)) for _nm, _p, sc in case_infos]
                uniq = set(case_vals)
                for nm, val in zip([ci[0] for ci in case_infos], case_vals):
                    color = "#c62828" if val == "<MISSING>" else ("#2e7d32" if len(uniq) == 1 else "#222222")
                    add_line(f"  {nm}: {val}", color)
                add_line("", "#222222")

        self._case_info_case_name = f"Multi-case Diff ({len(case_infos)} cases)"
        self._case_info_file_path = "input.cgyro scalar diff"
        self._case_info_lines = lines
        self._case_info_line_colors = colors
        self._enable_case_info_scroll()

    def _estimate_case_info_lines_per_view(self):
        """Estimate visible monospace lines using current axes pixel height."""
        # Rough line height for fontsize=9 monospace.
        line_px = 15.0
        ax_h_px = None
        try:
            ax_h_px = float(getattr(self.ax, "bbox", None).height)
        except Exception:
            ax_h_px = None

        if not ax_h_px or ax_h_px <= 0:
            try:
                fig_h_px = float(self.fig.get_figheight() * self.fig.dpi)
            except Exception:
                fig_h_px = 800.0
            ax_h_px = fig_h_px * 0.72

        # Reserve space for title/margins to avoid clipping the last line.
        usable_h_px = max(120.0, ax_h_px - 52.0)
        n_lines = int(usable_h_px / line_px)
        return max(6, n_lines)

    def _render_case_info_window(self, top_line_idx):
        """Render one vertical window of the loaded out.cgyro.info lines."""
        lines = list(getattr(self, "_case_info_lines", []))
        line_colors = list(getattr(self, "_case_info_line_colors", []))
        case_name = str(getattr(self, "_case_info_case_name", ""))

        if len(lines) == 0:
            self.ax.clear()
            self.ax.axis('off')
            self.ax.set_title("Case Info")
            self.ax.text(
                0.01,
                0.99,
                "No out.cgyro.info content.",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            return

        lines_per_view = self._estimate_case_info_lines_per_view()
        max_top = max(0, len(lines) - lines_per_view)
        top = max(0, min(int(top_line_idx), max_top))
        end = min(len(lines), top + lines_per_view)

        # Keep manual pager in sync with current view bounds.
        if getattr(self, "_manual_pager_active", False):
            self.total_frames = max_top + 1
            self.current_frame = top
            self._update_manual_pager_status()

        self.ax.clear()
        self.ax.axis('off')
        self.ax.set_title(f"Case Info: {case_name}  [{top + 1}-{end}/{len(lines)}]")
        view_lines = lines[top:end]
        view_colors = line_colors[top:end] if len(line_colors) == len(lines) else ["#222222"] * len(view_lines)
        y0 = 0.99
        dy = 0.96 / max(1, lines_per_view)
        for idx, txt in enumerate(view_lines):
            y = y0 - idx * dy
            if y < -0.05:
                break
            clr = view_colors[idx] if idx < len(view_colors) else "#222222"
            self.ax.text(
                0.01,
                y,
                txt,
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=9,
                family='monospace',
                color=clr,
                clip_on=True,
            )

    def _on_case_info_scroll(self, event):
        """Mouse-wheel handler for smooth scrolling in out.cgyro.info mode."""
        if not getattr(self, "_case_info_scroll_active", False):
            return
        if not getattr(self, "_manual_pager_active", False):
            return

        step = getattr(event, "step", 0)
        if step == 0:
            button = str(getattr(event, "button", "")).lower()
            if button == "up":
                step = 1
            elif button == "down":
                step = -1

        if step == 0:
            return

        # Scroll multiple lines per wheel notch for smoother movement.
        delta_lines = max(1, int(abs(step))) * 3
        total = int(getattr(self, "total_frames", 0))
        if total <= 0:
            return

        current = int(getattr(self, "current_frame", 0))
        if step > 0:
            new_top = max(0, current - delta_lines)
        else:
            new_top = min(total - 1, current + delta_lines)

        if new_top == current:
            return

        self.current_frame = new_top
        if hasattr(self, "anim_update_func") and self.anim_update_func:
            self.anim_update_func(self.current_frame)
        self._update_manual_pager_status()
        try:
            self.canvas.draw_idle()
        except Exception:
            self.canvas.draw()

    def _enable_case_info_scroll(self):
        """Enable continuous scroll controls (wheel + Prev/Next) for case info text."""
        lines = list(getattr(self, "_case_info_lines", []))
        if len(lines) == 0:
            return

        old_cid = getattr(self, "_case_info_scroll_cid", None)
        if old_cid is not None:
            try:
                self.canvas.mpl_disconnect(old_cid)
            except Exception:
                pass

        self._case_info_scroll_cid = self.canvas.mpl_connect('scroll_event', self._on_case_info_scroll)
        self._case_info_scroll_active = True

        lines_per_view = self._estimate_case_info_lines_per_view()
        max_top = max(0, len(lines) - lines_per_view)

        def _update(top_idx):
            self._render_case_info_window(top_idx)

        self._enable_manual_pager(
            total_frames=max_top + 1,
            update_func=_update,
            start_frame=0,
            label="Scroll",
        )
        _update(0)

    # ------------------------------------------------------------------
    # Plot Backends
    # ------------------------------------------------------------------





    def _resolve_axis_selection_index(
        self, axis_values, raw_text, axis_name, label, prefer_value=False
    ):
        """
        Resolve user text into one axis index.

        Supports both:
        - integer index (e.g. `3` or `idx:3`)
        - physical value mapped to nearest axis point (e.g. `0.25` or `val:0.25`)

        When `prefer_value=True`, plain numeric text is interpreted as value first.
        """
        axis = np.asarray(axis_values, dtype=float).reshape(-1)
        n = int(axis.size)
        if n <= 0:
            return 0

        txt = str(raw_text).strip()
        if len(txt) == 0:
            return 0

        def _nearest_value_idx(val):
            """Return index of finite axis entry closest to the requested value."""
            finite_mask = np.isfinite(axis)
            if not np.any(finite_mask):
                return None
            axis_finite = axis[finite_mask]
            idx_finite = int(np.argmin(np.abs(axis_finite - float(val))))
            idx_all = np.where(finite_mask)[0][idx_finite]
            return int(idx_all)

        txt_low = txt.lower()
        if txt_low.startswith('idx:'):
            try:
                idx = int(txt.split(':', 1)[1].strip())
                idx = max(0, min(idx, n - 1))
                return idx
            except Exception:
                print(f"Warning: invalid {axis_name} index '{txt}' for {label}; fallback to 0.")
                return 0
        if txt_low.startswith('val:'):
            try:
                idx = _nearest_value_idx(float(txt.split(':', 1)[1].strip()))
                if idx is not None:
                    return idx
            except Exception:
                pass
            print(f"Warning: invalid {axis_name} value '{txt}' for {label}; fallback to 0.")
            return 0

        if prefer_value:
            try:
                idx = _nearest_value_idx(float(txt))
                if idx is not None:
                    return idx
            except Exception:
                pass
            if re.match(r'^[+-]?\d+$', txt):
                idx = int(txt)
                idx = max(0, min(idx, n - 1))
                return idx
            print(f"Warning: invalid {axis_name} selection '{txt}' for {label}; fallback to index 0.")
            return 0

        if re.match(r'^[+-]?\d+$', txt):
            idx = int(txt)
            idx = max(0, min(idx, n - 1))
            return idx

        try:
            idx = _nearest_value_idx(float(txt))
            if idx is not None:
                return idx
        except Exception:
            pass

        print(f"Warning: invalid {axis_name} selection '{txt}' for {label}; fallback to index 0.")
        return 0

    def _build_theta_over_pi_axis(self, data, n_theta):
        """Build `z/pi` axis for POD mode-structure plotting."""
        if n_theta <= 0:
            return np.array([])

        # 0) Best match for kxky_* output: pygacode's sampled theta grid
        #    (constructed from THETA_PLOT selection in out.cgyro.grids).
        theta_p = getattr(data, 'thetap', None)
        if theta_p is not None:
            arr = np.asarray(theta_p, dtype=float).reshape(-1)
            arr = arr[np.isfinite(arr)]
            if arr.size == n_theta:
                return arr / np.pi

        candidates = [
            ('theta_over_pi', 1.0),
            ('theta_b_over_pi', 1.0),
            ('theta', 1.0 / np.pi),
            ('theta_plot_grid', 1.0 / np.pi),
        ]
        candidate_arrays = []
        for attr, scale in candidates:
            vals = getattr(data, attr, None)
            if vals is None:
                continue
            arr = np.asarray(vals, dtype=float).reshape(-1)
            if arr.size <= 1:
                continue
            arr = arr[np.isfinite(arr)]
            if arr.size <= 1:
                continue
            arr = arr * scale
            candidate_arrays.append((attr, arr))

        # Prefer exact-size theta grid when available.
        for _attr, arr in candidate_arrays:
            if arr.size == n_theta:
                return arr

        # If only full theta grid is available, mimic CGYRO THETA_PLOT downsampling:
        # select every m-th point where m = n_full / n_theta.
        for attr_name in ('theta', 'theta_over_pi'):
            vals = getattr(data, attr_name, None)
            if vals is None:
                continue
            arr = np.asarray(vals, dtype=float).reshape(-1)
            arr = arr[np.isfinite(arr)]
            if arr.size <= n_theta or (arr.size % n_theta) != 0:
                continue
            m = int(arr.size // n_theta)
            if m <= 0:
                continue
            arr_ds = arr[::m][:n_theta]
            if arr_ds.size != n_theta:
                continue
            if attr_name == 'theta':
                arr_ds = arr_ds / np.pi
            print(
                f"Info: downsampled {attr_name} ({arr.size} -> {n_theta}) "
                f"with stride {m} for POD parity z/pi axis."
            )
            return np.asarray(arr_ds, dtype=float).reshape(-1)

        # Fallback: resample the closest candidate to n_theta over the FULL range.
        # This avoids accidental truncation (e.g. old behavior arr[:n_theta]).
        if len(candidate_arrays) > 0:
            best_attr, best_arr = min(
                candidate_arrays, key=lambda item: abs(int(item[1].size) - int(n_theta))
            )
            if best_arr.size != n_theta:
                src_x = np.linspace(0.0, 1.0, int(best_arr.size))
                dst_x = np.linspace(0.0, 1.0, int(n_theta))
                best_arr = np.interp(dst_x, src_x, best_arr)
                print(
                    f"Info: resampled {best_attr} from {len(src_x)} to {n_theta} points "
                    f"for POD parity z/pi axis."
                )
            return np.asarray(best_arr, dtype=float).reshape(-1)

        # Fallback: normalized symmetric coordinate.
        return np.linspace(-1.0, 1.0, n_theta)

    @staticmethod
    def _build_extended_z_over_pi_axis(theta_over_pi, chain_field_idx, chain_kx=None, p_axis=None):
        """
        Build extended ballooning coordinate for connected-kx POD modes.

        Uses:
            z/pi = theta/pi + 2*p_rel
        where `p_rel` is centered at the kx~0 location. If a physical `p` axis
        is available from `out.cgyro.grids`, it is preferred; otherwise index-based
        connected order is used.
        """
        theta = np.asarray(theta_over_pi, dtype=float).reshape(-1)
        chain_idx = np.asarray(chain_field_idx, dtype=int).reshape(-1)
        n_theta = int(theta.size)
        n_chain = int(chain_idx.size)
        if n_theta <= 0 or n_chain <= 0:
            return np.array([])

        # Center by the mode closest to kx=0 when possible.
        i_center = n_chain // 2
        ckx = np.asarray(chain_kx if chain_kx is not None else [], dtype=float).reshape(-1)
        if ckx.size == n_chain and np.any(np.isfinite(ckx)):
            finite = np.where(np.isfinite(ckx))[0]
            if finite.size > 0:
                i_center = int(finite[np.argmin(np.abs(ckx[finite]))])

        p_rel = None
        p_arr = np.asarray(p_axis if p_axis is not None else [], dtype=float).reshape(-1)
        if p_arr.size > 0:
            if np.max(chain_idx) < p_arr.size:
                p_sel = p_arr[chain_idx]
                if np.all(np.isfinite(p_sel)):
                    p_rel = p_sel - float(p_sel[i_center])

        if p_rel is None:
            p_rel = np.arange(n_chain, dtype=float) - float(i_center)

        z_blocks = [theta + 2.0 * float(pv) for pv in p_rel]
        return np.concatenate(z_blocks)

    @staticmethod
    def _parity_factor(profile, z_over_pi):
        """
        Compute parity factor:
            P = abs(integral(A dz)) / integral(abs(A) dz)
        using z/pi as integration coordinate (constant scaling cancels).
        """
        z = np.asarray(z_over_pi, dtype=float).reshape(-1)
        a = np.asarray(profile, dtype=complex).reshape(-1)
        n = min(z.size, a.size)
        if n <= 1:
            return np.nan
        z = z[:n]
        a = a[:n]
        finite = np.isfinite(z) & np.isfinite(np.real(a)) & np.isfinite(np.imag(a))
        if np.sum(finite) <= 1:
            return np.nan
        z = z[finite]
        a = a[finite]
        # Enforce monotonic integration axis and align origin to physical center.
        order = np.argsort(z)
        z = z[order]
        a = a[order]
        z0 = z[int(np.argmin(np.abs(z)))]
        z = z - z0
        # NumPy 2.x prefers trapezoid; keep backward compatibility.
        integrate = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
        denom = integrate(np.abs(a), x=z)
        if not np.isfinite(denom) or denom <= 1.0e-14:
            return np.nan
        num = np.abs(integrate(a, x=z))
        return float(num / denom)

    @staticmethod
    def _prepare_pod_display_curve(z_axis, profile, z_window_pi):
        """
        Prepare one POD profile on displayed z/pi domain.

        Returns sorted and optionally window-cropped arrays `(z_plot, a_plot)`.
        """
        z = np.asarray(z_axis, dtype=float).reshape(-1)
        a = np.asarray(profile, dtype=complex).reshape(-1)
        n = min(z.size, a.size)
        if n <= 1:
            return np.array([]), np.array([], dtype=complex)

        z = z[:n]
        a = a[:n]
        finite = np.isfinite(z) & np.isfinite(np.real(a)) & np.isfinite(np.imag(a))
        if np.sum(finite) <= 1:
            return np.array([]), np.array([], dtype=complex)
        z = z[finite]
        a = a[finite]

        order = np.argsort(z)
        z = z[order]
        a = a[order]

        zmax = float(np.nanmax(np.abs(z))) if z.size > 0 else 0.0
        zwin = float(z_window_pi)
        if np.isfinite(zmax) and np.isfinite(zwin) and zwin > 0 and zmax > zwin + 1e-12:
            mask = np.abs(z) <= zwin
            if np.count_nonzero(mask) >= 2:
                z = z[mask]
                a = a[mask]

        return z, a

    @staticmethod
    def _parity_even_odd_ratios(profile, z_axis):
        """
        Compute mirror even/odd parity ratios on displayed domain:

            P_even = ∫ |A(z)+A(-z)|^2 dz / (4 ∫ |A(z)|^2 dz)
            P_odd  = ∫ |A(z)-A(-z)|^2 dz / (4 ∫ |A(z)|^2 dz)

        Returns `(P_even, P_odd)`.
        """
        z = np.asarray(z_axis, dtype=float).reshape(-1)
        a = np.asarray(profile, dtype=complex).reshape(-1)
        n = min(z.size, a.size)
        if n <= 2:
            return np.nan, np.nan
        z = z[:n]
        a = a[:n]
        finite = np.isfinite(z) & np.isfinite(np.real(a)) & np.isfinite(np.imag(a))
        if np.sum(finite) <= 2:
            return np.nan, np.nan
        z = z[finite]
        a = a[finite]

        order = np.argsort(z)
        z = z[order]
        a = a[order]

        zsym = min(abs(float(np.nanmin(z))), abs(float(np.nanmax(z))))
        if (not np.isfinite(zsym)) or zsym <= 1.0e-12:
            return np.nan, np.nan

        n_eval = max(256, int(len(z)))
        z_eval = np.linspace(-zsym, zsym, n_eval)

        a_re = np.interp(z_eval, z, np.real(a))
        a_im = np.interp(z_eval, z, np.imag(a))
        a_eval = a_re + 1j * a_im

        a_re_m = np.interp(-z_eval, z, np.real(a))
        a_im_m = np.interp(-z_eval, z, np.imag(a))
        a_m_eval = a_re_m + 1j * a_im_m

        integrate = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
        denom = 4.0 * integrate(np.abs(a_eval) ** 2, x=z_eval)
        if (not np.isfinite(denom)) or denom <= 1.0e-14:
            return np.nan, np.nan

        p_even = integrate(np.abs(a_eval + a_m_eval) ** 2, x=z_eval) / denom
        p_odd = integrate(np.abs(a_eval - a_m_eval) ** 2, x=z_eval) / denom

        if not np.isfinite(p_even):
            p_even = np.nan
        if not np.isfinite(p_odd):
            p_odd = np.nan
        return float(p_even), float(p_odd)

    @staticmethod
    def _connected_kx_smoothness(block):
        """Return a normalized adjacent-difference score (smaller is smoother)."""
        arr = np.asarray(block, dtype=complex)
        if arr.ndim < 1 or arr.shape[0] <= 1:
            return 0.0
        diff = arr[1:, ...] - arr[:-1, ...]
        num = np.nansum(np.abs(diff) ** 2)
        den = np.nansum(np.abs(arr) ** 2)
        if not np.isfinite(num):
            return np.inf
        if not np.isfinite(den) or den <= 1.0e-30:
            den = 1.0e-30
        return float(num / den)

    @staticmethod
    def _enforce_connected_kx_phase_continuity(block):
        """
        Align neighboring connected-kx slices by overlap phase.

        For adjacent slices v_{k-1}, v_k, apply v_k <- v_k * exp(-i*arg(<v_{k-1},v_k>))
        so their overlap is real-positive.
        """
        arr = np.asarray(block, dtype=complex).copy()
        if arr.ndim < 1 or arr.shape[0] <= 1:
            return arr

        eps = 1.0e-30
        for i in range(1, arr.shape[0]):
            prev = np.asarray(arr[i - 1], dtype=complex).reshape(-1)
            curr = np.asarray(arr[i], dtype=complex).reshape(-1)
            finite = (
                np.isfinite(np.real(prev)) & np.isfinite(np.imag(prev)) &
                np.isfinite(np.real(curr)) & np.isfinite(np.imag(curr))
            )
            if np.count_nonzero(finite) == 0:
                continue
            ov = np.vdot(prev[finite], curr[finite])
            if np.isfinite(np.real(ov)) and np.isfinite(np.imag(ov)) and abs(ov) > eps:
                arr[i] *= (np.conj(ov) / abs(ov))
        return arr

    def _estimate_connected_kx_balloon_phase(self, data, chain_field_idx, ky_val):
        """
        Estimate CGYRO-like radial phase factor from source mapping:
            exp(-2*pi*i*(px+px0)*k_theta_base*n_tor*rmin*sign_qs)
        using available metadata (p, q, ky0). Global phase is removed.
        """
        pvec = np.asarray(getattr(data, 'p', []), dtype=float).reshape(-1)
        idx = np.asarray(chain_field_idx, dtype=int).reshape(-1)
        if pvec.size <= 0 or idx.size <= 0:
            return None, None
        if np.any((idx < 0) | (idx >= pvec.size)):
            return None, None

        q = float(getattr(data, 'q', np.nan))
        ky0 = float(getattr(data, 'ky0', np.nan))
        if (not np.isfinite(q)) or (not np.isfinite(ky0)) or abs(ky0) <= 1.0e-12:
            return None, None

        n_tor = float(ky_val / ky0)
        if not np.isfinite(n_tor):
            return None, None
        n_tor_round = float(np.round(n_tor))
        if abs(n_tor - n_tor_round) < 0.15:
            n_tor = n_tor_round

        px = pvec[idx]
        i_ref = int(np.argmin(np.abs(px)))
        px_rel = px - float(px[i_ref])
        phase = np.exp(-2.0 * np.pi * 1j * px_rel * q * n_tor)
        if not np.all(np.isfinite(np.real(phase)) & np.isfinite(np.imag(phase))):
            return None, None
        return np.asarray(phase, dtype=complex), n_tor

    def _repair_connected_kx_stitching(self, data, label, block, chain_field_idx, chain_kx, ky_val):
        """
        Enforce alternating-phase unwrap for connected-kx stitching.

        Mandatory fix for CGYRO Eulerian alternating phase:
        apply (-1)^p relative to kx~0 center before SVD, then phase-continuize.
        """
        arr = np.asarray(block, dtype=complex)
        n_chain = int(arr.shape[0]) if arr.ndim > 0 else 0
        if n_chain <= 1:
            return arr, ""

        raw_score = self._connected_kx_smoothness(arr)

        if chain_kx.size == n_chain and np.all(np.isfinite(chain_kx)):
            i_center = int(np.argmin(np.abs(chain_kx)))
        else:
            i_center = n_chain // 2
        p_rel = np.arange(n_chain, dtype=int) - i_center
        alt_pattern = np.power(-1.0, p_rel).astype(float)

        # 1) Forced alternating-phase unwrap.
        arr_alt = arr * alt_pattern[:, None, None]
        alt_score = self._connected_kx_smoothness(arr_alt)

        # 2) Neighbor continuity clean-up after unwrap.
        arr_fix = self._enforce_connected_kx_phase_continuity(arr_alt)
        fix_score = self._connected_kx_smoothness(arr_fix)

        print(
            f"INFO: connected-kx stitching for {label}: "
            f"raw={raw_score:.3e} -> alt(-1)^p={alt_score:.3e} -> "
            f"alt(-1)^p+continuity={fix_score:.3e}"
        )
        return arr_fix, "forced-alt(-1)^p+continuity"








            # self.ax.set_yscale('log')

    def _get_zf_exb_phi_kx_t(self, data, label):
        """
        Extract collect-compatible zonal potential `phi(kx,t)` for ZF diagnostics.

        This follows `collect.py:get_zfshear` indexing logic as closely as possible:
        - use toroidal mode index `n=0` (array index 0, not nearest ky search)
        - use outboard-midplane index `i_theta_plot = 4*n_theta_plot//8`
          (effectively middle theta index)
        - use radial slice `[1:, ...]` (skip the first radial slot)
        """
        phi_raw = None
        if hasattr(data, 'kxky_phi'):
            phi_raw = np.asarray(data.kxky_phi)
        else:
            print(f"Loading big field data for {label}...")
            if self._ensure_bigfield_loaded(data, label) and hasattr(data, 'kxky_phi'):
                phi_raw = np.asarray(data.kxky_phi)

        # Fallback: infer nt from file size when metadata time length mismatches.
        if phi_raw is None and hasattr(data, 'extract'):
            try:
                _, fmt, raw = data.extract('.cgyro.kxky_phi', cmplx=True)
                if fmt != 'null':
                    n_radial = int(getattr(data, 'n_radial', 0))
                    theta_plot = int(getattr(data, 'theta_plot', 0))
                    n_n = int(getattr(data, 'n_n', 0))
                    spatial = n_radial * theta_plot * n_n
                    if spatial > 0:
                        nt = raw.size // spatial
                        rem = raw.size % spatial
                        if nt > 0:
                            if rem != 0:
                                print(
                                    f"Warning: truncating {rem} trailing entries in "
                                    f"bin.cgyro.kxky_phi for {label}."
                                )
                                raw = raw[:spatial * nt]
                            phi_raw = np.reshape(raw, (n_radial, theta_plot, n_n, nt), order='F')
                            data.kxky_phi = phi_raw
                            print(
                                f"Recovered kxky_phi for {label} via file-size inference "
                                f"(nt={nt})."
                            )
            except Exception as e:
                print(f"Fallback read failed for {label}: {e}")

        if phi_raw is None:
            print(f"No phi field data available for {label}")
            return None, None, None

        if np.iscomplexobj(phi_raw):
            if phi_raw.ndim == 4:
                phi = phi_raw
            elif phi_raw.ndim == 5:
                # Possible [nr,theta,species,nn,nt]
                phi = phi_raw[:, :, 0, :, :]
            else:
                print(f"Unsupported kxky_phi shape for {label}: {phi_raw.shape}")
                return None, None, None
        elif phi_raw.ndim == 5 and phi_raw.shape[0] == 2:
            # [2,nr,theta,nn,nt]
            phi = phi_raw[0] + 1j * phi_raw[1]
        elif phi_raw.ndim == 6 and phi_raw.shape[0] == 2:
            # [2,nr,theta,species,nn,nt]
            phi = phi_raw[0, :, :, 0, :, :] + 1j * phi_raw[1, :, :, 0, :, :]
        else:
            print(f"Unsupported kxky_phi shape for {label}: {phi_raw.shape}")
            return None, None, None

        # collect:get_zfshear uses n-index 0 directly.
        n_ky = int(phi.shape[2]) if phi.ndim >= 4 else 0
        if n_ky <= 0:
            print(f"No ky mode dimension available for {label}")
            return None, None, None
        ky_idx_0 = 0

        # collect:get_zfshear uses i_theta_plot = 4*n_theta_plot//8.
        n_theta = int(phi.shape[1]) if phi.ndim >= 4 else 0
        if n_theta <= 0:
            print(f"No theta dimension available for {label}")
            return None, None, None
        itheta0 = int((4 * n_theta) // 8)
        if itheta0 < 0:
            itheta0 = 0
        if itheta0 >= n_theta:
            itheta0 = n_theta - 1

        # collect:get_zfshear uses phi_cmplx[1:, i_theta_plot, 0, ind_t_ave].
        phi_zonal = phi[1:, itheta0, ky_idx_0, :]  # [n_radial-1, n_time]

        kx = np.asarray(getattr(data, 'kx', []), dtype=float).reshape(-1)
        n_radial = phi_zonal.shape[0]
        if kx.size == n_radial:
            phi_kx_t = phi_zonal
        elif kx.size == n_radial + 1:
            # If kx still contains the skipped first radial slot, drop it to align collect indexing.
            kx = kx[1:]
            phi_kx_t = phi_zonal
        elif kx.size > 0:
            n = min(kx.size, n_radial)
            print(
                f"Warning: mismatched kx/phi radial sizes for {label} "
                f"(kx={kx.size}, radial={n_radial}); truncating to {n}."
            )
            kx = kx[:n]
            phi_kx_t = phi_zonal[:n, :]
        else:
            length = float(getattr(data, 'length', 0.0))
            if length <= 0:
                print(f"No usable kx axis for {label}")
                return None, None, None
            p_index = np.arange(n_radial) - (n_radial // 2)
            kx = 2.0 * np.pi * p_index / length
            phi_kx_t = phi_zonal

        # Align with collect.get_zfshear: keep raw zonal phi here
        # (collect applies rho normalization in the final gamma_zf/zf_shear formula path).

        t = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        if t.size == 0:
            t = np.arange(phi_kx_t.shape[-1], dtype=float)

        nt = min(t.size, phi_kx_t.shape[-1])
        if nt <= 0:
            print(f"No time samples available for {label}")
            return None, None, None

        return kx, phi_kx_t[:, :nt], t[:nt]

    def _plot_zf_exb_shearing_rate(self, data, label, t_indices, t_start, t_end):
        """
        Legacy wrapper for zonal ExB shearing-rate time trace.

        New/clear name: `_plot_zf_exb_shearing_rate_time_trace`.
        """
        return self._plot_zf_exb_shearing_rate_time_trace(
            data, label, t_indices, t_start, t_end
        )

    def _plot_zf_exb_shearing_rate_time_trace(self, data, label, t_indices, t_start, t_end):
        """
        Plot zonal ExB shearing-rate time trace from zonal potential spectrum.

        Output:
        - `x`: time axis.
        - `y`: scalar shearing-rate estimate from zonal `phi(kx,t)`.
        - Legend includes averaging window and `(mean, std)` when a window is valid.
        """
        x, y = self._compute_zf_exb_shearing_time_series(data, label)
        if x is None:
            return

        n_t = min(x.size, y.size)
        if n_t <= 0:
            print(f"No time samples available for {label}")
            return

        x = x[:n_t]
        y = y[:n_t]
        valid_t = np.asarray(t_indices, dtype=int)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]

        if valid_t.size > 0:
            y_subset = y[valid_t]
            mean_val = float(np.mean(y_subset))
            std_val = float(np.std(y_subset))
            avg_inner = self._format_avg_range_from_axis(x, valid_t, prefix="Avg").strip(" ()")
            line, = self.ax.plot(
                x,
                y,
                label=(
                    f"{label} ({avg_inner}, "
                    f"Mean: {mean_val:.2e}, Std: {std_val:.2e})"
                ),
            )
            self.ax.plot(
                [x[valid_t[0]], x[valid_t[-1]]],
                [mean_val, mean_val],
                linestyle='--',
                color=line.get_color(),
                linewidth=1.5,
            )
        else:
            self.ax.plot(x, y, label=label)

        self.ax.set_xlabel(r'$t (a/c_s)$')
        self.ax.set_ylabel(r'$\omega_{E\times B}^{ZF}\ (c_s/a)$')

    # ------------------------------------------------------------------
    # Data Loading / Shape Coercion Helpers
    # ------------------------------------------------------------------

    def _resolve_case_dir(self, data):
        """Resolve case directory path from loaded data-object attributes."""
        dir_path = getattr(data, 'dir', None)
        if not dir_path:
            dir_path = getattr(data, 'path', None)
        if not dir_path:
            return None

        dir_path = str(dir_path)
        if os.path.isfile(dir_path):
            return os.path.dirname(os.path.abspath(dir_path))
        return os.path.abspath(dir_path.rstrip('/\\'))

    def _coerce_kxky_complex(self, raw, label, tag, species_dependent=False):
        """Normalize diverse `kxky_*` payload shapes into canonical complex arrays."""
        arr = np.asarray(raw)
        if arr.size == 0:
            print(f"Empty {tag} array for {label}.")
            return None

        if np.iscomplexobj(arr):
            if species_dependent:
                if arr.ndim == 5:
                    return arr
                if arr.ndim == 4:
                    # No explicit species dimension; treat as ns=1.
                    return arr[:, :, None, :, :]
            else:
                if arr.ndim == 4:
                    return arr
                if arr.ndim == 5 and arr.shape[2] == 1:
                    return arr[:, :, 0, :, :]

        # Packed real/imag first axis.
        if arr.ndim > 0 and arr.shape[0] == 2:
            if species_dependent:
                if arr.ndim == 6:
                    return arr[0] + 1j * arr[1]  # [nr,theta,ns,ky,t]
                if arr.ndim == 5:
                    c = arr[0] + 1j * arr[1]      # [nr,theta,ky,t]
                    return c[:, :, None, :, :]
            else:
                if arr.ndim == 5:
                    return arr[0] + 1j * arr[1]   # [nr,theta,ky,t]
                if arr.ndim == 6 and arr.shape[3] == 1:
                    return arr[0, :, :, 0, :, :] + 1j * arr[1, :, :, 0, :, :]

        print(f"Unsupported {tag} shape for {label}: {arr.shape}")
        return None

    def _load_kxky_complex(self, data, label, attr_name, file_suffix, species_dependent=False):
        """Load one `kxky_*` complex field with fallback file-size inference."""
        raw = getattr(data, attr_name, None)

        if raw is None and not getattr(data, "_cmp_bigfield_attempted", False):
            print(f"Loading big field data for {label}...")
            self._ensure_bigfield_loaded(data, label)
            setattr(data, "_cmp_bigfield_attempted", True)
            raw = getattr(data, attr_name, None)

        if raw is None and hasattr(data, 'extract'):
            try:
                _, fmt, flat = data.extract(file_suffix, cmplx=True)
                if fmt != 'null':
                    n_radial = int(getattr(data, 'n_radial', 0))
                    theta_plot = int(getattr(data, 'theta_plot', 0))
                    n_n = int(getattr(data, 'n_n', 0))
                    n_species = int(getattr(data, 'n_species', 1))
                    if n_radial > 0 and theta_plot > 0 and n_n > 0:
                        spatial = n_radial * theta_plot * n_n
                        if species_dependent:
                            spatial *= max(1, n_species)
                        nt = flat.size // spatial
                        rem = flat.size % spatial
                        if nt > 0:
                            if rem != 0:
                                print(
                                    f"Warning: truncating {rem} trailing entries in "
                                    f"{file_suffix} for {label}."
                                )
                                flat = flat[:spatial * nt]
                            shape = (
                                (n_radial, theta_plot, n_species, n_n, nt)
                                if species_dependent else
                                (n_radial, theta_plot, n_n, nt)
                            )
                            raw = np.reshape(flat, shape, order='F')
                            setattr(data, attr_name, raw)
                            print(
                                f"Recovered {attr_name} for {label} via file-size "
                                f"inference (nt={nt})."
                            )
            except Exception as e:
                print(f"Fallback read failed for {label} ({attr_name}): {e}")

        if raw is None:
            return None

        return self._coerce_kxky_complex(raw, label, attr_name, species_dependent=species_dependent)

    @staticmethod
    def _field_attr_suffix_from_name(field_name):
        """Map field/moment short name to `(attr_name, file_suffix)`."""
        key = str(field_name).strip().lower()
        mapping = {
            'phi': ('kxky_phi', '.cgyro.kxky_phi'),
            'apar': ('kxky_apar', '.cgyro.kxky_apar'),
            'bpar': ('kxky_bpar', '.cgyro.kxky_bpar'),
            'density': ('kxky_n', '.cgyro.kxky_n'),
            'energy': ('kxky_e', '.cgyro.kxky_e'),
            'velocity': ('kxky_v', '.cgyro.kxky_v'),
        }
        return mapping.get(key, (None, None))

    def _load_named_kxky_complex(self, data, label, field_name, species_dependent=False):
        """Load named `kxky_*` data using canonical field-name mapping."""
        attr_name, file_suffix = self._field_attr_suffix_from_name(field_name)
        if not attr_name:
            return None
        return self._load_kxky_complex(
            data,
            label,
            attr_name,
            file_suffix,
            species_dependent=species_dependent,
        )

    @staticmethod
    def _midplane_theta_index(data, n_theta):
        """Return robust outboard-midplane theta index `4*n_theta/8`."""
        n_th = int(n_theta)
        if n_th <= 0:
            return 0
        i_theta = int((4 * n_th) // 8)
        if i_theta < 0:
            i_theta = 0
        if i_theta >= n_th:
            i_theta = n_th - 1
        return i_theta

    def _extract_midplane_kykxt(self, field_complex, data, label, drop_radial0=False, species_idx=0):
        """
        Extract complex `[kx, ky, t]` slice at outboard-midplane from canonical kxky field.

        Supported inputs:
        - field shape `[nr, theta, ky, t]` (species independent),
        - field shape `[nr, theta, species, ky, t]` (species dependent).
        """
        if field_complex is None:
            return None

        arr = np.asarray(field_complex)
        if arr.ndim == 5:
            n_species = int(arr.shape[2])
            s_idx = int(species_idx)
            if s_idx < 0 or s_idx >= n_species:
                s_idx = max(0, min(n_species - 1, s_idx))
            arr = arr[:, :, s_idx, :, :]

        if arr.ndim != 4:
            print(f"Unsupported field shape for {label}: {arr.shape}")
            return None

        n_r, n_th = int(arr.shape[0]), int(arr.shape[1])
        i_theta = self._midplane_theta_index(data, n_th)
        i_start = 1 if (drop_radial0 and n_r > 1) else 0
        return np.asarray(arr[i_start:, i_theta, :, :], dtype=complex)

    def _build_kx_axis(self, data, n_r, label):
        """Construct robust kx axis for a given radial dimension size."""
        if n_r <= 0:
            return None, None

        kx = np.asarray(getattr(data, 'kx', []), dtype=float).reshape(-1)
        radial_idx = np.arange(n_r, dtype=int)

        if kx.size == n_r:
            return kx, radial_idx

        if kx.size == n_r - 1:
            # Typical CGYRO layout: one leftmost special radial slot.
            return kx, np.arange(1, n_r, dtype=int)

        if kx.size > 0:
            n = min(kx.size, n_r)
            print(
                f"Warning: mismatched kx size for {label} (kx={kx.size}, radial={n_r}); "
                f"truncating to {n}."
            )
            return kx[:n], np.arange(n, dtype=int)

        length = float(getattr(data, 'length', 0.0))
        if np.isfinite(length) and abs(length) > 1e-12:
            dkx = 2.0 * np.pi / length
        else:
            dkx = 1.0
        p_index = np.arange(n_r, dtype=float) - (n_r // 2)
        return p_index * dkx, radial_idx

















    def _compute_case_time_window(self, data):
        """Compute plotting time-window indices for one case."""
        t_indices = []
        t_start, t_end = 0, 0
        if hasattr(data, 't'):
            t_indices, t_start, t_end = self._get_time_indices(data.t)
        return t_indices, t_start, t_end




    def _dispatch_single_case_plot(
        self,
        data,
        label,
        plot_type,
        t_indices,
        t_start,
        t_end,
        species_override_index=None,
    ):
        """
        Route one case to plotting backend according to `plot_type`.

        Dispatch strategy (ordered):
        1) dedicated simple families (frequency/growth, flux, FFT),
        2) special disambiguation branches (energy-balance-single),
        3) generic fluctuation 1D (`vs ky/kx/time`),
        4) exact-identifier handlers.
        """
        exact_handlers = self._build_exact_plot_handler_map(
            data, label, plot_type, t_indices, t_start, t_end
        )

        if self._dispatch_common_plot_families(
            data, label, plot_type, t_indices, t_start, t_end, species_override_index
        ):
            return

        # Energy-balance single mode also contains "vs Time"/"vs ky" tokens in
        # its plot_type string, so it must be dispatched before the generic
        # fluctuation-1D branch below.
        if plot_type.startswith("Energy Balance Single "):
            exact_handlers["Energy Balance Single"]()
            return

        if ("vs ky" in plot_type or "vs kx" in plot_type or "vs Time" in plot_type) and ("Flux" not in plot_type):
            if any(x in plot_type for x in ["Phi", "Apar", "Bpar"]):
                self._plot_fluctuation_1d(data, label, plot_type, t_indices, t_start, t_end)
            return

        handler = exact_handlers.get(plot_type, None)
        if handler is not None:
            handler()

    def _build_exact_plot_handler_map(self, data, label, plot_type, t_indices, t_start, t_end):
        """
        Build exact `plot_type -> callable` mapping for one case.

        This isolates a large dictionary from the dispatch flow so the
        dispatcher logic remains short and easy to read.
        """
        return {
            "ZF ExB Shearing Rate": lambda: self._plot_zf_exb_shearing_rate_time_trace(
                data, label, t_indices, t_start, t_end
            ),
            "ZF ExB Shearing Spectrum": lambda: self._plot_zf_exb_shearing_rate_kx_spectrum(
                data, label, t_indices, t_start, t_end
            ),
            "ZF Phi Spectrum (theta0)": lambda: self._plot_zf_phi_kx_theta0_spectrum(
                data, label, t_indices, t_start, t_end
            ),
            "ZF ExB vs gamma_lin (kx=ky)": lambda: self._plot_zf_exb_vs_gamma_lin_kx_equals_ky(
                data, label, t_indices, t_start, t_end
            ),
            # Backward compatibility for historical UI label.
            "ZF ExB Fig4 (kx=ky)": lambda: self._plot_zf_exb_vs_gamma_lin_kx_equals_ky(
                data, label, t_indices, t_start, t_end
            ),
            "Energy Balance Entropy": lambda: self._plot_energy_balance_entropy(
                data, label, t_indices, t_start, t_end
            ),
            "Energy Balance ZF": lambda: self._plot_energy_balance_zf(
                data, label, t_indices, t_start, t_end
            ),
            "Energy Balance Gamma Eff": lambda: self._plot_energy_balance_gamma_eff_v3(
                data, label, t_indices, t_start, t_end
            ),
            "Energy Balance Single": lambda: self._plot_energy_balance_single_mode(
                data, label, t_indices, t_start, t_end
            ),
            "Fluctuation 2D": lambda: self._plot_fluctuation_2d(
                data, label, plot_type, t_indices, t_start, t_end
            ),
            "Fluctuation 2D vs xt": lambda: self._plot_xt_fluctuation_contours(
                data, label, plot_type, t_indices, t_start, t_end
            ),
            "Integration Error": lambda: self._plot_other_error(data, label),
            "Radial Correlation (rcorr_phi)": lambda: self._plot_other_rcorr_phi(
                data, label, t_indices, t_start, t_end
            ),
            "POD Parity": lambda: self._plot_other_apar_pod_parity(
                data, label, t_indices, t_start, t_end
            ),
        }

    def _dispatch_common_plot_families(
        self, data, label, plot_type, t_indices, t_start, t_end, species_override_index
    ):
        """
        Handle high-level plot families that are matched by pattern.

        Returns:
        - `True` when `plot_type` has been handled.
        - `False` when caller should continue with exact handler routing.
        """
        if plot_type in ["Frequency", "Growth Rate"]:
            self._plot_frequency_growth(data, label, plot_type, t_indices, t_start, t_end)
            return True

        if "Flux" in plot_type:
            self._plot_flux_diagnostics(
                data, label, plot_type, t_indices, t_start, t_end, species_override_index
            )
            return True

        if "FFT" in plot_type:
            self._plot_fluctuation_fft(data, label, plot_type, t_indices)
            return True

        return False

    def _plot_single_case(self, data, label, plot_type, species_override_index=None):
        """
        Plot one selected case using the current plot-type selection.

        Workflow:
        - Resolve per-case time window once.
        - Delegate to centralized dispatch.
        - Catch/report per-case exceptions so one bad case does not abort others.
        """
        t_indices, t_start, t_end = self._compute_case_time_window(data)

        try:
            self._dispatch_single_case_plot(
                data, label, plot_type, t_indices, t_start, t_end, species_override_index
            )
        except Exception as e:
            print(f"Error processing data for {label}: {e}")
            traceback.print_exc()
    


