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

try:
    from cgyro_comparison_bootstrap import (
        DEFAULT_LINEAR_GAMMA_FILE,
        DEFAULT_EXPORT_DIRNAME,
    )
except ImportError:
    from .cgyro_comparison_bootstrap import (
        DEFAULT_LINEAR_GAMMA_FILE,
        DEFAULT_EXPORT_DIRNAME,
    )

try:
    import scipy.signal as sp_signal
    from scipy.optimize import curve_fit as sp_curve_fit
except Exception:
    sp_signal = None
    sp_curve_fit = None


class CgyroPlottingMixin:
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
        elif plot_type == "ZF ExB Fig4 (kx=ky)":
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
        self.ax = self.fig.add_subplot(111)

    def _plot_all_species_flux_first_case(self, selected_cases, plot_type):
        """
        Plot flux for all species, restricted to the first selected case.

        This path is used when `Plot All Species` is enabled in flux mode.
        """
        if len(selected_cases) == 0:
            return

        case_name = list(selected_cases)[0]
        data = self.cases[case_name]

        # For non-estimated flux plots, loading ky_flux can still help some cases.
        # Estimated kx flux does not require ky_flux.
        if ("vs kx" not in plot_type.lower()) and (not hasattr(data, 'ky_flux')):
            try:
                data.getflux()
            except Exception:
                pass

        specs = self._get_case_species(data)
        if not specs:
            print(f"No species found for {case_name}")
            return

        for i, _spec in enumerate(specs):
            try:
                self._plot_single_case(data, case_name, plot_type, species_override_index=i)
            except Exception as e:
                print(f"Error plotting species {i} for {case_name}: {e}")
                traceback.print_exc()

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
                        axi.legend(fontsize=8, loc='best')
                    axi.grid(True, alpha=0.3)
            self.canvas.draw()
            return

        if not self._is_standard_line_plot(plot_type):
            self.canvas.draw()
            return

        handles, _labels = self.ax.get_legend_handles_labels()
        if handles:
            self.ax.legend()
        self.ax.grid(True)

        title = f"Comparison: {display_plot_type}"
        if self.norm_ky_var.get() and plot_type_selection in ["Frequency", "Growth Rate"]:
            title = f"{title}/ky"
        self.ax.set_title(title)

        self._apply_axis_labels(plot_type)

        if self.log_x_var.get():
            self.ax.set_xscale('log')
        if self.log_y_var.get():
            self.ax.set_yscale('log')

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

    def _plot_frequency_growth(self, data, label, plot_type, t_indices, t_start, t_end):
        """Plot frequency or growth-rate spectra, with optional time-window averaging."""
        if not hasattr(data, 'freq') or not hasattr(data, 'ky'):
            return

        x = self._positive_ky_axis(getattr(data, 'ky', []))
        n_ky = x.size
        if n_ky == 0:
            return

        comp_idx = 0 if plot_type == "Frequency" else 1
        freq = np.asarray(data.freq)
        y = np.array([])
        used_avg = False

        if freq.ndim == 2:
            if freq.shape[0] > comp_idx:
                y = np.asarray(freq[comp_idx]).reshape(-1)
        elif freq.ndim == 3 and freq.shape[0] > comp_idx:
            comp = np.asarray(freq[comp_idx])

            # Detect which axis is ky so we can average along time correctly.
            if comp.shape[0] == n_ky:
                time_axis = 1
            elif comp.shape[1] == n_ky:
                time_axis = 0
            else:
                # Fallback to CGYRO default layout [ky, t].
                time_axis = 1

            valid_t = np.asarray(t_indices, dtype=int)
            if valid_t.size > 0:
                t_len = comp.shape[time_axis]
                valid_t = valid_t[(valid_t >= 0) & (valid_t < t_len)]

            if valid_t.size > 0:
                y = np.mean(np.take(comp, valid_t, axis=time_axis), axis=time_axis)
                used_avg = True
            else:
                y = np.take(comp, -1, axis=time_axis)

        if y.size == 0:
            return

        if used_avg:
            label = self._append_avg_suffix(label, t_start, t_end, prefix="Avg")

        # Normalize by ky if requested
        if self.norm_ky_var.get():
            y = np.asarray(y, dtype=float).reshape(-1)
            x_safe = np.asarray(x, dtype=float).reshape(-1)
            with np.errstate(divide='ignore', invalid='ignore'):
                y = np.divide(
                    y,
                    x_safe,
                    out=np.full_like(y, np.nan),
                    where=np.abs(x_safe) > 1e-12,
                )

        self._plot_1d(x, y, label, plot_type)

    def _plot_other_error(self, data, label):
        """Plot CGYRO integration error traces (`err1`, `err2`) vs time."""
        t = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        err1 = np.asarray(getattr(data, 'err1', []), dtype=float).reshape(-1)
        err2 = np.asarray(getattr(data, 'err2', []), dtype=float).reshape(-1)

        n = min(t.size, err1.size, err2.size)
        if n <= 0:
            print(f"No error data available for {label}")
            return

        t = t[:n]
        err1 = err1[:n]
        err2 = err2[:n]

        # Match cgyro_plot style: skip the first two points.
        i0 = 2 if n > 2 else 0
        t = t[i0:]
        err1 = err1[i0:]
        err2 = err2[i0:]

        # Guard log-scale plotting.
        err1 = np.where(err1 > 0.0, err1, np.nan)
        err2 = np.where(err2 > 0.0, err2, np.nan)

        self.ax.plot(t, err1, label=f"{label} (Total error)")
        self.ax.plot(t, err2, linestyle='--', label=f"{label} (RK error)")
        self.ax.set_xlabel(r'$t\ (a/c_s)$')
        self.ax.set_ylabel(r'$\mathrm{Integration~Error}$')
        self.ax.set_yscale('log')

    def _plot_other_case_info(self, data, label):
        """Load one case `out.cgyro.info` and prepare continuous-scroll rendering."""
        case_dir = self._resolve_case_dir(data)
        if not case_dir:
            self.ax.axis('off')
            self.ax.text(
                0.01,
                0.99,
                f"Cannot resolve case directory for: {label}",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            self._case_info_scroll_active = False
            return

        info_path = os.path.join(case_dir, "out.cgyro.info")
        if not os.path.isfile(info_path):
            self.ax.axis('off')
            self.ax.text(
                0.01,
                0.99,
                f"out.cgyro.info not found:\n{info_path}",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            self._case_info_scroll_active = False
            return

        try:
            with open(info_path, "r", encoding="utf-8", errors="replace") as fh:
                content = fh.read()
        except Exception as e:
            self.ax.axis('off')
            self.ax.text(
                0.01,
                0.99,
                f"Failed to read out.cgyro.info:\n{info_path}\n\n{e}",
                transform=self.ax.transAxes,
                va='top',
                ha='left',
                fontsize=10,
                family='monospace',
            )
            self._case_info_scroll_active = False
            return

        header_lines = [
            f"Case: {label}",
            f"File: {info_path}",
            "-" * 110,
            "",
        ]
        # Use split('\n') to preserve a possible trailing empty line at EOF.
        body_lines = content.replace('\r\n', '\n').replace('\r', '\n').split('\n')

        self._case_info_case_name = str(label)
        self._case_info_file_path = str(info_path)
        self._case_info_lines = header_lines + body_lines
        self._case_info_line_colors = ["#0d47a1"] * len(header_lines) + ["#222222"] * len(body_lines)

        self._enable_case_info_scroll()

    def _plot_other_rcorr_phi(self, data, label, t_indices, t_start, t_end):
        """Plot radial auto-correlation for selected field and theta index."""
        if not hasattr(data, 'kxky_select'):
            print(f"rcorr_phi is not supported for this case object: {label}")
            return

        try:
            theta_idx = int(str(self.others_rcorr_theta_var.get()).strip())
        except Exception:
            theta_idx = -1

        field_name = str(self.others_rcorr_field_var.get()).strip().lower()
        field_idx = 0
        if field_name == "apar":
            field_idx = 1
        elif field_name == "bpar":
            field_idx = 2

        if not hasattr(data, 'kxky_phi'):
            if not self._ensure_bigfield_loaded(data, label):
                return

        try:
            # f shape: [n_radial-1, n_ky, n_t]
            f, ft = data.kxky_select(theta_idx, field_idx, 'phi', 0)
            f = np.asarray(f)
        except Exception as e:
            print(f"Failed to evaluate rcorr_phi source for {label}: {e}")
            return

        if f.ndim != 3:
            print(f"Unexpected rcorr_phi source shape for {label}: {f.shape}")
            return

        t = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        nt = min(int(f.shape[2]), t.size if t.size > 0 else int(f.shape[2]))
        if nt <= 0:
            print(f"No time data for rcorr_phi in {label}")
            return
        f = f[:, :, :nt]

        valid_t = np.asarray(t_indices, dtype=int)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < nt)]
        if valid_t.size == 0:
            valid_t = np.arange(nt, dtype=int)
            win_tag = "full"
        else:
            win_tag = f"{t_start:.1f}-{t_end:.1f}"

        n_r = int(getattr(data, 'n_radial', f.shape[0] + 1))
        if n_r <= 1:
            print(f"Insufficient radial points for rcorr_phi in {label}")
            return

        y = np.zeros((n_r, nt), dtype=float)
        n_use = min(f.shape[0], n_r - 1)
        if f.shape[1] > 1:
            y[1:1 + n_use, :] = np.sum(np.abs(f[:n_use, 1:, :]), axis=1)
        else:
            # Single-mode fallback
            y[1:1 + n_use, :] = np.abs(f[:n_use, 0, :])

        ave = np.mean(y[:, valid_t], axis=1)
        nx = int(ave.size)
        if nx <= 1:
            print(f"Invalid rcorr_phi average size for {label}")
            return

        ave = np.roll(ave, -nx // 2)
        ave[0] = 0.0
        corr = np.fft.fftshift(np.fft.fft(ave, nx))
        cmax = float(np.max(np.abs(corr)))
        if cmax > 0.0:
            corr = corr / cmax
        corr = np.real(corr)

        delta_r = np.fft.fftshift(np.fft.fftfreq(nx))
        kx = np.asarray(getattr(data, 'kx', []), dtype=float).reshape(-1)
        if kx.size > 1:
            dk = float(kx[1] - kx[0])
        else:
            length = float(getattr(data, 'length', 0.0))
            dk = (2.0 * np.pi / length) if np.isfinite(length) and length > 0 else 1.0
        if abs(dk) > 1e-12:
            delta_r = delta_r * (2.0 * np.pi / dk)

        plot_label = f"{label} ({field_name}, t={win_tag})"
        self.ax.plot(delta_r, corr, label=plot_label)
        self.ax.plot(delta_r, 0.0 * delta_r, color='k', linestyle='--', linewidth=0.8)

        # Optional exponential-envelope fit, as in cgyro_plot.
        if sp_signal is not None and sp_curve_fit is not None:
            try:
                corr_env = np.abs(sp_signal.hilbert(corr))

                def _absexp(x, tau):
                    """Absolute-value exponential envelope model for correlation fitting."""
                    return np.exp(-np.abs(x) / tau)

                p0 = [max(1.0, 0.1 * (np.max(delta_r) - np.min(delta_r)))]
                l_corr, _pcov = sp_curve_fit(_absexp, delta_r, corr_env, p0=p0, maxfev=10000)
                tau = float(l_corr[0])
                if np.isfinite(tau) and tau > 0.0:
                    self.ax.plot(
                        delta_r,
                        _absexp(delta_r, tau),
                        linestyle='-.',
                        linewidth=1.0,
                        label=f"{label} fit (l_corr={tau:.3g}, t={win_tag})",
                    )
                    print(f"INFO: (rcorr_phi) {label} l_corr = {tau:.3f}")
            except Exception as e:
                print(f"Warning: rcorr_phi fit failed for {label}: {e}")

        ft_txt = str(ft).strip() if ft is not None else r"\phi"
        if not ft_txt:
            ft_txt = r"\phi"
        self.ax.set_xlabel(r'$\Delta r/\rho_s$')
        self.ax.set_ylabel(r'$C_{' + ft_txt + r'}(\Delta r)$')
        self.ax.set_ylim(-1.0, 1.0)

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

    def _plot_other_apar_pod_parity(self, data, label, t_indices, t_start, t_end, field_override=None):
        """
        POD decomposition for selected field (A// or Phi) at selected ky,kx.

        Procedure (following the requested workflow):
        1) Select `(kx, ky)` and optionally include all connected kx modes.
        2) Build POD (SVD) on complex field data over selected time window.
        3) Compute parity factor P for first two POD modes.
        4) Assign ballooning/tearing tags from parity factor P.
        5) Remaining POD modes are grouped as residual.
        """
        field_name = (field_override or "").strip().lower()
        if not field_name:
            try:
                field_name = str(self.others_pod_field_var.get()).strip().lower()
            except Exception:
                field_name = "apar"
        is_phi_mode = (field_name == "phi")

        if is_phi_mode:
            field_attr = 'kxky_phi'
            field_suffix = '.cgyro.kxky_phi'
            field_title = "Phi"
            field_tex = r'\phi'
            decomp_title = "Phi POD parity decomposition"
        else:
            field_attr = 'kxky_apar'
            field_suffix = '.cgyro.kxky_apar'
            field_title = "A//"
            field_tex = r'A_{||}'
            decomp_title = "A// POD parity decomposition"

        pod_field = self._load_kxky_complex(
            data, label, field_attr, field_suffix, species_dependent=False
        )
        if pod_field is None:
            print(f"{field_title} POD parity unavailable for {label}: no {field_title} data.")
            return
        if pod_field.ndim != 4:
            print(f"Unsupported {field_title} shape for {label}: {pod_field.shape} (expected [nr,theta,ky,t]).")
            return

        n_r, n_theta, n_ky, n_t = pod_field.shape
        if n_r <= 0 or n_theta <= 0 or n_ky <= 0 or n_t <= 1:
            print(f"Insufficient {field_title} dimensions for POD in {label}: {pod_field.shape}")
            return
        if n_theta <= 1:
            self._show_pod_theta_resolution_warning()
            return

        ky_axis = self._positive_ky_axis(getattr(data, 'ky', []))
        if ky_axis.size != n_ky:
            ky_axis = np.arange(n_ky, dtype=float)

        kx_axis, radial_idx = self._build_kx_axis(data, n_r, label)
        if kx_axis is None or radial_idx is None:
            print(f"Could not construct kx axis for {label}.")
            return

        ky_sel_idx = self._resolve_axis_selection_index(
            ky_axis, self.others_pod_ky_var.get(), "ky", label
        )
        kx_sel_axis_idx = self._resolve_axis_selection_index(
            kx_axis, self.others_pod_kx_var.get(), "kx", label, prefer_value=True
        )

        kx_sel_field_idx = int(radial_idx[kx_sel_axis_idx])
        # Always use connected kx set for POD parity.
        use_connected = True
        chain_field_idx = np.asarray(radial_idx, dtype=int)

        # Map field-index -> physical kx for proxy-axis construction.
        field_to_kx = {}
        for i_ax, fidx in enumerate(np.asarray(radial_idx, dtype=int)):
            if i_ax < len(kx_axis):
                field_to_kx[int(fidx)] = float(kx_axis[i_ax])
        chain_kx = np.asarray(
            [field_to_kx.get(int(fidx), np.nan) for fidx in np.asarray(chain_field_idx, dtype=int)],
            dtype=float,
        )

        valid_t = np.asarray(t_indices, dtype=int)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size == 0:
            valid_t = np.arange(n_t, dtype=int)

        block = np.asarray(pod_field[chain_field_idx, :, ky_sel_idx, :][:, :, valid_t], dtype=complex)
        n_chain = int(block.shape[0])
        n_t_sel = int(block.shape[2])
        if n_chain <= 0 or n_t_sel <= 1:
            print(f"Not enough samples for POD in {label}.")
            return

        # Critical fix: enforce fftshift on connected-kx chain before POD/SVD
        # so the kx dimension is centered consistently (negative -> 0 -> positive).
        if n_chain > 1:
            block = np.fft.fftshift(block, axes=0)
            chain_field_idx = np.fft.fftshift(np.asarray(chain_field_idx, dtype=int))
            if chain_kx.size == n_chain:
                chain_kx = np.fft.fftshift(chain_kx)
            print(f"INFO: applied fftshift on connected-kx chain before POD for {label}.")

        # Keep physical order monotonic when kx metadata is available.
        if chain_kx.size == n_chain and np.all(np.isfinite(chain_kx)) and n_chain > 1:
            kx_order = np.argsort(chain_kx)
            block = block[kx_order, :, :]
            chain_field_idx = chain_field_idx[kx_order]
            chain_kx = chain_kx[kx_order]

        stitch_desc = "raw"
        if n_chain > 1:
            ky_sel_val = float(ky_axis[ky_sel_idx]) if ky_sel_idx < ky_axis.size else float(ky_sel_idx)
            block, stitch_desc = self._repair_connected_kx_stitching(
                data, label, block, chain_field_idx, chain_kx, ky_sel_val
            )

        # Build POD matrix: spatial dof = (connected kx, theta), temporal dof = time.
        x_mat = block.reshape(n_chain * n_theta, n_t_sel)
        x_mat = x_mat - np.mean(x_mat, axis=1, keepdims=True)

        if min(x_mat.shape) < 2:
            print(f"POD requires at least 2 spatial/time samples for {label}.")
            return

        try:
            u, s, vh = np.linalg.svd(x_mat, full_matrices=False)
        except Exception as e:
            print(f"POD SVD failed for {label}: {e}")
            return

        n_modes = int(min(2, u.shape[1]))
        if n_modes <= 0:
            print(f"No POD modes available for {label}.")
            return

        parity_axis_base = self._build_theta_over_pi_axis(data, n_theta)
        p_axis = np.asarray(getattr(data, 'p', []), dtype=float).reshape(-1)
        parity_axis = self._build_extended_z_over_pi_axis(
            parity_axis_base, chain_field_idx, chain_kx=chain_kx, p_axis=p_axis
        )
        if parity_axis.size != (n_chain * n_theta):
            # Safety fallback: build index-centered extended axis.
            parity_axis = self._build_extended_z_over_pi_axis(
                parity_axis_base, chain_field_idx, chain_kx=chain_kx, p_axis=None
            )
        parity_axis_label = r'$z/\pi$'
        zwin = float(DEFAULT_POD_Z_WINDOW_PI)
        use_window_xlim = False
        try:
            if len(parity_axis) > 0 and np.isfinite(np.nanmax(np.abs(parity_axis))):
                use_window_xlim = (np.nanmax(np.abs(parity_axis)) > zwin + 1e-12)
        except Exception:
            use_window_xlim = False

        mode_profiles = []
        mode_axes = []
        p_even_vals = []
        p_odd_vals = []
        for i_mode in range(n_modes):
            u_mode = u[:, i_mode].reshape(n_chain, n_theta)
            profile = np.asarray(u_mode.reshape(-1), dtype=complex)
            z_plot_i, prof_plot_i = self._prepare_pod_display_curve(parity_axis, profile, zwin)
            mode_axes.append(z_plot_i)
            mode_profiles.append(prof_plot_i)
            p_even, p_odd = self._parity_even_odd_ratios(prof_plot_i, z_plot_i)
            p_even_vals.append(p_even)
            p_odd_vals.append(p_odd)

        # Mirror parity thresholds on first two POD modes.
        # P_even > 0.7 -> tearing ; P_odd > 0.7 -> ballooning ; otherwise mixed.
        role_by_mode = {}
        ball_idx = None
        tear_idx = None
        mixed_indices = []
        for i_mode in range(min(2, n_modes)):
            p_even = p_even_vals[i_mode] if i_mode < len(p_even_vals) else np.nan
            p_odd = p_odd_vals[i_mode] if i_mode < len(p_odd_vals) else np.nan
            if np.isfinite(p_even) and np.isfinite(p_odd):
                if p_even > 0.7:
                    role_by_mode[i_mode] = "tearing"
                    if tear_idx is None:
                        tear_idx = i_mode
                elif p_odd > 0.7:
                    role_by_mode[i_mode] = "ballooning"
                    if ball_idx is None:
                        ball_idx = i_mode
                else:
                    role_by_mode[i_mode] = "mixed"
                    mixed_indices.append(i_mode)
            else:
                role_by_mode[i_mode] = "mixed"
                mixed_indices.append(i_mode)

        if n_modes == 1 and 0 not in role_by_mode:
            role_by_mode[0] = "mixed"
            mixed_indices.append(0)

        # Energy fractions for decomposition A = A_ball + A_tear + A_res.
        s2 = np.asarray(s, dtype=float) ** 2
        e_total = float(np.sum(s2)) if s2.size > 0 else 0.0
        e_ball = float(s2[ball_idx]) if (ball_idx is not None and ball_idx < s2.size) else 0.0
        e_tear = float(s2[tear_idx]) if (tear_idx is not None and tear_idx < s2.size) else 0.0
        e_mix = 0.0
        for i_mix in mixed_indices:
            if 0 <= i_mix < s2.size:
                e_mix += float(s2[i_mix])
        res_start = 2 if s2.size >= 2 else 1
        e_res = float(np.sum(s2[res_start:])) if s2.size > res_start else 0.0
        if e_total > 1.0e-14:
            f_ball = 100.0 * e_ball / e_total
            f_tear = 100.0 * e_tear / e_total
            f_mix = 100.0 * e_mix / e_total
            f_res = 100.0 * e_res / e_total
        else:
            f_ball = f_tear = f_mix = f_res = 0.0

        self.fig.clear()
        axes = self.fig.subplots(1, 2, squeeze=False)[0]
        self.ax = axes[0]

        for i_panel in range(2):
            axp = axes[i_panel]
            if i_panel >= n_modes:
                axp.axis("off")
                axp.text(0.5, 0.5, "No mode", ha="center", va="center", transform=axp.transAxes)
                continue

            prof = mode_profiles[i_panel]
            z_plot = mode_axes[i_panel]
            if len(z_plot) <= 1 or len(prof) <= 1:
                axp.axis("off")
                axp.text(0.5, 0.5, "No valid axis", ha="center", va="center", transform=axp.transAxes)
                continue
            scale = float(np.nanmax(np.abs(prof))) if prof.size > 0 else 0.0
            if np.isfinite(scale) and scale > 1.0e-14:
                prof_plot = prof / scale
            else:
                prof_plot = prof
            axp.plot(z_plot, np.real(prof_plot), 'k-', linewidth=1.8, label=rf'n={i_panel+1} Re${field_tex}$')
            axp.plot(z_plot, np.imag(prof_plot), 'r--', linewidth=1.8, label=rf'n={i_panel+1} Im${field_tex}$')
            if use_window_xlim:
                axp.set_xlim([-zwin, zwin])

            p_even_txt = "nan"
            p_odd_txt = "nan"
            if i_panel < len(p_even_vals) and np.isfinite(p_even_vals[i_panel]):
                p_even_txt = f"{p_even_vals[i_panel]:.3f}"
            if i_panel < len(p_odd_vals) and np.isfinite(p_odd_vals[i_panel]):
                p_odd_txt = f"{p_odd_vals[i_panel]:.3f}"
            role = role_by_mode.get(i_panel, "mode")
            axp.set_title(f"n={i_panel+1} ({role}, Pe={p_even_txt}, Po={p_odd_txt})")
            axp.set_xlabel(parity_axis_label)
            axp.set_ylabel(r'$' + field_tex + r'$ (arb.)')
            axp.grid(True, alpha=0.25)
            axp.legend(loc='best', fontsize=10, frameon=False)

        ky_val = float(ky_axis[ky_sel_idx]) if ky_sel_idx < ky_axis.size else float(ky_sel_idx)
        kx_val = float(kx_axis[kx_sel_axis_idx]) if kx_sel_axis_idx < kx_axis.size else float(kx_sel_axis_idx)
        t_arr = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        if t_arr.size > 0 and valid_t.size > 0 and np.max(valid_t) < t_arr.size:
            t_info = f"{float(t_arr[valid_t[0]]):.2f}-{float(t_arr[valid_t[-1]]):.2f}"
        else:
            t_info = f"{t_start:.2f}-{t_end:.2f}" if t_end > t_start else "full"
        conn_txt = "connected-kx" if use_connected else "single-kx"

        self.fig.suptitle(
            (
                f"{decomp_title}: {label} | "
                f"ky={ky_val:.4g}, kx={kx_val:.4g}, t={t_info}, {conn_txt}\n"
                f"Energy split: ball={f_ball:.1f}%  tear={f_tear:.1f}%  mixed={f_mix:.1f}%  res={f_res:.1f}%"
                f" | stitch={stitch_desc}"
            ),
            fontsize=11,
        )
        try:
            self.fig.tight_layout(rect=(0.01, 0.03, 0.99, 0.90))
        except Exception:
            pass

    def _plot_other_phi_pod_parity(self, data, label, t_indices, t_start, t_end):
        """POD parity decomposition for electrostatic potential `phi`."""
        self._plot_other_apar_pod_parity(
            data, label, t_indices, t_start, t_end, field_override="phi"
        )

    def _plot_fluctuation_2d(self, data, label, plot_type, t_indices, t_start, t_end):
        """Render ky-kx fluctuation contour animation or snapshot from bigfield data."""
        if not hasattr(data, 'kxky_phi'):
            print(f"Loading big field data for {label}...")
            if not self._ensure_bigfield_loaded(data, label):
                return

        use_x_rhoe = self._use_fluc2d_x_rhoe()
        x_unit_factor = self._get_rhos_to_rhoe_factor(data, label) if use_x_rhoe else 1.0
        x_unit_label = r'$x/\rho_e$' if use_x_rhoe else r'$x/\rho_s$'
        y_unit_label = r'$y/\rho_e$' if use_x_rhoe else r'$y/\rho_s$'

        moment = self.moment_var.get()
        n_species = int(getattr(data, 'n_species', 0))
        if n_species <= 0:
            n_species = max(1, len(self._get_case_species(data)))

        if moment in ["Density", "Energy", "Temperature"]:
            # Keep compatibility with previous behavior:
            # Main Ion selection is reduced to one ion for this 2D view path.
            target_indices, _spec_label = self._resolve_species_indices(
                data,
                n_species,
                label,
                fallback_first=True,
                main_ion_policy="first",
                single_species_only=(moment == "Temperature"),
            )
            if not target_indices:
                return
        else:
            # For field moments (Phi, Apar, Bpar), species doesn't matter.
            target_indices = [0]
        
        # Select data based on moment
        c_field = None
        if moment == "Phi":
            if hasattr(data, 'kxky_phi'):
                    c_field = data.kxky_phi
        elif moment == "Apar":
            if hasattr(data, 'kxky_apar'):
                    c_field = data.kxky_apar
        elif moment == "Bpar":
            if hasattr(data, 'kxky_bpar'):
                    c_field = data.kxky_bpar
        elif moment == "Density":
            if hasattr(data, 'kxky_n'):
                    c_field = data.kxky_n
        elif moment == "Energy":
            if hasattr(data, 'kxky_e'):
                    c_field = data.kxky_e
        elif moment == "Temperature":
            if hasattr(data, 'kxky_n') and hasattr(data, 'kxky_e'):
                c_field = (data.kxky_n, data.kxky_e)
            else:
                print(f"Data for Temperature (n and e) not available in {label}")
                return

        if c_field is None and moment != "Temperature":
            print(f"Data for {moment} not available in {label}")
            return

        # Theta index: midplane
        itheta = data.theta_plot * 4 // 8
        if itheta >= data.theta_plot: itheta = 0

        # Define calculation function
        def get_frame_data(ti):
            """Reconstruct one real-space frame for contour plotting at time index `ti`."""
            if moment == "Temperature":
                # Only support single species for Temperature for now
                species_idx = target_indices[0]
                
                # c_field is tuple (kxky_n, kxky_e)
                n_data = c_field[0]
                e_data = c_field[1]
                
                # Reconstruct complex arrays
                # data shape is [2, nr, theta, species, nn, nt]
                # we want [2, nr, nn]
                ndim = data.kxky_n.ndim
                shape = data.kxky_n.shape
                # print(f"DEBUG: kxky_phi shape: {shape}")

                if ndim == 6:
                    # Standard: [2, nr, theta, species, nn, nt]
                    n_real = n_data[0, :, itheta, species_idx, :, ti]
                    n_imag = n_data[1, :, itheta, species_idx, :, ti]
                    e_real = e_data[0, :, itheta, species_idx, :, ti]
                    e_imag = e_data[1, :, itheta, species_idx, :, ti]
                
                    cn = n_real + 1j * n_imag
                    ce = e_real + 1j * e_imag
        
                elif ndim == 5:
                    # Standard: [nr, theta, species, nn, nt]
                    cn = n_data[:, itheta, species_idx, :, ti]
                    ce = e_data[:, itheta, species_idx, :, ti]
                else:
                    print(f"Error: Unsupported kxky_phi dimensions: {ndim} {shape}")
                    return

                t0 = data.temp[species_idx]
                n0 = data.dens[species_idx]
                c = ((2.0/3.0) * ce - t0 * cn) / n0
                
            elif moment in ["Density", "Energy"]:
                    # [2, nr, theta, species, nn, nt]
                    ndim = c_field.ndim
                    shape = c_field.shape
                    
                    c = 0
                    
                    for s_idx in target_indices:
                        if ndim == 6:
                        # Standard: [2, nr, theta, nn, nt]
                            c_real = c_field[0, :, itheta, s_idx, :, ti]
                            c_imag = c_field[1, :, itheta, s_idx, :, ti]
                            c += c_real + 1j * c_imag
            
                        elif ndim == 5:
                        # Standard: [nr, theta, species, nn, nt]
                            c += c_field[:, itheta, s_idx, :, ti]
                        else:
                            print(f"Error: Unsupported kxky_phi dimensions: {ndim} {shape}")
                            return

            else:
                    # Fields are species independent
                    ndim = c_field.ndim
                    shape = c_field.shape

                    if ndim == 5:
                        c_real = c_field[0, :, itheta, :, ti]
                        c_imag = c_field[1, :, itheta, :, ti]
                        c = c_real + 1j * c_imag
                    elif ndim == 4:
                    # Standard: [nr, theta, nn, nt]
                        c = c_field[:, itheta, :, ti]
                    else:
                        print(f"Error: Unsupported kxky_phi dimensions: {ndim} {shape}")
                        return
            
            # FFT params
            nr = data.n_radial
            nn = data.n_n
            nx = nr + 1
            ny = 2 * nn - 1
            
            f = self._maptoreal_fft(nr, nn, nx, ny, c)
            f = 2 * f # Correct for half-sum

            # Grid
            x = np.arange(nx) * 2 * np.pi / nx
            y = np.arange(ny) * 2 * np.pi / ny
            
            xmax = data.length
            ky_axis = self._positive_ky_axis(getattr(data, 'ky', []))
            if ky_axis.size > 1 and abs(ky_axis[1]) > 1e-9:
                    ymax = (2 * np.pi) / np.abs(ky_axis[1])
            else:
                    ymax = 100.0 # Fallback
            
            xp = x / (2 * np.pi) * xmax
            yp = y / (2 * np.pi) * ymax
            if use_x_rhoe:
                xp = xp * x_unit_factor
                xmax = xmax * x_unit_factor
                yp = yp * x_unit_factor
                ymax = ymax * x_unit_factor
            
            # Periodic extensions
            xp = np.append(xp, xmax)
            yp = np.append(yp, ymax)
            fp = np.zeros([nx + 1, ny + 1])
            fp[0:nx, 0:ny] = f[:, :]
            fp[-1, :] = fp[0, :]
            fp[:, -1] = fp[:, 0]
            
            return xp, yp, fp, xmax, ymax

        # Determine if animation is needed
        # Only animate if user explicitly provided a time range (start or end)
        user_specified_time = self.t_start_var.get().strip() or self.t_end_var.get().strip()
        should_animate = user_specified_time and len(t_indices) > 1

        if should_animate:
            # Animation Logic
            self.total_frames = len(t_indices)
            self.current_frame = 0

            def update(frame):
                """Animation callback: render the current frame index."""
                # Handle loop or bounds if called manually
                frame = frame % self.total_frames
                self.current_frame = frame
                
                ti = t_indices[frame]
                frame_data = get_frame_data(ti)
                if frame_data is None:
                    return
                xp, yp, fp, xmax, ymax = frame_data
                
                self.ax.clear()
                levels = 50
                self.ax.contourf(xp, yp, np.transpose(fp), levels, cmap='RdBu_r')
                # Note: Colorbar in animation might be tricky if range changes.
                # Ideally set vmin/vmax fixed based on global min/max, but here we let it scale.
                if hasattr(data, 't') and len(data.t) > 0:
                    t_text = f"{data.t[ti]:.2f}"
                else:
                    t_text = "N/A"
                
                self.ax.set_xlabel(x_unit_label)
                self.ax.set_ylabel(y_unit_label)
                self.ax.set_title(f'{moment} Fluctuation: {label} (t={t_text})')
                self.ax.set_aspect('equal')
                self.ax.set_xlim([0, xmax])
                self.ax.set_ylim([0, ymax])

            self.anim_update_func = update
            self.ani = animation.FuncAnimation(self.fig, update, frames=self.total_frames, interval=100, repeat=True)
            
            # Enable controls
            self.btn_pause.config(state="normal", text="Pause")
            self.btn_prev.config(state="normal")
            self.btn_next.config(state="normal")
            
            self.canvas.draw()
            
        else:
            # Static Plot (Last time point or end of range)
            ti = -1
            if len(t_indices) > 0:
                ti = t_indices[-1]

            frame_data = get_frame_data(ti)
            if frame_data is None:
                return
            xp, yp, fp, xmax, ymax = frame_data
            
            levels = 50
            cs = self.ax.contourf(xp, yp, np.transpose(fp), levels, cmap='RdBu_r')
            self.fig.colorbar(cs, ax=self.ax)
            if hasattr(data, 't') and len(data.t) > 0:
                t_text = f"{data.t[ti]:.2f}"
            else:
                t_text = "N/A"
            
            self.ax.set_xlabel(x_unit_label)
            self.ax.set_ylabel(y_unit_label)
            self.ax.set_title(f'{moment} Fluctuation: {label} (t={t_text})')
            self.ax.set_aspect('equal')
            self.ax.set_xlim([0, xmax])
            self.ax.set_ylim([0, ymax])

    def _plot_xt_fluctuation_contours(self, data, label, plot_type, t_indices, t_start, t_end):
        """Render x-t fluctuation contours reconstructed from spectral data."""
        use_x_rhoe = self._use_fluc2d_x_rhoe()
        x_unit_factor = self._get_rhos_to_rhoe_factor(data, label) if use_x_rhoe else 1.0
        x_ylabel = r'$x/\rho_e$' if use_x_rhoe else r'$x/\rho_s$'

        moment = self.moment_var.get()
        spec_label = ""

        field = None
        if moment == "Phi":
            field = self._load_kxky_complex(data, label, 'kxky_phi', '.cgyro.kxky_phi', species_dependent=False)
        elif moment == "Apar":
            field = self._load_kxky_complex(data, label, 'kxky_apar', '.cgyro.kxky_apar', species_dependent=False)
        elif moment == "Bpar":
            field = self._load_kxky_complex(data, label, 'kxky_bpar', '.cgyro.kxky_bpar', species_dependent=False)
        elif moment in ["Density", "Energy"]:
            attr = 'kxky_n' if moment == "Density" else 'kxky_e'
            suffix = '.cgyro.kxky_n' if moment == "Density" else '.cgyro.kxky_e'
            species_field = self._load_kxky_complex(data, label, attr, suffix, species_dependent=True)
            if species_field is None:
                print(f"No {moment} data available for {label}")
                return
            n_species = int(species_field.shape[2]) if species_field.ndim >= 5 else 1
            target_indices, spec_label = self._resolve_species_indices(
                data,
                n_species,
                label,
                fallback_first=True,
                main_ion_policy="all",
                single_species_only=False,
            )
            if not target_indices:
                return
            field = np.sum(species_field[:, :, target_indices, :, :], axis=2)
        elif moment == "Temperature":
            mom_n = self._load_kxky_complex(data, label, 'kxky_n', '.cgyro.kxky_n', species_dependent=True)
            mom_e = self._load_kxky_complex(data, label, 'kxky_e', '.cgyro.kxky_e', species_dependent=True)
            if mom_n is None or mom_e is None:
                print(f"Data for Temperature (n and e) not available in {label}")
                return

            n_species = min(int(mom_n.shape[2]), int(mom_e.shape[2]))
            target_indices, spec_label = self._resolve_species_indices(
                data,
                n_species,
                label,
                fallback_first=True,
                main_ion_policy="all",
                single_species_only=True,
            )
            if not target_indices:
                return
            s_idx = int(target_indices[0])

            temp = np.asarray(getattr(data, 'temp', []), dtype=float).reshape(-1)
            dens = np.asarray(getattr(data, 'dens', []), dtype=float).reshape(-1)
            t0 = float(temp[s_idx]) if s_idx < temp.size else 1.0
            n0 = float(dens[s_idx]) if s_idx < dens.size else 1.0
            if abs(n0) < 1e-12:
                n0 = 1.0

            field = ((2.0 / 3.0) * mom_e[:, :, s_idx, :, :] - t0 * mom_n[:, :, s_idx, :, :]) / n0
        else:
            print(f"Unsupported moment for x-t contour: {moment}")
            return

        if field is None:
            print(f"No data for {moment} in {label}")
            return

        if field.ndim != 4:
            print(f"Unsupported {moment} shape for x-t contour in {label}: {field.shape}")
            return

        n_r, n_th, n_ky, n_t = field.shape
        if n_r <= 0 or n_th <= 0 or n_ky <= 0 or n_t <= 0:
            print(f"Invalid dimensions for x-t contour in {label}.")
            return

        i_theta = int(getattr(data, 'theta_plot', n_th) * 4 // 8)
        if i_theta < 0 or i_theta >= n_th:
            i_theta = max(0, min(n_th - 1, n_th // 2))

        ky_axis = self._positive_ky_axis(getattr(data, 'ky', []))
        if ky_axis.size > 0:
            n_ky = min(n_ky, ky_axis.size)
            ky_use = ky_axis[:n_ky]
        else:
            ky_use = np.arange(n_ky, dtype=float)

        if n_ky <= 0:
            print(f"No ky axis available for x-t contour in {label}.")
            return

        ky_idx_0 = int(np.argmin(np.abs(ky_use)))
        ky0 = float(ky_use[ky_idx_0])
        if abs(ky0) > 1e-6:
            print(
                f"Warning: ky=0 not found for {label}; using closest "
                f"ky={ky0:.4g} for x-t contour."
            )

        c_kx_t = np.asarray(field[:, i_theta, ky_idx_0, :], dtype=complex)

        t_axis = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        if t_axis.size <= 0:
            t_axis = np.arange(c_kx_t.shape[1], dtype=float)

        n_t = min(c_kx_t.shape[1], t_axis.size)
        if n_t < 2:
            print(f"Not enough time points for x-t contour in {label}.")
            return

        c_kx_t = c_kx_t[:, :n_t]
        t_axis = t_axis[:n_t]

        user_specified_time = bool(self.t_start_var.get().strip() or self.t_end_var.get().strip())
        if user_specified_time:
            valid_t = np.asarray(t_indices, dtype=int)
            valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        else:
            # Default behavior for x-t contour: use the full available time history.
            valid_t = np.arange(n_t, dtype=int)

        if valid_t.size >= 2:
            c_kx_t = c_kx_t[:, valid_t]
            t_plot = t_axis[valid_t]
            if user_specified_time:
                time_tag = f"{t_start:.1f}-{t_end:.1f}"
            else:
                time_tag = f"{t_plot[0]:.1f}-{t_plot[-1]:.1f}"
        else:
            t_plot = t_axis
            time_tag = f"{t_axis[0]:.1f}-{t_axis[-1]:.1f}"

        x_phase, f_tx = self._reconstruct_x_from_kx(c_kx_t, nx=n_r + 1)
        if x_phase is None or f_tx is None:
            print(f"Failed to reconstruct x profiles for {label}.")
            return

        length = float(getattr(data, 'length', 0.0))
        if not np.isfinite(length) or abs(length) < 1e-12:
            length = 2.0 * np.pi
        x_plot = x_phase / (2.0 * np.pi) * length
        if use_x_rhoe:
            x_plot = x_plot * x_unit_factor

        f_abs = np.nanmax(np.abs(f_tx))
        if np.isfinite(f_abs) and f_abs > 0:
            levels = np.linspace(-f_abs, f_abs, 80)
        else:
            levels = 50

        cs = self.ax.contourf(t_plot, x_plot, np.transpose(f_tx), levels=levels, cmap='RdBu_r')
        self.fig.colorbar(cs, ax=self.ax)
        self.ax.set_xlabel(r'$t\ (a/c_s)$')
        self.ax.set_ylabel(x_ylabel)

        moment_tag = moment
        if spec_label:
            moment_tag = f"{moment} ({spec_label})"
        self.ax.set_title(
            f'{moment_tag} x-t Contour: {label} '
            f'(ky={ky0:.4g}, t={time_tag})'
        )

    def _plot_fluctuation_fft(self, data, label, plot_type, t_indices):
        """Render FFT spectra from fluctuation time traces in selected FFT mode."""
        # Determine field name
        field_name = plot_type.split()[0] # "Phi", "Apar", "Bpar"

        # Check/Load data
        c_field = None
        if field_name == "Phi":
            if hasattr(data, 'kxky_phi'): c_field = data.kxky_phi
        elif field_name == "Apar":
            if hasattr(data, 'kxky_apar'): c_field = data.kxky_apar
        elif field_name == "Bpar":
            if hasattr(data, 'kxky_bpar'): c_field = data.kxky_bpar
            
        if c_field is None:
             print(f"Loading big field data for {label}...")
             if self._ensure_bigfield_loaded(data, label):
                 if field_name == "Phi" and hasattr(data, 'kxky_phi'): c_field = data.kxky_phi
                 elif field_name == "Apar" and hasattr(data, 'kxky_apar'): c_field = data.kxky_apar
                 elif field_name == "Bpar" and hasattr(data, 'kxky_bpar'): c_field = data.kxky_bpar
             else:
                 return

        if c_field is None:
            print(f"No {field_name} data available for {label}")
            return

        # kxky_field: [2, n_radial, theta_plot, n_n, nt] or similar
        ndim = c_field.ndim
        
        # Helper to extract complex field at midplane
        i_theta = data.theta_plot * 4 // 8
        if i_theta >= data.theta_plot: i_theta = 0
        
        field_t_all = None
        
        if ndim == 5:
            # Standard: [2, n_radial, theta_plot, n_n, nt]
            field_t_all = c_field[0, :, i_theta, :, :] + 1j * c_field[1, :, i_theta, :, :]
        elif ndim == 4:
             # [n_radial, theta_plot, n_n, nt]
             field_t_all = c_field[:, i_theta, :, :]
        else:
             print(f"Error: Unsupported {field_name} dimensions: {ndim}")
             return
        # field_t_all shape: [n_radial, n_n, nt]
        
        # Slice time if needed
        if len(t_indices) > 0:
            field_t_all = field_t_all[..., t_indices]
            t_array = data.t[t_indices]
        else:
            t_array = data.t
        
        nt_slice = len(t_array)
        if nt_slice < 2:
            print("Not enough time points for FFT")
            return

        # Perform FFT along time axis
        # omega = 2*pi*f
        dt_arr = np.diff(np.asarray(t_array, dtype=float))
        valid_dt = dt_arr[np.abs(dt_arr) > 0]
        if valid_dt.size == 0:
            print("Invalid or zero time spacing for FFT")
            return
        dt = float(np.mean(np.abs(valid_dt)))
        # Parse out.cgyro.info for ion direction to ensure correctness as per user request
        # Default to -1.0 (Standard CGYRO e^-iwt vs FFT e^-i2pift)
        freq_mult = -1.0
        
        try:
            data_dir = getattr(data, 'dir', getattr(data, 'path', ''))
            if data_dir:
                info_path = os.path.join(data_dir, 'out.cgyro.info')
                if os.path.exists(info_path):
                    with open(info_path, 'r') as f:
                        info_content = f.read()
                        if "Ion direction: omega > 0" in info_content:
                            freq_mult = -1.0
                        elif "Ion direction: omega < 0" in info_content:
                            freq_mult = -1.0
        except Exception as e:
            print(f"Warning checking info file: {e}")

        # Standard FFT freq first
        omega = np.fft.fftfreq(nt_slice, d=dt) * 2 * np.pi
        field_omega_all = np.fft.fft(field_t_all, axis=-1)
        
        # Shift zero frequency to center
        omega_shifted = np.fft.fftshift(omega)
        field_omega_all_shifted = np.fft.fftshift(field_omega_all, axes=-1)

        # Apply sign correction
        omega_shifted *= freq_mult
        
        # Flip to ensure ascending axis for contourf
        # Since freq_mult is -1.0, the axis is now descending. Flip it.
        if freq_mult < 0:
            omega_shifted = np.flip(omega_shifted)
            field_omega_all_shifted = np.flip(field_omega_all_shifted, axis=-1)

        try:
            use_power = self.fft_spectrum_var.get().strip().lower() == "power"
        except Exception:
            use_power = False
        spectrum_label = "Power" if use_power else "Amplitude"

        def spectrum_metric(arr):
            """Map complex FFT coefficients to selected spectrum metric."""
            amp = np.abs(arr)
            return amp * amp if use_power else amp
        
        # 1. kx=0 component
        i_r_0 = data.n_radial // 2
        field_kx0 = field_omega_all_shifted[i_r_0, :, :]
        mag_kx0 = spectrum_metric(field_kx0)

        # 3. ky=0 component (omega vs kx)
        # field_omega_all_shifted: [n_radial, n_ky, n_omega]
        # We need [n_radial, n_omega] for ky=0
        # ky=0 is usually at index 0, but let's verify
        ky = self._positive_ky_axis(getattr(data, 'ky', []))
        if ky.size == 0:
            print("No ky data available for FFT")
            return

        ky_idx_0 = 0
        if abs(ky[0]) > 1e-6:
                # Search for 0
                idx = np.where(np.abs(ky) < 1e-6)[0]
                if len(idx) > 0:
                    ky_idx_0 = idx[0]
                else:
                    print("Warning: ky=0 not found. Using first ky index.")
                    ky_idx_0 = 0
        
        field_ky0 = field_omega_all_shifted[:, ky_idx_0, :] # [n_radial, n_omega]
        
        # 2. Calculate quantities based on View Mode
        view_mode = self.fft_view_var.get()
        analysis_mode = self.fft_mode_var.get()
        
        # Pre-calculate kx array
        # data.kx usually contains only positive modes [0, 1, ..., n_radial-1] * dkx
        # But for FFT shifting we need the full negative to positive range.
        # Reconstruct full kx array based on data.kx[1] (dkx) and n_radial
        
        n_radial = field_t_all.shape[0]
        length = float(getattr(data, 'length', 0.0))
        
        # Check if data.kx is sufficient
        kx_data = np.asarray(getattr(data, 'kx', [])).reshape(-1)
        if kx_data.size > 1:
            dkx = float(kx_data[1] - kx_data[0]) # Assuming uniform spacing
            if abs(dkx) < 1e-12 and kx_data.size > 2:
                dkx = float(kx_data[2] - kx_data[1])
        else:
            # Fallback if kx is too short (e.g. n_radial=1)
            dkx = 2 * np.pi / length if length > 0 else 1.0
        
        # CGYRO radial index is centered: i = n_radial//2 corresponds to kx=0.
        # Build a kx axis directly from centered integer mode index p.
        p_index = np.arange(n_radial) - (n_radial // 2)
        kx_centered = p_index * dkx

        if view_mode == "Omega vs ky":
            # --- Logic for Omega vs ky ---
            # 1. kx=0 component
            # 2. Sum over all kx
            
            mag_plot1 = mag_kx0
            
            if analysis_mode == "Linear":
                # Coherent sum
                complex_sum = np.sum(field_omega_all_shifted, axis=0) # Sum over radial (kx)
                mag_plot2 = spectrum_metric(complex_sum)
                title2_suffix = "(Coherent sum over all kx)"
                
                # Per-ky normalization
                for i in range(mag_plot1.shape[0]):
                    if not np.all(np.isfinite(mag_plot1[i, :])):
                        mag_plot1[i, :] = 0.0
                    elif not use_power:
                        m = np.max(mag_plot1[i, :])
                        if m > 0:
                            mag_plot1[i, :] /= m
                        
                for i in range(mag_plot2.shape[0]):
                    if not np.all(np.isfinite(mag_plot2[i, :])):
                        mag_plot2[i, :] = 0.0
                    elif not use_power:
                        m = np.max(mag_plot2[i, :])
                        if m > 0:
                            mag_plot2[i, :] /= m
                        
            else: # Nonlinear
                # Incoherent sum
                mag_plot2 = np.sum(spectrum_metric(field_omega_all_shifted), axis=0)
                title2_suffix = "(Incoherent sum over all kx)"
                
                # Global normalization
                if not np.all(np.isfinite(mag_plot1)):
                    mag_plot1[~np.isfinite(mag_plot1)] = 0.0
                if not use_power:
                    m = np.max(mag_plot1)
                    if m > 0:
                        mag_plot1 /= m
                
                if not np.all(np.isfinite(mag_plot2)):
                    mag_plot2[~np.isfinite(mag_plot2)] = 0.0
                if not use_power:
                    m = np.max(mag_plot2)
                    if m > 0:
                        mag_plot2 /= m

            # Plotting
            ky_grid, omega_grid = np.meshgrid(ky, omega_shifted, indexing='ij')
            
            self.fig.clf()
            ax1 = self.fig.add_subplot(211)
            ax2 = self.fig.add_subplot(212)
            
            cs1 = ax1.contourf(ky_grid, omega_grid, mag_plot1, levels=50, cmap='jet')
            self.fig.colorbar(cs1, ax=ax1)
            ax1.set_title(f'{field_name} FFT {spectrum_label} (kx=0): {label}')
            ax1.set_ylabel(r'$\omega (c_s/a)$')
            
            cs2 = ax2.contourf(ky_grid, omega_grid, mag_plot2, levels=50, cmap='jet')
            self.fig.colorbar(cs2, ax=ax2)
            ax2.set_title(f'{field_name} FFT {spectrum_label} {title2_suffix}: {label}')
            ax2.set_ylabel(r'$\omega (c_s/a)$')
            ax2.set_xlabel(r'$k_y \rho_s$')
            
            # Overlay Linear Frequency from exported omega/gamma-vs-ky file (same as Zonal path logic).
            if self.fft_overlay_freq_var.get():
                file_path = self._resolve_linear_gamma_file_path(data)
                if file_path is None:
                    print(
                        f"Linear spectrum file not found for {label}. "
                        f"Set 'Linear File' to a valid omega/gamma-vs-ky file."
                    )
                else:
                    ky_lin, omega_lin, _gamma_lin = self._read_linear_omega_gamma_file(file_path, label)
                    if ky_lin is not None and omega_lin is not None:
                        finite = np.isfinite(ky_lin) & np.isfinite(omega_lin)
                        if np.any(finite):
                            freq_ky = np.asarray(ky_lin[finite], dtype=float).reshape(-1)
                            freq_omega = np.asarray(omega_lin[finite], dtype=float).reshape(-1)
                            order = np.argsort(freq_ky)
                            freq_ky = freq_ky[order]
                            freq_omega = freq_omega[order]
                            # Clip linear curve to the ky-range covered by current FFT view.
                            ky_fft = np.asarray(ky, dtype=float).reshape(-1)
                            ky_fft_finite = ky_fft[np.isfinite(ky_fft)]
                            if ky_fft_finite.size > 0:
                                ky_min = float(np.min(ky_fft_finite))
                                ky_max = float(np.max(ky_fft_finite))
                                if ky_min > ky_max:
                                    ky_min, ky_max = ky_max, ky_min
                                cov = (freq_ky >= ky_min - 1e-12) & (freq_ky <= ky_max + 1e-12)
                                if np.any(cov):
                                    freq_ky = freq_ky[cov]
                                    freq_omega = freq_omega[cov]
                                else:
                                    print(
                                        f"No overlap between linear ky-range and FFT ky-range for {label} "
                                        f"(FFT ky: [{ky_min:.4g}, {ky_max:.4g}])."
                                    )
                                    freq_ky = np.asarray([], dtype=float)
                                    freq_omega = np.asarray([], dtype=float)
                            # Exported linear omega already uses the desired sign convention; plot directly.
                            if freq_ky.size > 0:
                                freq_omega_plot = freq_omega
                                ax1.plot(freq_ky, freq_omega_plot, 'k--', linewidth=1.5, label='Linear Freq')
                                ax2.plot(freq_ky, freq_omega_plot, 'k--', linewidth=1.5, label='Linear Freq')
                                ax1.legend(loc='upper right')
                                ax2.legend(loc='upper right')
            
            self.ax = ax2

        elif view_mode == "Omega vs kx":
            # --- Logic for Omega vs kx ---
            # 1. ky=0 component
            # 2. Sum over all ky
            
            # field_omega_all_shifted radial axis already follows centered CGYRO p-index.
            # Do not fftshift radial again, otherwise kx and data become misaligned.
            field_radial = field_omega_all_shifted
            
            # ky=0 component
            field_ky0 = field_radial[:, ky_idx_0, :]
            mag_plot1 = spectrum_metric(field_ky0)
            
            # Sum over all ky
            if analysis_mode == "Linear":
                # Coherent sum
                complex_sum = np.sum(field_radial, axis=1) # Sum over ky
                mag_plot2 = spectrum_metric(complex_sum)
                title2_suffix = "(Coherent sum over all ky)"
                
                # Per-kx normalization
                for i in range(mag_plot1.shape[0]):
                    if not np.all(np.isfinite(mag_plot1[i, :])):
                        mag_plot1[i, :] = 0.0
                    elif not use_power:
                        m = np.max(mag_plot1[i, :])
                        if m > 0:
                            mag_plot1[i, :] /= m
                        
                for i in range(mag_plot2.shape[0]):
                    if not np.all(np.isfinite(mag_plot2[i, :])):
                        mag_plot2[i, :] = 0.0
                    elif not use_power:
                        m = np.max(mag_plot2[i, :])
                        if m > 0:
                            mag_plot2[i, :] /= m
                        
            else: # Nonlinear
                # Incoherent sum
                mag_plot2 = np.sum(spectrum_metric(field_radial), axis=1)
                title2_suffix = "(Incoherent sum over all ky)"
                
                # Global normalization
                if not np.all(np.isfinite(mag_plot1)):
                    mag_plot1[~np.isfinite(mag_plot1)] = 0.0
                if not use_power:
                    m = np.max(mag_plot1)
                    if m > 0:
                        mag_plot1 /= m
                
                if not np.all(np.isfinite(mag_plot2)):
                    mag_plot2[~np.isfinite(mag_plot2)] = 0.0
                if not use_power:
                    m = np.max(mag_plot2)
                    if m > 0:
                        mag_plot2 /= m

            # Plotting
            kx_plot = kx_centered
            if len(kx_plot) != mag_plot1.shape[0]:
                n_r = mag_plot1.shape[0]
                kx_plot = (np.arange(n_r) - (n_r // 2)) * dkx

            # Filter for positive kx only
            pos_kx_indices = np.where(kx_plot >= 0)[0]
            if len(pos_kx_indices) == 0:
                print("No non-negative kx points found; using full kx range.")
                pos_kx_indices = np.arange(len(kx_plot))
            kx_plot = kx_plot[pos_kx_indices]
            mag_plot1 = mag_plot1[pos_kx_indices, :]
            mag_plot2 = mag_plot2[pos_kx_indices, :]

            kx_grid, omega_grid_kx = np.meshgrid(kx_plot, omega_shifted, indexing='ij')
            
            self.fig.clf()
            ax1 = self.fig.add_subplot(211)
            ax2 = self.fig.add_subplot(212)
            
            cs1 = ax1.contourf(kx_grid, omega_grid_kx, mag_plot1, levels=50, cmap='jet')
            self.fig.colorbar(cs1, ax=ax1)
            ax1.set_title(f'{field_name} FFT {spectrum_label} (ky=0): {label}')
            ax1.set_ylabel(r'$\omega (c_s/a)$')
            
            cs2 = ax2.contourf(kx_grid, omega_grid_kx, mag_plot2, levels=50, cmap='jet')
            self.fig.colorbar(cs2, ax=ax2)
            ax2.set_title(f'{field_name} FFT {spectrum_label} {title2_suffix}: {label}')
            ax2.set_ylabel(r'$\omega (c_s/a)$')
            ax2.set_xlabel(r'$k_x \rho_s$')
            
            self.ax = ax2

        self.fig.tight_layout()

    def _plot_flux(self, data, label, plot_type, t_indices, t_start, t_end, species_override_index):
        """
        Legacy flux plot entrypoint.

        This wrapper is kept for backward compatibility with existing dispatch
        code and any external monkey-patching. New/clear name:
        `_plot_flux_diagnostics`.
        """
        return self._plot_flux_diagnostics(
            data, label, plot_type, t_indices, t_start, t_end, species_override_index
        )

    def _plot_flux_diagnostics(self, data, label, plot_type, t_indices, t_start, t_end, species_override_index):
        """
        Plot flux diagnostics versus `ky`, estimated `kx`, `time`, or `ky-time`.

        Responsibilities:
        - Resolve species selection and display suffixes.
        - Apply optional real-ion normalization (y-scale and x-axis conversion).
        - Route to dedicated sub-paths:
          1) estimated `kx` spectrum,
          2) `ky-time` contour panels,
          3) 1D `vs ky`,
          4) 1D `vs time`.
        - Ensure all time-averaged legends include explicit averaging windows.
        """
        moment_idx = 1 if "Energy" in plot_type else 0 # 0=particle, 1=energy
        is_decomp = "Decomp" in plot_type
        plot_type_lower = plot_type.lower()
        # Keep only estimated kx-spectrum path (lky_flux branch is removed).
        vs_kx_estimated = "vs kx" in plot_type_lower

        ky_flux = None

        if vs_kx_estimated:
            n_species = int(getattr(data, 'n_species', 0))
            if n_species <= 0:
                n_species = max(1, len(self._get_case_species(data)))
        else:
            # Ensure ky-flux data is loaded
            if not hasattr(data, 'ky_flux'):
                try:
                    print(f"Loading flux for {label}...")
                    data.getflux()
                except Exception as e:
                    print(f"Could not load flux for {label}: {e}")

            if not hasattr(data, 'ky_flux'):
                print(f"No flux data available for {label}")
                return

            ky_flux = np.asarray(data.ky_flux)
            if ky_flux.ndim != 5:
                print(f"Unsupported ky_flux shape for {label}: {ky_flux.shape}")
                return
            n_species = ky_flux.shape[0]

        target_indices, spec_label = self._resolve_species_indices(
            data,
            n_species,
            label,
            species_override_index=species_override_index,
            fallback_first=True,
            main_ion_policy="all",
            single_species_only=False,
        )
        if not target_indices:
            return

        if spec_label:
            label = f"{label} - {spec_label}"
        flux_norm_scale = self._get_flux_real_ion_norm_scale(data, moment_idx, label=label)
        x_ky_scale = 1.0
        x_t_scale = 1.0
        if self._use_flux_real_ion_norm():
            norm_ctx = self._get_flux_real_ion_norm_context(data, label=label)
            if norm_ctx.get('valid', False):
                try:
                    x_ky_scale = float(norm_ctx.get('rhoc', 1.0))
                except Exception:
                    x_ky_scale = 1.0
                try:
                    x_t_scale = float(norm_ctx.get('vc', 1.0))
                except Exception:
                    x_t_scale = 1.0
                if (not np.isfinite(x_ky_scale)) or x_ky_scale <= 0.0:
                    x_ky_scale = 1.0
                if (not np.isfinite(x_t_scale)) or x_t_scale <= 0.0:
                    x_t_scale = 1.0
        t_start_plot = float(t_start) * x_t_scale
        t_end_plot = float(t_end) * x_t_scale
        avg_suffix_plot = self._format_avg_suffix(t_start_plot, t_end_plot, prefix="Avg")

        if vs_kx_estimated:
            self._plot_flux_vs_kx_estimated(
                data,
                label,
                moment_idx,
                target_indices,
                t_indices,
                t_start_plot,
                t_end_plot,
                is_decomp,
                flux_norm_scale=flux_norm_scale,
            )
            return

        flux_sel = ky_flux[target_indices, moment_idx] # [species, field, ky, t]

        ky_axis = self._positive_ky_axis(getattr(data, 'ky', []))
        t_axis = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        t_valid = np.asarray(t_indices, dtype=int)
        t_valid = t_valid[(t_valid >= 0) & (t_valid < flux_sel.shape[-1])]

        if "vs ky_time" in plot_type:
            ky_plot = ky_axis if ky_axis.size > 0 else np.arange(flux_sel.shape[2], dtype=float)
            ky_plot = np.asarray(ky_plot, dtype=float) * float(x_ky_scale)

            if t_axis.size <= 0:
                t_axis = np.arange(flux_sel.shape[-1], dtype=float)
            t_plot = np.asarray(t_axis, dtype=float) * float(x_t_scale)

            n_ky = min(int(ky_plot.size), int(flux_sel.shape[2]))
            n_t = min(int(t_plot.size), int(flux_sel.shape[-1]))
            if n_ky <= 0 or n_t <= 0:
                print(f"No valid ky/time grid for Flux vs ky_time in {label}")
                return

            ky_plot = ky_plot[:n_ky]
            t_plot = t_plot[:n_t]
            flux_use = flux_sel[:, :, :n_ky, :n_t]

            # For ky-time evolution, default to full time unless user explicitly sets time window.
            user_specified_time = bool(self.t_start_var.get().strip() or self.t_end_var.get().strip())
            if user_specified_time and t_valid.size > 0:
                t_use = t_valid[t_valid < n_t]
                if t_use.size > 0:
                    flux_use = flux_use[:, :, :, t_use]
                    t_plot = t_plot[t_use]

            flux_use = np.asarray(flux_use, dtype=float) * float(flux_norm_scale)
            if flux_use.size <= 0 or t_plot.size <= 0:
                print(f"No valid flux data for Flux vs ky_time in {label}")
                return

            sub = self._get_flux_species_subscript()
            y_lbl = r"$t\ (a/c_s)$"
            x_lbl = r"$k_y \rho_s$"
            is_energy = (moment_idx == 1)
            field_names = ["Phi", "Apar", "Bpar"]
            cmap = 'RdBu_r'

            if is_decomp:
                n_fields = min(int(flux_use.shape[1]), 3)
                if n_fields <= 0:
                    print(f"No decomposed field channels for Flux vs ky_time in {label}")
                    return

                self.fig.clear()
                axes = self.fig.subplots(n_fields, 1, squeeze=False)[:, 0]
                self.ax = axes[0]

                for i_field in range(n_fields):
                    ax_i = axes[i_field]
                    z_ky_t = np.sum(flux_use[:, i_field, :, :], axis=0)  # [ky, t]
                    z_t_ky = np.transpose(z_ky_t)  # [t, ky]

                    cs = ax_i.contourf(ky_plot, t_plot, z_t_ky, levels=80, cmap=cmap)
                    cb = self.fig.colorbar(cs, ax=ax_i)
                    ax_i.set_xlabel(x_lbl)
                    ax_i.set_ylabel(y_lbl)

                    field_tag = field_names[i_field] if i_field < len(field_names) else f"Field{i_field}"
                    if is_energy:
                        ax_i.set_title(rf"$Q_{{{sub}}}^{{{field_tag}}}(k_y,t)$: {label}")
                        cb.set_label(rf"$Q_{{{sub}}}^{{{field_tag}}}(k_y,t)$")
                    else:
                        ax_i.set_title(rf"$\Gamma_{{{sub}}}^{{{field_tag}}}(k_y,t)$: {label}")
                        cb.set_label(rf"$\Gamma_{{{sub}}}^{{{field_tag}}}(k_y,t)$")

                try:
                    self.fig.tight_layout()
                except Exception:
                    pass
                return

            # Default single-panel mode focuses on electrostatic channel (phi).
            if flux_use.shape[1] > 0:
                z_ky_t = np.sum(flux_use[:, 0, :, :], axis=0)  # [ky, t], phi channel
            else:
                z_ky_t = np.sum(flux_use, axis=(0, 1))  # fallback
            z_t_ky = np.transpose(z_ky_t)  # [t, ky]
            cs = self.ax.contourf(ky_plot, t_plot, z_t_ky, levels=80, cmap=cmap)
            cb = self.fig.colorbar(cs, ax=self.ax)
            self.ax.set_xlabel(x_lbl)
            self.ax.set_ylabel(y_lbl)
            if is_energy:
                self.ax.set_title(rf"$Q_{{{sub}}}^{{\phi}}(k_y,t)$: {label}")
                cb.set_label(rf"$Q_{{{sub}}}^{{\phi}}(k_y,t)$")
            else:
                self.ax.set_title(rf"$\Gamma_{{{sub}}}^{{\phi}}(k_y,t)$: {label}")
                cb.set_label(rf"$\Gamma_{{{sub}}}^{{\phi}}(k_y,t)$")
            return

        if "vs ky" in plot_type:
            x = ky_axis if ky_axis.size > 0 else np.arange(flux_sel.shape[2])
            x = np.asarray(x, dtype=float) * float(x_ky_scale)
            flux_ky_t = np.sum(flux_sel, axis=(0, 1)) # [ky, t]

            if is_decomp:
                field_names = ["Phi", "Apar", "Bpar"]
                if len(t_valid) > 0:
                    y_total = np.mean(flux_ky_t[:, t_valid], axis=1)
                    label_total = f"{label} (Total){avg_suffix_plot}"
                else:
                    y_total = flux_ky_t[:, -1]
                    label_total = f"{label} (Total)"
                y_total = np.asarray(y_total, dtype=float) * float(flux_norm_scale)
                self._plot_1d(x, y_total, label_total, plot_type)

                n_fields = flux_sel.shape[1]
                for f_idx in range(n_fields):
                    flux_field = np.sum(flux_sel[:, f_idx, :, :], axis=0) # [ky, t]
                    if len(t_valid) > 0:
                        y_field = np.mean(flux_field[:, t_valid], axis=1)
                    else:
                        y_field = flux_field[:, -1]
                    y_field = np.asarray(y_field, dtype=float) * float(flux_norm_scale)

                    fname = field_names[f_idx] if f_idx < len(field_names) else f"Field{f_idx}"
                    if len(t_valid) > 0:
                        label_field = f"{label} ({fname}){avg_suffix_plot}"
                    else:
                        label_field = f"{label} ({fname})"
                    self.ax.plot(x[:len(y_field)], y_field[:len(x)], linestyle='--', marker='x', label=label_field)
                return

            if len(t_valid) > 0:
                y = np.mean(flux_ky_t[:, t_valid], axis=1)
                label = f"{label}{self._format_avg_suffix(t_start_plot, t_end_plot, prefix='Time Avg')}"
            else:
                y = flux_ky_t[:, -1]
            y = np.asarray(y, dtype=float) * float(flux_norm_scale)
            self._plot_1d(x, y, label, plot_type)
            return

        if "vs Time" in plot_type:
            x = t_axis if t_axis.size > 0 else np.arange(flux_sel.shape[-1])
            x = np.asarray(x, dtype=float) * float(x_t_scale)
            y_total = np.sum(flux_sel, axis=(0, 1, 2)) # [t]
            n_t = min(len(x), len(y_total))
            x = x[:n_t]
            y_total = y_total[:n_t]
            y_total = np.asarray(y_total, dtype=float) * float(flux_norm_scale)
            t_valid = t_valid[t_valid < n_t]

            if is_decomp:
                mean_str = ""
                range_str = ""
                if len(t_valid) > 0:
                    mean_str = f"(Mean: {np.mean(y_total[t_valid]):.2e})"
                    range_str = self._format_avg_range_from_axis(x, t_valid, prefix="Avg")

                label_total = f"{label} (Total) {range_str} {mean_str}".strip()
                self.ax.plot(x, y_total, label=label_total)

                field_names = ["Phi", "Apar", "Bpar"]
                n_fields = flux_sel.shape[1]
                for f_idx in range(n_fields):
                    flux_field = np.sum(flux_sel[:, f_idx, :, :], axis=(0, 1)) # [t]
                    flux_field = flux_field[:n_t]
                    flux_field = np.asarray(flux_field, dtype=float) * float(flux_norm_scale)
                    mean_str_f = ""
                    range_str_f = ""
                    if len(t_valid) > 0:
                        mean_str_f = f"(Mean: {np.mean(flux_field[t_valid]):.2e})"
                        range_str_f = self._format_avg_range_from_axis(x, t_valid, prefix="Avg")
                    fname = field_names[f_idx] if f_idx < len(field_names) else f"Field{f_idx}"
                    label_field = f"{label} ({fname}) {range_str_f} {mean_str_f}".strip()
                    self.ax.plot(x, flux_field, linestyle='--', linewidth=1.0, label=label_field)
                return

            y = y_total
            if len(t_valid) > 0:
                y_subset = y[t_valid]
                mean_val = np.mean(y_subset)
                std_val = np.std(y_subset)
                label = (
                    f"{label} ({self._format_avg_range_from_axis(x, t_valid, prefix='Avg').strip(' ()')}, "
                    f"Mean: {mean_val:.2e}, Std: {std_val:.2e})"
                )

                t_mean_start = x[t_valid[0]]
                t_mean_end = x[t_valid[-1]]
                line, = self.ax.plot(x, y, label=label)
                self.ax.plot(
                    [t_mean_start, t_mean_end],
                    [mean_val, mean_val],
                    linestyle='--',
                    color=line.get_color(),
                    linewidth=1.5,
                )
                return

            self._plot_1d(x, y, label, plot_type)

    def _plot_fluctuation_1d(self, data, label, plot_type, t_indices, t_start, t_end):
        """Plot 1D fluctuation cuts versus ky, kx, or time."""
        # Determine field name from plot_type e.g. "Phi vs ky", "Apar vs Time"
        field_name = plot_type.split()[0]
        
        c_field = None
        if field_name == "Phi":
            if hasattr(data, 'kxky_phi'): c_field = data.kxky_phi
        elif field_name == "Apar":
            if hasattr(data, 'kxky_apar'): c_field = data.kxky_apar
        elif field_name == "Bpar":
            if hasattr(data, 'kxky_bpar'): c_field = data.kxky_bpar
            
        if c_field is None:
             print(f"Loading big field data for {label}...")
             if self._ensure_bigfield_loaded(data, label):
                 if field_name == "Phi" and hasattr(data, 'kxky_phi'): c_field = data.kxky_phi
                 elif field_name == "Apar" and hasattr(data, 'kxky_apar'): c_field = data.kxky_apar
                 elif field_name == "Bpar" and hasattr(data, 'kxky_bpar'): c_field = data.kxky_bpar
             else:
                 return

        if c_field is None:
            return

        # Use helper to select field similar to data_plot.py kxky_select
        # We focus on theta_plot midplane (or whatever logic is standard)
        # Default midplane index
        i_theta = data.theta_plot * 4 // 8
        if i_theta >= data.theta_plot: i_theta = 0
        
        # Helper to extract complex field
        # Note: data.py kxky_select slices [1:] for radial dimension to remove p=0/kx=0 "special" element
        def get_field_complex(d):
            """Convert supported packed field layouts into complex `[kx, ky, t]` data."""
            ndim = d.ndim
            if ndim == 5: # [2, nr, theta, ny, nt]
                 return d[0, 1:, i_theta, :, :] + 1j * d[1, 1:, i_theta, :, :]
            elif ndim == 4: # [nr, theta, ny, nt] - unlikely for phi but possible
                 return d[1:, i_theta, :, :]
            elif ndim == 6: # [2, nr, theta, species, ny, nt]
                 # Not expected for phi, but maybe for moments
                 return d[0, 1:, i_theta, 0, :, :] + 1j * d[1, 1:, i_theta, 0, :, :]
            else:
                 print(f"Unknown field shape: {d.shape}")
                 return None
        
        field_data = get_field_complex(c_field)
        if field_data is None: return

        # Apply GyroBohm normalization if needed (data_plot.py does f/rhonorm)
        # cgyrodata usually stores raw. data_plot.py:getnorm('elec') sets rhonorm = rho.
        
        rho_norm = 1.0
        if hasattr(data, 'rho'):
            rho_norm = data.rho
        if abs(rho_norm) < 1e-12:
            rho_norm = 1.0
            
        # Normalize
        field_data = field_data / rho_norm
        
        # shape: [nr-1, ny, nt]
        # Calculate amplitude squared |field|^2
        field_sq = np.abs(field_data)**2
        ky = self._positive_ky_axis(getattr(data, 'ky', []))
        t_valid = np.asarray(t_indices, dtype=int)
        t_valid = t_valid[(t_valid >= 0) & (t_valid < field_sq.shape[-1])]
        
        if "vs ky" in plot_type:
            # Plot time-averaged RMS amplitude vs ky
            # Sum over kx (axis 0) -> [ny, nt]
            field_ky_t = np.sum(field_sq, axis=0)
            
            # Time average over selected window
            if len(t_valid) > 0:
                # Average over time indices
                y_vals = np.mean(field_ky_t[:, t_valid], axis=1)
                label = self._append_avg_suffix(label, t_start, t_end, prefix="Avg")
            else:
                # Fallback to last time point
                y_vals = field_ky_t[:, -1]
                
            # Plot sqrt(Sum |field|^2) -> RMS amplitude
            y = np.sqrt(y_vals)
            
            x = ky
            if x.size != y.size:
                n = min(x.size, y.size)
                x = x[:n]
                y = y[:n]
            
            self._plot_1d(x, y, label, plot_type)
            self.ax.set_xlabel(r'$k_y \rho_s$')
            self.ax.set_ylabel(fr'$\sqrt{{\sum | {field_name} |^2}}$')
            # self.ax.set_yscale('log')
            # self.ax.set_xscale('log') # Usually log-log for spectra

        elif "vs kx" in plot_type:
            # Plot time-averaged RMS amplitude vs kx
            # Sum over ky (axis 1) -> [kx, nt]
            field_kx_t = np.sum(field_sq, axis=1)

            if len(t_valid) > 0:
                y_vals = np.mean(field_kx_t[:, t_valid], axis=1)
                label = self._append_avg_suffix(label, t_start, t_end, prefix="Avg")
            else:
                y_vals = field_kx_t[:, -1]

            y = np.sqrt(y_vals)

            n_kx = y.size
            x = np.asarray(getattr(data, 'kx', [])).reshape(-1)
            if x.size == n_kx + 1:
                # Some datasets include one extra "special" leftmost element.
                x = x[1:]
            if x.size != n_kx:
                if x.size > 1:
                    dkx = float(x[1] - x[0])
                else:
                    length = float(getattr(data, 'length', 0.0))
                    dkx = 2 * np.pi / length if length > 0 else 1.0
                x = (np.arange(n_kx) - (n_kx // 2)) * dkx

            self._plot_1d(x, y, label, plot_type)
            self.ax.set_xlabel(r'$k_x \rho_s$')
            self.ax.set_ylabel(fr'$\sqrt{{\sum | {field_name} |^2}}$')

        
        elif "vs Time" in plot_type:
            # Plot total RMS amplitude vs Time
            # Matches data_plot.py plot_phi logic
            
            # Separate n=0 and n>0
            # field_sq: [nr-1, ny, nt]
            # ny dimension corresponds to ky.
            # ky=0 is usually index 0.
            
            ky_idx_0 = 0
            if ky.size > 0 and abs(ky[0]) > 1e-6:
                 # Search for 0
                 idx = np.where(np.abs(ky) < 1e-6)[0]
                 if len(idx) > 0: ky_idx_0 = idx[0]
            ky_idx_0 = min(ky_idx_0, field_sq.shape[1] - 1)
            
            # n=0 intensity: sum over kx (axis 0) at ky=0
            field_sq_n0 = np.sum(field_sq[:, ky_idx_0, :], axis=0) # [nt]
            
            # n>0 intensity: sum over kx (axis 0) and sum over ky!=0
            # Create mask for ky!=0
            mask_n = np.ones(field_sq.shape[1], dtype=bool)
            mask_n[ky_idx_0] = False
            
            field_sq_nn = np.sum(field_sq[:, mask_n, :], axis=(0, 1)) # [nt]
            
            # RMS
            y0 = np.sqrt(field_sq_n0)
            yn = np.sqrt(field_sq_nn)
            
            x = np.asarray(getattr(data, 't', [])).reshape(-1)
            if x.size == 0:
                return
            n_t = min(x.size, y0.size, yn.size)
            x = x[:n_t]
            y0 = y0[:n_t]
            yn = yn[:n_t]
            valid_t = np.asarray(t_indices, dtype=int)
            valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
            
            # Plot
            if len(valid_t) > 0:
                # n=0 stats
                y0_subset = y0[valid_t]
                mean_y0 = np.mean(y0_subset)
                std_y0 = np.std(y0_subset)
                t_mean_start = x[valid_t[0]]
                t_mean_end = x[valid_t[-1]]
                avg_inner = self._format_avg_suffix(t_mean_start, t_mean_end, prefix="Avg").strip(" ()")
                label_n0 = (
                    f"{label} (n=0) "
                    f"({avg_inner}, "
                    f"Mean: {mean_y0:.2e}, Std: {std_y0:.2e})"
                )

                # n>0 stats
                yn_subset = yn[valid_t]
                mean_yn = np.mean(yn_subset)
                std_yn = np.std(yn_subset)
                label_nn = (
                    f"{label} (n>0) "
                    f"({avg_inner}, "
                    f"Mean: {mean_yn:.2e}, Std: {std_yn:.2e})"
                )
                
                # Plot n=0
                line0, = self.ax.plot(x, y0, label=label_n0, linestyle='--')
                self.ax.plot([t_mean_start, t_mean_end], [mean_y0, mean_y0],
                             linestyle=':', color=line0.get_color(), linewidth=1.5)
                
                # Plot n>0
                linen, = self.ax.plot(x, yn, label=label_nn)
                self.ax.plot([t_mean_start, t_mean_end], [mean_yn, mean_yn],
                             linestyle='--', color=linen.get_color(), linewidth=1.5)
            else:
                self.ax.plot(x, y0, label=f"{label} (n=0)", linestyle='--')
                self.ax.plot(x, yn, label=f"{label} (n>0)")
            
            self.ax.set_xlabel(r'$t (a/c_s)$')
            self.ax.set_ylabel(fr'$\sqrt{{\sum | {field_name} |^2}}$')
            # self.ax.set_yscale('log')

    def _get_zf_exb_phi_kx_t(self, data, label):
        """Extract zonal potential spectrum `phi(kx,t)` and aligned `kx,t` axes."""
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

        ky = np.asarray(getattr(data, 'ky', []), dtype=float).reshape(-1)
        if ky.size == 0:
            print(f"No ky axis available for {label}")
            return None, None, None
        ky_idx_0 = int(np.argmin(np.abs(ky)))
        if abs(ky[ky_idx_0]) > 1e-6:
            print(
                f"Warning: ky=0 not found for {label}; "
                f"using closest ky={ky[ky_idx_0]:.4g}"
            )

        # Zonal component: ky~0 and average over theta.
        phi_zonal = np.mean(phi[:, :, ky_idx_0, :], axis=1)  # [n_radial, n_time]

        kx = np.asarray(getattr(data, 'kx', []), dtype=float).reshape(-1)
        n_radial = phi_zonal.shape[0]
        if kx.size == n_radial:
            phi_kx_t = phi_zonal
        elif kx.size == n_radial - 1:
            # Typical CGYRO non-ZF layout: kx excludes leftmost special slot.
            phi_kx_t = phi_zonal[1:, :]
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

        rho = float(getattr(data, 'rho', 1.0))
        if np.isfinite(rho) and abs(rho) > 1e-12:
            phi_kx_t = phi_kx_t / rho

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

    def _plot_flux_vs_kx_estimated(
        self,
        data,
        label,
        moment_idx,
        target_indices,
        t_indices,
        t_start,
        t_end,
        is_decomp,
        flux_norm_scale=1.0,
    ):
        """
        Legacy wrapper for estimated flux-vs-kx path.

        New/clear name: `_plot_flux_kx_estimated_spectrum`.
        """
        return self._plot_flux_kx_estimated_spectrum(
            data,
            label,
            moment_idx,
            target_indices,
            t_indices,
            t_start,
            t_end,
            is_decomp,
            flux_norm_scale=flux_norm_scale,
        )

    def _plot_flux_kx_estimated_spectrum(
        self,
        data,
        label,
        moment_idx,
        target_indices,
        t_indices,
        t_start,
        t_end,
        is_decomp,
        flux_norm_scale=1.0,
    ):
        """
        Estimate and plot flux-vs-kx from spectral ES/EM proxy channels.

        Method summary:
        - Build proxy channels from cross-phases:
          ES from `phi` with `n/e`, EM from `A_parallel` and `B_parallel`.
        - Collapse `[theta, ky]` to get `flux(kx, t)`.
        - Apply user-selected time average and optional decomposition overlay.
        - Keep the conventional non-negative `kx` branch for display.
        """
        try:
            flux_norm_scale = float(flux_norm_scale)
            if (not np.isfinite(flux_norm_scale)) or flux_norm_scale <= 0.0:
                flux_norm_scale = 1.0
        except Exception:
            flux_norm_scale = 1.0

        flux_name = "Energy" if moment_idx == 1 else "Particle"
        print(f"Computing Flux vs kx (estimated) for {label} [{flux_name}]...")

        # Required: ES field + corresponding moment.
        phi = self._load_kxky_complex(data, label, 'kxky_phi', '.cgyro.kxky_phi', species_dependent=False)
        mom_n = self._load_kxky_complex(data, label, 'kxky_n', '.cgyro.kxky_n', species_dependent=True)
        mom_e = self._load_kxky_complex(data, label, 'kxky_e', '.cgyro.kxky_e', species_dependent=True)
        mom_v = self._load_kxky_complex(data, label, 'kxky_v', '.cgyro.kxky_v', species_dependent=True)
        apar = self._load_kxky_complex(data, label, 'kxky_apar', '.cgyro.kxky_apar', species_dependent=False)
        bpar = self._load_kxky_complex(data, label, 'kxky_bpar', '.cgyro.kxky_bpar', species_dependent=False)

        primary = mom_e if moment_idx == 1 else mom_n
        if phi is None or primary is None:
            print(f"Insufficient data for estimated Flux vs kx in {label} (need phi and n/e moments).")
            return

        ns = int(primary.shape[2]) if primary.ndim >= 5 else 1
        target = [i for i in target_indices if 0 <= i < ns]
        if len(target) == 0:
            print(f"No valid species indices for estimated Flux vs kx in {label}.")
            return

        primary_sel = np.sum(primary[:, :, target, :, :], axis=2)  # [nr,theta,ky,t]
        mom_apar_sel = None
        if moment_idx == 0 and mom_v is not None:
            mom_apar_sel = np.sum(mom_v[:, :, target, :, :], axis=2)  # particle EM proxy
        elif moment_idx == 1:
            mom_apar_sel = primary_sel  # heat EM proxy (estimated)
        mom_bpar_sel = primary_sel

        ky = np.asarray(getattr(data, 'ky', []), dtype=float).reshape(-1)

        dims_r = [primary_sel.shape[0], phi.shape[0]]
        dims_th = [primary_sel.shape[1], phi.shape[1]]
        dims_ky = [primary_sel.shape[2], phi.shape[2]]
        dims_t = [primary_sel.shape[3], phi.shape[3]]

        if apar is not None:
            dims_r.append(apar.shape[0]); dims_th.append(apar.shape[1]); dims_ky.append(apar.shape[2]); dims_t.append(apar.shape[3])
        if bpar is not None:
            dims_r.append(bpar.shape[0]); dims_th.append(bpar.shape[1]); dims_ky.append(bpar.shape[2]); dims_t.append(bpar.shape[3])
        if mom_apar_sel is not None:
            dims_r.append(mom_apar_sel.shape[0]); dims_th.append(mom_apar_sel.shape[1]); dims_ky.append(mom_apar_sel.shape[2]); dims_t.append(mom_apar_sel.shape[3])

        n_r = int(min(dims_r))
        n_th = int(min(dims_th))
        n_ky = int(min(dims_ky))
        n_t = int(min(dims_t))
        if ky.size > 0:
            n_ky = min(n_ky, ky.size)
        if n_r <= 0 or n_th <= 0 or n_ky <= 0 or n_t <= 0:
            print(f"Invalid dimensions for estimated Flux vs kx in {label}.")
            return

        primary_sel = primary_sel[:n_r, :n_th, :n_ky, :n_t]
        phi = phi[:n_r, :n_th, :n_ky, :n_t]
        if apar is not None:
            apar = apar[:n_r, :n_th, :n_ky, :n_t]
        if bpar is not None:
            bpar = bpar[:n_r, :n_th, :n_ky, :n_t]
        if mom_apar_sel is not None:
            mom_apar_sel = mom_apar_sel[:n_r, :n_th, :n_ky, :n_t]

        if ky.size >= n_ky:
            ky_use = ky[:n_ky]
        else:
            ky_use = np.arange(n_ky, dtype=float)
        ky4 = ky_use.reshape(1, 1, n_ky, 1)

        def channel(moment_arr, field_arr, prefactor):
            """Compute one estimated ES/EM proxy channel as flux(kx,t)."""
            if moment_arr is None or field_arr is None:
                return None
            # Estimated cross-phase flux kernel in spectral form.
            kernel = np.real(np.conj(moment_arr) * (1j * ky4 * field_arr * prefactor))
            # Average over theta and sum over ky -> kx spectrum vs time.
            return np.mean(np.sum(kernel, axis=2), axis=1)  # [nr,t]

        flux_phi = channel(primary_sel, phi, -1.0)                    # ES proxy (CGYRO-normalized estimate)
        flux_apar = channel(mom_apar_sel, apar, 1.0)                  # EM flutter proxy (CGYRO-normalized estimate)
        flux_bpar = channel(mom_bpar_sel, bpar, 1.0)                  # EM compressional proxy

        channel_list = [f for f in [flux_phi, flux_apar, flux_bpar] if f is not None]
        if len(channel_list) == 0:
            print(f"No usable ES/EM channels for estimated Flux vs kx in {label}.")
            return

        flux_total = np.zeros_like(channel_list[0], dtype=float)
        for f in channel_list:
            flux_total += f

        kx_axis, radial_idx = self._build_kx_axis(data, n_r, label)
        if kx_axis is None or radial_idx is None:
            print(f"Could not construct kx axis for estimated Flux vs kx in {label}.")
            return

        flux_total = flux_total[radial_idx, :]
        if flux_phi is not None:
            flux_phi = flux_phi[radial_idx, :]
        if flux_apar is not None:
            flux_apar = flux_apar[radial_idx, :]
        if flux_bpar is not None:
            flux_bpar = flux_bpar[radial_idx, :]

        t_axis = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        if t_axis.size <= 0:
            t_axis = np.arange(flux_total.shape[-1], dtype=float)
        nt = min(t_axis.size, flux_total.shape[-1])
        if nt <= 0:
            print(f"No time samples for estimated Flux vs kx in {label}.")
            return

        flux_total = flux_total[:, :nt]
        if flux_phi is not None: flux_phi = flux_phi[:, :nt]
        if flux_apar is not None: flux_apar = flux_apar[:, :nt]
        if flux_bpar is not None: flux_bpar = flux_bpar[:, :nt]

        t_valid = np.asarray(t_indices, dtype=int)
        t_valid = t_valid[(t_valid >= 0) & (t_valid < nt)]
        avg_suffix = self._format_avg_suffix(t_start, t_end, prefix="Avg")

        if t_valid.size > 0:
            y_total = np.mean(flux_total[:, t_valid], axis=1)
            label_total = f"{label} (Estimated ES+EM){avg_suffix}"
        else:
            y_total = flux_total[:, -1]
            label_total = f"{label} (Estimated ES+EM)"
        y_total = np.asarray(y_total, dtype=float) * flux_norm_scale

        x = np.asarray(kx_axis, dtype=float).reshape(-1)
        y_total = np.asarray(y_total, dtype=float).reshape(-1)
        n = min(x.size, y_total.size)
        x = x[:n]
        y_total = y_total[:n]

        sort_idx = np.argsort(x)
        x = x[sort_idx]
        y_total = y_total[sort_idx]

        # Keep positive branch to match common kx-spectrum display convention.
        pos_mask = None
        if np.any(x < 0.0) and np.any(x > 0.0):
            pos_mask = x >= 0.0
            x = x[pos_mask]
            y_total = y_total[pos_mask]

        self._plot_1d(x, y_total, label_total, "Flux vs kx (estimated)")

        if is_decomp:
            def prep(y_arr):
                """Align one decomposition channel with plotted `kx` ordering/masking."""
                if y_arr is None:
                    return None
                y = np.mean(y_arr[:, t_valid], axis=1) if t_valid.size > 0 else y_arr[:, -1]
                y = np.asarray(y, dtype=float).reshape(-1)[:n] * flux_norm_scale
                y = y[sort_idx]
                if pos_mask is not None:
                    y = y[pos_mask]
                return y

            y_phi = prep(flux_phi)
            y_apar = prep(flux_apar)
            y_bpar = prep(flux_bpar)

            if y_phi is not None:
                lbl = f"{label} (ES: " + r"$\phi$" + ")"
                if t_valid.size > 0:
                    lbl += avg_suffix
                self.ax.plot(x[:len(y_phi)], y_phi[:len(x)], linestyle='--', marker='x', label=lbl)
            if y_apar is not None:
                lbl = f"{label} (EM: " + r"$A_\parallel$" + ")"
                if t_valid.size > 0:
                    lbl += avg_suffix
                self.ax.plot(x[:len(y_apar)], y_apar[:len(x)], linestyle='--', marker='x', label=lbl)
            if y_bpar is not None:
                lbl = f"{label} (EM: " + r"$B_\parallel$" + ")"
                if t_valid.size > 0:
                    lbl += avg_suffix
                self.ax.plot(x[:len(y_bpar)], y_bpar[:len(x)], linestyle='--', marker='x', label=lbl)

    def _resolve_linear_gamma_file_path(self, data):
        """Resolve linear omega/gamma reference-file path from GUI value and fallbacks."""
        user_text = ""
        try:
            user_text = self.linear_gamma_file_var.get().strip()
        except Exception:
            user_text = ""
        if not user_text:
            user_text = DEFAULT_LINEAR_GAMMA_FILE

        # 1) Absolute path from GUI.
        if os.path.isabs(user_text):
            return user_text if os.path.isfile(user_text) else None

        # 2) Relative to case directory (and one level above).
        case_dir = self._resolve_case_dir(data)
        if case_dir:
            candidates = [
                os.path.join(case_dir, user_text),
                os.path.join(os.path.dirname(case_dir), user_text),
            ]
            for c in candidates:
                if os.path.isfile(c):
                    return c

        # 3) Relative to current working directory.
        local = os.path.abspath(user_text)
        if os.path.isfile(local):
            return local

        # 4) Fallback export folder from CGYRO_vs_CGYRO exporter.
        fallback = os.path.join(os.getcwd(), DEFAULT_EXPORT_DIRNAME, user_text)
        if os.path.isfile(fallback):
            return os.path.abspath(fallback)

        return None

    def _read_linear_omega_gamma_file(self, file_path, label):
        """Read `ky, omega, gamma` columns from exported linear spectrum file."""
        rows = []
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    s = line.strip()
                    if not s or s.startswith('#'):
                        continue
                    parts = re.split(r'[\s,]+', s)
                    vals = []
                    for p in parts:
                        try:
                            vals.append(float(p))
                        except Exception:
                            pass
                    if len(vals) >= 3:
                        rows.append(vals[:3])  # ky, omega, gamma
        except Exception as e:
            print(f"Failed to read linear spectrum file for {label}: {file_path} ({e})")
            return None, None, None

        if len(rows) == 0:
            print(f"No usable data rows in linear spectrum file for {label}: {file_path}")
            return None, None, None

        arr = np.asarray(rows, dtype=float)
        ky = arr[:, 0]
        omega = arr[:, 1]
        gamma = arr[:, 2]

        order = np.argsort(ky)
        return ky[order], omega[order], gamma[order]

    def _compute_zf_exb_shearing_time_series(self, data, label):
        """Compute scalar zonal ExB shearing-rate time series from `phi(kx,t)`."""
        kx, phi_kx_t, t = self._get_zf_exb_phi_kx_t(data, label)
        if kx is None:
            return None, None

        # CGYRO does not provide a standalone C_xy metric output; use normalized estimate.
        omega_kx_t = -(kx[:, None] ** 2) * phi_kx_t
        shear_rate = np.sqrt(np.sum(np.abs(omega_kx_t) ** 2, axis=0))

        n_t = min(t.size, shear_rate.size)
        if n_t <= 0:
            print(f"No time samples available for {label}")
            return None, None

        return t[:n_t], shear_rate[:n_t]

    def _plot_zf_exb_shearing_spectrum(self, data, label, t_indices, t_start, t_end):
        """
        Legacy wrapper for zonal ExB shearing spectrum plot.

        New/clear name: `_plot_zf_exb_shearing_rate_kx_spectrum`.
        """
        return self._plot_zf_exb_shearing_rate_kx_spectrum(
            data, label, t_indices, t_start, t_end
        )

    def _plot_zf_exb_shearing_rate_kx_spectrum(self, data, label, t_indices, t_start, t_end):
        """
        Plot time-averaged zonal ExB shearing spectrum versus `kx`.

        This uses `| -kx^2 * phi_zonal(kx,t) |` as the spectral proxy and
        keeps the non-negative `kx` branch for conventional presentation.
        """
        kx, phi_kx_t, t = self._get_zf_exb_phi_kx_t(data, label)
        if kx is None:
            return

        # CGYRO does not provide a standalone C_xy metric output; use normalized estimate.
        omega_kx_t = np.abs(-(kx[:, None] ** 2) * phi_kx_t)
        n_t = min(t.size, omega_kx_t.shape[-1])
        if n_t <= 0:
            print(f"No time samples available for {label}")
            return

        valid_t = np.asarray(t_indices, dtype=int)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size > 0:
            y = np.mean(omega_kx_t[:, valid_t], axis=1)
            label = self._append_avg_suffix(label, t_start, t_end, prefix="Avg")
        else:
            y = omega_kx_t[:, -1]

        x = np.asarray(kx, dtype=float).reshape(-1)
        if x.size != y.size:
            n = min(x.size, y.size)
            x = x[:n]
            y = y[:n]

        # Match literature-style plotting (positive kx branch).
        if np.any(x < 0.0) and np.any(x > 0.0):
            mask = x >= 0.0
            x = x[mask]
            y = y[mask]

        order = np.argsort(x)
        x = x[order]
        y = y[order]

        self._plot_1d(x, y, label, "ZF ExB Shearing Spectrum")
        self.ax.set_xlabel(r'$k_x \rho_s$')
        self.ax.set_ylabel(r'$\omega_{E\times B}^{ZF}\ (c_s/a)$')

    def _plot_zf_exb_fig4_kx_eq_ky(self, data, label, t_indices, t_start, t_end):
        """
        Legacy wrapper for the FIG4-style ZF vs gamma_lin comparison.

        New/clear name: `_plot_zf_exb_fig4_kx_equals_ky_comparison`.
        """
        return self._plot_zf_exb_fig4_kx_equals_ky_comparison(
            data, label, t_indices, t_start, t_end
        )

    def _plot_zf_exb_fig4_kx_equals_ky_comparison(self, data, label, t_indices, t_start, t_end):
        """
        FIG4-like comparison on shared `k` axis under `kx = ky` convention.

        Curves:
        - `<omega_ZF>(kx)` from zonal potential.
        - `kx * <V_ZF>` from time-averaged zonal flow velocity proxy.
        - `gamma_lin(ky)` from external linear-spectrum file.

        Optional:
        - If user inputs one `ky`, draw a vertical marker and append
          `omegaZF/gamma_lin` and `kxVzf/gamma_lin` to legend labels.
        """
        # Optional ky marker requested by user for ratio outputs.
        ky_target = None
        try:
            ky_text = str(self.zf_gamma_lin_ky_var.get()).strip()
            if len(ky_text) > 0:
                ky_target = float(ky_text)
        except Exception:
            ky_target = None

        # 1) Nonlinear zonal shearing spectrum averaged over selected time window.
        kx, phi_kx_t, t = self._get_zf_exb_phi_kx_t(data, label)
        if kx is None:
            return

        omega_kx_t = np.abs(-(kx[:, None] ** 2) * phi_kx_t)
        # V_ZF(t) = 0.5 * sqrt(sum_kx |kx * rho_D * phi_bar(kx, ktheta=0, t)|^2).
        # In current normalization, rho_D is treated as 1.0.
        vzf_t = 0.5 * np.sqrt(np.sum(np.abs(kx[:, None] * phi_kx_t) ** 2, axis=0))
        n_t = min(t.size, omega_kx_t.shape[-1])
        if n_t <= 0:
            print(f"No time samples available for {label}")
            return

        valid_t = np.asarray(t_indices, dtype=int)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        avg_suffix = self._format_avg_suffix(t_start, t_end, prefix="Avg")
        if valid_t.size > 0:
            omega_kx = np.mean(omega_kx_t[:, valid_t], axis=1)
            vzf_mean = float(np.mean(vzf_t[valid_t]))
            zf_label = f"{label} " + r"$\langle\omega_{ZF}(k_x)\rangle$" + avg_suffix
            vzf_label = f"{label} " + r"$k_x\langle V_{ZF}\rangle$" + avg_suffix
        else:
            omega_kx = omega_kx_t[:, -1]
            vzf_mean = float(np.mean(vzf_t[:n_t]))
            zf_label = f"{label} " + r"$\langle\omega_{ZF}(k_x)\rangle$"
            vzf_label = f"{label} " + r"$k_x\langle V_{ZF}\rangle$"

        x_zf = np.asarray(kx, dtype=float).reshape(-1)
        y_zf = np.asarray(omega_kx, dtype=float).reshape(-1)
        y_kx_vzf = np.asarray(x_zf * vzf_mean, dtype=float).reshape(-1)
        n = min(x_zf.size, y_zf.size, y_kx_vzf.size)
        x_zf = x_zf[:n]
        y_zf = y_zf[:n]
        y_kx_vzf = y_kx_vzf[:n]

        # Use positive branch to match common k>=0 FIG4 presentation.
        if np.any(x_zf < 0.0) and np.any(x_zf > 0.0):
            mask_pos = x_zf >= 0.0
            x_zf = x_zf[mask_pos]
            y_zf = y_zf[mask_pos]
            y_kx_vzf = y_kx_vzf[mask_pos]
        else:
            mask_pos = x_zf >= 0.0
            if np.any(mask_pos):
                x_zf = x_zf[mask_pos]
                y_zf = y_zf[mask_pos]
                y_kx_vzf = y_kx_vzf[mask_pos]

        order = np.argsort(x_zf)
        x_zf = x_zf[order]
        y_zf = y_zf[order]
        y_kx_vzf = y_kx_vzf[order]
        if x_zf.size <= 0:
            print(f"No usable kx points for {label}")
            return

        # 2) Linear gamma(ky) from file.
        file_path = self._resolve_linear_gamma_file_path(data)
        if file_path is None:
            print(
                f"Linear spectrum file not found for {label}. "
                f"Set 'Linear File' to a valid omega/gamma-vs-ky file."
            )
            return

        ky_lin, _omega_lin, gamma_lin = self._read_linear_omega_gamma_file(file_path, label)
        if ky_lin is None:
            return

        finite = np.isfinite(ky_lin) & np.isfinite(gamma_lin)
        if not np.any(finite):
            print(f"No finite ky/gamma points in {file_path} for {label}")
            return

        x_lin = np.asarray(ky_lin[finite], dtype=float).reshape(-1)
        y_lin = np.asarray(gamma_lin[finite], dtype=float).reshape(-1)
        mask_nonneg = x_lin >= 0.0
        if np.any(mask_nonneg):
            x_lin = x_lin[mask_nonneg]
            y_lin = y_lin[mask_nonneg]
        order_lin = np.argsort(x_lin)
        x_lin = x_lin[order_lin]
        y_lin = y_lin[order_lin]
        if x_lin.size <= 0:
            print(f"No non-negative ky points in {file_path} for {label}")
            return

        # Restrict nonlinear kx curves only by the right edge of linear gamma ky range.
        x_max = float(np.max(x_lin))
        in_right_range = x_zf <= x_max
        if np.any(in_right_range):
            x_zf = x_zf[in_right_range]
            y_zf = y_zf[in_right_range]
            y_kx_vzf = y_kx_vzf[in_right_range]
        else:
            print(
                f"No overlapping k-range for {label}: "
                f"kx in [{np.min(x_zf):.4g}, {np.max(x_zf):.4g}] "
                f"vs ky<= {x_max:.4g}"
            )
            return

        ratio_wzf = None
        ratio_kxvzf = None

        # 4) Optional ratios at selected ky target (also appended to legend labels).
        if ky_target is not None:
            g_lin_sel = float(np.interp(ky_target, x_lin, y_lin, left=np.nan, right=np.nan)) if x_lin.size > 0 else np.nan
            wzf_sel = float(np.interp(ky_target, x_zf, y_zf, left=np.nan, right=np.nan)) if x_zf.size > 0 else np.nan
            kxvzf_sel = float(np.interp(ky_target, x_zf, y_kx_vzf, left=np.nan, right=np.nan)) if x_zf.size > 0 else np.nan

            if np.isfinite(g_lin_sel) and abs(g_lin_sel) > 1.0e-12:
                ratio_wzf = wzf_sel / g_lin_sel if np.isfinite(wzf_sel) else np.nan
                ratio_kxvzf = kxvzf_sel / g_lin_sel if np.isfinite(kxvzf_sel) else np.nan
                print(
                    f"{label}: ky={ky_target:.6g}, "
                    f"gamma_lin={g_lin_sel:.6g}, "
                    f"omegaZF={wzf_sel:.6g}, "
                    f"kxVzf={kxvzf_sel:.6g}, "
                    f"omegaZF/gamma_lin={ratio_wzf:.6g}, "
                    f"kxVzf/gamma_lin={ratio_kxvzf:.6g}"
                )
            else:
                print(
                    f"{label}: ky={ky_target:.6g}, "
                    "cannot compute omegaZF/gamma_lin and kxVzf/gamma_lin "
                    "because gamma_lin is missing/zero at selected ky."
                )

        zf_label_plot = zf_label + r"  $(k_x=k_y)$"
        vzf_label_plot = vzf_label + r"  $(k_x=k_y)$"
        if ratio_wzf is not None and np.isfinite(ratio_wzf):
            zf_label_plot += rf" $[\omega_{{ZF}}/\gamma_{{lin}}={ratio_wzf:.3g}]$"
        if ratio_kxvzf is not None and np.isfinite(ratio_kxvzf):
            vzf_label_plot += rf" $[k_xV_{{ZF}}/\gamma_{{lin}}={ratio_kxvzf:.3g}]$"

        # 3) Plot on shared k-axis under mapping kx=ky.
        self.ax.plot(
            x_lin,
            y_lin,
            marker='o',
            linestyle='-',
            label=f"{label} " + r"$\gamma_{lin}(k_y)$",
        )
        self.ax.plot(
            x_zf,
            y_zf,
            marker='x',
            linestyle='--',
            label=zf_label_plot,
        )
        self.ax.plot(
            x_zf,
            y_kx_vzf,
            marker='^',
            linestyle='-.',
            label=vzf_label_plot,
        )
        self.ax.set_xlim(right=x_max)
        self.ax.set_xlabel(r'$k\rho_s\ \ (k_x=k_y)$')
        self.ax.set_ylabel(r'$\langle\omega_{ZF}\rangle,\ k_x\langle V_{ZF}\rangle,\ \gamma_{lin}\ (c_s/a)$')

        # 5) Optional vertical marker at selected ky.
        if ky_target is not None:
            try:
                ky_mark = float(ky_target)
            except Exception:
                ky_mark = np.nan
            if np.isfinite(ky_mark):
                x_min_vis = float(min(np.min(x_lin), np.min(x_zf)))
                x_max_vis = float(max(np.max(x_lin), np.max(x_zf)))
                if x_min_vis <= ky_mark <= x_max_vis:
                    # Avoid drawing the same ky marker repeatedly when plotting multiple cases.
                    has_same_marker = False
                    for ln in self.ax.lines:
                        if not getattr(ln, "_zf_ky_marker", False):
                            continue
                        xdat = np.asarray(ln.get_xdata(), dtype=float).reshape(-1)
                        if xdat.size >= 1 and np.isfinite(xdat[0]) and abs(float(xdat[0]) - ky_mark) <= 1.0e-12:
                            has_same_marker = True
                            break
                    if not has_same_marker:
                        vln = self.ax.axvline(
                            ky_mark, color='0.35', linestyle='--', linewidth=1.2, alpha=0.9
                        )
                        try:
                            vln._zf_ky_marker = True
                        except Exception:
                            pass

    def _parse_energy_balance_species(self):
        """Map energy-balance species combo to triad_v2 species index."""
        try:
            txt = str(self.energy_balance_spec_var.get()).strip().lower()
        except Exception:
            txt = "total (-1)"
        if "electron" in txt:
            return 1
        if "main ion" in txt:
            return 0
        return -1

    def _parse_energy_balance_n_index(self):
        """Read selected toroidal-mode index n (integer)."""
        try:
            return int(float(str(self.energy_balance_n_var.get()).strip()))
        except Exception:
            return 0

    def _load_triad_if_needed(self, data, label):
        """Ensure triad and ky_flux arrays are loaded for energy-balance plotting."""
        if not hasattr(data, 'triad'):
            if not self._ensure_bigfield_loaded(data, label):
                # `getbigfield()` may fail on other large-field files before it reaches
                # triad loading. Try triad-only fallback before giving up.
                pass

        if not hasattr(data, 'triad'):
            triad_loaded = False
            try:
                if hasattr(data, 'extract'):
                    _tmsg, fmt, raw = data.extract('.cgyro.triad')
                    if fmt != 'null':
                        n_n = int(getattr(data, 'n_n', 0))
                        n_radial = int(getattr(data, 'n_radial', 0))
                        n_species = int(getattr(data, 'n_species', 0))
                        n_time = int(getattr(data, 'n_time', 0))
                        if n_n > 0 and n_radial > 0 and n_species > 0:
                            block = 2 * n_n * n_radial * n_species * 8
                            if n_time <= 0 and block > 0:
                                n_time = int(raw.size // block)
                                if n_time > 0:
                                    try:
                                        data.n_time = n_time
                                    except Exception:
                                        pass
                            nd = block * n_time
                            if raw.size >= nd:
                                if raw.size > nd:
                                    raw = raw[:nd]
                                data.triad = np.reshape(
                                    raw,
                                    (2, n_species, n_radial, 8, n_n, n_time),
                                    order='F',
                                )
                                triad_loaded = True
            except Exception as e:
                print(f"Triad fallback read failed for {label}: {e}")

            if (not triad_loaded) and (not hasattr(data, 'triad')):
                try:
                    case_dir = self._resolve_case_dir(data)
                except Exception:
                    case_dir = ""
                hint = f" (case_dir={case_dir})" if case_dir else ""
                print(f"No triad data available for {label}. Need bin/out.cgyro.triad.{hint}")
                return False
        if not hasattr(data, 'ky_flux'):
            try:
                data.getflux()
            except Exception as e:
                print(f"Could not load flux for {label}: {e}")
        if not hasattr(data, 'ky_flux'):
            print(f"No ky_flux data available for {label}.")
            return False
        return True

    def _load_triad_only_if_needed(self, data, label):
        """Ensure triad array is loaded (without requiring ky_flux)."""
        if not hasattr(data, 'triad'):
            if not self._ensure_bigfield_loaded(data, label):
                pass

        if not hasattr(data, 'triad'):
            triad_loaded = False
            try:
                if hasattr(data, 'extract'):
                    _tmsg, fmt, raw = data.extract('.cgyro.triad')
                    if fmt != 'null':
                        n_n = int(getattr(data, 'n_n', 0))
                        n_radial = int(getattr(data, 'n_radial', 0))
                        n_species = int(getattr(data, 'n_species', 0))
                        n_time = int(getattr(data, 'n_time', 0))
                        if n_n > 0 and n_radial > 0 and n_species > 0:
                            block = 2 * n_n * n_radial * n_species * 8
                            if n_time <= 0 and block > 0:
                                n_time = int(raw.size // block)
                                if n_time > 0:
                                    try:
                                        data.n_time = n_time
                                    except Exception:
                                        pass
                            nd = block * n_time
                            if raw.size >= nd:
                                if raw.size > nd:
                                    raw = raw[:nd]
                                data.triad = np.reshape(
                                    raw,
                                    (2, n_species, n_radial, 8, n_n, n_time),
                                    order='F',
                                )
                                triad_loaded = True
            except Exception as e:
                print(f"Triad-only fallback read failed for {label}: {e}")

            if (not triad_loaded) and (not hasattr(data, 'triad')):
                try:
                    case_dir = self._resolve_case_dir(data)
                except Exception:
                    case_dir = ""
                hint = f" (case_dir={case_dir})" if case_dir else ""
                print(f"No triad data available for {label}. Need bin/out.cgyro.triad.{hint}")
                return False
        return True

    def _plot_energy_balance_gamma_eff_v3(self, data, label, t_indices, t_start, t_end):
        """
        Integrate cgyro_plot `plot_triad_v3` behavior:
        - always plot gamma_eff^NZ and gamma_eff^Z
        - without linear file: also plot gamma_eff^NL = gamma_eff^NZ + gamma_eff^Z
        - with linear file: plot gamma_lin and gamma_eff^stable = gamma_lin - gamma_eff^NL
        """
        if not self._load_triad_only_if_needed(data, label):
            return

        triad = np.asarray(data.triad)
        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Unsupported triad shape for {label}: {triad.shape}")
            return
        # triad: [2, n_species, n_radial, 8, n_n, n_t]
        _ri, _n_species, n_radial, n_chan, n_n, n_t = triad.shape
        if n_chan < 8:
            print(f"Unsupported triad channels for {label}: {n_chan} (<8)")
            return
        if n_n <= 1:
            print(f"plot_triad_v3-style plot needs n_n > 1 for {label}.")
            return

        # Sum species first (matches upstream plot_triad_v3).
        f = triad[0].sum(axis=0) + 1j * triad[1].sum(axis=0)  # [n_radial, 8, n_n, n_t]

        # Time averaging window from selected indices.
        valid_t = np.asarray(t_indices, dtype=int)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size <= 0:
            valid_t = np.arange(n_t, dtype=int)
        if valid_t.size <= 0:
            print(f"No usable time window for {label}.")
            return
        avg_suffix = self._format_avg_suffix(float(t_start), float(t_end), prefix="Avg")

        # idx5 (1-based) -> channel 4: delta S_a,k.
        ent = np.real(f[:, 4, :, :])  # [radial, n_n, t]
        ent_avg = np.mean(ent[:, :, valid_t], axis=2)

        def _time_avg_channel(arr):
            return np.mean(arr[:, :, valid_t], axis=2)  # [radial, n_n]

        # TYPE 1: idx2 with ky=0 removed => gamma_eff^NZ
        a_nz = np.real(f[:, 1, :, :]).copy()   # [radial, n_n, t]
        if a_nz.shape[1] > 0:
            a_nz[:, 0, :] = 0.0
        a_nz_avg = _time_avg_channel(a_nz)
        gamma_eff_nz = np.divide(
            -a_nz_avg, 2.0 * ent_avg,
            out=np.full_like(a_nz_avg, np.nan, dtype=float),
            where=np.abs(ent_avg) > 0.0,
        )

        # TYPE 2: idx1 - (idx2 with ky=0 removed) => gamma_eff^Z
        a_total = np.real(f[:, 0, :, :])
        a_z = a_total - a_nz
        a_z_avg = _time_avg_channel(a_z)
        gamma_eff_z = np.divide(
            -a_z_avg, 2.0 * ent_avg,
            out=np.full_like(a_z_avg, np.nan, dtype=float),
            where=np.abs(ent_avg) > 0.0,
        )
        gamma_eff_nl = gamma_eff_nz + gamma_eff_z

        ky = np.asarray(getattr(data, 'kynorm', getattr(data, 'ky', [])), dtype=float).reshape(-1)
        if ky.size <= 0:
            print(f"No ky axis available for {label}")
            return
        n_ky_use = min(int(ky.size), int(n_n), int(gamma_eff_nz.shape[1]))
        if n_ky_use <= 0:
            print(f"No usable ky points for {label}")
            return

        ky = ky[:n_ky_use]
        if ky[-1] < 0.0:
            ky = -ky
        row = int(n_radial // 2)
        g_nz = np.asarray(gamma_eff_nz[row, :n_ky_use], dtype=float)
        g_z = np.asarray(gamma_eff_z[row, :n_ky_use], dtype=float)
        g_nl = np.asarray(gamma_eff_nl[row, :n_ky_use], dtype=float)

        # Keep nonnegative branch in ascending ky.
        mask_nonneg = np.isfinite(ky) & (ky >= 0.0)
        if np.any(mask_nonneg):
            ky = ky[mask_nonneg]
            g_nz = g_nz[mask_nonneg]
            g_z = g_z[mask_nonneg]
            g_nl = g_nl[mask_nonneg]
        order = np.argsort(ky)
        ky = ky[order]
        g_nz = g_nz[order]
        g_z = g_z[order]
        g_nl = g_nl[order]
        if ky.size <= 0:
            print(f"No non-negative ky points for {label}")
            return

        # Optional ky marker requested by user.
        ky_target = 0.0
        try:
            ky_target = float(str(self.energy_balance_n_var.get()).strip())
        except Exception:
            ky_target = 0.0

        def _extend_to_origin(xv, yv):
            """Force one exact origin point (0,0) for curve continuity."""
            x1 = np.asarray(xv, dtype=float).reshape(-1)
            y1 = np.asarray(yv, dtype=float).reshape(-1)
            m = np.isfinite(x1) & np.isfinite(y1)
            x1 = x1[m]
            y1 = y1[m]
            if x1.size <= 0:
                return np.asarray([0.0]), np.asarray([0.0])
            keep = np.abs(x1) > 1.0e-12
            x2 = np.concatenate(([0.0], x1[keep]))
            y2 = np.concatenate(([0.0], y1[keep]))
            order2 = np.argsort(x2)
            return x2[order2], y2[order2]

        # Optional linear gamma file overlay.
        file_path = self._resolve_linear_gamma_file_path(data)
        have_linear = bool(file_path is not None and os.path.isfile(file_path))
        ratio_nz = None
        ratio_z = None
        ratio_stable = None
        if have_linear:
            ky_lin, _omega_lin, gamma_lin = self._read_linear_omega_gamma_file(file_path, label)
            if ky_lin is not None:
                finite = np.isfinite(ky_lin) & np.isfinite(gamma_lin)
                ky_lin = np.asarray(ky_lin[finite], dtype=float).reshape(-1)
                gamma_lin = np.asarray(gamma_lin[finite], dtype=float).reshape(-1)
                pos = ky_lin >= 0.0
                if np.any(pos):
                    ky_lin = ky_lin[pos]
                    gamma_lin = gamma_lin[pos]
                if ky_lin.size > 0:
                    ord_lin = np.argsort(ky_lin)
                    ky_lin = ky_lin[ord_lin]
                    gamma_lin = gamma_lin[ord_lin]
                    # Formula from user:
                    # 0 = (gamma_lin - gamma_eff_stable) - gamma_eff_NL
                    # => gamma_eff_stable = gamma_lin - gamma_eff_NL
                    gamma_lin_on_ky = np.interp(ky, ky_lin, gamma_lin, left=np.nan, right=np.nan)
                    gamma_eff_stable = gamma_lin_on_ky - g_nl
                    ky_lin_p, gamma_lin_p = _extend_to_origin(ky_lin, gamma_lin)
                    ky_p, g_nz_p = _extend_to_origin(ky, g_nz)
                    _ky_p2, g_z_p = _extend_to_origin(ky, g_z)
                    _ky_p3, gamma_eff_stable_p = _extend_to_origin(ky, gamma_eff_stable)

                    if abs(ky_target) > 1.0e-12:
                        if ky_lin.size > 0:
                            g_lin_sel = float(np.interp(ky_target, ky_lin, gamma_lin, left=np.nan, right=np.nan))
                        else:
                            g_lin_sel = np.nan
                        g_nz_sel = float(np.interp(ky_target, ky, g_nz, left=np.nan, right=np.nan)) if ky.size > 0 else np.nan
                        g_z_sel = float(np.interp(ky_target, ky, g_z, left=np.nan, right=np.nan)) if ky.size > 0 else np.nan
                        g_st_sel = float(np.interp(ky_target, ky, gamma_eff_stable, left=np.nan, right=np.nan)) if ky.size > 0 else np.nan
                        if np.isfinite(g_lin_sel) and abs(g_lin_sel) > 1.0e-12:
                            ratio_nz = g_nz_sel / g_lin_sel if np.isfinite(g_nz_sel) else np.nan
                            ratio_z = g_z_sel / g_lin_sel if np.isfinite(g_z_sel) else np.nan
                            ratio_stable = g_st_sel / g_lin_sel if np.isfinite(g_st_sel) else np.nan
                            print(
                                f"{label}: ky={ky_target:.6g}, "
                                f"gamma_lin={g_lin_sel:.6g}, "
                                f"gamma_eff_NZ={g_nz_sel:.6g}, "
                                f"gamma_eff_Z={g_z_sel:.6g}, "
                                f"gamma_eff_stable={g_st_sel:.6g}, "
                                f"ratios: NZ/lin={ratio_nz:.6g}, Z/lin={ratio_z:.6g}, stable/lin={ratio_stable:.6g}"
                            )
                        else:
                            print(
                                f"{label}: ky={ky_target:.6g}, "
                                "cannot compute ratios because gamma_lin is missing/zero at selected ky."
                            )

                    lab_nz = rf'{label} ' + r'$\gamma_{eff}^{NZ}$' + avg_suffix
                    lab_z = rf'{label} ' + r'$\gamma_{eff}^{Z}$' + avg_suffix
                    lab_st = rf'{label} ' + r'$\gamma_{eff}^{stable}$' + avg_suffix
                    if ratio_nz is not None and np.isfinite(ratio_nz):
                        lab_nz += rf'$(/\gamma_{{lin}}={ratio_nz:.3g})$'
                    if ratio_z is not None and np.isfinite(ratio_z):
                        lab_z += rf'$(/\gamma_{{lin}}={ratio_z:.3g})$'
                    if ratio_stable is not None and np.isfinite(ratio_stable):
                        lab_st += rf'$(/\gamma_{{lin}}={ratio_stable:.3g})$'

                    self.ax.plot(
                        ky_lin_p, gamma_lin_p, '-ok', linewidth=2.0, markersize=5.0,
                        label=rf'{label} ' + r'$\gamma_{lin}$' + avg_suffix
                    )
                    self.ax.plot(ky_p, g_nz_p, '-sg', linewidth=2.0, markersize=5.0, label=lab_nz)
                    self.ax.plot(ky_p, g_z_p, '-^r', linewidth=2.0, markersize=5.0, label=lab_z)
                    self.ax.plot(ky_p, gamma_eff_stable_p, '--m', linewidth=2.0, label=lab_st)
                else:
                    have_linear = False
            else:
                have_linear = False

        if not have_linear:
            ky_p, g_nz_p = _extend_to_origin(ky, g_nz)
            _ky_p2, g_z_p = _extend_to_origin(ky, g_z)
            _ky_p3, g_nl_p = _extend_to_origin(ky, g_nl)
            self.ax.plot(
                ky_p, g_nz_p, '-sg', linewidth=2.0, markersize=5.0,
                label=rf'{label} ' + r'$\gamma_{eff}^{NZ}$' + avg_suffix
            )
            self.ax.plot(
                ky_p, g_z_p, '-^r', linewidth=2.0, markersize=5.0,
                label=rf'{label} ' + r'$\gamma_{eff}^{Z}$' + avg_suffix
            )
            self.ax.plot(
                ky_p, g_nl_p, '--k', linewidth=2.0,
                label=rf'{label} ' + r'$\gamma_{eff}^{NL}$' + avg_suffix
            )
            if abs(ky_target) > 1.0e-12:
                print(
                    f"{label}: ky={ky_target:.6g}, linear gamma file not available, "
                    "so ratios to gamma_lin are not computed."
                )

        if abs(ky_target) > 1.0e-12:
            self.ax.axvline(ky_target, linestyle='--', color='0.35', linewidth=1.5)

        self.ax.set_xlabel(r'$k_y \rho_s$')
        self.ax.set_ylabel(r'$\gamma\ (c_s/a)$')
        self.ax.set_title(r'Effective growth rate: $\gamma_{eff}^{Z}$ and $\gamma_{eff}^{NZ}$')

    def _compute_zonal_wes_time_series(self, data, label):
        """
        Compute zonal electrostatic potential energy W_es(t) from ky=0 phi(kx,t):

          W_es(t) = sum_a [ n_a z_a^2/(2 T_a) * sum_kx( |phi_kx(t)|^2 * (1-Gamma0_a(kx)) ) ]
          Gamma0_a = I0(b_a) * exp(-b_a),  b_a = (k_perp rho_a)^2

        with k_perp -> |kx| for the zonal (ky=0) branch.
        """
        kx, phi_kx_t, t = self._get_zf_exb_phi_kx_t(data, label)
        if kx is None or phi_kx_t is None:
            print(f"Cannot compute W_es for {label}: missing zonal phi(kx,t).")
            return None, None

        x = np.asarray(kx, dtype=float).reshape(-1)
        phi = np.asarray(phi_kx_t)
        if phi.ndim != 2:
            print(f"Cannot compute W_es for {label}: unsupported phi_kx_t shape {phi.shape}.")
            return None, None
        nk = min(x.size, phi.shape[0])
        nt = min(np.asarray(t).size, phi.shape[1])
        if nk <= 0 or nt <= 0:
            print(f"Cannot compute W_es for {label}: empty kx/time dimensions.")
            return None, None
        x = x[:nk]
        phi = phi[:nk, :nt]
        t = np.asarray(t, dtype=float).reshape(-1)[:nt]

        # Keep all kx components present in CGYRO output for the sum.
        phi2 = np.abs(phi) ** 2  # [kx, t]

        z = np.asarray(getattr(data, 'z', []), dtype=float).reshape(-1)
        dens = np.asarray(getattr(data, 'dens', []), dtype=float).reshape(-1)
        temp = np.asarray(getattr(data, 'temp', []), dtype=float).reshape(-1)
        mass = np.asarray(getattr(data, 'mass', []), dtype=float).reshape(-1)
        n_species = int(getattr(data, 'n_species', 0))
        ns = min(n_species, z.size, dens.size, temp.size, mass.size)
        if ns <= 0:
            print(f"Cannot compute W_es for {label}: missing species arrays (z/dens/temp/mass).")
            return None, None

        wes = np.zeros(nt, dtype=float)
        kperp = np.abs(x)
        for i in range(ns):
            zi = float(z[i])
            ni = float(dens[i])
            ti = float(temp[i])
            mi = float(mass[i])
            if (not np.isfinite(zi)) or (not np.isfinite(ni)) or (not np.isfinite(ti)) or (not np.isfinite(mi)):
                continue
            if abs(zi) <= 1.0e-12 or ti <= 1.0e-12 or ni == 0.0 or mi <= 0.0:
                continue

            # rho_a/rho_ref in current normalization (rho_ref ~ sqrt(Te*mD)/(eB)).
            rho_ratio = np.sqrt(abs(ti * mi)) / abs(zi)
            b = (kperp * rho_ratio) ** 2
            gamma0 = np.i0(b) * np.exp(-b)
            one_minus_gamma0 = 1.0 - gamma0
            coeff = ni * (zi ** 2) / (2.0 * ti)
            wes += coeff * np.sum(phi2 * one_minus_gamma0[:, np.newaxis], axis=0)

        return t, wes

    def _compute_triad_v2_terms(self, data, label):
        """
        Rebuild core triad_v2 terms used by cgyro_plot.

        Returns dict with time-series arrays or None on failure.
        """
        if not self._load_triad_if_needed(data, label):
            return None

        triad = np.asarray(data.triad)
        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Unsupported triad shape for {label}: {triad.shape}")
            return None
        # triad: [2, n_species, n_radial, 8, n_n, n_t]
        _ri, n_species, _n_radial, n_chan, n_n, n_t = triad.shape
        if n_chan < 8:
            print(f"Unsupported triad channels for {label}: {n_chan} (<8)")
            return None

        n_sel = self._parse_energy_balance_n_index()
        if n_sel < 0 or n_sel >= n_n:
            print(f"Energy balance n index out of range for {label}: n={n_sel}, valid=[0,{n_n-1}]")
            return None

        spec = self._parse_energy_balance_species()

        f_complex = triad[0] + 1j * triad[1]  # [species, radial, 8, n, t]
        # Per-species channel sums over kx(radial), used by custom Lz definitions.
        chan_srct = np.real(f_complex[:, :, :, n_sel, :])   # [species, radial, 8, t]
        chan_sct = np.sum(chan_srct, axis=1)                # [species, 8, t]

        Wes_total = None
        try:
            N_s = chan_sct[:, 1, :]
            Ent_s = chan_sct[:, 2, :]
            t_wes, Wes_total = self._compute_zonal_wes_time_series(data, label)
            if Wes_total is None:
                raise RuntimeError("W_es computation failed.")
            Wes_total = np.asarray(Wes_total, dtype=float).reshape(-1)
            t_wes = np.asarray(t_wes, dtype=float).reshape(-1)
            n_common_w = min(int(Wes_total.size), int(t_wes.size), int(N_s.shape[-1]), int(Ent_s.shape[-1]))
            if n_common_w <= 0:
                raise RuntimeError("Empty W_es time series.")
            N_s = N_s[:, :n_common_w]
            Ent_s = Ent_s[:, :n_common_w]
            Wes_total = Wes_total[:n_common_w]
            t_wes = t_wes[:n_common_w]
            # Keep per-species and total N available for total-Lz definitions.
            n_species_use = int(N_s.shape[0])
            n_total_series = np.sum(N_s[:n_species_use, :], axis=0)
            ent_total_series = np.sum(Ent_s[:n_species_use, :], axis=0)
        except Exception:
            n_total_series = None
            ent_total_series = None

        if (not isinstance(n_total_series, np.ndarray)) or (not isinstance(ent_total_series, np.ndarray)):
            print(f"Energy-balance total-Lz construction failed for {label}: cannot build total N/idx3 series.")
            return None

        if spec == -1:
            f_used = np.sum(f_complex, axis=0)  # [radial, 8, n, t]
        else:
            if spec < 0 or spec >= n_species:
                print(f"Species index {spec} is out of range for {label}; n_species={n_species}")
                return None
            f_used = f_complex[spec]  # [radial, 8, n, t]

        # bin.cgyro.triad channel convention (1-based idx in docs):
        # idx1->ch0: T_a
        # idx2->ch1: [ky=0] N_a  (for ky!=0 it corresponds to T_a^*)
        # idx3->ch2: d(T_a*delta S_a)/dt
        # idx4->ch3: dW_em/dt
        # idx5->ch4: delta S_a
        # idx6->ch5: D_r,a
        # idx7->ch6: D_theta,a
        # idx8->ch7: D_c,a
        T = np.real(f_used[:, 0, n_sel, :])      # [radial, t]
        N_raw = np.real(f_used[:, 1, n_sel, :])  # non-zonal pairs transfer
        Ent = np.real(f_used[:, 2, n_sel, :])    # [radial, t]
        Wkt = np.real(f_used[:, 3, n_sel, :])    # dW_em/dt, [radial, t]
        diss_r = np.real(f_used[:, 5, n_sel, :]) # [radial, t]
        diss_th = np.real(f_used[:, 6, n_sel, :])
        diss_c = np.real(f_used[:, 7, n_sel, :])

        T0 = np.sum(T, axis=0)
        N0 = np.sum(N_raw, axis=0)
        Ent0 = np.sum(Ent, axis=0)
        Wkt0 = np.sum(Wkt, axis=0)
        diss_r0 = np.sum(diss_r, axis=0)
        diss_th0 = np.sum(diss_th, axis=0)
        diss_c0 = np.sum(diss_c, axis=0)

        ky_flux = np.asarray(data.ky_flux)
        if ky_flux.ndim != 5:
            print(f"Unsupported ky_flux shape for {label}: {ky_flux.shape}")
            return None
        # ky_flux: [species, moment, field, n, t]
        ys = np.sum(ky_flux, axis=2)  # [species, moment, n, t]
        if ys.shape[2] <= n_sel:
            print(f"ky_flux n-index out of range for {label}: n={n_sel}, max={ys.shape[2]-1}")
            return None
        y = ys[:, 1, n_sel, :]   # Energy flux-like channel
        yn = ys[:, 0, n_sel, :]  # Particle flux-like channel

        dlntdr = np.asarray(getattr(data, 'dlntdr', []), dtype=float).reshape(-1)
        dlnndr = np.asarray(getattr(data, 'dlnndr', []), dtype=float).reshape(-1)
        temp = np.asarray(getattr(data, 'temp', []), dtype=float).reshape(-1)
        if dlntdr.size <= 1 or dlnndr.size <= 1 or temp.size <= 1:
            print(f"Missing profile arrays (dlntdr/dlnndr/temp) for {label}.")
            return None
        if abs(float(dlntdr[1])) < 1.0e-12:
            print(f"Invalid dlntdr[1] ~ 0 for {label}.")
            return None

        if spec == 1:
            Q = y[spec, :]
            G = (dlnndr[1] - 1.5 * dlntdr[1]) * yn[spec, :] / dlntdr[1]
        elif spec == 0:
            Q = y[spec, :] * (dlntdr[0] / dlntdr[1])
            G = (dlnndr[0] - 1.5 * dlntdr[0]) * yn[spec, :] / dlntdr[1] * temp[0]
        else:
            ns = min(y.shape[0], dlntdr.size, dlnndr.size, temp.size)
            Q = np.sum(y[:ns, :] * dlntdr[:ns, np.newaxis], axis=0) / dlntdr[1]
            G = np.sum(
                (dlnndr[:ns, np.newaxis] - 1.5 * dlntdr[:ns, np.newaxis])
                * yn[:ns, :]
                * temp[:ns, np.newaxis],
                axis=0
            ) / dlntdr[1]

        t = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        if t.size <= 0:
            t = np.arange(int(n_t), dtype=float)
        n_total_n = int(np.asarray(n_total_series).size)
        ent_total_n = int(np.asarray(ent_total_series).size)
        wes_n = int(np.asarray(Wes_total).size) if isinstance(Wes_total, np.ndarray) else int(1e30)
        n_time_use = min(
            t.size, T0.size, Ent0.size, Wkt0.size, diss_r0.size, diss_th0.size,
            diss_c0.size, Q.size, G.size, n_total_n, ent_total_n, wes_n
        )
        if n_time_use <= 0:
            return None

        t_use = np.asarray(t[:n_time_use], dtype=float)
        wes_use = np.asarray(Wes_total[:n_time_use], dtype=float)
        n0_use = np.asarray(N0[:n_time_use], dtype=float)
        ent0_use = np.asarray(Ent0[:n_time_use], dtype=float)
        n_total_use = np.asarray(n_total_series[:n_time_use], dtype=float)
        ent_total_use = np.asarray(ent_total_series[:n_time_use], dtype=float)
        if t_use.size >= 2 and np.any(np.diff(t_use) != 0.0):
            dWes_dt = np.gradient(wes_use, t_use)
        else:
            dWes_dt = np.zeros_like(wes_use)
        # User-requested strict definitions:
        # dSg/dt = sum(idx3) - dWes/dt,  Lz = dWes/dt - sum(idx2)
        dSg_dt_from_wes = ent_total_use - dWes_dt
        Lz_total_from_wes = dWes_dt - n_total_use

        return {
            't': t_use,
            'n_sel': int(n_sel),
            'spec': int(spec),
            'T0': np.asarray(T0[:n_time_use], dtype=float),
            'N0': n0_use,
            'Ent0': ent0_use,
            'Wkt0': np.asarray(Wkt0[:n_time_use], dtype=float),
            'diss_r0': np.asarray(diss_r0[:n_time_use], dtype=float),
            'diss_th0': np.asarray(diss_th0[:n_time_use], dtype=float),
            'diss_c0': np.asarray(diss_c0[:n_time_use], dtype=float),
            'G_plus_Q': np.asarray((G + Q)[:n_time_use], dtype=float),
            'Wes0': wes_use,
            'dWes_dt': np.asarray(dWes_dt, dtype=float),
            'dSg_dt': np.asarray(dSg_dt_from_wes, dtype=float),
            'Lz_total': np.asarray(Lz_total_from_wes, dtype=float),
            'N_total': n_total_use,
            'idx3_total': ent_total_use,
        }

    def _plot_energy_balance_entropy(self, data, label, t_indices, t_start, t_end):
        """Plot gyro-center entropy balance (Fig.5a style, triad_v2-derived)."""
        terms = self._compute_triad_v2_terms(data, label)
        if terms is None:
            return

        # Force paper-like white background regardless of global matplotlib style.
        try:
            self.ax.set_facecolor('white')
            self.fig.patch.set_facecolor('white')
        except Exception:
            pass

        t = terms['t']
        T0 = terms['T0']
        N0 = terms['N0']
        Ent0 = terms['Ent0']
        GQ = terms['G_plus_Q']
        d_r = terms['diss_r0']
        d_th = terms['diss_th0']
        d_c = terms['diss_c0']

        # Fig5(a)-like visual mapping:
        # (T-N)^{NZ->Z} -> (T0 - N0)
        # D_Z           -> D_r + D_theta + D_c
        # -L_Z          -> -[dWes/dt - sum_a idx2(a)]
        # dS_Z^{(g)}/dt -> sum(idx3) - dWes/dt
        transfer = np.asarray(T0 - N0, dtype=float)
        Dz = np.asarray(d_r + d_th + d_c, dtype=float)
        Lz_total_neg = -np.asarray(terms['Lz_total'], dtype=float)
        dSdt = np.asarray(terms['dSg_dt'], dtype=float)

        # Match paper line styles/colors.
        self.ax.plot(
            t, transfer, '--', color='k', linewidth=2.0,
            label=r'$(\mathcal{T}-\mathcal{N})^{\mathrm{NZ}\rightarrow \mathrm{Z}}$'
        )
        self.ax.plot(
            t, Dz, '-', color='0.45', linewidth=2.0,
            label=r'$D_Z$'
        )
        self.ax.plot(
            t, Lz_total_neg, '--', color='#2f62c9', linewidth=2.0,
            label=r'$-L_{Z,\mathrm{total}}$'
        )
        self.ax.plot(
            t, dSdt, '-', color='orange', linewidth=2.0,
            label=r'$dS_Z^{(g)}/dt$'
        )

        ne_txt = "?"
        try:
            dlnndr = np.asarray(getattr(data, 'dlnndr', []), dtype=float).reshape(-1)
            if dlnndr.size > 1 and np.isfinite(dlnndr[1]):
                ne_txt = f"{float(dlnndr[1]):.1f}"
        except Exception:
            pass

        self.ax.set_xlim(float(t[0]), float(t[-1]))
        self.ax.set_xlabel(r"$t\ (a/c_s)$")
        self.ax.set_ylabel(r"$[L_{T_e}Q_{GBD}]$")
        self.ax.set_title(rf"Eq. 12 - $\delta S_Z^{{(g)}}$ energy balance ($n_e={ne_txt}$)")

    def _plot_energy_balance_zf(self, data, label, t_indices, t_start, t_end):
        """Plot ZF energy-transfer balance (Fig.5b style, triad_v2-derived)."""
        terms = self._compute_triad_v2_terms(data, label)
        if terms is None:
            return

        # Force paper-like white background regardless of global matplotlib style.
        try:
            self.ax.set_facecolor('white')
            self.fig.patch.set_facecolor('white')
        except Exception:
            pass

        t = terms['t']
        T0 = terms['T0']
        N0 = terms['N0']
        GQ = terms['G_plus_Q']
        d_r = terms['diss_r0']
        d_th = terms['diss_th0']
        d_c = terms['diss_c0']

        # Fig5(b)-like visual mapping:
        # N^{NZ->Z}      -> N0
        # D'_Z           -> 0 (assumed negligible in current approximation)
        # L_{Z,total}    -> dW_es/dt - sum_a idx2(a)
        # L_{T_e} dW_es/dt -> d/dt of W_es computed from ky=0 phi(kx,t)
        transfer_n = np.asarray(N0, dtype=float)
        Dz_prime = np.zeros_like(transfer_n)
        Lz_total = np.asarray(terms['Lz_total'], dtype=float)
        dWes_dt = np.asarray(terms['dWes_dt'], dtype=float)

        # Match paper line styles/colors.
        self.ax.plot(
            t, transfer_n, '-', color='k', linewidth=2.0,
            label=r'$\mathcal{N}^{\mathrm{NZ}\rightarrow \mathrm{Z}}$'
        )
        self.ax.plot(
            t, Dz_prime, '-', color='0.45', linewidth=2.0,
            label=r"$D'_Z$"
        )
        self.ax.plot(
            t, Lz_total, '-', color='#2f62c9', linewidth=2.0,
            label=r'$L_{Z,\mathrm{total}}$'
        )
        self.ax.plot(
            t, dWes_dt, '-', color='red', linewidth=2.0,
            label=r'$L_{T_e}\,dW_{es}/dt$'
        )

        ne_txt = "?"
        try:
            dlnndr = np.asarray(getattr(data, 'dlnndr', []), dtype=float).reshape(-1)
            if dlnndr.size > 1 and np.isfinite(dlnndr[1]):
                ne_txt = f"{float(dlnndr[1]):.1f}"
        except Exception:
            pass

        self.ax.set_xlim(float(t[0]), float(t[-1]))
        self.ax.set_title(rf"Eq. 13 - $W_{{es}}$ energy balance ($n_e={ne_txt}$)")
        self.ax.set_xlabel(r"Simulation Time  $[c_s/a]$")
        self.ax.set_ylabel(r"$[L_{T_e}\,Q_{GBD}]$")

    def _compute_energy_balance_single_vs_ky(self, data, label, n_sel, quantity_key, t_indices):
        """
        Build single-quantity spectrum versus ky from triad channels.

        quantity_key:
          - 'T'      -> < sum_kx Re(idx1) >_t
          - 'N'      -> < sum_kx Re(idx2) >_t
          - 'T-N'    -> T - N
          - 'entropy'-> log( sum_kx <delta S_a(idx5)>_t ) for ion/electron species
        """
        if not self._load_triad_only_if_needed(data, label):
            return None, None

        triad = np.asarray(data.triad)
        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Unsupported triad shape for {label}: {triad.shape}")
            return None, None

        _ri, _n_species, n_radial, n_chan, n_n, n_t = triad.shape
        if n_chan < 8:
            print(f"Unsupported triad channels for {label}: {n_chan} (<8)")
            return None, None
        if n_n <= 1:
            print(f"Single-plot vs ky needs n_n > 1 for {label}.")
            return None, None

        n_use = int(n_sel)
        if n_use < 0 or n_use >= n_n:
            print(f"Energy balance n index out of range for {label}: n={n_use}, valid=[0,{n_n-1}]")
            return None, None

        f = triad[0].sum(axis=0) + 1j * triad[1].sum(axis=0)  # [radial, 8, n_n, t]
        valid_t = np.asarray(t_indices, dtype=int).reshape(-1)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size <= 0:
            valid_t = np.arange(n_t, dtype=int)
        if valid_t.size <= 0:
            return None, None

        ky = np.asarray(getattr(data, 'kynorm', getattr(data, 'ky', [])), dtype=float).reshape(-1)
        n_ky = min(int(ky.size), int(n_n))
        if n_ky <= 0:
            print(f"No ky axis available for {label}")
            return None, None
        ky = ky[:n_ky]
        if ky[-1] < 0.0:
            ky = -ky

        t_avg = np.mean(np.real(f[:, 0, :n_ky, :][:, :, valid_t]), axis=2)  # [radial, n_ky]
        n_avg = np.mean(np.real(f[:, 1, :n_ky, :][:, :, valid_t]), axis=2)  # [radial, n_ky]
        t_ky = np.sum(t_avg, axis=0)  # [n_ky]
        n_ky_vals = np.sum(n_avg, axis=0)  # [n_ky]
        tn_ky = t_ky - n_ky_vals

        if quantity_key == 'T':
            y = t_ky
        elif quantity_key == 'N':
            y = n_ky_vals
        elif quantity_key == 'T-N':
            y = tn_ky
        else:
            # Entropy spectrum should use idx5 (delta S), not idx3.
            # For multi-species view (ion/electron), handled by dedicated plot path.
            s_avg = np.mean(np.real(f[:, 4, :n_ky, :][:, :, valid_t]), axis=2)  # [radial, n_ky]
            s_ky = np.sum(s_avg, axis=0)  # [n_ky]
            y = s_ky

        ky = np.asarray(ky, dtype=float).reshape(-1)
        y = np.asarray(y, dtype=float).reshape(-1)
        mask = np.isfinite(ky) & np.isfinite(y) & (ky >= 0.0)
        if not np.any(mask):
            return None, None
        ky = ky[mask]
        y = y[mask]
        order = np.argsort(ky)
        return ky[order], y[order]

    def _plot_energy_balance_single_entropy_vs_ky(self, data, label, t_indices, avg_suffix=""):
        """
        Legacy wrapper for entropy single-plot (`vs ky`) renderer.

        New/clear name: `_plot_energy_balance_single_entropy_spectrum`.
        """
        return self._plot_energy_balance_single_entropy_spectrum(
            data, label, t_indices, avg_suffix=avg_suffix
        )

    def _plot_energy_balance_single_entropy_spectrum(self, data, label, t_indices, avg_suffix=""):
        """
        Plot entropy spectrum (`vs ky`) for single-plot energy-balance mode.

        Definition and conventions:
        - Entropy uses triad channel idx5 (0-based channel index 4), as requested.
        - Time average uses selected window; fallback is last 50% if needed.
        - Layout is fixed two panels:
          top = ion group (black), bottom = electron group (red).
        - Case-level line styles are cycled consistently across redraws.
        """
        if not self._load_triad_only_if_needed(data, label):
            return

        triad = np.asarray(data.triad)
        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Unsupported triad shape for {label}: {triad.shape}")
            return

        _ri, n_species, n_radial, n_chan, n_n, n_t = triad.shape
        if n_chan < 8:
            print(f"Unsupported triad channels for {label}: {n_chan} (<8)")
            return

        valid_t = np.asarray(t_indices, dtype=int).reshape(-1)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size <= 0:
            valid_t = np.arange(n_t // 2, n_t, dtype=int)
        if valid_t.size <= 0:
            return

        ky = np.asarray(getattr(data, 'kynorm', getattr(data, 'ky', [])), dtype=float).reshape(-1)
        n_ky = min(int(ky.size), int(n_n))
        if n_ky <= 0:
            print(f"No ky axis available for {label}")
            return
        ky = ky[:n_ky]
        if ky[-1] < 0.0:
            ky = -ky
        mask_k = np.isfinite(ky) & (ky >= 0.0)
        if not np.any(mask_k):
            return
        ky_plot = ky[mask_k]

        f = triad[0] + 1j * triad[1]  # [species, radial, 8, n, t]
        # idx5 -> channel 4: delta S
        ds = np.real(f[:, :, 4, :n_ky, :][:, :, :, valid_t])  # [species, radial, n_ky, t_sel]
        ds_tavg = np.mean(ds, axis=3)  # [species, radial, n_ky]
        ds_sum_kx = np.sum(ds_tavg, axis=1)  # [species, n_ky]

        z = np.asarray(getattr(data, 'z', []), dtype=float).reshape(-1)
        n_sp = min(int(n_species), int(ds_sum_kx.shape[0]), int(z.size))
        if n_sp <= 0:
            return

        ion_idx = [i for i in range(n_sp) if np.isfinite(z[i]) and z[i] > 0.0]
        ele_idx = [i for i in range(n_sp) if np.isfinite(z[i]) and z[i] < 0.0]
        if len(ele_idx) <= 0 and n_sp > 1:
            ele_idx = [1]
        if len(ion_idx) <= 0 and n_sp > 0:
            ion_idx = [0]

        eps = 1.0e-30
        ls_case = self._get_case_linestyle(
            case_label=label,
            map_attr="_energy_entropy_case_style",
            line_styles=['-', '--', '-.', ':'],
        )

        if not getattr(self, "_energy_entropy_axes_active", False):
            self.fig.clear()
            self._reset_figure_layout_defaults()
            ax_ion = self.fig.add_subplot(2, 1, 1)
            ax_ele = self.fig.add_subplot(2, 1, 2, sharex=ax_ion)
            self.ax = ax_ion
            self._energy_entropy_axes = (ax_ion, ax_ele)
            self._energy_entropy_axes_active = True

            ax_ion.set_title(r"Entropy of each species  $\log\sum_{k_x}\langle \delta S_a\rangle_t$")
            ax_ion.set_ylabel("Ion")
            ax_ele.set_ylabel("Electron")
            ax_ele.set_xlabel(r"$k_\theta \rho_D$")
            ax_ion.tick_params(labelbottom=False)
        else:
            ax_ion, ax_ele = self._energy_entropy_axes

        for j, i_sp in enumerate(ion_idx):
            y = np.asarray(ds_sum_kx[i_sp, :n_ky], dtype=float)[mask_k]
            y_log = np.log(np.maximum(np.abs(y), eps))
            ls = ls_case if j == 0 else '--'
            sp_tag = "" if len(ion_idx) == 1 else str(j + 1)
            ax_ion.plot(
                ky_plot, y_log, ls, color='k', linewidth=1.6,
                label=f"{label} ion{sp_tag}{avg_suffix}"
            )

        for j, i_sp in enumerate(ele_idx):
            y = np.asarray(ds_sum_kx[i_sp, :n_ky], dtype=float)[mask_k]
            y_log = np.log(np.maximum(np.abs(y), eps))
            ls = ls_case if j == 0 else '--'
            sp_tag = "" if len(ele_idx) == 1 else str(j + 1)
            ax_ele.plot(
                ky_plot, y_log, ls, color='r', linewidth=1.6,
                label=f"{label} elec{sp_tag}{avg_suffix}"
            )

    def _plot_energy_balance_single(self, data, label, t_indices, t_start, t_end):
        """
        Legacy wrapper for single-plot energy-balance mode.

        New/clear name: `_plot_energy_balance_single_mode`.
        """
        return self._plot_energy_balance_single_mode(
            data, label, t_indices, t_start, t_end
        )

    def _plot_energy_balance_single_mode(self, data, label, t_indices, t_start, t_end):
        """
        Render single-plot energy-balance mode with user-selected quantity and x-axis.

        Supported quantities:
        - `T`, `N`, `T-N`, `entropy`.
        Supported axes:
        - `vs time`: direct time-trace.
        - `vs ky`: time-averaged spectrum on chosen `n` index.

        Notes:
        - Any averaged result carries `(Avg: t0-t1)` in legend text.
        - `entropy vs ky` is delegated to the dedicated two-panel helper.
        """
        try:
            mode_txt = str(self.energy_balance_single_xaxis_var.get()).strip().lower()
        except Exception:
            mode_txt = "vs time"
        try:
            qty_txt = str(self.energy_balance_single_quantity_var.get()).strip().lower()
        except Exception:
            qty_txt = "t-n"
        if qty_txt not in ['t', 'n', 't-n', 'entropy']:
            qty_txt = 't-n'

        avg_suffix = ""
        try:
            if t_indices is not None and len(t_indices) > 0:
                avg_suffix = self._format_avg_suffix(float(t_start), float(t_end), prefix="Avg")
        except Exception:
            avg_suffix = ""

        terms = self._compute_triad_v2_terms(data, label)
        if terms is None:
            return

        t = np.asarray(terms['t'], dtype=float).reshape(-1)
        y_time = None
        y_label = ""
        qty_name = ""
        if qty_txt == 't':
            y_time = np.asarray(terms['T0'], dtype=float)
            y_label = r"$\mathcal{T}^{NZ\rightarrow Z}$"
            qty_name = "T"
        elif qty_txt == 'n':
            y_time = np.asarray(terms['N0'], dtype=float)
            y_label = r"$\mathcal{N}^{NZ\rightarrow Z}$"
            qty_name = "N"
        elif qty_txt == 'entropy':
            y_time = np.asarray(terms['dSg_dt'], dtype=float)
            y_label = r"$dS_g/dt$"
            qty_name = "entropy"
        else:
            y_time = np.asarray(terms['T0'], dtype=float) - np.asarray(terms['N0'], dtype=float)
            y_label = r"$(\mathcal{T}-\mathcal{N})^{NZ\rightarrow Z}$"
            qty_name = "T-N"

        if mode_txt == 'vs ky':
            n_sel = self._parse_energy_balance_n_index()
            if qty_txt == 'entropy':
                self._plot_energy_balance_single_entropy_spectrum(
                    data, label, t_indices, avg_suffix=avg_suffix
                )
            else:
                ky, y_ky = self._compute_energy_balance_single_vs_ky(data, label, n_sel, qty_name, t_indices)
                if ky is None or y_ky is None:
                    return
                self.ax.plot(
                    ky, y_ky, '-o', linewidth=2.0, markersize=5.0,
                    label=f"{label}{avg_suffix}"
                )
                self.ax.set_xlabel(r"$k_y \rho_s$")
                self.ax.set_ylabel(y_label)
                self.ax.set_title(
                    f"Energy balance single plot ({qty_name}, vs ky, n={n_sel})"
                )
            return

        n = min(t.size, y_time.size)
        if n <= 0:
            return
        t = t[:n]
        y_time = y_time[:n]
        self.ax.plot(t, y_time, '-', linewidth=2.0, label=f"{label}{avg_suffix}")
        self.ax.set_xlabel(r"$t\ (a/c_s)$")
        self.ax.set_ylabel(y_label)
        self.ax.set_title(
            f"Energy balance single plot ({qty_name}, vs time)"
        )

    def _compute_case_time_window(self, data):
        """Compute plotting time-window indices for one case."""
        t_indices = []
        t_start, t_end = 0, 0
        if hasattr(data, 't'):
            t_indices, t_start, t_end = self._get_time_indices(data.t)
        return t_indices, t_start, t_end

    def _resolve_flux_2d_time_window_indices(self, t_axis):
        """
        Resolve Flux-vs-2D averaging window.

        Rules:
        - If time start/end is empty, invalid, or fully out of range -> use last 50%.
        - Otherwise use clipped [t_start, t_end] intersection with available axis.
        """
        t = np.asarray(t_axis, dtype=float).reshape(-1)
        if t.size <= 0:
            return np.array([], dtype=int), 0.0, 0.0

        t_min = float(np.min(t))
        t_max = float(np.max(t))
        t_start_txt = ""
        t_end_txt = ""
        try:
            t_start_txt = str(self.t_start_var.get()).strip()
            t_end_txt = str(self.t_end_var.get()).strip()
        except Exception:
            t_start_txt = ""
            t_end_txt = ""

        use_last_half = False
        try:
            if len(t_start_txt) == 0 or len(t_end_txt) == 0:
                use_last_half = True
            else:
                t_start = float(t_start_txt)
                t_end = float(t_end_txt)
                if (not np.isfinite(t_start)) or (not np.isfinite(t_end)):
                    use_last_half = True
                elif t_start < t_min or t_end > t_max:
                    # User specified range exceeds available data range.
                    use_last_half = True
                else:
                    if t_start >= t_end:
                        use_last_half = True
        except Exception:
            use_last_half = True

        if use_last_half:
            mid = t_min + 0.5 * (t_max - t_min)
            idx = np.where(t >= mid)[0]
            if idx.size <= 0:
                idx = np.arange(max(0, t.size // 2), t.size, dtype=int)
            if idx.size <= 0:
                idx = np.arange(t.size, dtype=int)
            return idx, float(t[idx[0]]), float(t[idx[-1]])

        idx = np.where((t >= t_start) & (t <= t_end))[0]
        if idx.size <= 0:
            mid = t_min + 0.5 * (t_max - t_min)
            idx = np.where(t >= mid)[0]
            if idx.size <= 0:
                idx = np.arange(max(0, t.size // 2), t.size, dtype=int)
            if idx.size <= 0:
                idx = np.arange(t.size, dtype=int)
        return idx, float(t[idx[0]]), float(t[idx[-1]])

    def _plot_flux_vs_2d_selected_cases(self, selected_cases, plot_type):
        """Plot Flux time-averaged trends across cases vs selected input parameter."""
        try:
            x_param = str(self.flux_scan_xparam_var.get()).strip()
        except Exception:
            x_param = ""

        varying_params = []
        try:
            varying_params = list(self._get_flux_2d_varying_params(selected_cases))
        except Exception:
            varying_params = []

        if len(x_param) == 0 and len(varying_params) > 0:
            x_param = str(varying_params[0])
            try:
                self.flux_scan_xparam_var.set(x_param)
            except Exception:
                pass

        if len(x_param) == 0:
            messagebox.showwarning(
                "Warning",
                "Flux vs 2D requires at least one varying parameter from selected cases."
            )
            return

        # Build per-case scalar maps from input.cgyro.
        case_scalar_items = self._collect_flux_scan_case_scalars(selected_cases)
        case_to_scalars = {name: scalars for name, scalars in case_scalar_items}
        if len(case_to_scalars) <= 0:
            messagebox.showwarning(
                "Warning",
                "No readable input.cgyro scalar parameters found for selected cases."
            )
            return

        moment_idx = 1 if "Energy" in plot_type else 0
        flux_is_energy = (moment_idx == 1)
        other_params = [k for k in varying_params if str(k) != x_param]
        # Quasi-neutral scan convenience:
        # when x-axis is one DLNNDR channel, the coupled DLNNDR channels move
        # together and should not split legend groups. Keep only non-DLNNDR
        # parameters (e.g. MASS) for grouping/line connection.
        x_upper = str(x_param).strip().upper()
        if x_upper.startswith("DLNNDR_"):
            other_params = [k for k in other_params if not str(k).strip().upper().startswith("DLNNDR_")]

        grouped = {}
        skipped = 0
        for case_name in selected_cases:
            data = self.cases.get(case_name, None)
            scalars = case_to_scalars.get(case_name, None)
            if data is None or not hasattr(scalars, 'get'):
                skipped += 1
                continue

            x_raw = scalars.get(x_param, None)
            try:
                x_val = float(x_raw)
            except Exception:
                skipped += 1
                continue

            if not hasattr(data, 'ky_flux'):
                try:
                    data.getflux()
                except Exception as e:
                    print(f"Could not load flux for {case_name}: {e}")
            if not hasattr(data, 'ky_flux'):
                skipped += 1
                continue

            ky_flux = np.asarray(data.ky_flux)
            if ky_flux.ndim != 5:
                print(f"Unsupported ky_flux shape for {case_name}: {ky_flux.shape}")
                skipped += 1
                continue

            n_species = int(ky_flux.shape[0])
            target_indices, spec_label = self._resolve_species_indices(
                data,
                n_species,
                case_name,
                species_override_index=None,
                fallback_first=True,
                main_ion_policy="all",
                single_species_only=False,
            )
            if not target_indices:
                skipped += 1
                continue

            flux_sel = ky_flux[target_indices, moment_idx]  # [species, field, ky, t]
            y_t = np.sum(flux_sel, axis=(0, 1, 2)).reshape(-1)
            if y_t.size <= 0:
                skipped += 1
                continue

            t_axis = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
            if t_axis.size <= 0:
                t_axis = np.arange(y_t.size, dtype=float)
            n_t = min(t_axis.size, y_t.size)
            if n_t <= 0:
                skipped += 1
                continue
            t_axis = t_axis[:n_t]
            y_t = y_t[:n_t]

            t_idx, t0_used, t1_used = self._resolve_flux_2d_time_window_indices(t_axis)
            t_idx = t_idx[(t_idx >= 0) & (t_idx < n_t)]
            if t_idx.size <= 0:
                t_idx = np.arange(max(0, n_t // 2), n_t, dtype=int)
                if t_idx.size <= 0:
                    t_idx = np.arange(n_t, dtype=int)
            y_val = float(np.mean(y_t[t_idx])) if t_idx.size > 0 else float(y_t[-1])

            flux_norm_scale = self._get_flux_real_ion_norm_scale(
                data, moment_idx, label=f"{case_name}{' - ' + spec_label if spec_label else ''}"
            )
            y_val = float(y_val) * float(flux_norm_scale)

            # Group by all other varying parameters.
            grp_parts = []
            for k in other_params:
                try:
                    v = float(scalars.get(k))
                    grp_parts.append((str(k), v))
                except Exception:
                    grp_parts.append((str(k), np.nan))
            grp_key = tuple(grp_parts)
            grouped.setdefault(grp_key, []).append((float(x_val), float(y_val), case_name, t0_used, t1_used))

        if len(grouped) <= 0:
            messagebox.showwarning(
                "Warning",
                "No valid Flux-vs-2D points could be built from selected cases."
            )
            return

        # Plot one curve per group.
        for grp_key, pts in grouped.items():
            if len(pts) <= 0:
                continue

            x_vals = np.asarray([p[0] for p in pts], dtype=float)
            y_vals = np.asarray([p[1] for p in pts], dtype=float)
            finite = np.isfinite(x_vals) & np.isfinite(y_vals)
            if not np.any(finite):
                continue
            x_vals = x_vals[finite]
            y_vals = y_vals[finite]

            order = np.argsort(x_vals)
            x_vals = x_vals[order]
            y_vals = y_vals[order]

            # Merge repeated x by averaging.
            x_unique = []
            y_unique = []
            i0 = 0
            while i0 < len(x_vals):
                x0 = x_vals[i0]
                i1 = i0 + 1
                while i1 < len(x_vals) and abs(x_vals[i1] - x0) <= 1.0e-12:
                    i1 += 1
                x_unique.append(x0)
                y_unique.append(float(np.mean(y_vals[i0:i1])))
                i0 = i1

            if len(grp_key) <= 0:
                curve_label = "all other params fixed"
            else:
                parts = []
                for k, v in grp_key:
                    if np.isfinite(v):
                        parts.append(f"{k}={v:.4g}")
                    else:
                        parts.append(f"{k}=NA")
                curve_label = ", ".join(parts)

            self.ax.plot(
                np.asarray(x_unique, dtype=float),
                np.asarray(y_unique, dtype=float),
                linestyle='-',
                marker='o',
                markersize=5.0,
                linewidth=1.8,
                label=curve_label,
            )

        self.ax.set_xlabel(x_param)
        sub = self._get_flux_species_subscript()
        use_real_ion_norm = self._use_flux_real_ion_norm()
        if flux_is_energy:
            denom = r"Q_{GB,\mathrm{ri}}" if use_real_ion_norm else r"Q_{GB}"
            self.ax.set_ylabel(rf"$Q_{{{sub}}}/{denom}$")
        else:
            if use_real_ion_norm:
                self.ax.set_ylabel(rf"$\Gamma_{{{sub}}}/\Gamma_{{GB,\mathrm{{ri}}}}$")
            else:
                self.ax.set_ylabel(rf"$\Gamma_{{{sub}}}/Q_{{GB}}$")

        if skipped > 0:
            print(f"Flux vs 2D: skipped {skipped} case(s) due to missing parameter/flux/time data.")

    def _plot_energy_balance_vs_2d_selected_cases(self, selected_cases):
        """Plot energy-balance transfer rates vs scanned input parameter (2D mode style)."""
        try:
            x_param = str(self.flux_scan_xparam_var.get()).strip()
        except Exception:
            x_param = ""

        varying_params = []
        try:
            varying_params = list(self._get_flux_2d_varying_params(selected_cases))
        except Exception:
            varying_params = []

        if len(x_param) == 0 and len(varying_params) > 0:
            x_param = str(varying_params[0])
            try:
                self.flux_scan_xparam_var.set(x_param)
            except Exception:
                pass

        if len(x_param) == 0:
            messagebox.showwarning(
                "Warning",
                "Energy balance vs 2D requires at least one varying parameter from selected cases."
            )
            return

        case_scalar_items = self._collect_flux_scan_case_scalars(selected_cases)
        case_to_scalars = {name: scalars for name, scalars in case_scalar_items}
        if len(case_to_scalars) <= 0:
            messagebox.showwarning(
                "Warning",
                "No readable input.cgyro scalar parameters found for selected cases."
            )
            return

        n_sel = self._parse_energy_balance_n_index()
        other_params = [k for k in varying_params if str(k) != x_param]
        x_upper = str(x_param).strip().upper()
        if x_upper.startswith("DLNNDR_"):
            other_params = [k for k in other_params if not str(k).strip().upper().startswith("DLNNDR_")]

        grouped = {}
        skipped = 0

        for case_name in selected_cases:
            data = self.cases.get(case_name, None)
            scalars = case_to_scalars.get(case_name, None)
            if data is None or not hasattr(scalars, 'get'):
                skipped += 1
                continue

            try:
                x_val = float(scalars.get(x_param, None))
            except Exception:
                skipped += 1
                continue

            if not self._load_triad_only_if_needed(data, case_name):
                skipped += 1
                continue

            triad = np.asarray(data.triad)
            if triad.ndim != 6 or triad.shape[0] != 2:
                skipped += 1
                continue
            _ri, n_species, _n_radial, n_chan, n_n, n_t = triad.shape
            if n_chan < 8 or n_species < 2:
                skipped += 1
                continue
            if n_sel < 0 or n_sel >= n_n:
                skipped += 1
                continue

            t_axis = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
            if t_axis.size <= 0:
                t_axis = np.arange(n_t, dtype=float)
            n_t_use = min(t_axis.size, n_t)
            if n_t_use <= 0:
                skipped += 1
                continue
            t_axis = t_axis[:n_t_use]

            t_idx, _t0_used, _t1_used = self._resolve_flux_2d_time_window_indices(t_axis)
            t_idx = t_idx[(t_idx >= 0) & (t_idx < n_t_use)]
            if t_idx.size <= 0:
                t_idx = np.arange(max(0, n_t_use // 2), n_t_use, dtype=int)
                if t_idx.size <= 0:
                    t_idx = np.arange(n_t_use, dtype=int)
            if t_idx.size <= 0:
                skipped += 1
                continue

            f_complex = triad[0] + 1j * triad[1]  # [species, radial, 8, n, t]
            chan = np.real(f_complex[:, :, :, n_sel, :n_t_use])  # [species, radial, 8, t]

            # D: species index 0, e: species index 1 (as requested in current workflow)
            d_t = np.sum(chan[0, :, 0, :], axis=0)  # idx1 -> T
            d_n = np.sum(chan[0, :, 1, :], axis=0)  # idx2 -> N
            e_t = np.sum(chan[1, :, 0, :], axis=0)
            e_n = np.sum(chan[1, :, 1, :], axis=0)

            # S denominator definition (strict, from idx3):
            #   idx3 = d(T_a * deltaS_{a,k})/dt
            #   S(t) = sum_{a,p,n!=0} T_a*deltaS_{a,k}(p,n,t)
            # So we first sum idx3 over (a,p,n!=0), then integrate in time.
            dts_all = np.real(f_complex[:, :, 2, :, :n_t_use])  # [species, radial, n, t]
            if dts_all.shape[2] > 1:
                dts_nonzonal = dts_all[:, :, 1:, :]  # n != 0
            else:
                dts_nonzonal = dts_all[:, :, 0:0, :]

            if dts_nonzonal.shape[2] > 0:
                dS_dt_total = np.sum(dts_nonzonal, axis=(0, 1, 2))  # [t]
            else:
                # Fallback for degenerate grids without n>0 entries.
                dS_dt_total = np.sum(dts_all, axis=(0, 1, 2))  # [t]

            # Time integral with non-uniform dt support.
            s_total = np.zeros_like(dS_dt_total, dtype=float)
            if s_total.size > 1:
                dt = np.diff(t_axis)
                dt = np.where(np.isfinite(dt), dt, 0.0)
                dt = np.where(dt > 0.0, dt, 0.0)
                # Trapezoidal cumulative integral; integration constant set to 0.
                incr = 0.5 * (dS_dt_total[:-1] + dS_dt_total[1:]) * dt
                s_total[1:] = np.cumsum(incr)

            d_t_avg = float(np.mean(d_t[t_idx]))
            d_n_avg = float(np.mean(d_n[t_idx]))
            e_t_avg = float(np.mean(e_t[t_idx]))
            e_n_avg = float(np.mean(e_n[t_idx]))
            s_avg = float(np.mean(s_total[t_idx]))
            if (not np.isfinite(s_avg)) or abs(s_avg) <= 1.0e-12:
                skipped += 1
                continue

            grp_parts = []
            for k in other_params:
                try:
                    grp_parts.append((str(k), float(scalars.get(k))))
                except Exception:
                    grp_parts.append((str(k), np.nan))
            grp_key = tuple(grp_parts)
            grouped.setdefault(grp_key, []).append({
                'x': float(x_val),
                'N_D_S': float(d_n_avg / s_avg),
                'T_D_S': float(d_t_avg / s_avg),
                'N_e_S': float(e_n_avg / s_avg),
                'T_e_S': float(e_t_avg / s_avg),
                't0': float(t_axis[t_idx[0]]),
                't1': float(t_axis[t_idx[-1]]),
            })

        if len(grouped) <= 0:
            messagebox.showwarning(
                "Warning",
                "No valid Energy-balance-vs-2D points could be built from selected cases."
            )
            return

        for grp_key, rows in grouped.items():
            if len(rows) <= 0:
                continue
            x_vals = np.asarray([r['x'] for r in rows], dtype=float)
            nd_vals = np.asarray([r['N_D_S'] for r in rows], dtype=float)
            td_vals = np.asarray([r['T_D_S'] for r in rows], dtype=float)
            ne_vals = np.asarray([r['N_e_S'] for r in rows], dtype=float)
            te_vals = np.asarray([r['T_e_S'] for r in rows], dtype=float)

            finite = np.isfinite(x_vals) & np.isfinite(nd_vals) & np.isfinite(td_vals) & np.isfinite(ne_vals) & np.isfinite(te_vals)
            if not np.any(finite):
                continue
            x_vals = x_vals[finite]
            nd_vals = nd_vals[finite]
            td_vals = td_vals[finite]
            ne_vals = ne_vals[finite]
            te_vals = te_vals[finite]

            order = np.argsort(x_vals)
            x_vals = x_vals[order]
            nd_vals = nd_vals[order]
            td_vals = td_vals[order]
            ne_vals = ne_vals[order]
            te_vals = te_vals[order]

            # Merge duplicate x points by averaging.
            x_u, nd_u, td_u, ne_u, te_u = [], [], [], [], []
            i0 = 0
            while i0 < len(x_vals):
                x0 = x_vals[i0]
                i1 = i0 + 1
                while i1 < len(x_vals) and abs(x_vals[i1] - x0) <= 1.0e-12:
                    i1 += 1
                x_u.append(x0)
                nd_u.append(float(np.mean(nd_vals[i0:i1])))
                td_u.append(float(np.mean(td_vals[i0:i1])))
                ne_u.append(float(np.mean(ne_vals[i0:i1])))
                te_u.append(float(np.mean(te_vals[i0:i1])))
                i0 = i1

            if len(grp_key) <= 0:
                curve_suffix = "all other params fixed"
            else:
                parts = []
                for k, v in grp_key:
                    parts.append(f"{k}={v:.4g}" if np.isfinite(v) else f"{k}=NA")
                curve_suffix = ", ".join(parts)
            t0_vals = [float(r['t0']) for r in rows if np.isfinite(r.get('t0', np.nan))]
            t1_vals = [float(r['t1']) for r in rows if np.isfinite(r.get('t1', np.nan))]
            if len(t0_vals) > 0 and len(t1_vals) > 0:
                curve_suffix = f"{curve_suffix}, Avg: {min(t0_vals):.1f}-{max(t1_vals):.1f}"

            xu = np.asarray(x_u, dtype=float)
            self.ax.plot(xu, np.asarray(nd_u, dtype=float), '-s', color='k', linewidth=1.8, markersize=6.0, label=rf'$\mathcal{{N}}_D/S$ ({curve_suffix})')
            self.ax.plot(xu, np.asarray(td_u, dtype=float), '--s', color='k', linewidth=1.8, markersize=6.0, label=rf'$\mathcal{{T}}_D/S$ ({curve_suffix})')
            self.ax.plot(xu, np.asarray(ne_u, dtype=float), '-^', color='r', linewidth=1.8, markersize=6.0, label=rf'$\mathcal{{N}}_e/S$ ({curve_suffix})')
            self.ax.plot(xu, np.asarray(te_u, dtype=float), '--^', color='r', linewidth=1.8, markersize=6.0, label=rf'$\mathcal{{T}}_e/S$ ({curve_suffix})')

        self.ax.axhline(0.0, color='k', linewidth=1.2)
        self.ax.set_xlabel(x_param)
        self.ax.set_ylabel(r'$\mathcal{T}^{NZ\rightarrow Z}/S,\ \mathcal{N}^{NZ\rightarrow Z}/S$')
        self.ax.set_title('Free energy transferred from NZ to Z')

        if skipped > 0:
            print(f"Energy balance vs 2D: skipped {skipped} case(s) due to missing/invalid triad or scan data.")

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
            "ZF ExB Fig4 (kx=ky)": lambda: self._plot_zf_exb_fig4_kx_equals_ky_comparison(
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
    


