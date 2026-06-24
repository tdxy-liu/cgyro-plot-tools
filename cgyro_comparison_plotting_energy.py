"""
Option-specific plotting mixin for CGYRO comparison GUI.
Auto-extracted from cgyro_comparison_plotting.py during refactor.
"""

import numpy as np
import os
from tkinter import messagebox
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm, TwoSlopeNorm

FULLT_DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    'fullt_blue_white_red_suppressed_weak',
    [
        (0.00, '#053061'),
        (0.18, '#2166ac'),
        (0.34, '#92c5de'),
        (0.46, '#f7fbff'),
        (0.50, '#ffffff'),
        (0.54, '#fff5f0'),
        (0.66, '#f4a582'),
        (0.82, '#b2182b'),
        (1.00, '#67001f'),
    ],
    N=512,
)

class EnergyPlotting:
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

    def _parse_energy_balance_ky_scan(self):
        """Read selected physical ky value for single-plot scans."""
        try:
            return float(str(self.energy_balance_n_var.get()).strip())
        except Exception:
            return 0.0

    def _load_triad_common(self, data, label, require_flux):
        """Load triad data, optionally with ky_flux."""
        if not hasattr(data, 'triad'):
            triad_loaded = False
            try:
                triad = self._open_triad_memmap(data, label)
                if triad is not None:
                    data.triad = triad
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
        if not require_flux:
            return True

        if not hasattr(data, 'ky_flux'):
            try:
                data.getflux()
            except Exception as e:
                print(f"Could not load flux for {label}: {e}")
        if not hasattr(data, 'ky_flux'):
            print(f"No ky_flux data available for {label}.")
            return False
        return True

    def _open_triad_memmap(self, data, label):
        """Open bin/out.cgyro.triad as a Fortran-ordered memmap view."""
        try:
            case_dir = self._resolve_case_dir(data)
        except Exception:
            case_dir = getattr(data, "dir", None) or getattr(data, "path", None) or ""
        case_dir = str(case_dir or "")
        if not case_dir:
            return None

        paths = [
            os.path.join(case_dir, "bin.cgyro.triad"),
            os.path.join(case_dir, "bin", "bin.cgyro.triad"),
        ]
        path = None
        for candidate in paths:
            if os.path.isfile(candidate):
                path = candidate
                break
        if path is None:
            return None

        try:
            dtype = np.dtype(getattr(data, "BYTE", "float64"))
        except Exception:
            dtype = np.dtype("float64")

        n_n = int(getattr(data, "n_n", 0))
        n_radial = int(getattr(data, "n_radial", 0))
        n_species = int(getattr(data, "n_species", 0))
        if n_n <= 0 or n_radial <= 0 or n_species <= 0:
            return None

        try:
            n_elem = int(os.path.getsize(path) // dtype.itemsize)
        except Exception:
            return None
        block = 2 * n_species * n_radial * 8 * n_n
        if block <= 0 or n_elem % block != 0:
            print(
                f"Cannot infer triad shape for {label}: raw_size={n_elem}, "
                f"n_species={n_species}, n_radial={n_radial}, n_n={n_n}."
            )
            return None
        n_time = int(n_elem // block)
        if n_time <= 0:
            return None

        try:
            raw = np.memmap(path, dtype=dtype, mode="r", shape=(n_elem,))
            triad = np.ndarray(
                shape=(2, n_species, n_radial, 8, n_n, n_time),
                dtype=dtype,
                buffer=raw,
                order="F",
            )
        except Exception as e:
            print(f"Triad memmap open failed for {label}: {e}")
            return None

        try:
            setattr(data, "_triad_memmap_raw", raw)
            setattr(data, "_triad_memmap_path", path)
        except Exception:
            pass
        print(f"Triad for {label}: opened memmap {path} (n_time={n_time}).")
        return triad

    def _load_triad_if_needed(self, data, label):
        """Ensure triad and ky_flux arrays are loaded for energy-balance plotting."""
        return self._load_triad_common(data, label, require_flux=True)

    def _load_triad_only_if_needed(self, data, label):
        """Ensure triad array is loaded (without requiring ky_flux)."""
        return self._load_triad_common(data, label, require_flux=False)

    def _triad_species_view(self, triad, spec, label):
        """Return complex triad data for one species or total species sum."""
        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Unsupported triad shape for {label}: {triad.shape}")
            return None, None

        n_species = int(triad.shape[1])
        if spec < 0:
            return triad[0].sum(axis=0) + 1j * triad[1].sum(axis=0), "total species"
        if spec >= n_species:
            print(f"Energy balance species index out of range for {label}: {spec}, n_species={n_species}")
            return None, None
        return triad[0, spec] + 1j * triad[1, spec], f"species {spec}"

    def _triad_real_channel(
        self,
        triad,
        spec,
        channel,
        radial_sel=slice(None),
        ky_sel=slice(None),
        time_sel=slice(None),
    ):
        """Return real triad data without materializing a full complex array."""
        if triad.ndim != 6 or triad.shape[0] < 1:
            return None
        channel = int(channel)
        if channel < 0 or channel >= int(triad.shape[3]):
            return None

        def _select_axis(arr, axis, selector):
            if isinstance(selector, slice):
                index = [slice(None)] * arr.ndim
                index[int(axis)] = selector
                return arr[tuple(index)]
            return np.take(arr, selector, axis=int(axis))

        if int(spec) < 0:
            arr = triad[0, :, :, channel, :, :]  # [species, radial, ky, time]
            arr = _select_axis(arr, 3, time_sel)
            arr = _select_axis(arr, 2, ky_sel)
            arr = _select_axis(arr, 1, radial_sel)
            return np.sum(arr, axis=0)

        spec = int(spec)
        if spec < 0 or spec >= int(triad.shape[1]):
            spec = 0
        arr = triad[0, spec, :, channel, :, :]  # [radial, ky, time]
        arr = _select_axis(arr, 2, time_sel)
        arr = _select_axis(arr, 1, ky_sel)
        arr = _select_axis(arr, 0, radial_sel)
        return arr

    @staticmethod
    def _nearest_finite_index(axis, value):
        """Return the index of the finite axis value closest to `value`."""
        axis = np.asarray(axis, dtype=float).reshape(-1)
        try:
            value = float(value)
        except Exception:
            return None
        if not np.isfinite(value):
            return None
        finite_mask = np.isfinite(axis)
        if not np.any(finite_mask):
            return None
        axis_finite = axis[finite_mask]
        idx_local = int(np.argmin(np.abs(axis_finite - value)))
        return int(np.where(finite_mask)[0][idx_local])

    @staticmethod
    def _nearest_ky_index(axis, value):
        """Return nearest ky index using the physical signed ky value."""
        return EnergyPlotting._nearest_finite_index(axis, value)

    def _single_plot_ky_match(self, data, n_n, label, ky_scan=None):
        """Map a requested physical ky value to the nearest stored triad ky index."""
        if ky_scan is None:
            ky_req = self._parse_energy_balance_ky_scan()
        else:
            try:
                ky_req = float(ky_scan)
            except Exception:
                ky_req = 0.0

        ky_axis = np.asarray(getattr(data, 'kynorm', getattr(data, 'ky', [])), dtype=float).reshape(-1)
        if ky_axis.size >= n_n and np.any(np.isfinite(ky_axis[:n_n])):
            ky_use = self._nearest_ky_index(ky_axis[:n_n], ky_req)
        else:
            ky_use = int(round(ky_req))
        if ky_use is None:
            print(f"Cannot map ky scan={ky_req:.6g} for {label}")
            return None

        ky_use = int(ky_use)
        if ky_use < 0 or ky_use >= n_n:
            print(
                f"Energy balance ky scan out of range for {label}: "
                f"ky_scan={ky_req:.6g}, mapped index={ky_use}, valid=[0,{n_n-1}]"
            )
            return None
        ky_value = (
            float(ky_axis[ky_use])
            if ky_axis.size > ky_use and np.isfinite(ky_axis[ky_use])
            else float(ky_use)
        )
        return ky_use, ky_value, ky_req

    def _single_norm_mode(self):
        """Return the selected single-plot normalization mode."""
        try:
            text = str(self.energy_balance_single_norm_var.get()).strip().lower()
        except Exception:
            text = ""
        if "min" in text:
            return "min"
        if "max" in text:
            return "max"
        if text:
            return "none"
        # Backward compatibility for workspaces saved before the selector
        # existed.  The old checkbox meant "normalize by |min(T)|".
        try:
            if bool(self.energy_balance_single_norm_entropy_var.get()):
                return "min"
        except Exception:
            pass
        return "none"

    @staticmethod
    def _normalize_by_t_scale(values, t_values, mode, eps=1.0e-300):
        """Normalize `values` by |min(T)| or |max(T)| for one case/window."""
        arr = np.asarray(values, dtype=float)
        mode = str(mode).strip().lower()
        if mode in ("", "none"):
            return arr
        t_arr = np.asarray(t_values, dtype=float)
        finite_t = t_arr[np.isfinite(t_arr)]
        if finite_t.size <= 0:
            return np.full(arr.shape, np.nan, dtype=float)
        if mode == "min":
            ref = np.nanmin(finite_t)
        elif mode == "max":
            ref = np.nanmax(finite_t)
        else:
            return arr
        scale = abs(float(ref))
        if (not np.isfinite(scale)) or scale <= eps:
            return np.full(arr.shape, np.nan, dtype=float)
        return arr / scale

    @staticmethod
    def _single_norm_display(mode):
        """Return y-label, legend, and title suffixes for a norm mode."""
        mode = str(mode).strip().lower()
        if mode == "min":
            return r"$\,/\,|\min(T)|$", " / |min(T)|", " normalized by |min(T)|"
        if mode == "max":
            return r"$\,/\,|\max(T)|$", " / |max(T)|", " normalized by |max(T)|"
        return "", "", ""

    @staticmethod
    def _axis_cell_edges_from_centers(axis):
        """
        Convert monotonically increasing cell centers into pcolormesh edges.

        FULLT axes are physical spectral centers, not image pixel indices.  The
        half-step extrapolation below makes each color cell occupy the actual
        local dkx/dky interval instead of being stretched uniformly by imshow.
        """
        axis = np.asarray(axis, dtype=float).reshape(-1)
        if axis.size <= 0 or not np.all(np.isfinite(axis)):
            return None
        if axis.size == 1:
            v = float(axis[0])
            pad = 0.5 if abs(v) < 1.0e-12 else abs(v) * 0.05
            return np.asarray([v - pad, v + pad], dtype=float)

        diffs = np.diff(axis)
        if not np.all(np.isfinite(diffs)) or np.any(diffs <= 0.0):
            return None

        mid = 0.5 * (axis[:-1] + axis[1:])
        first = axis[0] - 0.5 * diffs[0]
        last = axis[-1] + 0.5 * diffs[-1]
        return np.concatenate(([first], mid, [last])).astype(float)

    @staticmethod
    def _diamond_lattice_from_centers(x_axis, y_axis, x_edges, y_edges, z_plot):
        """
        Build a node lattice that renders each FULLT sample as a diamond lobe.

        Real FULLT samples sit at center-center nodes; the surrounding edge
        nodes are zero so the color fades back to the neutral background before
        the next spectral sample.  Gouraud shading on this lattice gives the
        compact diamond blocks used for visual inspection without inventing new
        extrema between samples.
        """
        x_axis = np.asarray(x_axis, dtype=float).reshape(-1)
        y_axis = np.asarray(y_axis, dtype=float).reshape(-1)
        x_edges = np.asarray(x_edges, dtype=float).reshape(-1)
        y_edges = np.asarray(y_edges, dtype=float).reshape(-1)
        z_plot = np.asarray(z_plot, dtype=float)
        if z_plot.shape != (y_axis.size, x_axis.size):
            return None, None, None
        if x_edges.size != x_axis.size + 1 or y_edges.size != y_axis.size + 1:
            return None, None, None

        x_nodes = np.empty(2 * x_axis.size + 1, dtype=float)
        y_nodes = np.empty(2 * y_axis.size + 1, dtype=float)
        x_nodes[0::2] = x_edges
        x_nodes[1::2] = x_axis
        y_nodes[0::2] = y_edges
        y_nodes[1::2] = y_axis

        z_nodes = np.zeros((y_nodes.size, x_nodes.size), dtype=float)
        z_mask = np.zeros_like(z_nodes, dtype=bool)
        finite = np.isfinite(z_plot)
        z_nodes[1::2, 1::2] = np.where(finite, z_plot, 0.0)
        z_mask[1::2, 1::2] = ~finite

        x_grid, y_grid = np.meshgrid(x_nodes, y_nodes)
        return x_grid, y_grid, np.ma.array(z_nodes, mask=z_mask)

    @staticmethod
    def _sparse_axis_tick_values(axis, max_labels=11, include_zero=True):
        """Return a readable subset of axis-center ticks for dense spectral maps."""
        axis = np.asarray(axis, dtype=float).reshape(-1)
        axis = axis[np.isfinite(axis)]
        if axis.size <= 0:
            return np.asarray([], dtype=float)
        if axis.size <= max_labels:
            return axis

        idx = np.linspace(0, axis.size - 1, max_labels)
        idx = np.unique(np.rint(idx).astype(int))
        if include_zero:
            idx0 = int(np.argmin(np.abs(axis)))
            idx = np.unique(np.concatenate((idx, [idx0])))
        idx = idx[(idx >= 0) & (idx < axis.size)]
        return axis[np.sort(idx)]

    @staticmethod
    def _regular_axis_tick_values(vmin, vmax, step=0.1):
        """Return regular ticks covering [vmin, vmax], rounded for stable labels."""
        try:
            vmin = float(vmin)
            vmax = float(vmax)
            step = float(step)
        except Exception:
            return np.asarray([], dtype=float)
        if (not np.isfinite(vmin)) or (not np.isfinite(vmax)) or step <= 0.0:
            return np.asarray([], dtype=float)
        if vmax < vmin:
            vmin, vmax = vmax, vmin

        i0 = int(np.ceil((vmin - 1.0e-12) / step))
        i1 = int(np.floor((vmax + 1.0e-12) / step))
        if i1 < i0:
            return np.asarray([], dtype=float)
        ticks = np.arange(i0, i1 + 1, dtype=float) * step
        ticks[np.abs(ticks) < 1.0e-12] = 0.0
        return np.round(ticks, 12)

    @staticmethod
    def _fullt_x_tick_step_for_span(span):
        """Choose readable FULLT x-axis tick spacing from the current zoom span."""
        try:
            span = abs(float(span))
        except Exception:
            return 1.0
        if not np.isfinite(span):
            return 1.0
        if span > 4.0:
            return 1.0
        if span > 1.5:
            return 0.5
        if span > 0.8:
            return 0.2
        return 0.1

    def _apply_fullt_dynamic_x_ticks(self, ax):
        """Update FULLT x-axis labels after toolbar zoom/pan without changing grid lines."""
        if ax is None:
            return
        try:
            xmin, xmax = ax.get_xlim()
        except Exception:
            return
        if (not np.isfinite(xmin)) or (not np.isfinite(xmax)):
            return

        step = self._fullt_x_tick_step_for_span(xmax - xmin)
        ticks = self._regular_axis_tick_values(xmin, xmax, step=step)
        if ticks.size <= 0:
            return
        minor_ticks = self._regular_axis_tick_values(xmin, xmax, step=0.1)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{v:.3g}" for v in ticks])
        if minor_ticks.size > 0:
            ax.set_xticks(minor_ticks, minor=True)

    def _apply_fullt_dynamic_y_ticks(self, ax):
        """Use regular physical ky ticks for FULLT maps."""
        if ax is None:
            return
        try:
            ymin, ymax = ax.get_ylim()
        except Exception:
            return
        if (not np.isfinite(ymin)) or (not np.isfinite(ymax)):
            return

        step = self._fullt_x_tick_step_for_span(ymax - ymin)
        ticks = self._regular_axis_tick_values(ymin, ymax, step=step)
        if ticks.size <= 0:
            return
        minor_ticks = self._regular_axis_tick_values(ymin, ymax, step=0.1)
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{v:.3g}" for v in ticks])
        if minor_ticks.size > 0:
            ax.set_yticks(minor_ticks, minor=True)

    def _format_fullt_coord(self, x, y, kx_axis, ky_axis, z_plot):
        """Show cursor coordinates plus the nearest FULLT spectral sample."""
        try:
            x = float(x)
            y = float(y)
        except Exception:
            return ""
        if (not np.isfinite(x)) or (not np.isfinite(y)):
            return ""

        ix = self._nearest_finite_index(kx_axis, x)
        iy = self._nearest_finite_index(ky_axis, y)
        if ix is None or iy is None:
            return f"x={x:.4g} y={y:.4g}"

        try:
            kx = float(np.asarray(kx_axis, dtype=float).reshape(-1)[ix])
            ky = float(np.asarray(ky_axis, dtype=float).reshape(-1)[iy])
            val = float(np.asarray(z_plot, dtype=float)[iy, ix])
        except Exception:
            return f"x={x:.4g} y={y:.4g}"

        if np.isfinite(val):
            return f"x={x:.4g} y={y:.4g} | kx={kx:.4g} ky={ky:.4g} | T={val:.3e}"
        return f"x={x:.4g} y={y:.4g} | kx={kx:.4g} ky={ky:.4g}"

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
                        ky_lin_p, gamma_lin_p, '-o', linewidth=2.0, markersize=5.0,
                        color=self._get_gamma_lin_color(),
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
        S_raw = np.real(f_used[:, 4, n_sel, :])  # idx5: entropy-like state quantity
        Wkt = np.real(f_used[:, 3, n_sel, :])    # dW_em/dt, [radial, t]
        diss_r = np.real(f_used[:, 5, n_sel, :]) # [radial, t]
        diss_th = np.real(f_used[:, 6, n_sel, :])
        diss_c = np.real(f_used[:, 7, n_sel, :])

        T0 = np.sum(T, axis=0)
        N0 = np.sum(N_raw, axis=0)
        Ent0 = np.sum(Ent, axis=0)
        S0 = np.sum(S_raw, axis=0)
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
            'S0': np.asarray(S0[:n_time_use], dtype=float),
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

        colors = self._get_line_color_palette()
        component_styles = [
            dict(color=colors[0], linestyle='-', linewidth=2.0),
            dict(color=colors[1], linestyle='-', linewidth=2.0),
            dict(color=colors[2], linestyle='--', linewidth=2.0),
            dict(color=colors[3], linestyle='-', linewidth=2.0),
        ]
        self.ax.plot(
            t, transfer, **component_styles[0],
            label=r'$(\mathcal{T}-\mathcal{N})^{\mathrm{NZ}\rightarrow \mathrm{Z}}$'
        )
        self.ax.plot(
            t, Dz, **component_styles[1],
            label=r'$D_Z$'
        )
        self.ax.plot(
            t, Lz_total_neg, **component_styles[2],
            label=r'$-L_{Z,\mathrm{total}}$'
        )
        self.ax.plot(
            t, dSdt, **component_styles[3],
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

        colors = self._get_line_color_palette()
        component_styles = [
            dict(color=colors[0], linestyle='-', linewidth=2.0),
            dict(color=colors[1], linestyle='-', linewidth=2.0),
            dict(color=colors[2], linestyle='--', linewidth=2.0),
            dict(color=colors[3], linestyle='-', linewidth=2.0),
        ]
        self.ax.plot(
            t, transfer_n, **component_styles[0],
            label=r'$\mathcal{N}^{\mathrm{NZ}\rightarrow \mathrm{Z}}$'
        )
        self.ax.plot(
            t, Dz_prime, **component_styles[1],
            label=r"$D'_Z$"
        )
        self.ax.plot(
            t, Lz_total, **component_styles[2],
            label=r'$L_{Z,\mathrm{total}}$'
        )
        self.ax.plot(
            t, dWes_dt, **component_styles[3],
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

    def _compute_energy_balance_single_vs_ky(
        self, data, label, n_sel, quantity_key, t_indices, normalize_mode="none"
    ):
        """
        Build single-quantity spectrum versus ky from triad channels.

        quantity_key:
          - 'T'      -> < sum_kx Re(idx1) >_t
          - 'N'      -> < sum_kx Re(idx2) >_t
          - 'T-N'    -> T - N
          - 'Dr'     -> < sum_kx Re(idx6) >_t
          - 'Dtheta' -> < sum_kx Re(idx7) >_t
          - 'Dc'     -> < sum_kx Re(idx8) >_t
          - 'entropy'-> log( sum_kx <delta S_a(idx5)>_t ) for ion/electron species
        """
        del n_sel
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

        spec = self._parse_energy_balance_species()
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

        t_avg = np.mean(self._triad_real_channel(triad, spec, 0, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        n_avg = np.mean(self._triad_real_channel(triad, spec, 1, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        dr_avg = np.mean(self._triad_real_channel(triad, spec, 5, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        dtheta_avg = np.mean(self._triad_real_channel(triad, spec, 6, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        dc_avg = np.mean(self._triad_real_channel(triad, spec, 7, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        dz_avg = (
            dr_avg
            + dtheta_avg
            + dc_avg
        )
        t_ky = np.sum(t_avg, axis=0)  # [n_ky]
        n_ky_vals = np.sum(n_avg, axis=0)  # [n_ky]
        dr_ky = np.sum(dr_avg, axis=0)  # [n_ky]
        dtheta_ky = np.sum(dtheta_avg, axis=0)  # [n_ky]
        dc_ky = np.sum(dc_avg, axis=0)  # [n_ky]
        dz_ky = np.sum(dz_avg, axis=0)  # [n_ky]
        tn_ky = t_ky - n_ky_vals

        if quantity_key == 'T':
            y = t_ky
        elif quantity_key == 'N':
            y = n_ky_vals
        elif quantity_key == 'T-N':
            y = tn_ky
        elif quantity_key == 'Dr':
            y = dr_ky
        elif quantity_key == 'Dtheta':
            y = dtheta_ky
        elif quantity_key == 'Dc':
            y = dc_ky
        elif quantity_key == 'DZ':
            y = dz_ky
        else:
            # Entropy spectrum should use idx5 (delta S), not idx3.
            # For multi-species view (ion/electron), handled by dedicated plot path.
            s_avg = np.mean(self._triad_real_channel(triad, spec, 4, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
            s_ky = np.sum(s_avg, axis=0)  # [n_ky]
            y = s_ky
        if normalize_mode != "none" and quantity_key not in ['Dr', 'Dtheta', 'Dc', 'DZ', 'entropy']:
            y = self._normalize_by_t_scale(y, t_ky, normalize_mode)

        ky = np.asarray(ky, dtype=float).reshape(-1)
        y = np.asarray(y, dtype=float).reshape(-1)
        mask = np.isfinite(ky) & np.isfinite(y) & (ky >= 0.0)
        if not np.any(mask):
            return None, None
        ky = ky[mask]
        y = y[mask]
        order = np.argsort(ky)
        return ky[order], y[order]

    def _energy_balance_kx_selection(self, data, n_radial, label):
        """Return native cgyrodata kx axis plus matching triad radial indices."""
        kx = np.asarray(getattr(data, 'kx', []), dtype=float).reshape(-1)
        if kx.size == n_radial:
            return kx, np.arange(n_radial, dtype=int)

        if kx.size == n_radial - 1:
            # Native cgyrodata convention for non-zonal kxky arrays is:
            # data[1:, ...] is plotted against data.kx[:].  Keep the same
            # convention for triad single-plot kx diagnostics.
            return kx, np.arange(1, n_radial, dtype=int)

        print(
            f"Cannot map native cgyrodata kx axis for {label}: "
            f"kx length={kx.size}, radial length={n_radial}."
        )
        return None, None

    def _compute_energy_balance_single_vs_kx(
        self, data, label, ky_scan, quantity_key, t_indices, normalize_mode="none"
    ):
        """
        Build single-quantity spectrum versus kx at the nearest stored ky.

        quantity_key:
          - 'T'       -> < Re(idx1) >_t at selected ky
          - 'N'       -> < Re(idx2) >_t at selected ky
          - 'T-N'     -> T - N at selected ky
          - 'Dr'      -> < Re(idx6) >_t at selected ky
          - 'Dtheta'  -> < Re(idx7) >_t at selected ky
          - 'Dc'      -> < Re(idx8) >_t at selected ky
          - 'entropy' -> < Re(idx5) >_t at selected ky
        """
        if not self._load_triad_only_if_needed(data, label):
            return None, None, None

        triad = np.asarray(data.triad)
        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Unsupported triad shape for {label}: {triad.shape}")
            return None, None, None

        _ri, _n_species, n_radial, n_chan, n_n, n_t = triad.shape
        if n_chan < 8:
            print(f"Unsupported triad channels for {label}: {n_chan} (<8)")
            return None, None, None

        ky_match = self._single_plot_ky_match(data, n_n, label, ky_scan=ky_scan)
        if ky_match is None:
            return None, None, None
        ky_use, ky_value, ky_req = ky_match

        spec = self._parse_energy_balance_species()
        valid_t = np.asarray(t_indices, dtype=int).reshape(-1)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size <= 0:
            valid_t = np.arange(n_t, dtype=int)
        if valid_t.size <= 0:
            return None, None, None

        t_avg = np.mean(self._triad_real_channel(triad, spec, 0, ky_sel=ky_use, time_sel=valid_t), axis=1)
        n_avg = np.mean(self._triad_real_channel(triad, spec, 1, ky_sel=ky_use, time_sel=valid_t), axis=1)
        s_avg = np.mean(self._triad_real_channel(triad, spec, 4, ky_sel=ky_use, time_sel=valid_t), axis=1)
        dr_avg = np.mean(self._triad_real_channel(triad, spec, 5, ky_sel=ky_use, time_sel=valid_t), axis=1)
        dtheta_avg = np.mean(self._triad_real_channel(triad, spec, 6, ky_sel=ky_use, time_sel=valid_t), axis=1)
        dc_avg = np.mean(self._triad_real_channel(triad, spec, 7, ky_sel=ky_use, time_sel=valid_t), axis=1)
        dz_avg = (
            dr_avg
            + dtheta_avg
            + dc_avg
        )
        if quantity_key == 'T':
            y = t_avg
        elif quantity_key == 'N':
            y = n_avg
        elif quantity_key == 'T-N':
            y = t_avg - n_avg
        elif quantity_key == 'Dr':
            y = dr_avg
        elif quantity_key == 'Dtheta':
            y = dtheta_avg
        elif quantity_key == 'Dc':
            y = dc_avg
        elif quantity_key == 'DZ':
            y = dz_avg
        else:
            y = s_avg

        kx, radial_idx = self._energy_balance_kx_selection(data, n_radial, label)
        if kx is None or radial_idx is None:
            return None, None, None
        n_use = min(int(np.asarray(kx).size), int(np.asarray(radial_idx).size))
        if n_use <= 0:
            print(f"No usable kx points for {label}")
            return None, None, None

        kx = np.asarray(kx[:n_use], dtype=float).reshape(-1)
        radial_idx = np.asarray(radial_idx[:n_use], dtype=int).reshape(-1)
        valid_radial = (radial_idx >= 0) & (radial_idx < np.asarray(y).size)
        kx = kx[valid_radial]
        radial_idx = radial_idx[valid_radial]
        n_use = min(int(kx.size), int(radial_idx.size))
        kx = kx[:n_use]
        radial_idx = radial_idx[:n_use]
        y = np.asarray(y, dtype=float).reshape(-1)[radial_idx]
        t_ref = np.asarray(t_avg, dtype=float).reshape(-1)[radial_idx]
        mask = np.isfinite(kx) & np.isfinite(y)
        if not np.any(mask):
            print(f"No finite T-vs-kx points for {label}")
            return None, None, None

        kx = kx[mask]
        y = y[mask]
        t_ref = t_ref[mask]
        if normalize_mode != "none" and quantity_key not in ['Dr', 'Dtheta', 'Dc', 'DZ', 'entropy']:
            y = self._normalize_by_t_scale(y, t_ref, normalize_mode)
            mask = np.isfinite(kx) & np.isfinite(y)
            if not np.any(mask):
                print(f"No finite normalized T-vs-kx points for {label}")
                return None, None, None
            kx = kx[mask]
            y = y[mask]
        order = np.argsort(kx)

        return kx[order], y[order], (ky_value, ky_use, ky_req)

    def _compute_energy_balance_single_vs_kxky(
        self, data, label, n_sel, quantity_key, t_indices, normalize_mode="none"
    ):
        """
        Build single-quantity 2D map versus (kx, ky) from triad channels.

        quantity_key:
          - 'T'   -> < Re(idx1) >_t  on (kx, ky) grid
          - 'N'   -> < Re(idx2) >_t  on (kx, ky) grid
          - 'T-N' -> T - N on (kx, ky) grid
          - 'Dr'  -> < Re(idx6) >_t on (kx, ky) grid
          - 'Dtheta' -> < Re(idx7) >_t on (kx, ky) grid
          - 'Dc' -> < Re(idx8) >_t on (kx, ky) grid

        Returns (kx_grid, ky_grid, z_data) or None on failure.
        z_data shape: [n_kx, n_ky] (ready for contourf with transpose).
        """
        del n_sel
        if not self._load_triad_only_if_needed(data, label):
            return None

        triad = np.asarray(data.triad)
        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Unsupported triad shape for {label}: {triad.shape}")
            return None

        _ri, _n_species, n_radial, n_chan, n_n, n_t = triad.shape
        if n_chan < 8:
            print(f"Unsupported triad channels for {label}: {n_chan} (<8)")
            return None
        if n_n <= 1:
            print(f"Single-plot vs kxky needs n_n > 1 for {label}.")
            return None

        spec = self._parse_energy_balance_species()
        valid_t = np.asarray(t_indices, dtype=int).reshape(-1)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size <= 0:
            valid_t = np.arange(n_t, dtype=int)
        if valid_t.size <= 0:
            return None

        # Build ky axis
        ky = np.asarray(getattr(data, 'kynorm', getattr(data, 'ky', [])), dtype=float).reshape(-1)
        n_ky = min(int(ky.size), int(n_n))
        if n_ky <= 0:
            print(f"No ky axis available for {label}")
            return None
        ky = ky[:n_ky]
        if ky[-1] < 0.0:
            ky = -ky
        ky = np.asarray(ky, dtype=float).reshape(-1)

        kx, radial_idx = self._energy_balance_kx_selection(data, n_radial, label)
        if kx is None or radial_idx is None:
            return None

        # Time average over selected window
        if quantity_key == 'T':
            z = np.mean(self._triad_real_channel(triad, spec, 0, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        elif quantity_key == 'N':
            z = np.mean(self._triad_real_channel(triad, spec, 1, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        elif quantity_key == 'T-N':
            t_avg = np.mean(self._triad_real_channel(triad, spec, 0, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
            n_avg = np.mean(self._triad_real_channel(triad, spec, 1, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
            z = t_avg - n_avg
        elif quantity_key == 'Dr':
            z = np.mean(self._triad_real_channel(triad, spec, 5, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        elif quantity_key == 'Dtheta':
            z = np.mean(self._triad_real_channel(triad, spec, 6, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        elif quantity_key == 'Dc':
            z = np.mean(self._triad_real_channel(triad, spec, 7, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        else:
            # Fallback to T for unknown quantities
            z = np.mean(self._triad_real_channel(triad, spec, 0, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)
        t_ref = None
        if normalize_mode != "none" and quantity_key not in ['Dr', 'Dtheta', 'Dc']:
            t_ref = np.mean(self._triad_real_channel(triad, spec, 0, ky_sel=slice(0, n_ky), time_sel=valid_t), axis=2)

        # Align native cgyrodata kx length with the corresponding triad radial rows.
        n_kx_use = min(int(kx.size), int(radial_idx.size))
        if n_kx_use <= 0:
            print(f"No usable kx points for {label}")
            return None
        n_ky_use = min(int(ky.size), int(z.shape[1]))
        if n_ky_use <= 0:
            print(f"No usable ky points for {label}")
            return None
        kx = np.asarray(kx[:n_kx_use], dtype=float).reshape(-1)
        radial_idx = np.asarray(radial_idx[:n_kx_use], dtype=int).reshape(-1)
        valid_radial = (radial_idx >= 0) & (radial_idx < z.shape[0])
        kx = kx[valid_radial]
        radial_idx = radial_idx[valid_radial]
        n_kx_use = min(int(kx.size), int(radial_idx.size))
        kx = kx[:n_kx_use]
        radial_idx = radial_idx[:n_kx_use]
        ky = np.asarray(ky[:n_ky_use], dtype=float).reshape(-1)
        z = np.asarray(z[radial_idx, :n_ky_use], dtype=float)
        if t_ref is not None:
            t_ref = np.asarray(t_ref[radial_idx, :n_ky_use], dtype=float)

        finite_kx = np.isfinite(kx)
        finite_ky = np.isfinite(ky)
        if not np.any(finite_kx) or not np.any(finite_ky):
            print(f"No finite kx/ky points for {label}")
            return None
        kx = kx[finite_kx]
        z = z[finite_kx, :]
        if t_ref is not None:
            t_ref = t_ref[finite_kx, :]
        ky = ky[finite_ky]
        z = z[:, finite_ky]
        if t_ref is not None:
            t_ref = t_ref[:, finite_ky]

        kx_order = np.argsort(kx)
        ky_order = np.argsort(ky)
        kx = kx[kx_order]
        ky = ky[ky_order]
        z = z[np.ix_(kx_order, ky_order)]
        if t_ref is not None:
            t_ref = t_ref[np.ix_(kx_order, ky_order)]
            z = self._normalize_by_t_scale(z, t_ref, normalize_mode)

        # Build meshgrid for contourf
        kx_grid, ky_grid = np.meshgrid(kx, ky, indexing='ij')
        return kx_grid, ky_grid, np.asarray(z, dtype=float)

    def _compute_energy_balance_single_vs_time(
        self, data, label, ky_scan, quantity_key, normalize_mode="none"
    ):
        """
        Build a single-quantity time trace.

        This intentionally mirrors the manual triad sum used by the kx/ky views:
        select species, select ky, select channel, then sum over kx.
        """
        if not self._load_triad_only_if_needed(data, label):
            return None, None, None

        triad = np.asarray(data.triad)
        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Unsupported triad shape for {label}: {triad.shape}")
            return None, None, None

        _ri, _n_species, _n_radial, n_chan, n_n, n_t = triad.shape
        if n_chan < 8:
            print(f"Unsupported triad channels for {label}: {n_chan} (<8)")
            return None, None, None

        ky_match = self._single_plot_ky_match(data, n_n, label, ky_scan=ky_scan)
        if ky_match is None:
            return None, None, None
        ky_use, ky_value, ky_req = ky_match

        spec = self._parse_energy_balance_species()
        t_series = np.sum(self._triad_real_channel(triad, spec, 0, ky_sel=ky_use), axis=0)
        n_series = np.sum(self._triad_real_channel(triad, spec, 1, ky_sel=ky_use), axis=0)
        s_series = np.sum(self._triad_real_channel(triad, spec, 4, ky_sel=ky_use), axis=0)
        dr_series = np.sum(self._triad_real_channel(triad, spec, 5, ky_sel=ky_use), axis=0)
        dtheta_series = np.sum(self._triad_real_channel(triad, spec, 6, ky_sel=ky_use), axis=0)
        dc_series = np.sum(self._triad_real_channel(triad, spec, 7, ky_sel=ky_use), axis=0)
        dz_series = (
            dr_series
            + dtheta_series
            + dc_series
        )

        if quantity_key == 'T':
            y = t_series
        elif quantity_key == 'N':
            y = n_series
        elif quantity_key == 'T-N':
            y = t_series - n_series
        elif quantity_key == 'Dr':
            y = dr_series
        elif quantity_key == 'Dtheta':
            y = dtheta_series
        elif quantity_key == 'Dc':
            y = dc_series
        elif quantity_key == 'DZ':
            y = dz_series
        else:
            y = s_series
        if normalize_mode != "none" and quantity_key not in ['Dr', 'Dtheta', 'Dc', 'entropy']:
            y = self._normalize_by_t_scale(y, t_series, normalize_mode)

        t = np.asarray(getattr(data, 't', []), dtype=float).reshape(-1)
        if t.size <= 0:
            t = np.arange(int(n_t), dtype=float)
        n_use = min(int(t.size), int(np.asarray(y).size))
        if n_use <= 0:
            return None, None, None
        t = np.asarray(t[:n_use], dtype=float)
        y = np.asarray(y[:n_use], dtype=float)
        finite = np.isfinite(t) & np.isfinite(y)
        if not np.any(finite):
            return None, None, None
        return t[finite], y[finite], (ky_value, ky_use, ky_req)

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

        # idx5 -> channel 4: delta S
        ds = triad[0, :, :, 4, :n_ky, :][:, :, :, valid_t]  # [species, radial, n_ky, t_sel]
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
        - `T`, `N`, `T-N`, `Dr`, `Dtheta`, `Dc`, `DZ`, `entropy`.
        Supported axes:
        - `vs time`: direct time-trace.
        - `vs ky`: time-averaged spectrum summed over kx.
        - `vs kx`: time-averaged spectrum at the nearest stored `ky scan`.

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
        if qty_txt not in ['t', 'n', 't-n', 'dr', 'dtheta', 'dc', 'dz', 'entropy']:
            qty_txt = 't-n'
        normalize_mode = self._single_norm_mode() if qty_txt not in ['dr', 'dtheta', 'dc', 'dz', 'entropy'] else "none"

        avg_suffix = ""
        try:
            if t_indices is not None and len(t_indices) > 0:
                avg_suffix = self._format_avg_suffix(float(t_start), float(t_end), prefix="Avg")
        except Exception:
            avg_suffix = ""

        y_time = None
        y_label = ""
        qty_name = ""
        if qty_txt == 't':
            y_label = r"$\mathcal{T}^{NZ\rightarrow Z}$"
            qty_name = "T"
        elif qty_txt == 'n':
            y_label = r"$\mathcal{N}^{NZ\rightarrow Z}$"
            qty_name = "N"
        elif qty_txt == 'entropy':
            y_label = r"$S_g$"
            qty_name = "entropy"
        elif qty_txt == 'dr':
            y_label = r"$D_r$"
            qty_name = "Dr"
        elif qty_txt == 'dtheta':
            y_label = r"$D_\theta$"
            qty_name = "Dtheta"
        elif qty_txt == 'dc':
            y_label = r"$D_c$"
            qty_name = "Dc"
        elif qty_txt == 'dz':
            y_label = r"$D_Z$"
            qty_name = "DZ"
        else:
            y_label = r"$(\mathcal{T}-\mathcal{N})^{NZ\rightarrow Z}$"
            qty_name = "T-N"
        y_label_suffix, norm_suffix, norm_label = self._single_norm_display(normalize_mode)
        if y_label_suffix:
            y_label = y_label + y_label_suffix

        if mode_txt == 'vs ky':
            if qty_txt == 'entropy':
                self._plot_energy_balance_single_entropy_spectrum(
                    data, label, t_indices, avg_suffix=avg_suffix
                )
            else:
                ky, y_ky = self._compute_energy_balance_single_vs_ky(
                    data, label, None, qty_name, t_indices,
                    normalize_mode=normalize_mode,
                )
                if ky is None or y_ky is None:
                    return
                self.ax.plot(
                    ky, y_ky, '-o', linewidth=2.0, markersize=5.0,
                    label=f"{label}{norm_suffix}{avg_suffix}"
                )
                self.ax.set_xlabel(r"$k_y \rho_s$")
                self.ax.set_ylabel(y_label)
                self.ax.set_title(
                    f"Energy balance single plot ({qty_name}{norm_label}, vs ky){avg_suffix}"
                )
            return

        if mode_txt == 'vs kx':
            ky_scan = self._parse_energy_balance_ky_scan()
            kx, y_kx, ky_match = self._compute_energy_balance_single_vs_kx(
                data, label, ky_scan, qty_name, t_indices,
                normalize_mode=normalize_mode,
            )
            if kx is None or y_kx is None:
                return
            try:
                ky_value, ky_index, ky_request = ky_match
            except Exception:
                ky_value, ky_index, ky_request = np.nan, -1, ky_scan
            self.ax.plot(
                kx, y_kx, '-o', linewidth=2.0, markersize=5.0,
                label=f"{label}{norm_suffix} ky scan={ky_request:.6g} -> ky={ky_value:.6g} (idx={int(ky_index)}){avg_suffix}"
            )
            self.ax.axhline(0.0, color='0.45', linestyle=':', linewidth=0.9)
            self.ax.set_xlabel(r"$k_x \rho_s$")
            self.ax.set_ylabel(y_label)
            self.ax.set_title(
                f"Energy balance single plot ({qty_name}{norm_label}, vs kx, ky scan={ky_request:.6g}){avg_suffix}"
            )
            return

        if mode_txt == 'vs kxky':
            if qty_txt == 'dz':
                print(f"DZ is only supported for vs time/ky/kx in Single plot; skipped {label}.")
                return
            if qty_txt == 'entropy':
                print(f"entropy vs kxky not supported for {label}; use T, N, T-N, Dr, Dtheta or Dc.")
                return
            kxky_data = self._compute_energy_balance_single_vs_kxky(
                data, label, None, qty_name, t_indices,
                normalize_mode=normalize_mode,
            )
            if kxky_data is None:
                return
            kx_grid, ky_grid, z_data = kxky_data
            kx_axis = np.asarray(kx_grid[:, 0], dtype=float).reshape(-1)
            ky_axis = np.asarray(ky_grid[0, :], dtype=float).reshape(-1)
            z_plot = np.asarray(z_data, dtype=float)
            finite = z_plot[np.isfinite(z_plot)]
            vmax = 1.0
            if finite.size > 0:
                clip = float(np.nanpercentile(np.abs(finite), 99.5))
                if np.isfinite(clip) and clip > 0.0:
                    vmax = clip
                else:
                    raw_max = float(np.nanmax(np.abs(finite)))
                    if np.isfinite(raw_max) and raw_max > 0.0:
                        vmax = raw_max
            norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
            kx_edges = self._axis_cell_edges_from_centers(kx_axis)
            ky_edges = self._axis_cell_edges_from_centers(ky_axis)
            if kx_edges is None or ky_edges is None:
                print(f"Invalid kx/ky grid for {label} {qty_name} map.")
                return
            x_grid, y_grid = np.meshgrid(kx_edges, ky_edges)
            cs = self.ax.pcolormesh(
                x_grid,
                y_grid,
                np.ma.masked_invalid(z_plot.T),
                shading='flat',
                cmap=FULLT_DIVERGING_CMAP,
                norm=norm,
                rasterized=True,
            )
            cbar = self.fig.colorbar(cs, ax=self.ax, extend='both')
            cbar.set_label(y_label)
            self.ax.set_xlabel(r"$k_x \rho_s$")
            self.ax.set_ylabel(r"$k_y \rho_s$")
            self.ax.set_xlim(float(kx_edges[0]), float(kx_edges[-1]))
            self.ax.set_ylim(float(ky_edges[0]), float(ky_edges[-1]))
            self.ax.set_title(
                f"{qty_name}{norm_label} map{avg_suffix}"
            )
            return

        ky_scan = self._parse_energy_balance_ky_scan()
        t, y_time, ky_match = self._compute_energy_balance_single_vs_time(
            data, label, ky_scan, qty_name, normalize_mode=normalize_mode
        )
        if t is None or y_time is None:
            return
        try:
            ky_value, ky_index, ky_request = ky_match
        except Exception:
            ky_value, ky_index, ky_request = np.nan, -1, ky_scan

        n = min(t.size, y_time.size)
        if n <= 0:
            return
        t = t[:n]
        y_time = y_time[:n]
        valid_t = np.asarray(t_indices, dtype=int)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n)]
        valid_t = valid_t[np.isfinite(y_time[valid_t])]
        mean_suffix = ""
        if valid_t.size > 0:
            mean_val = float(np.mean(y_time[valid_t]))
            avg_inner = self._format_avg_range_from_axis(t, valid_t, prefix="Avg").strip(" ()")
            mean_suffix = f" ({avg_inner}, Mean: {mean_val:.2e})"
        line, = self.ax.plot(
            t,
            y_time,
            '-',
            linewidth=2.0,
            label=f"{label}{norm_suffix}{'' if qty_name == 'DZ' else ' ky='+format(ky_value, '.6g')}{mean_suffix}",
        )
        if valid_t.size > 0:
            self.ax.plot(
                [float(t[valid_t[0]]), float(t[valid_t[-1]])],
                [mean_val, mean_val],
                linestyle='--',
                color=line.get_color(),
                linewidth=1.5,
                label='_nolegend_',
            )
        self.ax.set_xlabel(r"$t\ (a/c_s)$")
        self.ax.set_ylabel(y_label)
        if qty_name == 'DZ':
            self.ax.set_title(f"Energy balance single plot ({qty_name}, vs time)")
        else:
            self.ax.set_title(
                f"Energy balance single plot ({qty_name}{norm_label}, vs time, ky scan={ky_request:.6g})"
            )

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

            # S denominator (direct from idx5):
            #   S(t) = sum_{a,p,n!=0} Re[idx5](a,p,n,t)
            s_all = np.real(f_complex[:, :, 4, :, :n_t_use])  # [species, radial, n, t]
            if s_all.shape[2] > 1:
                s_nonzonal = s_all[:, :, 1:, :]  # n != 0
            else:
                s_nonzonal = s_all[:, :, 0:0, :]
            if s_nonzonal.shape[2] > 0:
                s_total_idx5 = np.sum(s_nonzonal, axis=(0, 1, 2))  # [t]
            else:
                # Fallback for degenerate grids without n>0 entries.
                s_total_idx5 = np.sum(s_all, axis=(0, 1, 2))  # [t]

            d_t_avg = float(np.mean(d_t[t_idx]))
            d_n_avg = float(np.mean(d_n[t_idx]))
            e_t_avg = float(np.mean(e_t[t_idx]))
            e_n_avg = float(np.mean(e_n[t_idx]))
            s_avg = float(np.mean(s_total_idx5[t_idx]))
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
            # Keep legend concise in vs-2D mode: no time-window text in curve labels.

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

    def _fullt_integer_p_axis(self, data, n_radial):
        """Return integer radial mode indices matching FULLT's target-kx axis."""
        p = np.asarray(getattr(data, 'p', []), dtype=float).reshape(-1)
        if p.size == n_radial and np.all(np.isfinite(p)):
            p_int = np.rint(p).astype(int)
            if np.allclose(p, p_int, rtol=0.0, atol=1.0e-8):
                return p_int
        return np.arange(n_radial, dtype=int) - (n_radial // 2)

    def _construct_fullt_asym_from_fullt(self, full_t, data, label):
        """Construct FULLT_ASYM(k') = FULLT(k') + FULLT(k-k') in plotting layout."""
        full_t = np.asarray(full_t)
        if full_t.ndim != 6 or full_t.shape[0] != 2:
            print(f"Cannot construct FULLT_ASYM for {label}: unsupported FULLT shape {full_t.shape}.")
            return None, None

        _ri, n_source, n_radial, n_signed, _n_channel, _n_time = full_t.shape
        if n_signed != 2 * n_source - 1:
            print(
                f"Cannot construct FULLT_ASYM for {label}: signed ky dimension "
                f"{n_signed} does not match source ky dimension {n_source}."
            )
            return None, None

        p_axis = self._fullt_integer_p_axis(data, n_radial)
        p_to_ir = {int(p): int(i) for i, p in enumerate(p_axis)}
        ky_shift = n_source - 1
        ky_min = -ky_shift
        ky_max = ky_shift
        asym = np.zeros_like(full_t)
        stats = {
            'constructed': 0,
            'skipped_ky': 0,
            'skipped_radial_partner': 0,
            'skipped_radial_nyquist': 0,
        }

        p_min = int(np.min(p_axis)) if p_axis.size else 0
        implicit_p_nyquist = -p_min

        for source_ky in range(n_source):
            for iky in range(n_signed):
                ky_prime = iky - ky_shift
                ky_second = source_ky - ky_prime
                if ky_second < ky_min or ky_second > ky_max:
                    stats['skipped_ky'] += n_radial
                    continue
                iky_second = ky_second + ky_shift

                for ir, p_prime in enumerate(p_axis):
                    p_second = -int(p_prime)
                    ir_second = p_to_ir.get(p_second)
                    if ir_second is None:
                        if p_second == implicit_p_nyquist:
                            stats['skipped_radial_nyquist'] += 1
                        else:
                            stats['skipped_radial_partner'] += 1
                        continue

                    asym[:, source_ky, ir, iky, :, :] = (
                        full_t[:, source_ky, ir, iky, :, :]
                        + full_t[:, source_ky, ir_second, iky_second, :, :]
                    )
                    stats['constructed'] += 1

        return asym, stats

    def _load_fullt_if_needed(
        self,
        data,
        label,
        asym=False,
        source_ky_value=None,
        source_kx_value=None,
        time_indices=None,
    ):
        """Load and normalize FULLT data to [ri,source_ky,target_kx,target_ky,channel,time]."""
        diag_name = "FULLT_ASYM" if asym else "FULLT"
        file_suffix = ".cgyro.fullt_asym" if asym else ".cgyro.fullt"
        cache_attr = "fullt_asym" if asym else "fullt"
        reader_attr = "full_t_asym" if asym else "full_t"
        channel_attr = cache_attr + "_n_channel"
        source_attr = cache_attr + "_source"
        case_token_attr = cache_attr + "_case_token"
        cache_version_attr = cache_attr + "_cache_version"
        cache_version = 4
        try:
            case_token = self._resolve_case_dir(data)
        except Exception:
            case_token = getattr(data, 'dir', None) or getattr(data, 'path', None) or ""
        case_token = str(case_token)
        # Current CGYRO writes real FULLT channels:
        #   full_t(target_kx, target_ky, channel, source_ky_local, source_kx)
        # channel 0 stores Re(T), channel 1 stores Im(T) when enabled.  Older
        # files used a complex binary container for the same logical channels.
        # The GUI keeps its historical 6D plotting layout by selecting the
        # requested source_kx slice when a full source-kx file is present.
        # Note: cgyro_write_distributed_breal appends local source-ky blocks by
        # toroidal rank.  The global binary order used by pygacode/verifiers is
        # therefore (..., source_kx, global_source_ky, time).
        n_n = int(getattr(data, 'n_n', 0))
        n_radial = int(getattr(data, 'n_radial', 0))

        def _clear_fullt_cache(reason):
            attrs = [
                cache_attr,
                channel_attr,
                source_attr,
                case_token_attr,
                cache_version_attr,
                cache_attr + "_source_ky_axis",
                cache_attr + "_source_ky_full_index",
                cache_attr + "_n_source_ky_full",
                cache_attr + "_source_kx",
                cache_attr + "_source_kx_index",
                cache_attr + "_n_source_kx",
                cache_attr + "_time_indices",
            ]
            for attr in attrs:
                try:
                    if hasattr(data, attr):
                        delattr(data, attr)
                except Exception:
                    pass
            if reason:
                print(f"{diag_name} cache for {label}: {reason}; reloading.")

        def _real_dtype():
            try:
                return np.dtype(getattr(data, "BYTE", "float64"))
            except Exception:
                return np.dtype("float64")

        def _case_file_path():
            try:
                case_dir = self._resolve_case_dir(data)
            except Exception:
                case_dir = getattr(data, "dir", None) or getattr(data, "path", None) or ""
            case_dir = str(case_dir or "")
            if not case_dir:
                return None
            candidates = [
                os.path.join(case_dir, "bin" + file_suffix),
                os.path.join(case_dir, "bin", "bin" + file_suffix),
            ]
            for path in candidates:
                if os.path.isfile(path):
                    return path
            return None

        def _case_time_count_from_file():
            """Read current case time count from out.cgyro.time on disk."""
            try:
                case_dir = self._resolve_case_dir(data)
            except Exception:
                case_dir = getattr(data, "dir", None) or getattr(data, "path", None) or ""
            case_dir = str(case_dir or "")
            if not case_dir:
                return 0
            path = os.path.join(case_dir, "out.cgyro.time")
            if not os.path.isfile(path):
                return 0
            try:
                t_file = np.fromfile(path, dtype=float, sep=" ")
                return int(np.asarray(t_file).size)
            except Exception:
                return 0

        def _input_fullt_hints():
            """
            Build FULLT shape candidates from current-case input.cgyro.

            Do not use data.t or data.n_time here: those can be stale after
            switching cases in a long-running GUI session.
            """
            try:
                case_dir = self._resolve_case_dir(data)
            except Exception:
                case_dir = getattr(data, "dir", None) or getattr(data, "path", None) or ""
            try:
                scalars = self._read_input_cgyro_scalars(case_dir)
            except Exception:
                scalars = {}

            try:
                real_only = int(round(float(scalars.get("FULL_T_REAL_ONLY", 0))))
            except Exception:
                real_only = 0
            channel_hint = 1 if real_only == 1 else 2

            try:
                kx0 = int(round(float(scalars.get("FULL_T_KX0", 0))))
            except Exception:
                kx0 = 0
            source_kx_hint = 1 if kx0 == 1 else n_radial

            channel_candidates = [channel_hint, 2, 1]
            source_kx_candidates = [source_kx_hint]
            if n_radial > 1:
                source_kx_candidates.extend([n_radial, 1])
            else:
                source_kx_candidates.append(1)

            channel_candidates = list(dict.fromkeys([
                int(x) for x in channel_candidates if int(x) in (1, 2)
            ]))
            source_kx_candidates = list(dict.fromkeys([
                int(x) for x in source_kx_candidates if int(x) > 0
            ]))
            return channel_candidates, source_kx_candidates

        def _source_ky_axis(n_source):
            axis = np.asarray(getattr(data, 'full_t_source_ky', []), dtype=float).reshape(-1)
            if axis.size >= n_source and np.all(np.isfinite(axis[:n_source])):
                return axis[:n_source]
            axis = np.asarray(getattr(data, 'ky', []), dtype=float).reshape(-1)
            if axis.size >= n_source and np.all(np.isfinite(axis[:n_source])):
                return axis[:n_source]
            return np.arange(n_source, dtype=float)

        def _select_source_ky_index(n_source):
            if n_source <= 1:
                return 0
            if source_ky_value is None:
                return 0
            axis = _source_ky_axis(n_source)
            try:
                return int(np.nanargmin(np.abs(axis - float(source_ky_value))))
            except Exception:
                return 0

        def _source_kx_axis(n_source_kx):
            if n_source_kx == 1:
                return np.asarray([0.0], dtype=float)

            p = np.asarray(getattr(data, 'p', []), dtype=float).reshape(-1)
            length = float(getattr(data, 'length', 0.0) or 0.0)
            if p.size >= n_source_kx and np.all(np.isfinite(p[:n_source_kx])):
                if np.isfinite(length) and abs(length) > 1.0e-12:
                    return 2.0 * np.pi * p[:n_source_kx] / length

            axis = np.asarray(getattr(data, 'kxnorm', getattr(data, 'kx', [])), dtype=float).reshape(-1)
            if axis.size >= n_source_kx and np.all(np.isfinite(axis[:n_source_kx])):
                return axis[:n_source_kx]

            axis = np.asarray(getattr(data, 'full_t_source_kx', []), dtype=float).reshape(-1)
            if axis.size >= n_source_kx and np.all(np.isfinite(axis[:n_source_kx])):
                return axis[:n_source_kx]

            p = np.asarray(getattr(data, 'p', []), dtype=float).reshape(-1)
            if p.size >= n_source_kx and np.all(np.isfinite(p[:n_source_kx])):
                return p[:n_source_kx]
            return np.arange(n_source_kx, dtype=float) - (n_source_kx // 2)

        def _select_source_kx_index(n_source_kx):
            if n_source_kx <= 1:
                return 0
            axis = _source_kx_axis(n_source_kx)
            try:
                target = 0.0 if source_kx_value is None else float(source_kx_value)
                return int(np.nanargmin(np.abs(axis - target)))
            except Exception:
                return n_source_kx // 2

        def _record_source_kx(n_source_kx, source_kx_idx):
            axis = _source_kx_axis(n_source_kx)
            try:
                setattr(data, cache_attr + "_source_kx", float(axis[source_kx_idx]))
                setattr(data, cache_attr + "_source_kx_index", int(source_kx_idx))
                setattr(data, cache_attr + "_n_source_kx", int(n_source_kx))
            except Exception:
                pass

        def _requested_time_indices():
            if time_indices is None:
                return None
            try:
                idx = np.asarray(time_indices, dtype=int).reshape(-1)
            except Exception:
                return None
            idx = idx[idx >= 0]
            if idx.size <= 0:
                return None
            return idx

        def _infer_file_shape(n_elem):
            block_no_time = n_radial * (2 * n_n - 1) * n_n
            if block_no_time <= 0:
                return 0, 0, 0, False
            channel_candidates, source_kx_candidates = _input_fullt_hints()
            n_time_file = _case_time_count_from_file()

            if n_time_file > 0:
                for scale in (1, 2):
                    for n_source_try in source_kx_candidates:
                        for n_channel_try in channel_candidates:
                            denom = scale * block_no_time * n_source_try * n_channel_try
                            if denom > 0 and n_elem == denom * n_time_file:
                                return int(n_channel_try), int(n_time_file), int(n_source_try), bool(scale == 2)

            for scale in (1, 2):
                for n_source_try in source_kx_candidates:
                    for n_channel_try in channel_candidates:
                        denom = scale * block_no_time * n_source_try * n_channel_try
                        if denom > 0 and n_elem % denom == 0:
                            n_time_try = int(n_elem // denom)
                            if n_time_try > 0:
                                return int(n_channel_try), int(n_time_try), int(n_source_try), bool(scale == 2)
            return 0, 0, 0, False

        def _load_direct_source_slice():
            path = _case_file_path()
            if path is None or n_n <= 0 or n_radial <= 0:
                return None
            dtype = _real_dtype()
            try:
                n_elem = int(os.path.getsize(path) // dtype.itemsize)
            except Exception:
                return None
            n_channel_file, n_time, n_source_kx, legacy_complex = _infer_file_shape(n_elem)
            if n_channel_file not in (1, 2) or n_time <= 0 or n_source_kx <= 0:
                return None

            source_ky_idx_full = _select_source_ky_index(n_n)
            source_kx_idx = _select_source_kx_index(n_source_kx)
            _record_source_kx(n_source_kx, source_kx_idx)
            source_ky_axis_full = _source_ky_axis(n_n)
            source_ky_used = float(source_ky_axis_full[source_ky_idx_full]) if source_ky_axis_full.size > source_ky_idx_full else float(source_ky_idx_full)
            n_signed = 2 * n_n - 1
            time_idx = _requested_time_indices()
            if time_idx is None:
                time_idx = np.arange(n_time, dtype=int)
            else:
                time_idx = time_idx[time_idx < n_time]
                if time_idx.size <= 0:
                    time_idx = np.arange(n_time, dtype=int)
            n_time_read = int(time_idx.size)
            try:
                mm = np.memmap(path, dtype=dtype, mode="r", shape=(n_elem,))
                arr = np.zeros((2, 1, n_radial, n_signed, 1, n_time_read), dtype=float)
                if legacy_complex:
                    raw = np.reshape(
                        mm,
                        (2, n_radial, n_signed, n_channel_file, n_source_kx, n_n, n_time),
                        order="F",
                    )
                    src = raw[:, :, :, 0, source_kx_idx, source_ky_idx_full, :]
                    arr[:, 0, :, :, 0, :] = np.asarray(src[:, :, :, time_idx], dtype=float)
                else:
                    raw = np.reshape(
                        mm,
                        (n_radial, n_signed, n_channel_file, n_source_kx, n_n, n_time),
                        order="F",
                    )
                    src = raw[:, :, 0, source_kx_idx, source_ky_idx_full, :]
                    arr[0, 0, :, :, 0, :] = np.asarray(src[:, :, time_idx], dtype=float)
            except Exception as e:
                print(f"{diag_name} direct slice read failed for {label}: {e}")
                return None

            setattr(data, cache_attr, arr)
            setattr(data, channel_attr, 1)
            setattr(data, case_token_attr, case_token)
            setattr(data, cache_version_attr, cache_version)
            setattr(data, cache_attr + "_source_ky_axis", np.asarray([source_ky_used], dtype=float))
            setattr(data, cache_attr + "_source_ky_full_index", int(source_ky_idx_full))
            setattr(data, cache_attr + "_n_source_ky_full", int(n_n))
            setattr(data, cache_attr + "_time_indices", np.asarray(time_idx, dtype=int))
            storage = "legacy complex" if legacy_complex else "real"
            source_kx_msg = ""
            try:
                source_kx_msg = f", source_kx={getattr(data, cache_attr + '_source_kx'):.6g}"
            except Exception:
                pass
            setattr(
                data,
                source_attr,
                f"direct slice {os.path.basename(path)} ({storage}, source_ky={source_ky_used:.6g}, "
                f"source_kx_count={n_source_kx}{source_kx_msg}, channel=0, n_time={n_time_read})",
            )
            print(
                f"{diag_name} for {label}: read direct source slice from {path} "
                f"({storage}, source_ky={source_ky_used:.6g}, source_kx_count={n_source_kx}{source_kx_msg}, "
                f"channel=0, n_time={n_time_read}/{n_time})."
            )
            return arr

        def _accept_existing(arr):
            arr = np.asarray(arr)
            if arr.ndim not in (6, 7) or arr.shape[0] != 2:
                return None
            if n_n > 0 and n_radial > 0:
                n_signed = 2 * n_n - 1
                # New pygacode layout.
                if arr.ndim == 6 and arr.shape[1] == n_n and arr.shape[2] == n_radial and arr.shape[3] == n_signed:
                    return arr
                # New pygacode layout with source_kx:
                # [ri,source_ky,source_kx,target_kx,target_ky,channel,time].
                if arr.ndim == 7 and arr.shape[1] == n_n and arr.shape[3] == n_radial and arr.shape[4] == n_signed:
                    source_kx_idx = _select_source_kx_index(int(arr.shape[2]))
                    _record_source_kx(int(arr.shape[2]), source_kx_idx)
                    return arr[:, :, source_kx_idx, :, :, :, :]
                # Raw Fortran reshape layout.
                if arr.ndim == 6 and arr.shape[1] == n_radial and arr.shape[2] == n_signed and arr.shape[4] == n_n:
                    # Convert [ri,target_kx,target_ky,channel,source_ky,time]
                    # into [ri,source_ky,target_kx,target_ky,channel,time].
                    return np.transpose(arr, (0, 4, 1, 2, 3, 5))
                if arr.ndim == 7 and arr.shape[1] == n_radial and arr.shape[2] == n_signed and arr.shape[5] == n_n:
                    source_kx_idx = _select_source_kx_index(int(arr.shape[4]))
                    _record_source_kx(int(arr.shape[4]), source_kx_idx)
                    return np.transpose(arr[:, :, :, :, source_kx_idx, :, :], (0, 4, 1, 2, 3, 5))
            return None

        if hasattr(data, cache_attr):
            arr = _accept_existing(getattr(data, cache_attr))
            if arr is not None:
                cached_token = str(getattr(data, case_token_attr, case_token))
                cached_version = int(getattr(data, cache_version_attr, 0))
                if cached_version != cache_version:
                    _clear_fullt_cache(
                        f"stale cache version {cached_version}, expected {cache_version}"
                    )
                elif cached_token and case_token and cached_token != case_token:
                    _clear_fullt_cache("was tied to another case")
                elif source_ky_value is not None and hasattr(data, cache_attr + "_source_ky_axis"):
                    ky_cached = np.asarray(getattr(data, cache_attr + "_source_ky_axis"), dtype=float).reshape(-1)
                    kx_cached = np.asarray([getattr(data, cache_attr + "_source_kx", np.nan)], dtype=float).reshape(-1)
                    time_cached = np.asarray(getattr(data, cache_attr + "_time_indices", []), dtype=int).reshape(-1)
                    try:
                        ky_full = _source_ky_axis(n_n)
                        ky_target_idx = _select_source_ky_index(n_n)
                        ky_target = float(ky_full[ky_target_idx]) if ky_full.size > ky_target_idx else float(source_ky_value)
                        ky_delta = abs(float(ky_cached[0]) - ky_target)
                    except Exception:
                        ky_delta = 0.0
                    try:
                        n_source_kx_cached = int(getattr(data, cache_attr + "_n_source_kx", 1))
                        kx_axis = _source_kx_axis(n_source_kx_cached)
                        kx_target = 0.0 if source_kx_value is None else float(source_kx_value)
                        kx_idx_target = int(np.nanargmin(np.abs(kx_axis - kx_target))) if kx_axis.size > 0 else 0
                        kx_used_target = float(kx_axis[kx_idx_target]) if kx_axis.size > kx_idx_target else kx_target
                        kx_delta = abs(float(kx_cached[0]) - kx_used_target)
                    except Exception:
                        kx_delta = 0.0
                    time_requested = _requested_time_indices()
                    time_delta = False
                    if time_requested is not None and time_cached.size > 0:
                        time_delta = (
                            time_cached.size != time_requested.size
                            or np.any(time_cached != time_requested)
                        )
                    if ky_cached.size > 0 and (ky_delta > 1.0e-12 or kx_delta > 1.0e-12 or time_delta):
                        _clear_fullt_cache(
                            f"used source ky={float(ky_cached[0]):.6g}, "
                            f"source kx={float(kx_cached[0]):.6g}; reloading for "
                            f"source ky={ky_target:.6g}, source kx={kx_used_target:.6g}"
                        )
                    else:
                        arr = np.array(arr, copy=False)
                        setattr(data, cache_attr, arr)
                        setattr(data, channel_attr, int(arr.shape[4]))
                        setattr(data, case_token_attr, case_token)
                        setattr(data, cache_version_attr, cache_version)
                        source = getattr(data, source_attr, "cache")
                        print(f"{diag_name} for {label}: using {source}.")
                        return True
                else:
                    arr = np.array(arr, copy=False)
                    setattr(data, cache_attr, arr)
                    setattr(data, channel_attr, int(arr.shape[4]))
                    setattr(data, case_token_attr, case_token)
                    setattr(data, cache_version_attr, cache_version)
                    source = getattr(data, source_attr, "cache")
                    print(f"{diag_name} for {label}: using {source}.")
                    return True

        if n_n <= 0 or n_radial <= 0:
            print(f"Cannot load {diag_name} for {label}: missing n_n/n_radial metadata.")
            return False

        # Avoid `data.extract()` for FULLT whenever possible: it reads the
        # complete dense file into memory.  FULLT transfer maps only need one
        # source-ky/source-kx slice, so direct memmap keeps large all-kx cases
        # from exhausting RAM.
        if source_ky_value is not None:
            arr = _load_direct_source_slice()
            if arr is not None:
                return True

        if hasattr(data, reader_attr):
            arr = _accept_existing(getattr(data, reader_attr))
            if arr is not None:
                # Copy the reader attribute into a GUI-owned cache.  This avoids
                # subtle aliasing if a reader implementation reuses internal
                # buffers while several cases are inspected in one GUI session.
                arr = np.array(arr, copy=True)
                setattr(data, cache_attr, arr)
                setattr(data, channel_attr, int(arr.shape[4]))
                setattr(data, case_token_attr, case_token)
                setattr(data, cache_version_attr, cache_version)
                setattr(data, source_attr, f"reader attribute {reader_attr}")
                print(f"{diag_name} for {label}: read from reader attribute {reader_attr}.")
                return True

        if not hasattr(data, 'extract'):
            print(f"Cannot load {diag_name} for {label}: data object has no extract() method.")
            return False

        try:
            _tmsg, fmt, raw = data.extract(file_suffix)
        except Exception as e:
            print(f"{diag_name} read failed for {label}: {e}")
            return False

        if fmt == 'null':
            if asym:
                print(f"No native FULLT_ASYM data for {label}; constructing from FULLT.")
                if not self._load_fullt_if_needed(data, label, asym=False, source_ky_value=source_ky_value):
                    print(f"Cannot construct FULLT_ASYM for {label}: FULLT is unavailable.")
                    return False
                base = _accept_existing(getattr(data, 'fullt', None))
                if base is None:
                    print(f"Cannot construct FULLT_ASYM for {label}: cached FULLT has invalid shape.")
                    return False
                arr, stats = self._construct_fullt_asym_from_fullt(base, data, label)
                if arr is None:
                    return False
                arr = np.array(arr, copy=True)
                setattr(data, cache_attr, arr)
                setattr(data, channel_attr, int(arr.shape[4]))
                setattr(data, case_token_attr, case_token)
                setattr(data, cache_version_attr, cache_version)
                setattr(data, source_attr, "constructed from FULLT")
                skipped = (
                    stats['skipped_ky']
                    + stats['skipped_radial_partner']
                    + stats['skipped_radial_nyquist']
                )
                print(
                    f"FULLT_ASYM for {label}: constructed from FULLT "
                    f"(constructed cells={stats['constructed']}, skipped cells={skipped}, "
                    f"Nyquist skipped={stats['skipped_radial_nyquist']})."
                )
                return True
            print(f"No {diag_name} data available for {label}. Need bin/out{file_suffix}.")
            return False

        raw = np.asarray(raw)
        if raw.size <= 0:
            print(f"Empty {diag_name} data for {label}.")
            return False
        if np.iscomplexobj(raw):
            raw = np.column_stack((np.real(raw).reshape(-1), np.imag(raw).reshape(-1))).reshape(-1)
        raw = np.asarray(raw, dtype=float).reshape(-1)

        block_no_time = n_radial * (2 * n_n - 1) * n_n
        if block_no_time <= 0:
            print(f"Invalid {diag_name} dimensions for {label}.")
            return False

        n_channel = 0
        n_time = 0
        n_source_kx = 0
        legacy_complex = False
        channel_candidates, source_kx_candidates = _input_fullt_hints()
        n_time_file = _case_time_count_from_file()

        if n_time_file > 0:
            for scale in (1, 2):
                if n_channel > 0:
                    break
                for n_source_try in source_kx_candidates:
                    if n_channel > 0:
                        break
                    for n_channel_try in channel_candidates:
                        denom = scale * block_no_time * n_source_try * n_channel_try
                        if denom > 0 and raw.size == denom * n_time_file:
                            n_channel = int(n_channel_try)
                            n_time = int(n_time_file)
                            n_source_kx = int(n_source_try)
                            legacy_complex = (scale == 2)
                            break

        if n_channel <= 0:
            for scale in (1, 2):
                if n_channel > 0:
                    break
                for n_source_try in source_kx_candidates:
                    if n_channel > 0:
                        break
                    for n_channel_try in channel_candidates:
                        denom = scale * block_no_time * n_source_try * n_channel_try
                        if denom > 0 and raw.size % denom == 0:
                            n_time_try = int(raw.size // denom)
                            if n_time_try > 0:
                                n_channel = int(n_channel_try)
                                n_time = int(n_time_try)
                                n_source_kx = int(n_source_try)
                                legacy_complex = (scale == 2)
                                break

        if n_channel not in (1, 2) or n_time <= 0 or n_source_kx <= 0:
            print(
                f"Cannot infer {diag_name} shape for {label}: raw_size={raw.size}, "
                f"n_radial={n_radial}, n_n={n_n}."
            )
            return False

        scale = 2 if legacy_complex else 1
        nd = scale * block_no_time * n_source_kx * n_channel * n_time
        if raw.size < nd:
            print(f"{diag_name} data is smaller than inferred shape for {label}.")
            return False
        source_kx_idx = _select_source_kx_index(n_source_kx)
        _record_source_kx(n_source_kx, source_kx_idx)
        try:
            # FULLT is written in Fortran order with target_kx fastest among
            # the physical dimensions.  Reshaping in C order would swap axes and
            # make the target/source-ky verification fail even if the binary is
            # correct.
            if legacy_complex:
                raw7 = np.reshape(
                    raw[:nd],
                    (2, n_radial, 2 * n_n - 1, n_channel, n_source_kx, n_n, n_time),
                    order='F',
                )
                arr = np.array(np.transpose(raw7[:, :, :, :, source_kx_idx, :, :], (0, 4, 1, 2, 3, 5)), copy=True)
            else:
                raw6 = np.reshape(
                    raw[:nd],
                    (n_radial, 2 * n_n - 1, n_channel, n_source_kx, n_n, n_time),
                    order='F',
                )
                arr = np.zeros((2, n_n, n_radial, 2 * n_n - 1, n_channel, n_time), dtype=float)
                arr[0, :, :, :, :, :] = np.transpose(raw6[:, :, :, source_kx_idx, :, :], (3, 0, 1, 2, 4))
        except Exception as e:
            print(f"{diag_name} reshape failed for {label}: {e}")
            return False

        setattr(data, cache_attr, arr)
        setattr(data, channel_attr, int(n_channel))
        setattr(data, case_token_attr, case_token)
        setattr(data, cache_version_attr, cache_version)
        storage = "legacy complex" if legacy_complex else "real"
        source_kx_msg = ""
        try:
            source_kx_msg = f", source_kx={getattr(data, cache_attr + '_source_kx'):.6g}"
        except Exception:
            pass
        setattr(data, source_attr, f"{fmt}{file_suffix} ({storage}, source_kx_count={n_source_kx}{source_kx_msg})")
        print(f"{diag_name} for {label}: read from {fmt}{file_suffix} ({storage}, source_kx_count={n_source_kx}{source_kx_msg}).")
        return True

    def _energy_fullt_target_kx_axis(self, data, n_radial, label):
        """Build the full target-kx axis, including the leftmost Nyquist storage bin."""
        p = np.asarray(getattr(data, 'p', []), dtype=float).reshape(-1)
        length = float(getattr(data, 'length', 0.0))
        if p.size == n_radial and np.isfinite(length) and abs(length) > 1.0e-12:
            return 2.0 * np.pi * p / length, np.arange(n_radial, dtype=int)

        kx = np.asarray(getattr(data, 'kx', []), dtype=float).reshape(-1)
        if kx.size == n_radial:
            return kx, np.arange(n_radial, dtype=int)

        if np.isfinite(length) and abs(length) > 1.0e-12:
            dkx = 2.0 * np.pi / length
        else:
            dkx = 1.0
        p_index = np.arange(n_radial, dtype=float) - (n_radial // 2)
        if kx.size == n_radial - 1 and n_radial > 1:
            test = p_index[1:] * dkx
            if np.allclose(test[:kx.size], kx, rtol=1.0e-5, atol=1.0e-12):
                return p_index * dkx, np.arange(n_radial, dtype=int)
        if kx.size > 0:
            print(
                f"Warning: FULLT target-kx axis for {label} inferred from n_radial; "
                f"metadata kx length is {kx.size}, FULLT radial length is {n_radial}."
            )
        return p_index * dkx, np.arange(n_radial, dtype=int)

    def _fullt_target_ky_axis(self, data, n_n, n_signed):
        """Build signed target-ky axis matching FULLT target_ky storage order."""
        ky_native = np.asarray(getattr(data, 'full_t_target_ky', []), dtype=float).reshape(-1)
        if ky_native.size == n_signed and np.all(np.isfinite(ky_native)):
            return ky_native

        # Source ky follows the physical sign stored in CGYRO grids.  Target ky
        # scans signed mode indices, so multiply those indices by the signed
        # source-ky step rather than forcing ky positive.
        ky_source = np.asarray(getattr(data, 'kynorm', getattr(data, 'ky', [])), dtype=float).reshape(-1)
        if ky_source.size >= n_n and n_n > 0:
            ky_source = ky_source[:n_n]
            dky = float(ky_source[1]) if n_n > 1 else 1.0
            ky_axis = np.arange(-(n_n - 1), n_n, dtype=float) * dky
        elif n_n > 1:
            ky_axis = np.arange(-(n_n - 1), n_n, dtype=float)
        else:
            ky_axis = np.asarray([0.0], dtype=float)

        if ky_axis.size != n_signed:
            ky_axis = np.arange(-(n_signed // 2), n_signed - (n_signed // 2), dtype=float)
        return ky_axis

    def _fullt_source_ky_axis(self, data, n_source, preferred_attr=None):
        """Build signed source-ky axis for FULLT source_ky index selection."""
        attrs = []
        if preferred_attr:
            attrs.append(str(preferred_attr))
        attrs.extend(["fullt_source_ky_axis", "fullt_asym_source_ky_axis"])
        seen = set()
        for attr in attrs:
            if attr in seen:
                continue
            seen.add(attr)
            ky_slice = np.asarray(getattr(data, attr, []), dtype=float).reshape(-1)
            if ky_slice.size >= n_source and np.all(np.isfinite(ky_slice[:n_source])):
                return ky_slice[:n_source]
        ky_native = np.asarray(getattr(data, 'full_t_source_ky', []), dtype=float).reshape(-1)
        if ky_native.size >= n_source and np.all(np.isfinite(ky_native[:n_source])):
            return ky_native[:n_source]

        ky = np.asarray(getattr(data, 'kynorm', getattr(data, 'ky', [])), dtype=float).reshape(-1)
        ky = np.asarray(ky, dtype=float).reshape(-1)
        if ky.size >= n_source:
            return ky[:n_source]
        return np.arange(n_source, dtype=float)

    def _fullt_asym_pair_symmetry_error(self, full_t_asym, data, source_ky_idx, valid_t):
        """Return pair-symmetry diagnostics for FULLT_ASYM at one source ky."""
        arr = np.asarray(full_t_asym)
        if arr.ndim != 6 or arr.shape[0] != 2:
            return None

        _ri, n_source, n_radial, n_signed, n_channel, n_t = arr.shape
        if (
            source_ky_idx < 0
            or source_ky_idx >= n_source
            or n_signed != 2 * n_source - 1
            or n_channel <= 0
            or n_t <= 0
        ):
            return None

        valid_t = np.asarray(valid_t, dtype=int).reshape(-1)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size <= 0:
            valid_t = np.arange(n_t, dtype=int)

        p_axis = self._fullt_integer_p_axis(data, n_radial)
        p_to_ir = {int(p): int(i) for i, p in enumerate(p_axis)}
        ky_shift = n_source - 1
        ky_min = -ky_shift
        ky_max = ky_shift

        sum_abs2 = 0.0
        sum_ref2 = 0.0
        max_abs = 0.0
        n_value = 0
        n_pair = 0
        for iky in range(n_signed):
            ky_prime = iky - ky_shift
            ky_second = int(source_ky_idx) - ky_prime
            if ky_second < ky_min or ky_second > ky_max:
                continue
            iky_second = ky_second + ky_shift
            for ir, p_prime in enumerate(p_axis):
                ir_second = p_to_ir.get(-int(p_prime))
                if ir_second is None:
                    continue
                a = np.asarray(arr[0, source_ky_idx, ir, iky, 0, valid_t], dtype=float)
                b = np.asarray(arr[0, source_ky_idx, ir_second, iky_second, 0, valid_t], dtype=float)
                mask = np.isfinite(a) & np.isfinite(b)
                if not np.any(mask):
                    continue
                diff = a[mask] - b[mask]
                ref = 0.5 * (np.abs(a[mask]) + np.abs(b[mask]))
                sum_abs2 += float(np.sum(diff * diff))
                sum_ref2 += float(np.sum(ref * ref))
                max_abs = max(max_abs, float(np.max(np.abs(diff))))
                n_value += int(diff.size)
                n_pair += 1

        if n_value <= 0:
            return None
        rms_abs = float(np.sqrt(sum_abs2 / max(n_value, 1)))
        rms_ref = float(np.sqrt(sum_ref2 / max(n_value, 1)))
        rms_rel = rms_abs / max(rms_ref, 1.0e-300)
        return {
            "rms_abs": rms_abs,
            "rms_rel": rms_rel,
            "max_abs": max_abs,
            "n_pair": n_pair,
            "n_value": n_value,
        }

    def _fallback_construct_fullt_asym_for_plot(self, data, label):
        """Replace cached FULLT_ASYM with a plot-local construction from FULLT."""
        base_existing = np.asarray(getattr(data, "fullt", None))
        if base_existing.ndim != 6 or base_existing.shape[0] != 2 or base_existing.shape[1] <= 1:
            print(
                f"FULLT_ASYM fallback for {label}: no complete cached FULLT is available; "
                "skipping construction to avoid loading the full dense FULLT file."
            )
            return None
        base = base_existing
        arr, stats = self._construct_fullt_asym_from_fullt(base, data, label)
        if arr is None:
            return None
        arr = np.array(arr, copy=True)
        setattr(data, "fullt_asym", arr)
        setattr(data, "fullt_asym_n_channel", int(arr.shape[4]))
        setattr(data, "fullt_asym_source", "constructed from FULLT after native pair-symmetry check")
        try:
            setattr(data, "fullt_asym_case_token", str(self._resolve_case_dir(data)))
        except Exception:
            pass
        skipped = (
            stats["skipped_ky"]
            + stats["skipped_radial_partner"]
            + stats["skipped_radial_nyquist"]
        )
        print(
            f"FULLT_ASYM for {label}: plotting constructed data from FULLT "
            f"(constructed cells={stats['constructed']}, skipped cells={skipped})."
        )
        return arr

    def _plot_energy_balance_fullt(self, data, label, t_indices, t_start, t_end):
        """Plot FULLT for one fixed source ky over scanned target kx/target ky."""
        try:
            use_asym = bool(self.energy_balance_transfer_asym_var.get())
        except Exception:
            use_asym = False
        try:
            normalize_by_max_t = bool(self.energy_balance_transfer_norm_max_var.get())
        except Exception:
            normalize_by_max_t = False
        diag_name = "FULLT asym" if use_asym else "FULLT"
        data_attr = "fullt_asym" if use_asym else "fullt"
        transfer_symbol = r"T_{\mathrm{asym}}^{\Phi}" if use_asym else r"T^{\Phi}"
        try:
            view_txt = str(self.energy_balance_transfer_xaxis_var.get()).strip().lower()
        except Exception:
            view_txt = "vs kxky"
        view_mode = "vs kx" if view_txt == "vs kx" else "vs kxky"

        try:
            ky_sel = float(str(self.energy_balance_transfer_ky_var.get()).strip())
        except Exception:
            print(f"Invalid {diag_name} fixed source-ky input for {label}")
            return
        try:
            kx_sel = float(str(self.energy_balance_transfer_kx_var.get()).strip())
        except Exception:
            print(f"Invalid {diag_name} fixed source-kx input for {label}")
            return

        valid_t_request = np.asarray(t_indices, dtype=int).reshape(-1)
        if not self._load_fullt_if_needed(
            data,
            label,
            asym=use_asym,
            source_ky_value=ky_sel,
            source_kx_value=kx_sel,
            time_indices=valid_t_request,
        ):
            return

        fullt = np.asarray(getattr(data, data_attr))
        if fullt.ndim != 6 or fullt.shape[0] != 2:
            print(f"Unsupported {diag_name} shape for {label}: {fullt.shape}")
            return

        _ri, n_source, n_radial, n_signed, n_channel, n_t = fullt.shape
        if n_source <= 0 or n_radial <= 0 or n_signed <= 0 or n_channel <= 0 or n_t <= 0:
            print(f"Unsupported {diag_name} shape for {label}: {fullt.shape}")
            return

        kx_axis, radial_idx = self._energy_fullt_target_kx_axis(data, n_radial, label)
        n_source_for_target_axis = int(getattr(data, data_attr + "_n_source_ky_full", n_source))
        ky_target_axis = self._fullt_target_ky_axis(data, n_source_for_target_axis, n_signed)
        # Keep this call compatible with older tool deployments where
        # _fullt_source_ky_axis only accepts (data, n_source).
        ky_source_axis = self._fullt_source_ky_axis(data, n_source)
        ky_source_preferred = np.asarray(getattr(data, data_attr + "_source_ky_axis", []), dtype=float).reshape(-1)
        if ky_source_preferred.size >= n_source and np.all(np.isfinite(ky_source_preferred[:n_source])):
            ky_source_axis = ky_source_preferred[:n_source]
        source_ky_idx_cached = getattr(data, data_attr + "_source_ky_full_index", None)
        if kx_axis is None or radial_idx is None or ky_target_axis is None or ky_source_axis is None:
            print(f"Missing {diag_name} axes for {label}")
            return

        kx_axis = np.asarray(kx_axis, dtype=float).reshape(-1)
        radial_idx = np.asarray(radial_idx, dtype=int).reshape(-1)
        ky_target_axis = np.asarray(ky_target_axis, dtype=float).reshape(-1)
        ky_source_axis = np.asarray(ky_source_axis, dtype=float).reshape(-1)

        n_x = min(kx_axis.size, radial_idx.size, n_radial)
        n_y = min(ky_target_axis.size, n_signed)
        if n_x <= 0 or n_y <= 0:
            print(f"Empty {diag_name} axes for {label}")
            return
        kx_axis = kx_axis[:n_x]
        radial_idx = radial_idx[:n_x]
        ky_target_axis = ky_target_axis[:n_y]
        ky_idx = np.arange(n_y, dtype=int)

        finite_x = np.isfinite(kx_axis)
        finite_y = np.isfinite(ky_target_axis)
        if not np.any(finite_x) or not np.any(finite_y):
            print(f"No finite {diag_name} axes for {label}")
            return
        if not np.all(finite_x):
            kx_axis = kx_axis[finite_x]
            radial_idx = radial_idx[finite_x]
            n_x = int(kx_axis.size)
        if not np.all(finite_y):
            ky_target_axis = ky_target_axis[finite_y]
            ky_idx = ky_idx[finite_y]
            n_y = int(ky_target_axis.size)
        if n_x <= 0 or n_y <= 0:
            print(f"No finite {diag_name} axis points for {label}")
            return

        if ky_source_preferred.size == 1 and n_source == 1:
            source_ky_idx = 0
        else:
            source_ky_idx = self._nearest_ky_index(ky_source_axis[:n_source], ky_sel)
        if source_ky_idx is None:
            print(f"Could not map selected fixed source ky for {label}")
            return
        fixed_ky = float(ky_source_axis[source_ky_idx])
        if ky_source_preferred.size > 0 and np.isfinite(ky_source_preferred[0]):
            fixed_ky = float(ky_source_preferred[0])
        try:
            source_ky_full_idx_msg = int(source_ky_idx_cached)
        except Exception:
            source_ky_full_idx_msg = int(source_ky_idx)
        try:
            fixed_kx = float(getattr(data, data_attr + "_source_kx"))
        except Exception:
            fixed_kx = float(kx_sel)

        cached_time_idx = np.asarray(getattr(data, data_attr + "_time_indices", []), dtype=int).reshape(-1)
        if cached_time_idx.size == n_t:
            valid_t = np.arange(n_t, dtype=int)
        else:
            valid_t = np.asarray(t_indices, dtype=int).reshape(-1)
            valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
            if valid_t.size <= 0:
                valid_t = np.arange(n_t, dtype=int)
        if valid_t.size <= 0:
            print(f"No time samples for {label}")
            return

        if use_asym:
            pair_diag = self._fullt_asym_pair_symmetry_error(
                fullt, data, source_ky_idx, valid_t
            )
            if pair_diag is not None:
                print(
                    f"FULLT_ASYM pair symmetry for {label}: "
                    f"rms_rel={pair_diag['rms_rel']:.3e}, "
                    f"rms_abs={pair_diag['rms_abs']:.3e}, "
                    f"max_abs={pair_diag['max_abs']:.3e}, "
                    f"pairs={pair_diag['n_pair']}"
                )
                if (
                    np.isfinite(pair_diag["rms_rel"])
                    and pair_diag["rms_rel"] > 1.0e-8
                    and pair_diag["rms_abs"] > 1.0e-14
                ):
                    print(
                        f"FULLT_ASYM native data for {label} failed pair-symmetry check; "
                        "using FULLT-constructed ASYM for this plot."
                    )
                    fallback = self._fallback_construct_fullt_asym_for_plot(data, label)
                    if fallback is None:
                        print(f"Cannot construct FULLT_ASYM fallback for {label}.")
                        return
                    fullt = np.asarray(fallback)
                    _ri, n_source, n_radial, n_signed, n_channel, n_t = fullt.shape
                    if source_ky_idx >= n_source:
                        print(f"FULLT_ASYM fallback source index out of range for {label}.")
                        return
                    valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
                    if valid_t.size <= 0:
                        valid_t = np.arange(n_t, dtype=int)
                    if valid_t.size <= 0:
                        print(f"No fallback time samples for {label}.")
                        return

        z_map = np.asarray(fullt[0, source_ky_idx, radial_idx[:, None], ky_idx[None, :], 0, :], dtype=float)
        quantity_label = "Re"
        use_symmetric_scale = True

        # z_map is indexed explicitly as [target_kx, target_ky, time].
        # The plotted array must be [target_ky, target_kx].  Do not infer this
        # from shape: square grids make (n_x,n_y) indistinguishable from
        # (n_y,n_x) and silently swap kx/ky on the figure.
        z_avg = np.mean(z_map[:, :, valid_t], axis=2)
        if z_avg.shape != (n_x, n_y):
            print(
                f"Unsupported {diag_name} plot shape for {label}: z_avg={z_avg.shape}, "
                f"expected ({n_x},{n_y}) from [target_kx,target_ky]."
            )
            return
        z_plot = z_avg.T
        if np.all(~np.isfinite(z_avg)):
            print(f"No valid {diag_name} data for {label}")
            return

        kx_order = np.argsort(kx_axis)
        ky_order = np.argsort(ky_target_axis)
        kx_axis = kx_axis[kx_order]
        radial_idx = radial_idx[kx_order]
        ky_target_axis = ky_target_axis[ky_order]
        z_plot = np.asarray(z_plot, dtype=float)[np.ix_(ky_order, kx_order)]

        norm_pos = np.nan
        norm_neg = np.nan
        if normalize_by_max_t:
            finite_norm = z_plot[np.isfinite(z_plot)]
            if finite_norm.size > 0:
                positive = finite_norm[finite_norm > 0.0]
                negative = finite_norm[finite_norm < 0.0]
                norm_pos = float(np.nanmax(positive)) if positive.size > 0 else np.nan
                norm_neg = float(abs(np.nanmin(negative))) if negative.size > 0 else np.nan
                z_norm = np.zeros_like(z_plot, dtype=float)
                valid_pos = np.isfinite(norm_pos) and norm_pos > 0.0
                valid_neg = np.isfinite(norm_neg) and norm_neg > 0.0
                if valid_pos:
                    pos_mask = np.isfinite(z_plot) & (z_plot > 0.0)
                    z_norm[pos_mask] = z_plot[pos_mask] / norm_pos
                if valid_neg:
                    neg_mask = np.isfinite(z_plot) & (z_plot < 0.0)
                    z_norm[neg_mask] = z_plot[neg_mask] / norm_neg
                z_norm[~np.isfinite(z_plot)] = np.nan
                z_plot = z_norm
                print(
                    f"{diag_name} normalization for {label}: "
                    f"positive max T={norm_pos:.6e}, negative max |T|={norm_neg:.6e}."
                )
            else:
                print(f"{diag_name} normalization for {label}: no finite T values.")

        try:
            abs_z = np.abs(np.asarray(z_plot, dtype=float))
            flat = abs_z.reshape(-1)
            finite_flat = np.isfinite(flat)
            if np.any(finite_flat):
                finite_indices = np.where(finite_flat)[0]
                top_flat = finite_indices[np.argsort(flat[finite_flat])[-5:]][::-1]
                top_parts = []
                for flat_idx in top_flat:
                    iy_top, ix_top = np.unravel_index(int(flat_idx), z_plot.shape)
                    top_parts.append(
                        f"(kx={kx_axis[ix_top]:.4g}, ky={ky_target_axis[iy_top]:.4g}, "
                        f"{quantity_label}={z_plot[iy_top, ix_top]:.3e})"
                    )
                try:
                    case_dir_dbg = self._resolve_case_dir(data)
                except Exception:
                    case_dir_dbg = getattr(data, 'dir', None) or getattr(data, 'path', None) or ""
                source_dbg = getattr(data, data_attr + "_source", "unknown")
                first_t = int(valid_t[0]) if valid_t.size > 0 else -1
                last_t = int(valid_t[-1]) if valid_t.size > 0 else -1
                print(
                    f"{diag_name} {view_mode} debug for {label}: selected source ky input={ky_sel:.6g}, "
                    f"selected source kx input={kx_sel:.6g}, "
                    f"source_idx={source_ky_idx}, source_ky={fixed_ky:.6g}, source_kx={fixed_kx:.6g}, "
                    f"shape={fullt.shape}, data_id={id(data)}, array_id={id(fullt)}, "
                    f"source={source_dbg}, time_idx={first_t}:{last_t}/{n_t}, "
                    f"normalized_by_max_T={normalize_by_max_t}, "
                    f"case_dir={case_dir_dbg}, "
                    f"target_ky_range=[{np.nanmin(ky_target_axis):.6g},{np.nanmax(ky_target_axis):.6g}], "
                    f"top |T|: " + "; ".join(top_parts)
                )
        except Exception as e:
            print(f"{diag_name} {view_mode} debug failed for {label}: {e}")

        target_ky_mark_idx = self._nearest_finite_index(ky_target_axis, fixed_ky)
        target_kx_mark_idx = self._nearest_finite_index(kx_axis, fixed_kx)
        if view_mode == "vs kx":
            if target_ky_mark_idx is None:
                print(f"Could not map {diag_name} vs kx slice ky for {label}")
                return
            y_line = np.asarray(z_plot[target_ky_mark_idx, :], dtype=float).reshape(-1)
            x_line = np.asarray(kx_axis, dtype=float).reshape(-1)
            n_line = min(x_line.size, y_line.size)
            if n_line <= 0:
                print(f"No {diag_name} vs kx points for {label}")
                return
            x_line = x_line[:n_line]
            y_line = y_line[:n_line]
            finite_line = np.isfinite(x_line) & np.isfinite(y_line)
            if not np.any(finite_line):
                print(f"No finite {diag_name} vs kx points for {label}")
                return
            x_line = x_line[finite_line]
            y_line = y_line[finite_line]
            line_order = np.argsort(x_line)
            x_line = x_line[line_order]
            y_line = y_line[line_order]
            slice_ky = float(ky_target_axis[target_ky_mark_idx])

            self.ax.plot(
                x_line,
                y_line,
                '-o',
                linewidth=2.0,
                markersize=5.0,
                label=f"{label} source ky={fixed_ky:.6g}, slice ky'={slice_ky:.6g}",
            )
            if not getattr(self.ax, "_cgyro_fullt_vs_kx_reference_drawn", False):
                self.ax.axhline(0.0, color='0.45', linestyle=':', linewidth=0.9)
                self.ax.axvline(0.0, color='0.45', linestyle=':', linewidth=0.9)
                self.ax._cgyro_fullt_vs_kx_reference_drawn = True
            self.ax.set_xlabel(r"$k_x'\rho_s$")
            if normalize_by_max_t:
                self.ax.set_ylabel(rf"$\langle \Re\,{transfer_symbol}\rangle_t/\max |T|_{{\pm}}$")
            else:
                self.ax.set_ylabel(rf"$\langle \Re\,{transfer_symbol}\rangle_t$")
            norm_title = " normalized by max T" if normalize_by_max_t else ""
            self.ax.set_title(
                f"{diag_name} vs kx{norm_title}: "
                rf"$k_x\rho_s={fixed_kx:.3g},\ k_y\rho_s={fixed_ky:.3g},\ k_y'\rho_s={slice_ky:.3g}$"
            )
            self._apply_fullt_dynamic_x_ticks(self.ax)
            self.ax.grid(True, linestyle=':', linewidth=0.7, alpha=0.55)
            return

        kx_edges = self._axis_cell_edges_from_centers(kx_axis)
        ky_edges = self._axis_cell_edges_from_centers(ky_target_axis)
        if kx_edges is None or ky_edges is None:
            print(f"Invalid {diag_name} plot grid for {label}")
            return
        # Mark the diagonal target-ky == source-ky location as a visual anchor.
        diamond_x, diamond_y, z_diamond = self._diamond_lattice_from_centers(
            kx_axis,
            ky_target_axis,
            kx_edges,
            ky_edges,
            z_plot,
        )
        if z_diamond is None:
            print(f"Invalid {diag_name} diamond plot grid for {label}")
            return

        self.fig.clear()
        self._reset_figure_layout_defaults()
        try:
            toolbar = getattr(self, "toolbar", None)
            if toolbar is not None and hasattr(toolbar, "_nav_stack"):
                toolbar._nav_stack.clear()
                toolbar.set_history_buttons()
        except Exception:
            pass
        self.ax = self.fig.add_subplot(111)
        try:
            self.ax._cgyro_custom_grid = True
        except Exception:
            pass
        try:
            self.ax.set_facecolor('white')
        except Exception:
            pass

        finite = z_plot[np.isfinite(z_plot)]
        if finite.size > 0:
            vmax = float(np.nanmax(np.abs(finite)))
            zmax = float(np.nanmax(finite))
            zmin = float(np.nanmin(finite))
        else:
            vmax, zmin, zmax = 1.0, 0.0, 1.0
        if (not np.isfinite(vmax)) or vmax <= 0.0:
            vmax = 1.0

        if normalize_by_max_t:
            vmax = 1.0
            zmin = -1.0 if use_symmetric_scale else 0.0
            zmax = 1.0
        elif finite.size > 0:
            clip = float(np.nanpercentile(np.abs(finite), 99.7))
            if np.isfinite(clip) and clip > 0.0:
                vmax = clip
        if use_symmetric_scale and finite.size > 0 and np.isfinite(vmax) and vmax > 0.0:
            weak_cut = 0.18 * vmax
            z_diamond = np.ma.masked_array(
                np.where(np.abs(np.asarray(z_diamond, dtype=float)) < weak_cut, 0.0, z_diamond),
                mask=np.ma.getmaskarray(z_diamond),
            )

        use_symlog_color = False
        try:
            use_symlog_color = bool(self.log_y_var.get())
        except Exception:
            use_symlog_color = False

        if use_symmetric_scale:
            norm = None
            color_limits = dict(vmin=-vmax, vmax=vmax)
            if use_symlog_color:
                linthresh = max(vmax * 2.0e-2, np.finfo(float).tiny)
                norm = SymLogNorm(
                    linthresh=linthresh,
                    linscale=1.0,
                    vmin=-vmax,
                    vmax=vmax,
                    base=10,
                )
                color_limits = dict(norm=norm)
            pcm = self.ax.pcolormesh(
                diamond_x,
                diamond_y,
                z_diamond,
                cmap=FULLT_DIVERGING_CMAP,
                shading='gouraud',
                edgecolors='none',
                linewidth=0.0,
                antialiased=True,
                rasterized=True,
                **color_limits,
            )
            colorbar_extend = 'both'
        else:
            if (not np.isfinite(zmin)) or (not np.isfinite(zmax)) or zmax <= zmin:
                zmin, zmax = 0.0, vmax
            pcm = self.ax.pcolormesh(
                diamond_x,
                diamond_y,
                z_diamond,
                cmap='viridis',
                shading='gouraud',
                edgecolors='none',
                linewidth=0.0,
                antialiased=True,
                rasterized=True,
                vmin=zmin,
                vmax=zmax,
            )
            colorbar_extend = 'max'

        cbar = self.fig.colorbar(pcm, ax=self.ax, extend=colorbar_extend)
        if normalize_by_max_t:
            cbar.set_label(rf"$\langle {quantity_label}\,{transfer_symbol}\rangle_t/\max |T|_{{\pm}}$")
        else:
            cbar.set_label(rf"$\langle {quantity_label}\,{transfer_symbol}\rangle_t$")
        if use_symmetric_scale and np.isfinite(vmax) and vmax > 0.0:
            cbar.set_ticks(np.linspace(-vmax, vmax, 9))

        kx_mark_idx = target_kx_mark_idx
        if kx_mark_idx is not None and target_ky_mark_idx is not None:
            kx_mark = float(kx_axis[kx_mark_idx])
            ky_mark = float(ky_target_axis[target_ky_mark_idx])
            if np.isfinite(kx_mark) and np.isfinite(ky_mark):
                self.ax.scatter(
                    [kx_mark],
                    [ky_mark],
                    marker='D',
                    s=150,
                    facecolor='k',
                    edgecolor='white',
                    linewidth=1.1,
                    zorder=6,
                )
            else:
                print(f"Could not map finite {diag_name} marker point for {label}")
        else:
            print(f"Could not map {diag_name} marker point for {label}")

        self.ax.set_xlim(float(kx_edges[0]), float(kx_edges[-1]))
        self.ax.set_ylim(float(ky_edges[0]), float(ky_edges[-1]))
        self._apply_fullt_dynamic_x_ticks(self.ax)
        self._apply_fullt_dynamic_y_ticks(self.ax)
        self.ax.tick_params(axis='x', labelrotation=0, labelsize=8)
        self.ax.tick_params(axis='y', labelsize=8)
        self.ax.tick_params(axis='x', which='minor', length=3.0, width=0.75)
        self.ax.grid(False)
        self.ax.vlines(
            kx_axis,
            float(ky_edges[0]),
            float(ky_edges[-1]),
            color='0.26',
            linewidth=0.48,
            alpha=0.42,
            linestyles=':',
            zorder=4,
        )
        self.ax.hlines(
            ky_target_axis,
            float(kx_edges[0]),
            float(kx_edges[-1]),
            color='0.26',
            linewidth=0.48,
            alpha=0.32,
            linestyles=':',
            zorder=4,
        )
        self.ax.axhline(0.0, color='0.16', linewidth=0.9, alpha=0.46, zorder=5)
        self.ax.axvline(0.0, color='0.16', linewidth=0.9, alpha=0.46, zorder=5)
        self.ax.set_xlabel(r"$k_x'\rho_s$")
        self.ax.set_ylabel(r"$k_y'\rho_s$")
        if use_asym:
            title_math = (
                rf"$T^{{\Phi}}_{{\mathrm{{asym}},\,"
                rf"k_x\rho_s={fixed_kx:.3g},\ k_y\rho_s={fixed_ky:.3g}}}(k')$"
            )
        else:
            title_math = (
                rf"$T^{{\Phi}}_{{k_x\rho_s={fixed_kx:.3g},\ "
                rf"k_y\rho_s={fixed_ky:.3g}}}(k')$"
            )
        norm_title = " normalized by max T" if normalize_by_max_t else ""
        self.ax.set_title(f"{diag_name}{norm_title}: {title_math}")
        case_label = str(label)
        if len(case_label) > 56:
            case_label = "..." + case_label[-53:]
        self.ax.text(
            0.015,
            0.985,
            f"case: {case_label}\nsource ky index: {source_ky_full_idx_msg}",
            transform=self.ax.transAxes,
            ha='left',
            va='top',
            fontsize=8,
            color='0.25',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.65, pad=2.0),
        )
        coord_kx = np.asarray(kx_axis, dtype=float).copy()
        coord_ky = np.asarray(ky_target_axis, dtype=float).copy()
        coord_z = np.asarray(z_plot, dtype=float).copy()
        if hasattr(self, "_record_current_plot_xyz_dataset"):
            dataset_norm = " normalized_by_max_T" if normalize_by_max_t else ""
            self._record_current_plot_xyz_dataset(
                f"{diag_name} {quantity_label}{dataset_norm} {label}",
                coord_kx,
                coord_ky,
                coord_z,
            )
        self.ax.format_coord = lambda x, y: self._format_fullt_coord(
            x,
            y,
            coord_kx,
            coord_ky,
            coord_z,
        )
        try:
            self.ax.callbacks.connect(
                "xlim_changed",
                lambda ax: self._apply_fullt_dynamic_x_ticks(ax),
            )
        except Exception:
            pass
        try:
            self.ax.callbacks.connect(
                "ylim_changed",
                lambda ax: self._apply_fullt_dynamic_y_ticks(ax),
            )
        except Exception:
            pass
        try:
            t_start_txt = float(t_start)
            t_end_txt = float(t_end)
            if not (np.isfinite(t_start_txt) and np.isfinite(t_end_txt)):
                raise ValueError("non-finite time window")
            self.ax.text(
                0.985,
                0.985,
                rf"$\langle\cdot\rangle_t:\ {t_start_txt:.1f}-{t_end_txt:.1f}$",
                transform=self.ax.transAxes,
                ha='right',
                va='top',
                fontsize=8,
                color='0.25',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.65, pad=2.0),
            )
        except Exception:
            pass


