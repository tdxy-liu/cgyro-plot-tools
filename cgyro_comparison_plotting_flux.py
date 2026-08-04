"""
Option-specific plotting mixin for CGYRO comparison GUI.
Auto-extracted from cgyro_comparison_plotting.py during refactor.
"""

import numpy as np
import traceback
from tkinter import messagebox

class FluxPlotting:
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
            # Estimated-kx spectra are rebuilt from other case metadata, so
            # they do not require `ky_flux` to be present in memory.
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
                # Real-ion normalization changes both the plotted flux scale
                # and the natural axes: ky*rho_s and t*c_s/a must be converted
                # with the same reference-ion factors used for Q/Gamma.
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
            cmap = self._get_2d_contour_cmap()

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
                    if hasattr(self, "_record_current_plot_xyz_dataset"):
                        # The plotted contour has physical axes (ky, time).
                        # Cache them before Matplotlib converts the contour to
                        # polygon artists, which are hard to export back to XYZ.
                        self._record_current_plot_xyz_dataset(
                            f"{field_tag} flux ky-time: {label}",
                            ky_plot,
                            t_plot,
                            z_t_ky,
                        )

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
            if hasattr(self, "_record_current_plot_xyz_dataset"):
                # Keep Origin export consistent with the decomposed branch.
                self._record_current_plot_xyz_dataset(
                    f"Phi flux ky-time: {label}",
                    ky_plot,
                    t_plot,
                    z_t_ky,
                )
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

    def _plot_flux_vs_2d_selected_cases(self, selected_cases, plot_type):
        """Plot Flux time-averaged trends across cases vs selected input parameter."""
        try:
            x_param = str(self.flux_scan_xparam_var.get()).strip()
        except Exception:
            x_param = ""

        try:
            show_errorbar = bool(self.flux_2d_errorbar_var.get())
        except Exception:
            show_errorbar = False

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
            y_window = np.asarray(y_t[t_idx], dtype=float)
            y_val = float(np.mean(y_window)) if y_window.size > 0 else float(y_t[-1])
            y_err = float(np.std(y_window)) if y_window.size > 1 else 0.0

            flux_norm_scale = self._get_flux_real_ion_norm_scale(
                data, moment_idx, label=f"{case_name}{' - ' + spec_label if spec_label else ''}"
            )
            y_val = float(y_val) * float(flux_norm_scale)
            y_err = float(y_err) * abs(float(flux_norm_scale))

            # Group by all other varying parameters.
            grp_parts = []
            for k in other_params:
                try:
                    v = float(scalars.get(k))
                    grp_parts.append((str(k), v))
                except Exception:
                    grp_parts.append((str(k), np.nan))
            grp_key = tuple(grp_parts)
            grouped.setdefault(grp_key, []).append((float(x_val), float(y_val), float(y_err), case_name, t0_used, t1_used))

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
            y_err_vals = np.asarray([p[2] for p in pts], dtype=float)
            finite = np.isfinite(x_vals) & np.isfinite(y_vals)
            if show_errorbar:
                finite = finite & np.isfinite(y_err_vals)
            if not np.any(finite):
                continue
            x_vals = x_vals[finite]
            y_vals = y_vals[finite]
            y_err_vals = y_err_vals[finite]

            order = np.argsort(x_vals)
            x_vals = x_vals[order]
            y_vals = y_vals[order]
            y_err_vals = y_err_vals[order]

            # Merge repeated x by averaging.
            x_unique = []
            y_unique = []
            yerr_unique = []
            i0 = 0
            while i0 < len(x_vals):
                x0 = x_vals[i0]
                i1 = i0 + 1
                while i1 < len(x_vals) and abs(x_vals[i1] - x0) <= 1.0e-12:
                    i1 += 1
                x_unique.append(x0)
                y_slice = y_vals[i0:i1]
                y_mean = float(np.mean(y_slice))
                y_unique.append(y_mean)
                if show_errorbar:
                    e_slice = np.maximum(y_err_vals[i0:i1], 0.0)
                    scatter = y_slice - y_mean
                    yerr_unique.append(float(np.sqrt(np.mean(e_slice**2 + scatter**2))))
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

            x_plot = np.asarray(x_unique, dtype=float)
            y_plot = np.asarray(y_unique, dtype=float)
            if show_errorbar:
                yerr_plot = np.asarray(yerr_unique, dtype=float)
                err_container = self.ax.errorbar(
                    x_plot,
                    y_plot,
                    yerr=yerr_plot,
                    linestyle='-',
                    marker='o',
                    markersize=5.0,
                    linewidth=1.8,
                    capsize=3.0,
                    capthick=1.8,
                    elinewidth=1.8,
                    barsabove=True,
                    zorder=10,
                    label=curve_label,
                )
                errorbar_alpha = 0.60
                for capline in err_container[1]:
                    capline.set_alpha(errorbar_alpha)
                for barlinecol in err_container[2]:
                    barlinecol.set_alpha(errorbar_alpha)
            else:
                self.ax.plot(
                    x_plot,
                    y_plot,
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
            denom = r"Q_{GB,\mathrm{ri}}" if use_real_ion_norm else r"Q_{GB,\mathrm{D}}"
            self.ax.set_ylabel(rf"$Q_{{{sub}}}/{denom}$")
        else:
            if use_real_ion_norm:
                self.ax.set_ylabel(rf"$\Gamma_{{{sub}}}/\Gamma_{{GB,\mathrm{{ri}}}}$")
            else:
                self.ax.set_ylabel(rf"$\Gamma_{{{sub}}}/\Gamma_{{GB,\mathrm{{D}}}}$")

        if skipped > 0:
            print(f"Flux vs 2D: skipped {skipped} case(s) due to missing parameter/flux/time data.")

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

