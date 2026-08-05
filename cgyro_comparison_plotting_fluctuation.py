"""
Option-specific plotting mixin for CGYRO comparison GUI.
Auto-extracted from cgyro_comparison_plotting.py during refactor.
"""

import numpy as np
import matplotlib.animation as animation
from matplotlib.colors import LogNorm

class FluctuationPlotting:
    def _get_fluctuation_case_color(self, case_label):
        """Return a stable line color for all fluctuation traces from one case."""
        color_map = getattr(self, "_fluctuation_case_color_map", None)
        if not isinstance(color_map, dict):
            color_map = {}
            self._fluctuation_case_color_map = color_map

        if case_label not in color_map:
            try:
                palette = self._get_line_color_palette()
            except Exception:
                palette = []
            if not palette:
                palette = [None]
            color_map[case_label] = palette[len(color_map) % len(palette)]
        return color_map[case_label]

    def _fluctuation_max_normalization_enabled(self):
        """Return whether 1D fluctuation profiles should be scaled by max(abs(y))."""
        var = getattr(self, "fluc_norm_max_var", None)
        try:
            return bool(var.get()) if var is not None else False
        except Exception:
            return False

    def _normalize_fluctuation_profile(self, values):
        """Normalize a plotted 1D profile by its finite maximum when requested."""
        profile = np.asarray(values, dtype=float)
        if not self._fluctuation_max_normalization_enabled():
            return profile, False

        finite = np.isfinite(profile)
        if not np.any(finite):
            return profile, False
        max_value = float(np.max(np.abs(profile[finite])))
        if not np.isfinite(max_value) or max_value <= 1.0e-15:
            return profile, False
        return profile / max_value, True

    def _fluc_case_axis_text_for_plot(self, case_label, axis_name):
        """Resolve a shared or per-case kx/ky selector for plotting."""
        getter = getattr(self, "_get_fluc_case_axis_text", None)
        if callable(getter):
            return getter(case_label, axis_name)

        var_name = "fluc_theta_kx_var" if axis_name == "kx" else "fluc_theta_ky_var"
        var = getattr(self, var_name, None)
        try:
            return str(var.get()).strip() if var is not None else ""
        except Exception:
            return ""

    def _load_fluctuation_moment_field(self, data, label, moment, main_ion_policy="all"):
        """
        Load fluctuation moment into unified complex shape `[nr, theta, ky, t]`.

        Returns:
        - field_4d: complex ndarray `[nr, theta, ky, t]` or None
        - spec_label: species label suffix for title/legend
        """
        spec_label = ""
        moment_key = str(moment).strip()

        if moment_key in ["Phi", "Apar", "Bpar"]:
            field = self._load_named_kxky_complex(
                data, label, moment_key, species_dependent=False
            )
            return field, spec_label

        if moment_key in ["Density", "Energy"]:
            species_field = self._load_named_kxky_complex(
                data, label, moment_key, species_dependent=True
            )
            if species_field is None:
                print(f"No {moment_key} data available for {label}")
                return None, spec_label
            n_species = int(species_field.shape[2]) if species_field.ndim >= 5 else 1
            target_indices, spec_label = self._resolve_species_indices(
                data,
                n_species,
                label,
                fallback_first=True,
                main_ion_policy=main_ion_policy,
                single_species_only=False,
            )
            if not target_indices:
                return None, spec_label
            field = np.sum(species_field[:, :, target_indices, :, :], axis=2)
            return field, spec_label

        if moment_key == "Temperature":
            mom_n = self._load_named_kxky_complex(
                data, label, "Density", species_dependent=True
            )
            mom_e = self._load_named_kxky_complex(
                data, label, "Energy", species_dependent=True
            )
            if mom_n is None or mom_e is None:
                print(f"Data for Temperature (n and e) not available in {label}")
                return None, spec_label

            n_species = min(int(mom_n.shape[2]), int(mom_e.shape[2]))
            target_indices, spec_label = self._resolve_species_indices(
                data,
                n_species,
                label,
                fallback_first=True,
                main_ion_policy=main_ion_policy,
                single_species_only=True,
            )
            if not target_indices:
                return None, spec_label
            s_idx = int(target_indices[0])

            temp = np.asarray(getattr(data, 'temp', []), dtype=float).reshape(-1)
            dens = np.asarray(getattr(data, 'dens', []), dtype=float).reshape(-1)
            t0 = float(temp[s_idx]) if s_idx < temp.size else 1.0
            n0 = float(dens[s_idx]) if s_idx < dens.size else 1.0
            if abs(n0) < 1e-12:
                n0 = 1.0
            field = ((2.0 / 3.0) * mom_e[:, :, s_idx, :, :] - t0 * mom_n[:, :, s_idx, :, :]) / n0
            return field, spec_label

        print(f"Unsupported fluctuation moment: {moment_key}")
        return None, spec_label

    def _resolve_optional_axis_indices(self, axis_values, raw_text, axis_name, label):
        """Return all indices for blank input or the nearest physical value."""
        axis = np.asarray(axis_values, dtype=float).reshape(-1)
        if axis.size <= 0:
            return np.array([], dtype=int), "none"

        text = str(raw_text).strip()
        if not text:
            return np.arange(axis.size, dtype=int), "average"

        index = self._resolve_axis_selection_index(
            axis,
            text,
            axis_name,
            label,
            # kx/ky entry fields represent physical quantities.  Use idx:n
            # when an explicit array index is needed.
            prefer_value=True,
        )
        index = max(0, min(int(index), axis.size - 1))
        return np.array([index], dtype=int), f"{axis[index]:.4g}"

    def _plot_fluctuation_theta(self, data, label, plot_type, t_indices, t_start, t_end):
        """Plot a time-averaged Phi/Apar/Bpar profile versus parallel theta."""
        field_name = plot_type.split()[0]
        field_complex = self._load_named_kxky_complex(
            data, label, field_name, species_dependent=False
        )
        if field_complex is None:
            print(f"No {field_name} data available for {label}")
            return

        field = np.asarray(field_complex, dtype=complex)
        if field.ndim != 4:
            print(f"Unsupported {field_name} shape for theta plot in {label}: {field.shape}")
            return

        n_r, n_theta, n_ky, n_t = field.shape
        if n_r <= 0 or n_theta <= 0 or n_ky <= 0 or n_t <= 0:
            print(f"Invalid {field_name} dimensions for theta plot in {label}: {field.shape}")
            return

        kx_axis, radial_idx = self._build_kx_axis(data, n_r, label)
        if kx_axis is None or radial_idx is None:
            print(f"No usable kx axis for theta plot in {label}")
            return
        kx_axis = np.asarray(kx_axis, dtype=float).reshape(-1)
        radial_idx = np.asarray(radial_idx, dtype=int).reshape(-1)
        n_kx = min(kx_axis.size, radial_idx.size)
        if n_kx <= 0:
            print(f"No usable kx points for theta plot in {label}")
            return
        kx_axis = kx_axis[:n_kx]
        radial_idx = radial_idx[:n_kx]
        valid_radial = (radial_idx >= 0) & (radial_idx < n_r)
        kx_axis = kx_axis[valid_radial]
        radial_idx = radial_idx[valid_radial]

        ky_axis = self._positive_ky_axis(getattr(data, "ky", []))
        if ky_axis.size != n_ky:
            ky_axis = np.arange(n_ky, dtype=float)
        else:
            ky_axis = ky_axis[:n_ky]

        kx_text = self._fluc_case_axis_text_for_plot(label, "kx")
        ky_text = self._fluc_case_axis_text_for_plot(label, "ky")
        kx_indices, kx_selection = self._resolve_optional_axis_indices(
            kx_axis, kx_text, "kx", label
        )
        ky_indices, ky_selection = self._resolve_optional_axis_indices(
            ky_axis, ky_text, "ky", label
        )
        if kx_indices.size <= 0 or ky_indices.size <= 0:
            print(f"No usable kx/ky selection for theta plot in {label}")
            return

        valid_t = np.asarray(t_indices, dtype=int).reshape(-1)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size <= 0:
            valid_t = np.arange(n_t, dtype=int)
            time_suffix = " (full time)"
        else:
            time_suffix = self._format_avg_suffix(
                t_start, t_end, prefix=f"Avg-{self._average_mode_short_tag()}"
            )

        selected = field[
            np.ix_(radial_idx[kx_indices], np.arange(n_theta), ky_indices, valid_t)
        ]
        rho_norm = float(getattr(data, "rho", 1.0))
        if not np.isfinite(rho_norm) or abs(rho_norm) < 1.0e-12:
            rho_norm = 1.0
        selected = selected / rho_norm

        amplitude = np.abs(selected)
        average_axes = (0, 2, 3)
        if self._use_mean_absolute_average():
            profile = np.mean(amplitude, axis=average_axes)
            y_label = rf"$\langle |{field_name}/\rho_s| \rangle_{{k_x,k_y,t}}$"
        else:
            profile = np.sqrt(np.mean(amplitude ** 2, axis=average_axes))
            y_label = rf"$\sqrt{{\langle |{field_name}/\rho_s|^2 \rangle_{{k_x,k_y,t}}}}$"

        profile, max_normalized = self._normalize_fluctuation_profile(profile)
        if max_normalized:
            y_label = f"{y_label} (max=1)"

        theta_axis = self._build_theta_over_pi_axis(data, n_theta)
        theta_axis = np.asarray(theta_axis, dtype=float).reshape(-1)
        if theta_axis.size != n_theta:
            theta_axis = np.linspace(-1.0, 1.0, n_theta)

        normalization_suffix = ", max-normalized" if max_normalized else ""
        plot_label = (
            f"{label} (kx={kx_selection}, ky={ky_selection}"
            f"{time_suffix}, {self._average_mode_name()}{normalization_suffix})"
        )
        self._plot_1d(theta_axis, profile, plot_label, plot_type)
        self.ax.set_xlabel(r"$\theta/\pi$")
        self.ax.set_ylabel(y_label)

    def _plot_fluctuation_1d(self, data, label, plot_type, t_indices, t_start, t_end):
        """Plot 1D fluctuation cuts versus ky, kx, time, or theta."""
        if "vs theta" in plot_type:
            self._plot_fluctuation_theta(data, label, plot_type, t_indices, t_start, t_end)
            return

        # Determine field name from plot_type e.g. "Phi vs ky", "Apar vs Time"
        field_name = plot_type.split()[0]
        field_complex = self._load_named_kxky_complex(
            data, label, field_name, species_dependent=False
        )
        if field_complex is None:
            print(f"No {field_name} data available for {label}")
            return
        # Match previous behavior: drop leftmost radial slot and use midplane theta.
        field_data = self._extract_midplane_kykxt(
            field_complex,
            data,
            label,
            drop_radial0=True,
            species_idx=0,
        )
        if field_data is None:
            return

        # Apply GyroBohm normalization if needed (data_plot.py does f/rhonorm)
        # cgyrodata usually stores raw. data_plot.py:getnorm('elec') sets rhonorm = rho.
        
        rho_norm = 1.0
        if hasattr(data, 'rho'):
            rho_norm = data.rho
        if abs(rho_norm) < 1e-12:
            rho_norm = 1.0
            
        # Normalize
        field_data = field_data / rho_norm
        avg_mode = self._average_mode_name()
        avg_tag = self._average_mode_short_tag()
        if self._use_mean_absolute_average():
            y_label_rho_norm = fr'$\sum \langle | {field_name}/\rho_s | \rangle_t$'
        else:
            y_label_rho_norm = fr'$\sqrt{{\sum \langle | {field_name}/\rho_s |^2 \rangle_t}}$'
        
        # shape: [nr-1, ny, nt]
        # Calculate amplitude squared |field|^2
        field_sq = np.abs(field_data)**2
        ky = self._positive_ky_axis(getattr(data, 'ky', []))
        if ky.size != field_data.shape[1]:
            ky = np.arange(field_data.shape[1], dtype=float)
        else:
            ky = ky[:field_data.shape[1]]
        t_valid = np.asarray(t_indices, dtype=int)
        t_valid = t_valid[(t_valid >= 0) & (t_valid < field_sq.shape[-1])]

        if "vs ky" in plot_type:
            # Plot time-averaged amplitude vs ky.
            # Blank kx input averages over all kx modes.
            kx_axis, radial_idx = self._build_kx_axis(data, field_data.shape[0], label)
            if kx_axis is None or radial_idx is None:
                print(f"No usable kx axis for {label}")
                return
            kx_axis = np.asarray(kx_axis, dtype=float).reshape(-1)
            radial_idx = np.asarray(radial_idx, dtype=int).reshape(-1)
            n_kx = min(kx_axis.size, radial_idx.size)
            kx_axis = kx_axis[:n_kx]
            radial_idx = radial_idx[:n_kx]
            valid_kx = (
                (radial_idx >= 0)
                & (radial_idx < field_data.shape[0])
                & np.isfinite(kx_axis)
            )
            kx_axis = kx_axis[valid_kx]
            radial_idx = radial_idx[valid_kx]
            kx_text = self._fluc_case_axis_text_for_plot(label, "kx")
            kx_indices, kx_selection = self._resolve_optional_axis_indices(
                kx_axis, kx_text, "kx", label
            )
            if kx_indices.size <= 0:
                print(f"No usable kx selection for {label}")
                return
            selected = field_data[radial_idx[kx_indices], :, :]
            if self._use_mean_absolute_average():
                field_ky_t = np.mean(np.abs(selected), axis=0)
            else:
                field_ky_t = np.mean(np.abs(selected) ** 2, axis=0)
            
            if self._use_mean_absolute_average():
                y_label = fr'$\langle | {field_name}/\rho_s | \rangle_{{k_x,t}}$'
            else:
                y_label = fr'$\sqrt{{\langle | {field_name}/\rho_s |^2 \rangle_{{k_x,t}}}}$'

            # Time average over selected window
            plot_label = f"{label} (kx={kx_selection}, ky=all, {avg_mode})"
            if t_valid.size > 0:
                if self._use_mean_absolute_average():
                    y_vals = np.mean(field_ky_t[:, t_valid], axis=1)
                else:
                    y_vals = np.sqrt(np.mean(field_ky_t[:, t_valid], axis=1))
                plot_label = self._append_avg_suffix(
                    plot_label, t_start, t_end, prefix=f"Avg-{avg_tag}"
                )
            else:
                # Fallback to last time point
                if self._use_mean_absolute_average():
                    y_vals = np.abs(field_ky_t[:, -1])
                else:
                    y_vals = np.sqrt(np.maximum(field_ky_t[:, -1], 0.0))

            y, max_normalized = self._normalize_fluctuation_profile(y_vals)
            if max_normalized:
                y_label = f"{y_label} (max=1)"
                plot_label += " [max-normalized]"
            x = ky
            n = min(x.size, y.size)
            self._plot_1d(x[:n], y[:n], plot_label, plot_type)
            self.ax.set_xlabel(r'$k_y \rho_s$')
            self.ax.set_ylabel(y_label)
            # self.ax.set_yscale('log')
            # self.ax.set_xscale('log') # Usually log-log for spectra

        elif "vs kx" in plot_type:
            # The x-axis is kx; blank fixed-ky input averages over all ky modes.
            ky_text = self._fluc_case_axis_text_for_plot(label, "ky")
            ky_indices, ky_selection = self._resolve_optional_axis_indices(
                ky, ky_text, "ky", label
            )
            if ky_indices.size <= 0:
                print(f"No usable ky selection for {label}")
                return
            kx_axis, radial_idx = self._build_kx_axis(data, field_data.shape[0], label)
            if kx_axis is None or radial_idx is None:
                print(f"No usable kx axis for {label}")
                return
            kx_axis = np.asarray(kx_axis, dtype=float).reshape(-1)
            radial_idx = np.asarray(radial_idx, dtype=int).reshape(-1)
            n_kx = min(kx_axis.size, radial_idx.size)
            kx_axis = kx_axis[:n_kx]
            radial_idx = radial_idx[:n_kx]
            valid_kx = (
                (radial_idx >= 0)
                & (radial_idx < field_data.shape[0])
                & np.isfinite(kx_axis)
            )
            kx_axis = kx_axis[valid_kx]
            radial_idx = radial_idx[valid_kx]
            selected = field_data[radial_idx, :, :][:, ky_indices, :]
            if self._use_mean_absolute_average():
                field_kx_t = np.mean(np.abs(selected), axis=1)
            else:
                field_kx_t = np.mean(np.abs(selected) ** 2, axis=1)

            if self._use_mean_absolute_average():
                y_label = fr'$\langle | {field_name}/\rho_s | \rangle_{{k_y,t}}$'
            else:
                y_label = fr'$\sqrt{{\langle | {field_name}/\rho_s |^2 \rangle_{{k_y,t}}}}$'

            plot_label = f"{label} (kx=all, ky={ky_selection}, {avg_mode})"
            if t_valid.size > 0:
                if self._use_mean_absolute_average():
                    y_vals = np.mean(field_kx_t[:, t_valid], axis=1)
                else:
                    y_vals = np.sqrt(np.mean(field_kx_t[:, t_valid], axis=1))
                plot_label = self._append_avg_suffix(
                    plot_label, t_start, t_end, prefix=f"Avg-{avg_tag}"
                )
            else:
                if self._use_mean_absolute_average():
                    y_vals = np.abs(field_kx_t[:, -1])
                else:
                    y_vals = np.sqrt(np.maximum(field_kx_t[:, -1], 0.0))

            y, max_normalized = self._normalize_fluctuation_profile(y_vals)
            if max_normalized:
                y_label = f"{y_label} (max=1)"
                plot_label += " [max-normalized]"

            x = kx_axis
            n = min(x.size, y.size)
            self._plot_1d(x[:n], y[:n], plot_label, plot_type)
            self.ax.set_xlabel(r'$k_x \rho_s$')
            self.ax.set_ylabel(y_label)

        
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
            if self._use_mean_absolute_average():
                field_abs_n0 = np.sum(np.abs(field_data[:, ky_idx_0, :]), axis=0)  # [nt]
            else:
                field_sq_n0 = np.sum(field_sq[:, ky_idx_0, :], axis=0)  # [nt]
            
            # n>0 intensity: sum over kx (axis 0) and sum over ky!=0
            # Create mask for ky!=0
            mask_n = np.ones(field_sq.shape[1], dtype=bool)
            mask_n[ky_idx_0] = False
            
            if self._use_mean_absolute_average():
                field_abs_nn = np.sum(np.abs(field_data[:, mask_n, :]), axis=(0, 1))  # [nt]
                y0 = np.abs(field_abs_n0)
                yn = np.abs(field_abs_nn)
            else:
                field_sq_nn = np.sum(field_sq[:, mask_n, :], axis=(0, 1))  # [nt]
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
            case_color = self._get_fluctuation_case_color(label)
            
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
                    f"({avg_inner}, {avg_mode}, "
                    f"Mean: {mean_y0:.2e}, Std: {std_y0:.2e})"
                )

                # n>0 stats
                yn_subset = yn[valid_t]
                mean_yn = np.mean(yn_subset)
                std_yn = np.std(yn_subset)
                label_nn = (
                    f"{label} (n>0) "
                    f"({avg_inner}, {avg_mode}, "
                    f"Mean: {mean_yn:.2e}, Std: {std_yn:.2e})"
                )
                
                # Plot n=0
                self.ax.plot(x, y0, label=label_n0, linestyle='--', color=case_color)
                self.ax.plot([t_mean_start, t_mean_end], [mean_y0, mean_y0],
                             linestyle=':', color=case_color, linewidth=1.5)
                
                # Plot n>0
                self.ax.plot(x, yn, label=label_nn, linestyle='-', color=case_color)
                self.ax.plot([t_mean_start, t_mean_end], [mean_yn, mean_yn],
                             linestyle='-.', color=case_color, linewidth=1.5)
            else:
                self.ax.plot(x, y0, label=f"{label} (n=0)", linestyle='--', color=case_color)
                self.ax.plot(x, yn, label=f"{label} (n>0)", linestyle='-', color=case_color)
            
            self.ax.set_xlabel(r'$t (a/c_s)$')
            self.ax.set_ylabel(y_label_rho_norm)

    def _plot_fluctuation_kxky_map_from_2d(self, data, label, plot_type, t_indices, t_start, t_end):
        """Render a spectral kx-ky map from the Fluctuation 2D Moment selector."""
        # This path intentionally hangs off the Fluctuation-2D controls: the
        # user chooses a physical moment there, then the view selector decides
        # whether it is rendered in real space (xy/xt) or spectral space (kxky).
        moment = str(self.moment_var.get()).strip() if hasattr(self, "moment_var") else "Phi"
        field_4d, spec_label = self._load_fluctuation_moment_field(
            data,
            label,
            moment,
            main_ion_policy="all",
        )
        if field_4d is None:
            print(f"No {moment} data available for {label}")
            return

        field_data = self._extract_midplane_kykxt(
            field_4d,
            data,
            label,
            drop_radial0=True,
            species_idx=0,
        )
        if field_data is None:
            return

        normalized_by_rho = moment in ["Phi", "Apar", "Bpar"]
        if normalized_by_rho:
            # Field amplitudes are displayed in the same rho_s-normalized units
            # as the existing 1D fluctuation plots, so comparing Phi vs ky and
            # Phi vs kxky does not silently change units.
            rho_norm = float(getattr(data, 'rho', 1.0))
            if abs(rho_norm) < 1e-12:
                rho_norm = 1.0
            field_data = field_data / rho_norm

        map_label = f"{label}{spec_label}" if spec_label else label
        self._plot_fluctuation_kxky_map(
            data,
            map_label,
            moment,
            field_data,
            t_indices,
            t_start,
            t_end,
            normalized_by_rho=normalized_by_rho,
        )

    def _plot_fluctuation_kxky_map(
        self,
        data,
        label,
        field_name,
        field_data,
        t_valid,
        t_start,
        t_end,
        normalized_by_rho=True,
    ):
        """Plot a time-averaged midplane fluctuation amplitude map on the kx-ky plane."""
        arr = np.asarray(field_data, dtype=complex)
        if arr.ndim != 3:
            print(f"Unsupported {field_name} kxky data shape for {label}: {arr.shape}")
            return

        n_kx, n_ky, n_t = arr.shape
        if n_kx <= 0 or n_ky <= 0 or n_t <= 0:
            print(f"Empty {field_name} kxky data for {label}")
            return

        valid_t = np.asarray(t_valid, dtype=int).reshape(-1)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        avg_mode = self._average_mode_name()
        avg_tag = self._average_mode_short_tag()

        if valid_t.size > 0:
            arr_t = arr[:, :, valid_t]
            avg_suffix = self._format_avg_suffix(t_start, t_end, prefix=f"Avg-{avg_tag}")
        else:
            # If no valid window is selected, show the last available sample
            # rather than failing; this matches the behavior of other snapshot
            # style plots in the comparison tool.
            arr_t = arr[:, :, -1:]
            avg_suffix = " (last time)"

        if self._use_mean_absolute_average():
            z = np.mean(np.abs(arr_t), axis=2)
            value_label = rf'{field_name}/\rho_s' if normalized_by_rho else str(field_name)
            cbar_label = rf'$\langle |{value_label}| \rangle_t$'
        else:
            z = np.sqrt(np.mean(np.abs(arr_t) ** 2, axis=2))
            value_label = rf'{field_name}/\rho_s' if normalized_by_rho else str(field_name)
            cbar_label = rf'$\sqrt{{\langle |{value_label}|^2 \rangle_t}}$'

        kx_axis, radial_idx = self._build_kx_axis(data, n_kx, label)
        if kx_axis is None:
            print(f"No usable kx axis for {label}")
            return
        kx_axis = np.asarray(kx_axis, dtype=float).reshape(-1)
        radial_idx = np.asarray(radial_idx, dtype=int).reshape(-1)
        n_x = min(kx_axis.size, radial_idx.size, z.shape[0])
        if n_x <= 0:
            print(f"No usable kx points for {label}")
            return
        kx_axis = kx_axis[:n_x]
        # `radial_idx` maps the plotted kx axis back into the CGYRO radial
        # storage dimension.  This is needed because some outputs omit the
        # leftmost special/Nyquist slot while FULLT-like outputs keep it.
        z = z[radial_idx[:n_x], :]

        ky_axis = self._positive_ky_axis(getattr(data, 'ky', []))
        ky_axis = np.asarray(ky_axis, dtype=float).reshape(-1)
        if ky_axis.size <= 0:
            ky_axis = np.arange(n_ky, dtype=float)
        n_y = min(ky_axis.size, z.shape[1])
        if n_y <= 0:
            print(f"No usable ky points for {label}")
            return
        ky_axis = ky_axis[:n_y]
        z = np.asarray(z[:, :n_y], dtype=float)

        finite_kx = np.isfinite(kx_axis)
        finite_ky = np.isfinite(ky_axis)
        if not np.any(finite_kx) or not np.any(finite_ky):
            print(f"No finite kx/ky axis values for {label}")
            return
        kx_axis = kx_axis[finite_kx]
        z = z[finite_kx, :]
        ky_axis = ky_axis[finite_ky]
        z = z[:, finite_ky]

        kx_order = np.argsort(kx_axis)
        ky_order = np.argsort(ky_axis)
        # Sort before imshow so the image extent is monotonic.  Without this,
        # signed radial storage order can produce a visually blank or flipped
        # kxky map even when the underlying array is valid.
        kx_axis = kx_axis[kx_order]
        ky_axis = ky_axis[ky_order]
        z = z[np.ix_(kx_order, ky_order)]

        z_plot = np.ma.masked_invalid(z.T)
        use_log_z = self._use_fluc2d_log_z()
        image_norm = None
        title_scale_suffix = ""
        if use_log_z:
            positive_values = z[np.isfinite(z) & (z > 0.0)]
            if positive_values.size > 0:
                z_min = float(np.min(positive_values))
                z_max = float(np.max(positive_values))
                if z_max > z_min:
                    z_plot = np.ma.masked_less_equal(z_plot, 0.0)
                    image_norm = LogNorm(vmin=z_min, vmax=z_max)
                    cbar_label += " (log scale)"
                    title_scale_suffix = " [Log z]"
                else:
                    print(f"Log z unavailable for {label}: kxky map has only one positive value.")
            else:
                print(f"Log z unavailable for {label}: kxky map has no positive values.")

        def _extent_bounds(axis):
            axis = np.asarray(axis, dtype=float).reshape(-1)
            if axis.size == 1:
                v = float(axis[0])
                pad = 0.5 if abs(v) < 1.0e-12 else abs(v) * 0.05
                return v - pad, v + pad
            return float(np.nanmin(axis)), float(np.nanmax(axis))

        x0, x1 = _extent_bounds(kx_axis)
        y0, y1 = _extent_bounds(ky_axis)
        pcm = self.ax.imshow(
            z_plot,
            extent=[x0, x1, y0, y1],
            origin='lower',
            aspect='auto',
            interpolation='bicubic',
            cmap=self._get_2d_contour_cmap(),
            norm=image_norm,
        )
        cbar = self.fig.colorbar(pcm, ax=self.ax)
        cbar.set_label(cbar_label)
        self.ax.set_xlabel(r'$k_x \rho_s$')
        self.ax.set_ylabel(r'$k_y \rho_s$')
        self.ax.set_title(
            f'{field_name} vs kxky: {label}{avg_suffix} ({avg_mode}){title_scale_suffix}'
        )

        if hasattr(self, "_record_current_plot_xyz_dataset"):
            # Store the physical axes and z values before Matplotlib turns them
            # into pixels.  Data -> Save current plot data can then export an
            # Origin-ready XYZ table matching the displayed map.
            self._record_current_plot_xyz_dataset(
                f"{field_name} vs kxky: {label}",
                kx_axis,
                ky_axis,
                z_plot,
            )

    def _plot_fluctuation_2d(self, data, label, plot_type, t_indices, t_start, t_end):
        """Render ky-kx fluctuation contour animation or snapshot from bigfield data."""
        use_x_rhoe = self._use_fluc2d_x_rhoe()
        x_unit_factor = self._get_rhos_to_rhoe_factor(data, label) if use_x_rhoe else 1.0
        x_unit_label = r'$x/\rho_e$' if use_x_rhoe else r'$x/\rho_s$'
        y_unit_label = r'$y/\rho_e$' if use_x_rhoe else r'$y/\rho_s$'

        moment = self.moment_var.get()
        # Keep compatibility with previous behavior in 2D:
        # species moments select one ion by default.
        field_4d, _spec_label = self._load_fluctuation_moment_field(
            data, label, moment, main_ion_policy="first"
        )
        if field_4d is None:
            print(f"Data for {moment} not available in {label}")
            return

        c_midplane = self._extract_midplane_kykxt(
            field_4d,
            data,
            label,
            drop_radial0=False,
            species_idx=0,
        )
        if c_midplane is None:
            return

        # Define calculation function
        def get_frame_data(ti):
            """Reconstruct one real-space frame for contour plotting at time index `ti`."""
            if ti < 0 or ti >= c_midplane.shape[-1]:
                return
            c = c_midplane[:, :, ti]
            
            # FFT params
            nr = int(c_midplane.shape[0])
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
                self.ax.contourf(xp, yp, np.transpose(fp), levels, cmap=self._get_2d_contour_cmap())
                if hasattr(self, "_clear_current_plot_data"):
                    self._clear_current_plot_data()
                if hasattr(self, "_record_current_plot_xyz_dataset"):
                    self._record_current_plot_xyz_dataset(
                        f"{moment} fluctuation xy: {label}",
                        xp,
                        yp,
                        np.transpose(fp),
                    )
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
            cs = self.ax.contourf(xp, yp, np.transpose(fp), levels, cmap=self._get_2d_contour_cmap())
            self.fig.colorbar(cs, ax=self.ax)
            if hasattr(self, "_record_current_plot_xyz_dataset"):
                self._record_current_plot_xyz_dataset(
                    f"{moment} fluctuation xy: {label}",
                    xp,
                    yp,
                    np.transpose(fp),
                )
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
        field, spec_label = self._load_fluctuation_moment_field(
            data, label, moment, main_ion_policy="all"
        )

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

        cs = self.ax.contourf(t_plot, x_plot, np.transpose(f_tx), levels=levels, cmap=self._get_2d_contour_cmap())
        self.fig.colorbar(cs, ax=self.ax)
        if hasattr(self, "_record_current_plot_xyz_dataset"):
            self._record_current_plot_xyz_dataset(
                f"{moment} fluctuation x-t: {label}",
                t_plot,
                x_plot,
                np.transpose(f_tx),
            )
        self.ax.set_xlabel(r'$t\ (a/c_s)$')
        self.ax.set_ylabel(x_ylabel)

        moment_tag = moment
        if spec_label:
            moment_tag = f"{moment} ({spec_label})"
        self.ax.set_title(
            f'{moment_tag} x-t Contour: {label} '
            f'(ky={ky0:.4g}, t={time_tag})'
        )

