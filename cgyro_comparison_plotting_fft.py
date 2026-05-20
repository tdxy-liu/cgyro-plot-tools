"""
Option-specific plotting mixin for CGYRO comparison GUI.
Auto-extracted from cgyro_comparison_plotting.py during refactor.
"""

import numpy as np
import os

try:
    import scipy.signal as sp_signal
    from scipy.optimize import curve_fit as sp_curve_fit
except Exception:
    sp_signal = None
    sp_curve_fit = None

DEFAULT_POD_Z_WINDOW_PI = 8.0

class FftPlotting:
    def _plot_fluctuation_fft(self, data, label, plot_type, t_indices):
        """Render FFT spectra from fluctuation time traces in selected FFT mode."""
        # Determine field name
        field_name = plot_type.split()[0] # "Phi", "Apar", "Bpar"
        field_complex = self._load_named_kxky_complex(
            data, label, field_name, species_dependent=False
        )
        if field_complex is None:
            print(f"No {field_name} data available for {label}")
            return
        field_t_all = self._extract_midplane_kykxt(
            field_complex,
            data,
            label,
            drop_radial0=False,
            species_idx=0,
        )
        if field_t_all is None:
            return
        # field_t_all shape: [n_radial, n_ky, nt]
        
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
            
            cs1 = ax1.contourf(ky_grid, omega_grid, mag_plot1, levels=50, cmap=self._get_2d_contour_cmap())
            self.fig.colorbar(cs1, ax=ax1)
            ax1.set_title(f'{field_name} FFT {spectrum_label} (kx=0): {label}')
            ax1.set_ylabel(r'$\omega (c_s/a)$')
            
            cs2 = ax2.contourf(ky_grid, omega_grid, mag_plot2, levels=50, cmap=self._get_2d_contour_cmap())
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
                                ax1.plot(
                                    freq_ky, freq_omega_plot, '--',
                                    color=self._get_gamma_lin_color(), linewidth=1.5, label='Linear Freq'
                                )
                                ax2.plot(
                                    freq_ky, freq_omega_plot, '--',
                                    color=self._get_gamma_lin_color(), linewidth=1.5, label='Linear Freq'
                                )
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
            
            cs1 = ax1.contourf(kx_grid, omega_grid_kx, mag_plot1, levels=50, cmap=self._get_2d_contour_cmap())
            self.fig.colorbar(cs1, ax=ax1)
            ax1.set_title(f'{field_name} FFT {spectrum_label} (ky=0): {label}')
            ax1.set_ylabel(r'$\omega (c_s/a)$')
            
            cs2 = ax2.contourf(kx_grid, omega_grid_kx, mag_plot2, levels=50, cmap=self._get_2d_contour_cmap())
            self.fig.colorbar(cs2, ax=ax2)
            ax2.set_title(f'{field_name} FFT {spectrum_label} {title2_suffix}: {label}')
            ax2.set_ylabel(r'$\omega (c_s/a)$')
            ax2.set_xlabel(r'$k_x \rho_s$')
            
            self.ax = ax2

        self.fig.tight_layout()

