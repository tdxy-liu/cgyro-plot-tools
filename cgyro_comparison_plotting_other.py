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

class OtherPlotting:
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

