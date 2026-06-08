"""
Zonal-flow plotting mixin for CGYRO comparison GUI.

This module isolates ZF-related plotting and linear-gamma comparison logic from
the main plotting mixin to keep file size manageable and responsibilities clear.
"""

import os
import re

import numpy as np

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


class ZfPlotting:
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
        # Linear scans are often stored next to a nonlinear case rather than
        # inside it, so checking the parent directory saves repeated browsing.
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
            with open(file_path, "r", encoding="utf-8") as f:
                for line in f:
                    s = line.strip()
                    if not s or s.startswith("#"):
                        continue
                    parts = re.split(r"[\s,]+", s)
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

        # zf_kx(t) contribution uses |phi(kx,t)| and shearing scalar is
        # sum_kx kx^2 * |phi(kx,t)| / rho.
        rho = self._rho_scalar_for_norm(data, label)
        shear_rate = np.sum((kx[:, None] ** 2) * np.abs(phi_kx_t), axis=0) / rho

        n_t = min(t.size, shear_rate.size)
        if n_t <= 0:
            print(f"No time samples available for {label}")
            return None, None

        return t[:n_t], shear_rate[:n_t]

    def _plot_zf_exb_shearing_spectrum(self, data, label, t_indices, t_start, t_end):
        """Legacy wrapper for zonal ExB shearing spectrum plot."""
        return self._plot_zf_exb_shearing_rate_kx_spectrum(
            data, label, t_indices, t_start, t_end
        )

    def _compute_zf_phi_theta0_kx_profile(self, data, label, t_indices, t_start, t_end):
        """
        Compute time-averaged |phi_ZF|/rho_s profile vs kx at theta=0.

        Returns dict with:
        - x: sorted kx axis
        - y_phi: sorted <|phi_ZF|/rho_s>_t
        - label: label with avg suffix when a valid time window exists
        - kx_raw, phi_kx_t, valid_t, n_t: raw arrays for downstream reuse
        """
        kx, phi_kx_t, t = self._get_zf_exb_phi_kx_t(data, label)
        if kx is None:
            return None

        # ZF plots use the raw zonal phi extracted by `_get_zf_exb_phi_kx_t`.
        # Normalize only at the plotting metric level, matching the existing
        # collect/cgyro_plot shearing-rate convention.
        rho = self._rho_scalar_for_norm(data, label)

        phi_abs_kx_t = np.abs(phi_kx_t) / rho
        n_t = min(t.size, phi_abs_kx_t.shape[-1])
        if n_t <= 0:
            print(f"No time samples available for {label}")
            return None

        valid_t = np.asarray(t_indices, dtype=int)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size > 0:
            y_phi = np.mean(phi_abs_kx_t[:, valid_t], axis=1)
            label_out = self._append_avg_suffix(label, t_start, t_end, prefix="Avg")
        else:
            y_phi = phi_abs_kx_t[:, -1]
            label_out = label

        x = np.asarray(kx, dtype=float).reshape(-1)
        y_phi = np.asarray(y_phi, dtype=float).reshape(-1)
        if x.size != y_phi.size:
            n = min(x.size, y_phi.size)
            x = x[:n]
            y_phi = y_phi[:n]

        order = np.argsort(x)
        # Sort signed kx for line plotting; storage order can be centered or
        # shifted depending on whether the original binary kept the special
        # radial slot.
        x = x[order]
        y_phi = y_phi[order]

        return {
            "x": x,
            "y_phi": y_phi,
            "label": label_out,
            "kx_raw": np.asarray(kx, dtype=float).reshape(-1),
            "phi_kx_t": np.asarray(phi_kx_t),
            "valid_t": valid_t,
            "n_t": int(n_t),
        }

    def _plot_zf_exb_shearing_rate_kx_spectrum(self, data, label, t_indices, t_start, t_end):
        """Plot time-averaged zonal ExB shearing spectrum versus `kx`."""
        prof = self._compute_zf_phi_theta0_kx_profile(data, label, t_indices, t_start, t_end)
        if prof is None:
            return

        x = prof["x"]
        y_phi = prof["y_phi"]
        label_plot = prof["label"]
        y = (x ** 2) * y_phi

        self._plot_1d(x, y, label_plot, "ZF ExB Shearing Spectrum")
        self.ax.set_xlabel(r"$k_x \rho_s$")
        self.ax.set_ylabel(r"$\omega_{E\times B}^{ZF}\ (c_s/a)$")

    def _plot_zf_phi_kx_theta0_spectrum(self, data, label, t_indices, t_start, t_end):
        """Plot zonal potential spectrum |phi_ZF| vs kx at theta=0."""
        prof = self._compute_zf_phi_theta0_kx_profile(data, label, t_indices, t_start, t_end)
        if prof is None:
            return

        self._plot_1d(prof["x"], prof["y_phi"], prof["label"], "ZF Phi Spectrum (theta0)")
        self.ax.set_xlabel(r"$k_x \rho_s$")
        self.ax.set_ylabel(r"$\langle |\phi_{ZF}(k_x,\theta=0)|/\rho_s \rangle_t$")

    def _plot_zf_exb_fig4_kx_eq_ky(self, data, label, t_indices, t_start, t_end):
        """Legacy wrapper kept for backward compatibility."""
        return self._plot_zf_exb_vs_gamma_lin_kx_equals_ky(
            data, label, t_indices, t_start, t_end
        )

    def _plot_zf_exb_fig4_kx_equals_ky_comparison(self, data, label, t_indices, t_start, t_end):
        """Backward-compatible wrapper for renamed vs-gamma-linear plotter."""
        return self._plot_zf_exb_vs_gamma_lin_kx_equals_ky(
            data, label, t_indices, t_start, t_end
        )

    def _plot_zf_exb_vs_gamma_lin_kx_equals_ky(self, data, label, t_indices, t_start, t_end):
        """
        Compare zonal ExB metrics against linear growth rate on shared `k` axis
        under `kx = ky` convention.
        """
        ky_target = None
        try:
            ky_text = str(self.zf_gamma_lin_ky_var.get()).strip()
            if len(ky_text) > 0:
                ky_target = float(ky_text)
        except Exception:
            ky_target = None

        prof = self._compute_zf_phi_theta0_kx_profile(data, label, t_indices, t_start, t_end)
        if prof is None:
            return

        x_zf = np.asarray(prof["x"], dtype=float).reshape(-1)
        y_phi = np.asarray(prof["y_phi"], dtype=float).reshape(-1)
        y_zf = (x_zf ** 2) * y_phi

        kx_raw = np.asarray(prof["kx_raw"], dtype=float).reshape(-1)
        phi_kx_t = np.asarray(prof["phi_kx_t"])
        valid_t = np.asarray(prof["valid_t"], dtype=int).reshape(-1)
        n_t = int(prof["n_t"])

        rho = self._rho_scalar_for_norm(data, label)

        phi_abs_kx_t = np.abs(np.asarray(phi_kx_t)[:, :n_t]) / rho
        if valid_t.size > 0:
            zf_kx_norm = np.mean(phi_abs_kx_t[:, valid_t], axis=1)
        else:
            zf_kx_norm = phi_abs_kx_t[:, -1]
        n_k = min(kx_raw.size, zf_kx_norm.size)
        kx_for_vzf = np.asarray(kx_raw[:n_k], dtype=float).reshape(-1)
        zf_kx_for_vzf = np.asarray(zf_kx_norm[:n_k], dtype=float).reshape(-1)
        vzf_mean = float(
            0.5 * np.sqrt(np.sum(np.abs(kx_for_vzf * zf_kx_for_vzf) ** 2))
        )
        # Fig4-style comparison plots ky*<V_ZF> against gamma_lin by identifying
        # kx and ky numerically.  This is a diagnostic convention rather than a
        # new CGYRO output quantity.

        avg_suffix = self._format_avg_suffix(t_start, t_end, prefix="Avg")
        if valid_t.size > 0:
            zf_label = f"{label} " + r"$\langle\omega_{ZF}(k_x)\rangle$" + avg_suffix
            vzf_label = f"{label} " + r"$k_y\langle V_{ZF}\rangle$" + avg_suffix
        else:
            zf_label = f"{label} " + r"$\langle\omega_{ZF}(k_x)\rangle$"
            vzf_label = f"{label} " + r"$k_y\langle V_{ZF}\rangle$"

        y_kx_vzf = np.asarray(x_zf * vzf_mean, dtype=float).reshape(-1)
        n = min(x_zf.size, y_zf.size, y_kx_vzf.size)
        x_zf = x_zf[:n]
        y_zf = y_zf[:n]
        y_kx_vzf = y_kx_vzf[:n]

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
                    f"kyVzf={kxvzf_sel:.6g}, "
                    f"omegaZF/gamma_lin={ratio_wzf:.6g}, "
                    f"kyVzf/gamma_lin={ratio_kxvzf:.6g}"
                )
            else:
                print(
                    f"{label}: ky={ky_target:.6g}, "
                    "cannot compute omegaZF/gamma_lin and kyVzf/gamma_lin "
                    "because gamma_lin is missing/zero at selected ky."
                )

        zf_label_plot = zf_label + r"  $(k_x=k_y)$"
        vzf_label_plot = vzf_label + r"  $(k_x=k_y)$"
        if ratio_wzf is not None and np.isfinite(ratio_wzf):
            zf_label_plot += rf" $[\omega_{{ZF}}/\gamma_{{lin}}={ratio_wzf:.3g}]$"
        if ratio_kxvzf is not None and np.isfinite(ratio_kxvzf):
            vzf_label_plot += rf" $[k_yV_{{ZF}}/\gamma_{{lin}}={ratio_kxvzf:.3g}]$"

        self.ax.plot(
            x_lin,
            y_lin,
            marker="o",
            linestyle="-",
            color=self._get_gamma_lin_color(),
            label=f"{label} " + r"$\gamma_{lin}(k_y)$",
        )
        self.ax.plot(
            x_zf,
            y_zf,
            marker="x",
            linestyle="--",
            label=zf_label_plot,
        )
        self.ax.plot(
            x_zf,
            y_kx_vzf,
            marker="^",
            linestyle="-.",
            label=vzf_label_plot,
        )
        self.ax.set_xlim(right=x_max)
        self.ax.set_xlabel(r"$k\rho_s\ \ (k_x=k_y)$")
        self.ax.set_ylabel(r"$\langle\omega_{ZF}\rangle,\ k_y\langle V_{ZF}\rangle,\ \gamma_{lin}\ (c_s/a)$")

        if ky_target is not None:
            try:
                ky_mark = float(ky_target)
            except Exception:
                ky_mark = np.nan
            if np.isfinite(ky_mark):
                x_min_vis = float(min(np.min(x_lin), np.min(x_zf)))
                x_max_vis = float(max(np.max(x_lin), np.max(x_zf)))
                if x_min_vis <= ky_mark <= x_max_vis:
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
                            ky_mark, color="0.35", linestyle="--", linewidth=1.2, alpha=0.9
                        )
                        try:
                            vln._zf_ky_marker = True
                        except Exception:
                            pass
