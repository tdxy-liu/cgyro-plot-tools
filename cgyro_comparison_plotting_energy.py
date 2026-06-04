"""
Option-specific plotting mixin for CGYRO comparison GUI.
Auto-extracted from cgyro_comparison_plotting.py during refactor.
"""

import numpy as np
import os
from tkinter import messagebox
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm, TwoSlopeNorm

try:
    import scipy.signal as sp_signal
    from scipy.optimize import curve_fit as sp_curve_fit
except Exception:
    sp_signal = None
    sp_curve_fit = None

DEFAULT_POD_Z_WINDOW_PI = 8.0
FULLT_DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    'fullt_blue_white_red',
    [
        (0.00, '#2166ac'),
        (0.42, '#d5e8f3'),
        (0.50, '#ffffff'),
        (0.58, '#f7d7cd'),
        (1.00, '#b2182b'),
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

    def _load_triad_common(self, data, label, require_flux):
        """Load triad data, optionally with ky_flux."""
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

    @staticmethod
    def _nearest_finite_index(axis, value):
        """Return the index of the finite axis value closest to `value`."""
        axis = np.asarray(axis, dtype=float).reshape(-1)
        finite_mask = np.isfinite(axis)
        if not np.any(finite_mask):
            return None
        axis_finite = axis[finite_mask]
        idx_local = int(np.argmin(np.abs(axis_finite - float(value))))
        return int(np.where(finite_mask)[0][idx_local])

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

    def _compute_energy_balance_single_vs_kxky(self, data, label, n_sel, quantity_key, t_indices):
        """
        Build single-quantity 2D map versus (kx, ky) from triad channels.

        quantity_key:
          - 'T'   -> < Re(idx1) >_t  on (kx, ky) grid
          - 'N'   -> < Re(idx2) >_t  on (kx, ky) grid
          - 'T-N' -> T - N on (kx, ky) grid

        Returns (kx_grid, ky_grid, z_data) or None on failure.
        z_data shape: [n_kx, n_ky] (ready for contourf with transpose).
        """
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
        f, _spec_label = self._triad_species_view(triad, spec, label)
        if f is None:
            return None
        # f: [radial, 8, n_n, t]
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

        # Build kx axis
        kx = np.asarray(getattr(data, 'kx', []), dtype=float).reshape(-1)
        if kx.size == n_radial:
            pass
        elif kx.size == n_radial - 1:
            # CGYRO often drops the first radial slot; pad with a dummy zero
            kx = np.concatenate(([0.0], kx))
        else:
            length = float(getattr(data, 'length', 0.0))
            if length > 1e-12:
                dkx = 2.0 * np.pi / length
            else:
                dkx = 1.0
            p_index = np.arange(n_radial, dtype=float) - (n_radial // 2)
            kx = p_index * dkx

        # Time average over selected window
        if quantity_key == 'T':
            z = np.mean(np.real(f[:, 0, :n_ky, :][:, :, valid_t]), axis=2)  # [radial, n_ky]
        elif quantity_key == 'N':
            z = np.mean(np.real(f[:, 1, :n_ky, :][:, :, valid_t]), axis=2)  # [radial, n_ky]
        elif quantity_key == 'T-N':
            t_avg = np.mean(np.real(f[:, 0, :n_ky, :][:, :, valid_t]), axis=2)
            n_avg = np.mean(np.real(f[:, 1, :n_ky, :][:, :, valid_t]), axis=2)
            z = t_avg - n_avg
        else:
            # Fallback to T for unknown quantities
            z = np.mean(np.real(f[:, 0, :n_ky, :][:, :, valid_t]), axis=2)

        # Align kx length with radial dimension.
        n_kx_use = min(int(kx.size), int(z.shape[0]))
        if n_kx_use <= 0:
            print(f"No usable kx points for {label}")
            return None
        n_ky_use = min(int(ky.size), int(z.shape[1]))
        if n_ky_use <= 0:
            print(f"No usable ky points for {label}")
            return None
        kx = np.asarray(kx[:n_kx_use], dtype=float).reshape(-1)
        ky = np.asarray(ky[:n_ky_use], dtype=float).reshape(-1)
        z = np.asarray(z[:n_kx_use, :n_ky_use], dtype=float)

        finite_kx = np.isfinite(kx)
        finite_ky = np.isfinite(ky)
        if not np.any(finite_kx) or not np.any(finite_ky):
            print(f"No finite kx/ky points for {label}")
            return None
        kx = kx[finite_kx]
        z = z[finite_kx, :]
        ky = ky[finite_ky]
        z = z[:, finite_ky]

        kx_order = np.argsort(kx)
        ky_order = np.argsort(ky)
        kx = kx[kx_order]
        ky = ky[ky_order]
        z = z[np.ix_(kx_order, ky_order)]

        # Build meshgrid for contourf
        kx_grid, ky_grid = np.meshgrid(kx, ky, indexing='ij')
        return kx_grid, ky_grid, np.asarray(z, dtype=float)

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
        else:
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
                    f"Energy balance single plot ({qty_name}, vs ky){avg_suffix}"
                )
            return

        if mode_txt == 'vs kxky':
            n_sel = self._parse_energy_balance_n_index()
            if qty_txt == 'entropy':
                print(f"entropy vs kxky not supported for {label}; use T, N or T-N.")
                return
            kxky_data = self._compute_energy_balance_single_vs_kxky(
                data, label, n_sel, qty_name, t_indices
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
            cs = self.ax.imshow(
                np.ma.masked_invalid(z_plot.T),
                extent=[
                    float(np.nanmin(kx_axis)), float(np.nanmax(kx_axis)),
                    float(np.nanmin(ky_axis)), float(np.nanmax(ky_axis)),
                ],
                origin='lower',
                aspect='auto',
                interpolation='bicubic',
                cmap=FULLT_DIVERGING_CMAP,
                norm=norm,
            )
            cbar = self.fig.colorbar(cs, ax=self.ax, extend='both')
            cbar.set_label(y_label)
            self.ax.set_xlabel(r"$k_x \rho_s$")
            self.ax.set_ylabel(r"$k_y \rho_s$")
            self.ax.set_title(
                f"{qty_name} map{avg_suffix}"
            )
            return

        terms = self._compute_triad_v2_terms(data, label)
        if terms is None:
            return
        t = np.asarray(terms['t'], dtype=float).reshape(-1)
        if qty_txt == 't':
            y_time = np.asarray(terms['T0'], dtype=float)
        elif qty_txt == 'n':
            y_time = np.asarray(terms['N0'], dtype=float)
        elif qty_txt == 'entropy':
            y_time = np.asarray(terms.get('S0', terms['dSg_dt']), dtype=float)
        else:
            y_time = np.asarray(terms['T0'], dtype=float) - np.asarray(terms['N0'], dtype=float)

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

    def _load_fullt_if_needed(self, data, label):
        """Load and normalize bin.cgyro.fullt to [ri,source_ky,target_kx,target_ky,channel,time]."""
        n_n = int(getattr(data, 'n_n', 0))
        n_radial = int(getattr(data, 'n_radial', 0))
        n_time_hint = int(getattr(data, 'n_time', 0))

        def _accept_existing(arr):
            arr = np.asarray(arr)
            if arr.ndim != 6 or arr.shape[0] != 2:
                return None
            if n_n > 0 and n_radial > 0:
                n_signed = 2 * n_n - 1
                # New pygacode layout.
                if arr.shape[1] == n_n and arr.shape[2] == n_radial and arr.shape[3] == n_signed:
                    return arr
                # Raw Fortran reshape layout.
                if arr.shape[1] == n_radial and arr.shape[2] == n_signed and arr.shape[4] == n_n:
                    return np.transpose(arr, (0, 4, 1, 2, 3, 5))
            return None

        if hasattr(data, 'fullt'):
            arr = _accept_existing(getattr(data, 'fullt'))
            if arr is not None:
                data.fullt = arr
                data.fullt_n_channel = int(arr.shape[4])
                return True

        if hasattr(data, 'full_t'):
            arr = _accept_existing(getattr(data, 'full_t'))
            if arr is not None:
                data.fullt = arr
                data.fullt_n_channel = int(arr.shape[4])
                return True

        if n_n <= 0 or n_radial <= 0:
            print(f"Cannot load FULLT for {label}: missing n_n/n_radial metadata.")
            return False
        if not hasattr(data, 'extract'):
            print(f"Cannot load FULLT for {label}: data object has no extract() method.")
            return False

        try:
            _tmsg, fmt, raw = data.extract('.cgyro.fullt')
        except TypeError:
            _tmsg, fmt, raw = data.extract('.cgyro.fullt', cmplx=False)
        except Exception as e:
            print(f"FULLT read failed for {label}: {e}")
            return False

        if fmt == 'null':
            print(f"No FULLT data available for {label}. Need bin/out.cgyro.fullt.")
            return False

        raw = np.asarray(raw)
        if raw.size <= 0:
            print(f"Empty FULLT data for {label}.")
            return False
        if np.iscomplexobj(raw):
            raw = np.column_stack((np.real(raw).reshape(-1), np.imag(raw).reshape(-1))).reshape(-1)
        raw = np.asarray(raw, dtype=float).reshape(-1)

        block_no_time = 2 * n_radial * (2 * n_n - 1) * n_n
        if block_no_time <= 0:
            print(f"Invalid FULLT dimensions for {label}.")
            return False

        n_channel = 0
        n_time = 0
        if n_time_hint > 0 and raw.size % (block_no_time * n_time_hint) == 0:
            n_channel_try = int(raw.size // (block_no_time * n_time_hint))
            if n_channel_try in (1, 2):
                n_channel = n_channel_try
                n_time = n_time_hint

        if n_channel <= 0:
            preferred = int(getattr(data, 'fullt_n_channel', 0))
            candidates = [preferred, 2, 1] if preferred in (1, 2) else [2, 1]
            for n_channel_try in candidates:
                denom = block_no_time * n_channel_try
                if denom > 0 and raw.size % denom == 0:
                    n_time_try = int(raw.size // denom)
                    if n_time_try > 0:
                        n_channel = int(n_channel_try)
                        n_time = int(n_time_try)
                        break

        if n_channel not in (1, 2) or n_time <= 0:
            print(
                f"Cannot infer FULLT shape for {label}: raw_size={raw.size}, "
                f"n_radial={n_radial}, n_n={n_n}, n_time_hint={n_time_hint}."
            )
            return False

        nd = block_no_time * n_channel * n_time
        if raw.size < nd:
            print(f"FULLT data is smaller than inferred shape for {label}.")
            return False
        try:
            raw6 = np.reshape(
                raw[:nd],
                (2, n_radial, 2 * n_n - 1, n_channel, n_n, n_time),
                order='F',
            )
        except Exception as e:
            print(f"FULLT reshape failed for {label}: {e}")
            return False

        data.fullt = np.transpose(raw6, (0, 4, 1, 2, 3, 5))
        data.fullt_n_channel = int(n_channel)
        try:
            data.n_time = int(n_time)
        except Exception:
            pass
        return True

    def _fullt_target_kx_axis(self, data, n_radial, label):
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
        ky_pos = self._positive_ky_axis(getattr(data, 'kynorm', getattr(data, 'ky', [])))
        ky_pos = np.asarray(ky_pos, dtype=float).reshape(-1)
        if ky_pos.size >= n_n and n_n > 0:
            ky_pos = ky_pos[:n_n]
            ky_axis = np.concatenate((-ky_pos[1:][::-1], ky_pos))
        elif n_n > 1:
            ky_axis = np.arange(-(n_n - 1), n_n, dtype=float)
        else:
            ky_axis = np.asarray([0.0], dtype=float)

        if ky_axis.size != n_signed:
            ky_axis = np.arange(-(n_signed // 2), n_signed - (n_signed // 2), dtype=float)
        return ky_axis

    def _fullt_source_ky_axis(self, data, n_source):
        """Build non-negative source-ky axis for FULLT source_ky index selection."""
        ky = self._positive_ky_axis(getattr(data, 'kynorm', getattr(data, 'ky', [])))
        ky = np.asarray(ky, dtype=float).reshape(-1)
        if ky.size >= n_source:
            return ky[:n_source]
        return np.arange(n_source, dtype=float)

    def _plot_energy_balance_fullt(self, data, label, t_indices, t_start, t_end):
        """Plot FULLT for one fixed source ky over scanned target kx/target ky."""
        if not self._load_fullt_if_needed(data, label):
            return

        fullt = np.asarray(data.fullt)
        if fullt.ndim != 6 or fullt.shape[0] != 2:
            print(f"Unsupported FULLT shape for {label}: {fullt.shape}")
            return

        _ri, n_source, n_radial, n_signed, n_channel, n_t = fullt.shape
        if n_source <= 0 or n_radial <= 0 or n_signed <= 0 or n_channel <= 0 or n_t <= 0:
            print(f"Unsupported FULLT shape for {label}: {fullt.shape}")
            return

        qty = str(self.energy_balance_transfer_quantity_var.get()).strip().lower()
        try:
            ky_sel = float(str(self.energy_balance_transfer_ky_var.get()).strip())
        except Exception:
            print(f"Invalid FULLT fixed-ky input for {label}")
            return

        kx_axis, radial_idx = self._fullt_target_kx_axis(data, n_radial, label)
        ky_target_axis = self._fullt_target_ky_axis(data, n_source, n_signed)
        ky_source_axis = self._fullt_source_ky_axis(data, n_source)
        if kx_axis is None or radial_idx is None or ky_target_axis is None or ky_source_axis is None:
            print(f"Missing FULLT axes for {label}")
            return

        kx_axis = np.asarray(kx_axis, dtype=float).reshape(-1)
        radial_idx = np.asarray(radial_idx, dtype=int).reshape(-1)
        ky_target_axis = np.asarray(ky_target_axis, dtype=float).reshape(-1)
        ky_source_axis = np.asarray(ky_source_axis, dtype=float).reshape(-1)

        n_x = min(kx_axis.size, radial_idx.size, n_radial)
        n_y = min(ky_target_axis.size, n_signed)
        if n_x <= 0 or n_y <= 0:
            print(f"Empty FULLT axes for {label}")
            return
        kx_axis = kx_axis[:n_x]
        radial_idx = radial_idx[:n_x]
        ky_target_axis = ky_target_axis[:n_y]

        source_ky_idx = self._nearest_finite_index(ky_source_axis[:n_source], ky_sel)
        if source_ky_idx is None:
            print(f"Could not map selected fixed ky for {label}")
            return
        fixed_ky = float(ky_source_axis[source_ky_idx])
        target_ky_mark_idx = self._nearest_finite_index(ky_target_axis, fixed_ky)

        valid_t = np.asarray(t_indices, dtype=int).reshape(-1)
        valid_t = valid_t[(valid_t >= 0) & (valid_t < n_t)]
        if valid_t.size <= 0:
            valid_t = np.arange(n_t, dtype=int)
        if valid_t.size <= 0:
            print(f"No time samples for {label}")
            return

        real_map = np.asarray(fullt[0, source_ky_idx, radial_idx, :n_y, 0, :], dtype=float)
        if n_channel >= 2:
            imag_map = np.asarray(fullt[0, source_ky_idx, radial_idx, :n_y, 1, :], dtype=float)
        else:
            imag_map = np.asarray(fullt[1, source_ky_idx, radial_idx, :n_y, 0, :], dtype=float)

        if qty == 'im':
            z_map = imag_map
            quantity_label = "Im"
            title_quantity = "imaginary diagnostic"
            use_symmetric_scale = True
        elif qty == 'abs':
            z_map = np.sqrt(real_map * real_map + imag_map * imag_map)
            quantity_label = "Abs"
            title_quantity = "complex magnitude"
            use_symmetric_scale = False
        else:
            z_map = real_map
            quantity_label = "Re"
            title_quantity = "real transfer"
            use_symmetric_scale = True

        z_avg = np.mean(z_map[:, :, valid_t], axis=2)
        if z_avg.shape == (n_y, n_x):
            z_plot = z_avg
        elif z_avg.shape == (n_x, n_y):
            z_plot = z_avg.T
        else:
            print(
                f"Unsupported FULLT plot shape for {label}: z_avg={z_avg.shape}, "
                f"expected ({n_x},{n_y}) or ({n_y},{n_x})."
            )
            return
        if np.all(~np.isfinite(z_avg)):
            print(f"No valid FULLT data for {label}")
            return

        self.fig.clear()
        self._reset_figure_layout_defaults()
        self.ax = self.fig.add_subplot(111)
        try:
            self.ax.set_facecolor('white')
        except Exception:
            pass

        finite = z_avg[np.isfinite(z_avg)]
        if finite.size > 0:
            vmax = float(np.nanmax(np.abs(finite)))
            zmax = float(np.nanmax(finite))
            zmin = float(np.nanmin(finite))
        else:
            vmax, zmin, zmax = 1.0, 0.0, 1.0
        if (not np.isfinite(vmax)) or vmax <= 0.0:
            vmax = 1.0

        if finite.size > 0:
            clip = float(np.nanpercentile(np.abs(finite), 99.8))
            if np.isfinite(clip) and clip > 0.0:
                vmax = clip

        use_symlog_color = False
        try:
            use_symlog_color = bool(self.log_y_var.get())
        except Exception:
            use_symlog_color = False

        if use_symmetric_scale:
            norm = None
            color_limits = dict(vmin=-vmax, vmax=vmax)
            if use_symlog_color:
                linthresh = max(vmax * 1.0e-3, np.finfo(float).tiny)
                norm = SymLogNorm(
                    linthresh=linthresh,
                    linscale=1.0,
                    vmin=-vmax,
                    vmax=vmax,
                    base=10,
                )
                color_limits = dict(norm=norm)
            pcm = self.ax.imshow(
                np.ma.masked_invalid(z_plot),
                extent=[float(np.nanmin(kx_axis)), float(np.nanmax(kx_axis)),
                        float(np.nanmin(ky_target_axis)), float(np.nanmax(ky_target_axis))],
                origin='lower',
                aspect='auto',
                interpolation='bicubic',
                cmap=FULLT_DIVERGING_CMAP,
                **color_limits,
            )
            colorbar_extend = 'both'
        else:
            if (not np.isfinite(zmin)) or (not np.isfinite(zmax)) or zmax <= zmin:
                zmin, zmax = 0.0, vmax
            pcm = self.ax.imshow(
                np.ma.masked_invalid(z_plot),
                extent=[float(np.nanmin(kx_axis)), float(np.nanmax(kx_axis)),
                        float(np.nanmin(ky_target_axis)), float(np.nanmax(ky_target_axis))],
                origin='lower',
                aspect='auto',
                interpolation='bicubic',
                cmap='viridis',
                vmin=zmin,
                vmax=zmax,
            )
            colorbar_extend = 'max'

        cbar = self.fig.colorbar(pcm, ax=self.ax, extend=colorbar_extend)
        cbar.set_label(rf"$\langle {quantity_label}\,T^{{\Phi}}\rangle_t$")

        kx_mark_idx = self._nearest_finite_index(kx_axis, 0.0)
        if kx_mark_idx is not None and target_ky_mark_idx is not None:
            kx_mark = float(kx_axis[kx_mark_idx])
            ky_mark = float(ky_target_axis[target_ky_mark_idx])
            self.ax.scatter(
                [kx_mark],
                [ky_mark],
                marker='D',
                s=115,
                facecolor='k',
                edgecolor='white',
                linewidth=1.0,
                zorder=5,
            )
        else:
            print(f"Could not map FULLT marker point for {label}")

        self.ax.grid(False)
        self.ax.axhline(0.0, color='0.40', linewidth=0.65, alpha=0.22)
        self.ax.axvline(0.0, color='0.40', linewidth=0.65, alpha=0.22)
        self.ax.set_xlabel(r"$k_x'\rho_s$")
        self.ax.set_ylabel(r"$k_y'\rho_s$")
        self.ax.set_title(
            rf"$T^{{\Phi}}_{{k_x\rho_s=0,\ k_y\rho_s={fixed_ky:.3g}}}(k')$"
        )
        try:
            self.ax.text(
                0.985,
                0.985,
                rf"$\langle\cdot\rangle_t:\ {float(t_start):.1f}-{float(t_end):.1f}$",
                transform=self.ax.transAxes,
                ha='right',
                va='top',
                fontsize=8,
                color='0.25',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.65, pad=2.0),
            )
        except Exception:
            pass


