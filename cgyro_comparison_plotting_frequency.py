"""
Option-specific plotting mixin for CGYRO comparison GUI.
Auto-extracted from cgyro_comparison_plotting.py during refactor.
"""

import numpy as np

class FrequencyPlotting:
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
            # Different pygacode versions have exposed `freq` as [component,ky,t]
            # or [component,t,ky].  Shape-based detection keeps the plot robust
            # without requiring a version check.
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
            # Avoid silently creating infinities at ky=0; those are exported as
            # NaN and skipped by Matplotlib/Origin rather than polluting scales.
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

