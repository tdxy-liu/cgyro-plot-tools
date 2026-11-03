"""
Option-specific plotting mixin for CGYRO comparison GUI.
Auto-extracted from cgyro_comparison_plotting.py during refactor.
"""

import numpy as np

class FrequencyPlotting:
    @staticmethod
    def _frequency_data(data, n_ky):
        """Load only the frequency diagnostic skipped by pygacode's fast=True."""
        cached = getattr(data, 'freq', None)
        if cached is not None:
            return np.asarray(cached)

        # getdata() also loads ballooning potentials and other large arrays.
        # extract() respects the case precision and binary/ASCII precedence.
        extract = getattr(data, 'extract', None)
        if not callable(extract):
            raise ValueError("Frequency data are unavailable: the case reader has no extract() method.")
        _elapsed, fmt, values = extract('.cgyro.freq')
        if fmt == 'null':
            raise ValueError(
                "Missing out.cgyro.freq / bin.cgyro.freq in " + str(getattr(data, 'dir', 'the case folder'))
                + ". This case has no saved frequency/growth-rate diagnostic."
            )
        values = np.asarray(values)
        n_time = np.asarray(getattr(data, 't', [])).size or int(getattr(data, 'n_time', 0))
        expected = 2 * n_ky * n_time
        if n_time <= 0 or values.size != expected:
            raise ValueError(
                f"Invalid {fmt}.cgyro.freq: found {values.size} values; expected {expected} "
                f"(2 components x {n_ky} ky x {n_time} times). Reload the case after its files finish writing."
            )
        freq = values.reshape((2, n_ky, n_time), order='F')
        data.freq = freq
        return freq

    def _plot_frequency_growth(self, data, label, plot_type, t_indices, t_start, t_end):
        """Plot frequency or growth-rate spectra, with optional time-window averaging."""
        x = self._positive_ky_axis(getattr(data, 'ky', []))
        n_ky = x.size
        if n_ky == 0 or not np.all(np.isfinite(x)):
            raise ValueError("Frequency plotting requires a nonempty, finite ky grid.")

        comp_idx = 0 if plot_type == "Frequency" else 1
        freq = self._frequency_data(data, n_ky)
        y = np.array([])
        used_avg = False

        if freq.shape == (2, n_ky):
            y = np.asarray(freq[comp_idx]).reshape(-1)
        elif freq.ndim == 3 and freq.shape[0] == 2:
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
                raise ValueError(f"Frequency shape {freq.shape} does not match the {n_ky}-point ky grid.")

            valid_t = np.asarray(t_indices, dtype=int)
            t_len = comp.shape[time_axis]
            if t_len == 0:
                raise ValueError("The frequency diagnostic has no time samples.")
            if np.any((valid_t < 0) | (valid_t >= t_len)):
                raise ValueError("The selected time window is not covered by the frequency diagnostic. Reload the case.")

            if valid_t.size > 0:
                y = np.mean(np.take(comp, valid_t, axis=time_axis), axis=time_axis)
                used_avg = True
            else:
                y = np.take(comp, -1, axis=time_axis)

        if y.size == 0:
            raise ValueError(f"Unsupported frequency shape {freq.shape}; expected (2, ky) or (2, ky, time).")

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

        if not np.any(np.isfinite(y)):
            raise ValueError("No finite frequency values in the selected time window (check the data and division by ky=0).")
        self._plot_1d(x, y, label, plot_type)

