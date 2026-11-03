"""Regression coverage for frequency diagnostics on fast-loaded CGYRO cases."""

from types import SimpleNamespace
import unittest
from unittest.mock import Mock

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy as np

from cgyro_comparison_plotting import Plotting


class FrequencyPlotTests(unittest.TestCase):
    def plot(self):
        plot = Plotting()
        plot.fig = Figure()
        FigureCanvasAgg(plot.fig)
        plot.ax = plot.fig.subplots()
        plot.norm_ky_var = SimpleNamespace(get=lambda: False)
        return plot

    def case(self, **kwargs):
        values = np.arange(24., dtype=float).reshape(2, 3, 4)
        data = SimpleNamespace(
            dir="/case/", ky=np.array([0., -.1, -.2]), t=np.arange(4.),
            extract=Mock(return_value=("elapsed", "out", values.ravel(order="F"))),
            getdata=Mock(side_effect=AssertionError("must not read unrelated diagnostics")),
        )
        data.__dict__.update(kwargs)
        return data, values

    def test_fast_loaded_case_reads_only_frequency_once_and_averages_selected_times(self):
        plot = self.plot()
        data, values = self.case()
        for component, kind in enumerate(("Frequency", "Growth Rate")):
            plot._plot_frequency_growth(data, "case", kind, [2, 3], 2., 3.)
            line = plot.ax.lines[-1]
            np.testing.assert_allclose(line.get_xdata(), [0., .1, .2])
            np.testing.assert_allclose(line.get_ydata(), values[component, :, 2:].mean(axis=1))
            self.assertIn("Avg: 2.0-3.0", line.get_label())
        data.extract.assert_called_once_with('.cgyro.freq')
        data.getdata.assert_not_called()
        np.testing.assert_array_equal(data.freq, values)

    def test_cached_transposed_and_snapshot_layouts_remain_supported(self):
        for layout in ('time-first', 'snapshot'):
            with self.subTest(layout=layout):
                plot = self.plot()
                data, values = self.case()
                data.freq = values.transpose(0, 2, 1) if layout == 'time-first' else values[:, :, -1]
                plot._plot_frequency_growth(data, 'case', 'Growth Rate', [1, 2], 1., 2.)
                expected = values[1, :, 1:3].mean(axis=1) if layout == 'time-first' else values[1, :, -1]
                np.testing.assert_allclose(plot.ax.lines[0].get_ydata(), expected)
                data.extract.assert_not_called()

    def test_missing_or_incomplete_file_is_reported_and_not_cached(self):
        for fmt, raw, message in (('null', [], 'Missing out.cgyro.freq'),
                                  ('out', [], 'expected 24'),
                                  ('bin', np.arange(23.), 'expected 24')):
            with self.subTest(fmt=fmt, length=len(raw)):
                plot = self.plot()
                data, _ = self.case(extract=Mock(return_value=('', fmt, raw)))
                with self.assertRaisesRegex(ValueError, message):
                    plot._plot_frequency_growth(data, 'case', 'Frequency', [2, 3], 2., 3.)
                self.assertFalse(hasattr(data, 'freq'))
                self.assertFalse(plot.ax.lines)

    def test_invalid_shape_and_uncovered_time_window_do_not_silently_return(self):
        for freq, indices, message in ((np.zeros((2, 2, 4)), [1], 'ky grid'),
                                      (np.zeros((2, 3, 0)), [], 'no time samples'),
                                      (np.zeros((2, 3, 2)), [3], 'time window'),
                                      (np.zeros((1, 3)), [0], 'Unsupported frequency shape')):
            with self.subTest(shape=freq.shape, indices=indices):
                data, _ = self.case(freq=freq)
                with self.assertRaisesRegex(ValueError, message):
                    self.plot()._plot_frequency_growth(data, 'case', 'Frequency', indices, 0., 3.)

    def test_division_by_zero_ky_keeps_other_modes_and_reports_all_invalid(self):
        plot = self.plot()
        plot.norm_ky_var = SimpleNamespace(get=lambda: True)
        data, values = self.case()
        plot._plot_frequency_growth(data, 'case', 'Frequency', [2, 3], 2., 3.)
        np.testing.assert_allclose(plot.ax.lines[0].get_ydata(), values[0, 1:, 2:].mean(axis=1) / [.1, .2])
        data.freq[:] = np.nan
        with self.assertRaisesRegex(ValueError, 'No finite frequency values'):
            plot._plot_frequency_growth(data, 'case', 'Frequency', [2, 3], 2., 3.)


if __name__ == '__main__':
    unittest.main()
