"""Positive one-sided ky and FULLT signed-partner plotting regressions."""

import io
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.collections import PathCollection
import numpy as np

from cgyro_comparison_plotting import Plotting
from cgyro_data_export import CgyroDataExportMixin
from cgyro_fullt_trace_core import Mode, TracePanel, TraceResult, TraceSettings
from cgyro_fullt_trace_plot import FigureSettings, build_figure
from cgyro_fullt_style import fullt_color_scale


def variable(value):
    return SimpleNamespace(get=lambda: value)


class PlotHarness(CgyroDataExportMixin, Plotting):
    """Use the same export/plotting MRO as the actual GUI without a Tk window."""

    def __init__(self, *, asym=False, normalize=False, view="vs kxky", ky="-.19"):
        self.fig = Figure(figsize=(9, 5))
        FigureCanvasAgg(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.energy_balance_transfer_asym_var = variable(asym)
        self.energy_balance_transfer_norm_max_var = variable(normalize)
        self.energy_balance_transfer_xaxis_var = variable(view)
        self.energy_balance_transfer_ky_var = variable(ky)
        self.energy_balance_transfer_kx_var = variable("0")
        self._load_fullt_if_needed = Mock(return_value=True)
        self._clear_current_plot_data()


def fullt_case(sign=-1, native=False):
    source, radial, partner, time = np.indices((3, 4, 5, 2))
    values = 100 * source + 10 * radial + partner + time - 220
    payload = np.zeros((2, 3, 4, 5, 1, 2))
    payload[0, :, :, :, 0, :] = values
    data = SimpleNamespace(
        ky=sign * np.arange(3) * .1, n_n=3, n_radial=4,
        p=np.arange(-2, 2), length=2 * np.pi / .1,
        t=np.array([1., 2.]), fullt=payload,
    )
    data.fullt_asym, _ = Plotting()._construct_fullt_asym_from_fullt(payload, data, "fixture")
    if native:
        data.full_t_source_ky = data.ky.copy()
        data.full_t_target_ky = sign * np.arange(-2, 3) * .1
    return data


class KyPlotConventionTests(unittest.TestCase):
    def test_positive_and_negative_requests_select_same_physical_mode(self):
        plot = Plotting()
        for sign in (-1, 1):
            axis = sign * np.array([0., .1, .2, 1.])
            for text, expected in (("-.19", 2), (".19", 2), ("val:-.19", 2),
                                   ("-1", 3), ("1", 3), ("idx:1", 1)):
                with self.subTest(sign=sign, text=text):
                    indices, label = plot._resolve_optional_axis_indices(axis, text, "ky", "case")
                    np.testing.assert_array_equal(indices, [expected])
                    self.assertEqual(float(label), abs(axis[expected]))
            indices, label = plot._resolve_optional_axis_indices(axis, "", "ky", "case")
            np.testing.assert_array_equal(indices, np.arange(4))
            self.assertEqual(label, "average")
        # Radial modes remain signed.
        self.assertEqual(plot._resolve_axis_selection_index([-.2, 0., .2], "-.2", "kx", "case", True), 0)

    def test_fullt_map_labels_markers_and_export_keep_coordinate_value_pairing(self):
        for sign in (-1, 1):
            for native in (False, True):
                for asym in (False, True):
                    with self.subTest(sign=sign, native=native, asym=asym):
                        plot = PlotHarness(asym=asym)
                        data = fullt_case(sign, native)
                        before = data.ky.copy()
                        plot._plot_energy_balance_fullt(data, "case-A", [0, 1], 1., 2.)
                        plot.fig.canvas.draw()
                        self.assertNotIn("FULLT", plot.ax.get_title())
                        self.assertNotIn("asym", plot.ax.get_title())
                        self.assertTrue(plot.ax.get_title().startswith(r"$T_{k_x\rho_s="))
                        self.assertIn(r"k_y\rho_s=0.2", plot.ax.get_title())
                        self.assertEqual(plot.fig.axes[1].get_ylabel(), "")
                        self.assertEqual(plot.fig.axes[1].get_title(), "")
                        texts = [text.get_text() for text in plot.ax.texts]
                        self.assertIn("case: case-A", texts)
                        self.assertFalse(any("source ky index" in text for text in texts))
                        marker = next(item for item in plot.ax.collections if isinstance(item, PathCollection))
                        np.testing.assert_allclose(marker.get_offsets(), [[0., .2]])
                        payload = data.fullt_asym if asym else data.fullt
                        expected_z = payload[0, 2, :, :, 0, :].mean(axis=-1).T
                        datasets = [item for item in plot._collect_current_plot_xyz_datasets()
                                    if item[0].startswith("FULLT")]
                        self.assertEqual(len(datasets), 1)
                        _, x, y, z = datasets[0]
                        expected_x, expected_y = np.meshgrid(data.p * .1, np.arange(-2, 3) * .1)
                        np.testing.assert_allclose(x, expected_x.ravel())
                        np.testing.assert_allclose(y, expected_y.ravel())
                        np.testing.assert_array_equal(z, expected_z.ravel())
                        np.testing.assert_array_equal(data.ky, before)
                        # Raw file export retains the original coordinate convention.
                        np.testing.assert_array_equal(plot._fullt_source_ky_axis(data, 3), before)
                        np.testing.assert_allclose(plot._fullt_target_ky_axis(data, 3, 5), sign * np.arange(-2, 3) * .1)

    def test_fullt_kx_cut_uses_selected_positive_branch(self):
        plot = PlotHarness(asym=True, view="vs kx")
        data = fullt_case()
        plot._plot_energy_balance_fullt(data, "cut", [0, 1], 1., 2.)
        plot.fig.canvas.draw()
        line = plot.ax.lines[0]
        np.testing.assert_allclose(line.get_xdata(), data.p * .1)
        np.testing.assert_array_equal(line.get_ydata(), data.fullt_asym[0, 2, :, 4, 0, :].mean(axis=-1))
        self.assertIn("ky=0.2", line.get_label())
        self.assertNotIn("FULLT", plot.ax.get_title())
        self.assertNotIn("asym", plot.ax.get_ylabel())

    def test_cached_selected_source_ky_is_positive_and_uses_full_partner_grid(self):
        plot = PlotHarness(asym=True)
        data = fullt_case()
        data.fullt_asym = data.fullt_asym[:, 2:3]
        data.fullt_asym_source_ky_axis = np.array([-.2])
        data.fullt_asym_n_source_ky_full = 3
        plot._plot_energy_balance_fullt(data, "cached", [0, 1], 1., 2.)
        self.assertIn(r"k_y\rho_s=0.2", plot.ax.get_title())
        _, _, y, z = plot._collect_current_plot_xyz_datasets()[0]
        np.testing.assert_allclose(np.unique(y), np.arange(-2, 3) * .1)
        np.testing.assert_array_equal(z, data.fullt_asym[0, 0, :, :, 0, :].mean(axis=-1).T.ravel())

    def test_fullt_trace_has_clean_title_annotation_and_colorbar(self):
        plot = PlotHarness(asym=True, view="trace pxqx")
        plot.energy_balance_transfer_kx_var = variable(".1")
        trace = np.zeros((4, 4, 5, 2))
        trace[:, :, 4, :] = np.arange(32).reshape(4, 4, 2)
        plot._load_fullt_trace_source_slice = Mock(return_value={
            "trace": trace, "source_kx_axis": np.arange(-2, 2) * .1,
            "target_kx_axis": np.arange(-2, 2) * .1,
            "target_ky_axis": np.arange(-2, 3) * .1,
            "source_ky": .2, "source_ky_idx": 2, "time_indices": [0, 1],
        })
        plot._plot_energy_balance_fullt(fullt_case(), "trace", [0, 1], 1., 2.)
        plot.fig.canvas.draw()
        self.assertEqual(plot.fig.axes[1].get_ylabel(), "")
        self.assertNotIn("FULLT", plot.ax.get_title())
        self.assertNotIn("asym", plot.ax.get_title())
        self.assertFalse(any("source ky index" in text.get_text() for text in plot.ax.texts))

    def test_normalized_fullt_map_keeps_signed_transfer_and_blank_colorbar_label(self):
        plot = PlotHarness(normalize=True)
        data = fullt_case()
        plot._plot_energy_balance_fullt(data, "normalized", [0, 1], 1., 2.)
        z = plot._collect_current_plot_xyz_datasets()[0][3]
        raw = data.fullt[0, 2, :, :, 0, :].mean(axis=-1).T
        np.testing.assert_allclose(z, (raw / np.max(np.abs(raw))).ravel())
        self.assertAlmostEqual(float(np.nanmin(z)), -1.)
        self.assertEqual(plot.fig.axes[1].get_ylabel(), "")
        self.assertNotIn("FULLT", plot.ax.get_title())

    def test_tools_share_fullt_colors_and_ticks_for_the_same_map(self):
        data = fullt_case()
        raw = data.fullt[0, 2, :, :, 0, :].mean(axis=-1).T
        result = TraceResult("fixture", "synthetic", TraceSettings(ky=.2, steps=1),
                             data.p * .1, np.arange(-2, 3) * .1, data.t, np.arange(2), {}, "fixture",
                             panels=[TracePanel(Mode(0., .2), raw, None)])
        for normalize, limit in ((False, None), (True, None), (False, 9.), (True, .5), (False, 2.4e-4)):
            with self.subTest(normalize=normalize, limit=limit):
                plot = PlotHarness(normalize=normalize)
                plot.energy_balance_transfer_color_max_var = variable("" if limit is None else str(limit))
                plot._plot_energy_balance_fullt(data, "fixture", [0, 1], 1., 2.)
                plot._apply_unified_visual_style_to_figure()
                traced = build_figure(result, FigureSettings(normalize=normalize, color_max=limit,
                                                             panel_scaling="none"))
                plot.fig.canvas.draw()
                FigureCanvasAgg(traced).draw()
                left = plot.ax.collections[0]
                right = traced._fullt_trace_axes[0].collections[0]
                trace_ax = traced._fullt_trace_axes[0]
                self.assertGreater(trace_ax.xaxis.get_zorder(), right.get_zorder())
                self.assertTrue(all(line.get_visible() for line in trace_ax.get_xgridlines()))
                self.assertTrue(all(line.get_visible() for line in trace_ax.get_ygridlines()))
                self.assertEqual(left.get_clim(), right.get_clim())
                np.testing.assert_allclose(left.get_array(), right.get_array())
                np.testing.assert_allclose(left.to_rgba(left.get_array()), right.to_rgba(right.get_array()))
                expected = raw / np.max(np.abs(raw)) if normalize else raw
                np.testing.assert_allclose(left.get_array().reshape(-1, 9)[:, 0], expected.ravel())
                np.testing.assert_allclose(plot.fig.axes[-1].get_yticks(), traced.axes[-1].get_yticks())
                self.assertEqual(plot.fig.axes[-1].yaxis.get_offset_text().get_text(),
                                 traced.axes[-1].yaxis.get_offset_text().get_text())
                if limit == 2.4e-4:
                    for cax in (plot.fig.axes[-1], traced.axes[-1]):
                        multiplier = cax.yaxis.get_offset_text()
                        self.assertIn("10^{-4}", multiplier.get_text().replace("−", "-"))
                        self.assertEqual(multiplier.get_ha(), "left")
                        self.assertEqual(multiplier.get_fontsize(), cax.get_yticklabels()[0].get_fontsize())

    def test_fullt_scale_handles_invalid_limits_and_zero_or_missing_maps(self):
        for limit in ("invalid", 0, -1, np.nan, np.inf):
            with self.subTest(limit=limit), self.assertRaisesRegex(ValueError, "finite and positive"):
                fullt_color_scale([1., -2.], color_max=limit)
        for values in (np.zeros((2, 3)), np.full((2, 3), np.nan)):
            divisor, norm, maximum = fullt_color_scale(values, normalize=True)
            self.assertEqual((divisor, norm.vmin, norm.vmax, maximum), (1., -1., 1., 0.))
        divisor, norm, maximum = fullt_color_scale([0., -2., 1., np.nan], normalize=True, logarithmic=True)
        self.assertEqual((divisor, norm.vmin, norm.vmax, maximum), (2., -1., 1., 2.))
        self.assertAlmostEqual(float(norm(0)), .5)

    def test_negative_linear_spectrum_sorts_coordinates_with_values(self):
        plot = Plotting()
        text = "# ky omega gamma\n-0.3 30 3\n0.0 0 0\n-0.1 10 1\n"
        with patch("builtins.open", return_value=io.StringIO(text)):
            ky, omega, gamma = plot._read_linear_omega_gamma_file("fixture", "case")
        np.testing.assert_allclose(ky, [0., .1, .3])
        np.testing.assert_array_equal(omega, [0., 10., 30.])
        np.testing.assert_array_equal(gamma, [0., 1., 3.])


if __name__ == "__main__":
    unittest.main()
