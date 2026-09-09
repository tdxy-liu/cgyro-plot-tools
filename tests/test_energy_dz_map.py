"""DZ map regression checks; run with python -m unittest discover -s tests."""

import unittest
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from cgyro_comparison_plotting import Plotting
from cgyro_data_export import CgyroDataExportMixin


def value_var(value):
    return SimpleNamespace(get=lambda: value)


class PlotHarness(CgyroDataExportMixin, Plotting):
    def __init__(self):
        self.energy_balance_spec_var = value_var("Total (-1)")
        self.energy_balance_single_quantity_var = value_var("DZ")
        self.energy_balance_single_xaxis_var = value_var("vs kxky")
        self.energy_balance_single_norm_var = value_var("Min T")
        self.fig = Figure(figsize=(6, 4))
        FigureCanvasAgg(self.fig)
        self.ax = self.fig.add_subplot(111)
        self._clear_current_plot_data()


class DzMapTests(unittest.TestCase):
    def setUp(self):
        self.plot = PlotHarness()
        # Distinct species, radial, ky, and time dependence exposes axis swaps.
        s, r, k, t = np.indices((2, 5, 3, 4))
        base = 10 * s + 2 * r - 3 * k + t
        triad = np.full((2, 2, 5, 8, 3, 4), 10000, dtype=np.float32)
        triad[0, :, :, 5, :, :] = -base
        triad[0, :, :, 6, :, :] = 2 * base + 4
        triad[0, :, :, 7, :, :] = -6
        self.data = SimpleNamespace(
            triad=triad,
            kx=np.array([0.2, -0.2, 0.1, 0.0, -0.1]),
            kynorm=np.array([0.3, 0.0, 0.1]),
            t=np.array([10.0, 20.0, 30.0, 40.0]),
        )
        self.time_indices = np.array([0, 2])
        # Dr + Dtheta + Dc = base - 2; mean selected time index is 1.
        self.expected = (10 * s + 2 * r - 3 * k - 1)[..., 0]

    def map(self, quantity="DZ", normalization="none"):
        return self.plot._compute_energy_balance_single_vs_kxky(
            self.data, "synthetic", None, quantity, self.time_indices,
            normalize_mode=normalization,
        )

    def test_signed_species_and_time_average(self):
        kx_order = np.argsort(self.data.kx)
        ky_order = np.argsort(self.data.kynorm)
        for text, expected in (
            ("Main ion (0)", self.expected[0]),
            ("Electron (1)", self.expected[1]),
            ("Total (-1)", self.expected.sum(axis=0)),
        ):
            with self.subTest(species=text):
                self.plot.energy_balance_spec_var = value_var(text)
                x, y, z = self.map()
                np.testing.assert_array_equal(x[:, 0], self.data.kx[kx_order])
                np.testing.assert_array_equal(y[0, :], self.data.kynorm[ky_order])
                np.testing.assert_array_equal(z, expected[np.ix_(kx_order, ky_order)])

    def test_native_kx_excludes_first_radial_row(self):
        self.data.kx = self.data.kx[1:]
        self.plot.energy_balance_spec_var = value_var("Main ion (0)")
        _, _, z = self.map()
        expected = self.expected[0, 1:, :]
        expected = expected[np.ix_(np.argsort(self.data.kx), np.argsort(self.data.kynorm))]
        np.testing.assert_array_equal(z, expected)

    def test_component_sum_and_1d_spectra_agree(self):
        x, y, z = self.map()
        components = sum(self.map(quantity)[2] for quantity in ("Dr", "Dtheta", "Dc"))
        np.testing.assert_array_equal(z, components)
        ky, spectrum = self.plot._compute_energy_balance_single_vs_ky(
            self.data, "synthetic", None, "DZ", self.time_indices,
        )
        np.testing.assert_array_equal(ky, y[0, :])
        np.testing.assert_allclose(spectrum, z.sum(axis=0))
        for j, ky_value in enumerate(ky):
            kx, spectrum, _ = self.plot._compute_energy_balance_single_vs_kx(
                self.data, "synthetic", ky_value, "DZ", self.time_indices,
            )
            np.testing.assert_array_equal(kx, x[:, 0])
            np.testing.assert_allclose(spectrum, z[:, j])

    def test_stale_t_normalization_does_not_scale_dz(self):
        for mode in ("min", "max"):
            with self.subTest(normalization=mode):
                np.testing.assert_array_equal(self.map()[2], self.map(normalization=mode)[2])

    def test_render_and_export_preserve_physical_axes_and_values(self):
        x, y, z = self.map()
        self.plot._plot_energy_balance_single_mode(
            self.data, "synthetic", self.time_indices, 10.0, 30.0,
        )
        self.plot.fig.canvas.draw()
        self.assertEqual(len(self.plot.ax.collections), 1)
        mesh = self.plot.ax.collections[0]
        np.testing.assert_array_equal(np.asarray(mesh.get_array()).reshape(z.T.shape), z.T)
        self.assertEqual(mesh.norm.vcenter, 0.0)
        self.assertAlmostEqual(mesh.norm.vmin, -mesh.norm.vmax)
        self.assertEqual(self.plot.fig.axes[1].get_ylabel(), r"$D_Z$")
        self.assertIn("DZ", self.plot.ax.get_title())
        self.assertIn("Avg", self.plot.ax.get_title())
        datasets = self.plot._collect_current_plot_xyz_datasets()
        self.assertEqual(len(datasets), 1)
        label, export_x, export_y, export_z = datasets[0]
        self.assertIn("synthetic DZ", label)
        np.testing.assert_array_equal(export_x, x.T.ravel())
        np.testing.assert_array_equal(export_y, y.T.ravel())
        np.testing.assert_array_equal(export_z, z.T.ravel())

    def test_zero_map_renders(self):
        self.data.triad.fill(0)
        self.plot._plot_energy_balance_single_mode(
            self.data, "zero", self.time_indices, 10.0, 30.0,
        )
        self.plot.fig.canvas.draw()
        self.assertEqual(len(self.plot.ax.collections), 1)
        np.testing.assert_array_equal(self.plot.ax.collections[0].get_array(), 0)


if __name__ == "__main__":
    unittest.main()
