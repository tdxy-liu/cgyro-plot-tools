"""Native-file integration and publication/export checks for FULLT Trace Studio."""

from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg

from cgyro_fullt_trace_core import (FullTSource, Mode, TraceCancelled, TraceError,
                                    TraceSettings, _enable_private_decoder, select_peak, trace_case)
from cgyro_fullt_trace_plot import FigureSettings, build_figure, export_bundle
from cgyro_fullt_reader import FTZWriter, Layout


def make_case(folder, *, compressed=True, sign=-1, ns=8):
    """Explicit SYNTHETIC fixture; never used as application fallback/demo data."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    m = Layout(8, 4, 1, ns, 4, .1, sign * .1)
    blocks = [np.zeros(m.shape, dtype=m.dtype) for _ in range(m.ny)]
    # Distinct positive and negative four-panel paths; include a negative-ky
    # partner that requires reversing BOTH kx and ky at the next step.
    paths = [((4, 1), (5, 4), 10), ((5, 1), (6, 5), 8),
             ((6, 2), (5, 0), 6), ((3, 3), (4, 3), 4),
             ((4, 1), (3, 4), -9), ((3, 1), (2, 5), -7),
             ((2, 2), (4, 6), -5), ((4, 3), (4, 3), -3)]
    for (source_x, source_y), (partner_x, partner_y), value in paths:
        if ns == 1 and source_x != 4:
            continue
        blocks[source_y][partner_x, partner_y, 0, 0 if ns == 1 else source_x] = value
    times = np.array([100., 200., 300., 400.])
    (folder / "input.cgyro").write_text(
        f"FULL_T_ASYM_COMPRESSION={int(compressed)}\nFULL_T_REAL_ONLY=1\n"
        f"FULL_T_KX0={int(ns == 1)}\nHIPREC_FLAG=0\nTOROIDALS_PER_PROC=2\n")
    grid = [4, 2, 4, 8, 2, 1, 1, 2, 2 * np.pi / m.dx, 0, 1]
    grid += list(np.arange(-4, 4)) + [0] * 4 + [0] * 8 + list(np.arange(4) * m.dy)
    np.savetxt(folder / "out.cgyro.grids", grid)
    np.savetxt(folder / "out.cgyro.time", np.c_[times, np.zeros((4, 3))])
    if compressed:
        _enable_private_decoder()
        with FTZWriter(folder / "bin.cgyro.fullt_asym.ftz", m) as writer:
            for i, time in enumerate(times, 1):
                writer.append([block * i for block in blocks], time)
    else:
        with (folder / "bin.cgyro.fullt_asym").open("wb") as handle:
            for i in range(1, 5):
                for rank in range(2):
                    handle.write(np.stack([b * i for b in blocks[rank * 2:rank * 2 + 2]], axis=3).tobytes(order="F"))
    return folder


class FullTTraceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.case = make_case(self.folder / "synthetic_fixture")
        self.settings = TraceSettings(ky=.1, steps=4, t_start=100., t_end=300.)

    def test_four_step_positive_chain_and_negative_partner(self):
        result = trace_case(self.case, self.settings)
        np.testing.assert_allclose([[p.mode.kx, p.mode.ky] for p in result.panels],
                                   [[0, .1], [.1, .1], [.2, .2], [-.1, .3]])
        np.testing.assert_allclose([p.peak.value for p in result.panels], [20, 16, 12, 8])
        self.assertEqual(result.stop_code, "completed")
        self.assertTrue(result.panels[2].peak.conjugated)
        np.testing.assert_allclose([result.panels[2].peak.kx, result.panels[2].peak.ky], [.1, -.3])

    def test_negative_chain_preserves_signed_T(self):
        result = trace_case(self.case, replace(self.settings, direction="negative"))
        np.testing.assert_allclose([[p.mode.kx, p.mode.ky] for p in result.panels],
                                   [[0, .1], [-.1, .1], [-.2, .2], [0, .3]])
        np.testing.assert_allclose([p.peak.value for p in result.panels], [-18, -14, -10, -6])

    def test_raw_and_ftz_and_ky_conventions_agree(self):
        expected = trace_case(self.case, self.settings)
        for compressed, sign in ((False, -1), (False, 1), (True, 1)):
            case = make_case(self.folder / f"case_{compressed}_{sign}", compressed=compressed, sign=sign)
            result = trace_case(case, replace(self.settings, ky=-.11))
            for got, want in zip(result.panels, expected.panels):
                self.assertEqual(got.mode, want.mode)
                np.testing.assert_array_equal(got.values, want.values)

    def test_initial_snap_default_window_and_empty_window(self):
        result = trace_case(self.case, TraceSettings(kx=.012, ky=-.11, steps=1))
        self.assertEqual(result.panels[0].mode, Mode(0., .1))
        np.testing.assert_array_equal(result.times, [300., 400.])
        self.assertAlmostEqual(result.panels[0].peak.value, 35.)
        with self.assertRaisesRegex(TraceError, "No saved samples"):
            trace_case(self.case, replace(self.settings, t_start=500., t_end=600.))

    def test_shared_ky_time_means_are_reused(self):
        with FullTSource(self.case, self.settings) as source:
            with patch.object(source.reader, "block", wraps=source.reader.block) as block:
                first = source.mean_map((4, 1))
                second = source.mean_map((5, 1))
                self.assertEqual(block.call_count, 3)
                self.assertAlmostEqual(first[4, 5], 20.)
                self.assertAlmostEqual(second[5, 6], 16.)

    def test_repeated_mode_policies_and_missing_sign(self):
        for policy, count, code in (("stop", 1, "cycle"), ("allow", 4, "completed"), ("skip", 1, "no_extremum")):
            with self.subTest(policy=policy):
                settings = replace(self.settings, kx=0., ky=.3, direction="negative", cycles=policy)
                with patch.object(FullTSource, "mean_map") as read:
                    values = np.zeros((7, 8))
                    values[6, 4] = -1.
                    read.return_value = values
                    result = trace_case(self.case, settings)
                self.assertEqual(len(result.panels), count)
                self.assertEqual(result.stop_code, code)
        result = trace_case(self.case, replace(self.settings, ky=0.))
        self.assertEqual(result.stop_code, "no_extremum")

    def test_missing_full_kx_does_not_silently_snap_next_to_zero(self):
        case = make_case(self.folder / "kx0_only", ns=1)
        result = trace_case(case, self.settings)
        self.assertEqual(len(result.panels), 1)
        self.assertEqual(result.stop_code, "unavailable_mode")
        self.assertIn("FULL_T_KX0=1", result.stop_reason)
        with self.assertRaisesRegex(TraceError, "kx=0"):
            trace_case(case, replace(self.settings, kx=.1))

    def test_skip_excludes_both_half_plane_representations_of_visited_mode(self):
        kx = np.array([-.1, 0., .1])
        ky = np.array([-.1, 0., .1])
        values = np.array([[10., 0., 0.], [0., 0., 0.], [4., 0., 8.]])
        peak, count = select_peak(values, kx, ky, "positive", [Mode(.1, .1)])
        self.assertEqual(count, 2)
        self.assertEqual(peak.next_mode, Mode(-.1, .1))
        self.assertEqual(peak.value, 4.)

    def test_nonfinite_extrema_and_unavailable_nyquist_do_not_corrupt_chain(self):
        values = np.zeros((7, 8))
        values[0, 7] = np.nan
        values[1, 0] = 50.
        with patch.object(FullTSource, "mean_map", return_value=values):
            result = trace_case(self.case, self.settings)
        self.assertEqual(result.stop_code, "unavailable_mode")
        self.assertAlmostEqual(result.panels[0].peak.next_mode.kx, .4)
        self.assertEqual(result.metadata()["panels"][0]["nonfinite_points"], 1)

    def test_cancel_closes_reader_and_produces_no_result(self):
        cancel = threading.Event()
        cancel.set()
        with self.assertRaises(TraceCancelled):
            trace_case(self.case, self.settings, cancel=cancel)

    def test_vertical_panels_use_one_shared_norm_and_no_colorbar_title(self):
        result = trace_case(self.case, self.settings)
        for width in (89., 183.):
            fig = build_figure(result, FigureSettings(width_mm=width))
            FigureCanvasAgg(fig).draw()
            axes = fig._fullt_trace_axes
            self.assertEqual(len(axes), 4)
            self.assertTrue(all(a.collections[0].norm is axes[0].collections[0].norm for a in axes))
            self.assertAlmostEqual(axes[0].collections[0].norm.vmax, 20.)
            self.assertEqual(fig.axes[-1].get_ylabel(), "")
            self.assertTrue(all(axes[i].get_position().y0 > axes[i + 1].get_position().y0 for i in range(3)))
            for i, panel in enumerate(result.panels):
                np.testing.assert_array_equal(axes[i].collections[0].get_array(), panel.values.ravel())
            # All header text remains inside the physical figure.
            renderer = fig.canvas.get_renderer()
            for ax in fig.axes:
                boxes = []
                for text in ax.texts:
                    bounds = text.get_window_extent(renderer)
                    self.assertGreaterEqual(bounds.x0, 0)
                    self.assertLessEqual(bounds.x1, fig.bbox.width)
                    self.assertFalse(any(bounds.overlaps(old) for old in boxes))
                    boxes.append(bounds)
                for label in (ax.xaxis.label, ax.yaxis.label):
                    if label.get_text():
                        self.assertGreaterEqual(label.get_window_extent(renderer).x0, 0)

    @unittest.skipUnless(os.environ.get("CGYRO_TEST_TK_LAYOUT") == "1", "optional Tk integration")
    def test_gui_traces_and_restyles_without_rereading_and_scrolls_options(self):
        import tkinter as tk
        from cgyro_fullt_trace_ui import FullTTraceApp

        errors_patch = patch("cgyro_fullt_trace_ui.messagebox.showerror")
        errors = errors_patch.start()
        self.addCleanup(errors_patch.stop)
        warnings_patch = patch("cgyro_fullt_trace_ui.messagebox.showwarning")
        warnings = warnings_patch.start()
        self.addCleanup(warnings_patch.stop)
        root = tk.Tk()
        root.attributes("-alpha", 0)
        try:
            app = FullTTraceApp(root, [self.case], trace=self.settings)
            root.geometry("1040x720")
            root.update()
            app.run_trace()
            deadline = time.monotonic() + 20
            while app.busy and time.monotonic() < deadline:
                root.update()
                time.sleep(.02)
            self.assertFalse(app.busy)
            self.assertEqual(len(app.results[0].panels), 4)
            self.assertEqual(len(app.table.get_children()), 4)
            app.vars["normalize"].set(True)
            with patch("cgyro_fullt_trace_ui.trace_case", side_effect=AssertionError("style must use cached maps")):
                app.redraw()
            self.assertEqual(app.canvas.figure._fullt_trace_scale["limit"], 1.)
            errors.assert_not_called()
            warnings.assert_not_called()
            app.option_tabs.select(1)
            root.update()
            frame, canvas = app.option_canvases[1]
            self.assertGreater(frame.winfo_height(), canvas.winfo_height())
            canvas.yview_moveto(1.)
            root.update()
            self.assertAlmostEqual(canvas.yview()[1], 1.)
            before = canvas.yview()
            app.case_list.event_generate("<MouseWheel>", delta=120)
            root.update()
            self.assertEqual(canvas.yview(), before)
            malformed = self.folder / "bad_session.json"
            malformed.write_text("[]", encoding="utf-8")
            with patch("cgyro_fullt_trace_ui.filedialog.askopenfilename", return_value=str(malformed)), \
                    patch("cgyro_fullt_trace_ui.messagebox.showerror") as error:
                app.load_session()
                error.assert_called_once()
            self.assertEqual(app.cases, [str(self.case.absolute())])
        finally:
            for timer in root.tk.call("after", "info"):
                root.after_cancel(timer)
            root.destroy()

    def test_export_keeps_raw_values_and_never_overwrites(self):
        result = trace_case(self.case, self.settings)
        settings = FigureSettings(dpi=300, normalize=True)
        output = export_bundle(result, settings, self.folder / "out", formats=("svg", "png"))
        second = export_bundle(result, settings, self.folder / "out", formats=("svg",))
        self.assertNotEqual(output, second)
        self.assertIn("<text", (output / "trace.svg").read_text(encoding="utf-8"))
        metadata = json.loads((output / "trace.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["shared_color_scale"]["divisor"], 20.)
        with np.load(output / "maps.npz") as data:
            np.testing.assert_array_equal(data["step_01"], result.panels[0].values)
        table = np.loadtxt(output / "plot_data.csv", delimiter=",", skiprows=1)
        np.testing.assert_allclose(table[:, 4], table[:, 3] / 20.)

    def test_cli_batch_and_bad_case_failure(self):
        script = str(Path(__file__).resolve().parents[1] / "cgyro_fullt_trace.py")
        completed = subprocess.run([sys.executable, script, "--case", str(self.case), "--ky", ".1", "--steps", "3",
                                    "--output", str(self.folder / "cli"), "--formats", "svg"],
                                   capture_output=True, text=True, timeout=45)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Exported 3 panels", completed.stdout)
        failed = subprocess.run([sys.executable, script, "--case", str(self.folder / "missing"), "--output", str(self.folder / "cli")],
                                capture_output=True, text=True, timeout=30)
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("No native", failed.stderr)


if __name__ == "__main__":
    unittest.main()
