"""Opt-in Tk layout checks: set CGYRO_TEST_TK_LAYOUT=1 with a display available."""

import os
import unittest
from contextlib import contextmanager
from unittest.mock import patch


@unittest.skipUnless(os.environ.get("CGYRO_TEST_TK_LAYOUT") == "1", "optional Tk layout integration")
class OptionLayoutTests(unittest.TestCase):
    @contextmanager
    def application(self, scaling=1.5):
        import tkinter as tk
        from cgyro_comparison import CGYRO_Comparison

        root = tk.Tk()
        root.attributes("-alpha", 0)
        root.tk.call("tk", "scaling", scaling)
        try:
            with patch.dict(os.environ, {"CGYRO_AUTO_WORKSPACE": ""}), patch.object(
                CGYRO_Comparison, "_check_for_updates_silently", lambda self: None
            ):
                app = CGYRO_Comparison(root)
                root.geometry("1400x1000")
                app.plot_type_var.set("Fluctuation 1D")
                app.fluc_xaxis_var.set("v.s theta")
                app.update_options()
                root.update_idletasks()
                yield app
        finally:
            for timer in root.tk.call("after", "info"):
                root.after_cancel(timer)
            root.destroy()

    def assert_readable(self, app, widgets):
        from tkinter import font as tkfont

        rectangles = []
        for widget in widgets:
            self.assertTrue(widget.winfo_ismapped(), str(widget))
            left = widget.winfo_rootx() - app.options_frame.winfo_rootx()
            right = left + widget.winfo_width()
            self.assertGreaterEqual(left, 0, str(widget))
            self.assertLessEqual(right, app.options_frame.winfo_width(), str(widget))
            if widget.winfo_class() in ("TEntry", "TCombobox"):
                font = tkfont.Font(root=app.root, font=widget.cget("font") or "TkTextFont")
                chrome = widget.winfo_reqwidth() - int(widget.cget("width")) * font.measure("0")
                required = font.measure(widget.get()) + chrome
            else:
                required = widget.winfo_reqwidth()
            self.assertGreaterEqual(widget.winfo_width(), required, str(widget))
            x, y = widget.winfo_rootx(), widget.winfo_rooty()
            rectangles.append((x, y, x + widget.winfo_width(), y + widget.winfo_height()))
        for i, first in enumerate(rectangles):
            for second in rectangles[i + 1:]:
                overlap = min(first[2], second[2]) > max(first[0], second[0]) and min(first[3], second[3]) > max(first[1], second[1])
                self.assertFalse(overlap, (first, second))

    @staticmethod
    def time_controls(app):
        return [app.t_start_label, app.t_start_entry, app.t_end_label, app.t_end_entry,
                app.log_x_check, app.log_y_check, app.clear_time_button]

    @staticmethod
    def theta_controls(app):
        return [app.fluc_field_combo, app.fluc_xaxis_combo, app.fluc_advanced_check,
                app.fluc_theta_kx_label, app.fluc_theta_kx_entry,
                app.fluc_theta_ky_label, app.fluc_theta_ky_entry]

    def test_theta_is_readable_at_narrow_widths_and_dpi_scales(self):
        for scaling in (1.0, 1.5, 2.0):
            with self.subTest(scaling=scaling), self.application(scaling) as app:
                for width in (340, 400, 500):
                    with self.subTest(width=width):
                        app.main_pane.sashpos(0, width)
                        app.root.update_idletasks()
                        self.assert_readable(app, self.time_controls(app) + self.theta_controls(app))
                        self.assertEqual(app.fluc_theta_kx_label.cget("text"), "kx:")
                        self.assertEqual(app.fluc_theta_ky_label.cget("text"), "ky:")

    def test_mode_switches_keep_values_and_restore_theta_columns(self):
        with self.application() as app:
            app.main_pane.sashpos(0, 340)
            app.t_start_var.set("500")
            app.t_end_var.set("1200")
            app.fluc_theta_kx_var.set("0")
            app.fluc_theta_ky_var.set("0.3")
            for mode in ("v.s ky", "v.s kx", "v.s Time", "fft", "v.s theta"):
                app.fluc_xaxis_var.set(mode)
                app.update_options()
                app.root.update_idletasks()
                self.assert_readable(app, self.time_controls(app))
            self.assert_readable(app, self.theta_controls(app))
            self.assertEqual(app.fluc_theta_kx_entry.grid_info()["columnspan"], 1)
            self.assertEqual(app.fluc_theta_ky_entry.grid_info()["columnspan"], 1)
            self.assertEqual(app.fluc_theta_ky_entry.grid_info()["column"], 3)
            self.assertEqual((app.t_start_var.get(), app.t_end_var.get()), ("500", "1200"))
            self.assertEqual((app.fluc_theta_kx_var.get(), app.fluc_theta_ky_var.get()), ("0", "0.3"))
            app.fluc_advanced_var.set(True)
            app.update_options()
            app.root.update_idletasks()
            self.assertTrue(app.fluc_advanced_button.winfo_ismapped())
            self.assertFalse(app.fluc_axis_frame.winfo_ismapped())
            self.assert_readable(app, [app.fluc_advanced_check, app.fluc_advanced_button])
            app.fluc_advanced_var.set(False)
            app.update_options()
            app.root.update_idletasks()
            self.assertFalse(app.fluc_advanced_button.winfo_ismapped())
            self.assert_readable(app, self.theta_controls(app))

    def test_time_controls_survive_plot_type_switches_and_clear(self):
        with self.application() as app:
            app.main_pane.sashpos(0, 340)
            for plot_type in ("Flux", "Fluctuation 2D", "Energy balance", "Frequency", "Fluctuation 1D"):
                app.plot_type_var.set(plot_type)
                app.update_options()
                app.root.update_idletasks()
                self.assert_readable(app, self.time_controls(app))
            app.t_start_var.set("500")
            app.t_end_var.set("1200")
            app.clear_time_button.invoke()
            app.root.update_idletasks()
            self.assertEqual((app._get_entry_value("t_start_var"), app._get_entry_value("t_end_var")), ("", ""))
            self.assertEqual((app.t_start_entry.get(), app.t_end_entry.get()), ("50% End", "End"))
            self.assertEqual((app._get_entry_value("fluc_theta_kx_var"), app._get_entry_value("fluc_theta_ky_var")), ("", ""))
            self.assertEqual((app.fluc_theta_kx_entry.get(), app.fluc_theta_ky_entry.get()), ("Avg", "Avg"))
            self.assert_readable(app, self.time_controls(app) + self.theta_controls(app))
            app.fluc_theta_kx_entry.event_generate("<FocusIn>")
            self.assertEqual(app.fluc_theta_kx_entry.get(), "")
            app.fluc_theta_kx_entry.insert(0, "0.25")
            self.assertEqual(app._get_entry_value("fluc_theta_kx_var"), "0.25")


if __name__ == "__main__":
    unittest.main()
