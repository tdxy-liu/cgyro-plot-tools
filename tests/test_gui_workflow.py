"""Opt-in native Tk workflow and language regression checks."""

from contextlib import contextmanager
from pathlib import Path
import json
import os
import tempfile
import unittest
from unittest.mock import patch


@unittest.skipUnless(os.environ.get("CGYRO_TEST_TK_LAYOUT") == "1", "optional Tk workflow integration")
class GuiWorkflowTests(unittest.TestCase):
    def test_frequency_action_plots_fast_case_and_shows_missing_case_warning(self):
        from types import SimpleNamespace
        from unittest.mock import Mock
        import numpy as np
        with self.app() as app:
            data = SimpleNamespace(ky=np.array([0., -.1]), t=np.arange(4.),
                                   extract=Mock(return_value=('', 'out', np.arange(16.))))
            missing = SimpleNamespace(ky=np.array([0., -.1]), t=np.arange(4.), dir='/missing/',
                                      extract=Mock(return_value=('', 'null', [])))
            app.cases = {'good': data, 'missing': missing}
            for name in app.cases:
                app.case_listbox.insert('end', name)
            app.plot_type_var.set('Frequency')
            with patch('cgyro_comparison_plotting.messagebox.showwarning') as warning:
                app.plot_comparison()
                self.assertEqual(len(app.ax.lines), 1)
                warning.assert_called_once()
                self.assertIn('missing: Missing out.cgyro.freq', warning.call_args.args[1])
                self.assertNotIn('good:', warning.call_args.args[1])
            self.assertEqual(app.plot_button.cget('text'), 'Plot')
            app.case_listbox.selection_set(0)
            with patch('cgyro_comparison_plotting.messagebox.showwarning') as warning:
                app.plot_comparison()
                warning.assert_not_called()
            data.extract.assert_called_once()

    @contextmanager
    def app(self, kind="main", scaling=1.5):
        import tkinter as tk
        from cgyro_comparison import CGYRO_Comparison
        from cgyro_fullt_trace_ui import FullTTraceApp
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {
                "CGYRO_UI_LANGUAGE": "en", "CGYRO_AUTO_WORKSPACE": "",
                "CGYRO_UI_CONFIG": str(Path(directory) / "ui.json")}), \
                patch.object(CGYRO_Comparison, "_check_for_updates_silently", lambda self: None):
            root = tk.Tk()
            root.attributes("-alpha", 0)
            root.tk.call("tk", "scaling", scaling)
            try:
                app = CGYRO_Comparison(root) if kind == "main" else FullTTraceApp(root)
                root.geometry("1280x900")
                root.update_idletasks()
                yield app
            finally:
                # Finish Matplotlib's widget-owned idle draw before cancelling root timers.
                root.update_idletasks()
                for timer in root.tk.call("after", "info"):
                    root.after_cancel(timer)
                root.destroy()

    def assert_form_fits(self, form):
        parent = form.parent
        if not parent.winfo_ismapped():
            return
        boxes = []
        for row in form.rows:
            for widget in row:
                if not widget.winfo_ismapped():
                    continue
                x, y = widget.winfo_x(), widget.winfo_y()
                right, bottom = x + widget.winfo_width(), y + widget.winfo_height()
                self.assertGreaterEqual(x, 0, str(widget))
                self.assertLessEqual(right, parent.winfo_width(), str(widget))
                if widget.winfo_class() in ("TEntry", "TCombobox"):
                    self.assertGreater(widget.winfo_width(), 45, str(widget))
                elif widget.winfo_class() == "TLabel":
                    self.assertGreaterEqual(widget.winfo_width(), widget.winfo_reqwidth(), str(widget))
                box = (x, y, right, bottom)
                for previous in boxes:
                    self.assertFalse(min(right, previous[2]) > max(x, previous[0]) and
                                     min(bottom, previous[3]) > max(y, previous[1]), (widget, box, previous))
                boxes.append(box)

    def test_all_diagnostics_reflow_without_overlap_at_normal_and_narrow_sizes(self):
        for scaling in (1., 1.5, 2.):
            with self.subTest(scaling=scaling), self.app(scaling=scaling) as app:
                modes = [(name, None) for name in app.plot_type_combo._values]
                modes += [("Energy balance", name) for name in app._ENERGY_BALANCE_MODE_OPTIONS]
                modes += [("Fluctuation 1D", name) for name in app._FLUC_XAXIS_OPTIONS]
                for geometry, width in (("920x640", 320), ("1280x900", 410), ("1600x1000", 540)):
                    app.root.geometry(geometry)
                    app.root.update_idletasks()
                    app.main_pane.sashpos(0, width)
                    for mode, detail in modes:
                        app.plot_type_var.set(mode)
                        if detail:
                            (app.energy_balance_mode_var if mode == "Energy balance" else app.fluc_xaxis_var).set(detail)
                        app.update_options()
                        app.root.update_idletasks()
                        for form in app._option_forms.values():
                            self.assert_form_fits(form)
                        self.assert_form_fits(app._time_form)
                        self.assertLess(app.plot_button.winfo_rooty() + app.plot_button.winfo_height(), app.root.winfo_rooty() + app.root.winfo_height())

    def test_language_switch_keeps_model_values_selection_and_plot(self):
        from tkinter import font as tkfont
        from cgyro_ui_language import get_language
        for kind in ("main", "trace"):
            with self.subTest(kind=kind), self.app(kind) as app:
                lang = get_language(app.root)
                self.assertEqual(tkfont.Font(root=app.root, name="TkDefaultFont", exists=True).cget("size"), 11)
                if kind == "main":
                    app.plot_type_var.set("Fluctuation 1D")
                    app.fluc_xaxis_var.set("v.s theta")
                    app.t_start_var.set("300")
                    app.t_end_var.set("600")
                    app.lock_case_selection_var.set(True)
                    app.update_options()
                    original = app._get_workspace_state()
                    figure = app.fig
                else:
                    app.vars["panel_multipliers"].set("1, 1, 5")
                    app._add_path("case with spaces")
                    original = app._options()
                with patch.object(tkfont, "families", side_effect=AssertionError("font inventory must be cached")):
                    for code in ("zh", "en", "zh", "en"):
                        lang.variable.set(code)
                        lang.apply()
                        app.root.update_idletasks()
                        if kind == "main":
                            self.assertEqual(app._get_workspace_state(), original)
                            self.assertIs(app.fig, figure)
                            self.assertEqual(app.plot_type_combo.get(), lang.text("Fluctuation 1D"))
                            for artist, chinese, english in app._empty_plot_labels:
                                self.assertEqual(artist.get_text(), chinese if lang.chinese else english)
                        else:
                            self.assertEqual(app._options(), original)
                            self.assertEqual(app.case_list.curselection(), (0,))
                            self.assertEqual(app.vars["panel_scaling"].get(), "Custom multipliers")
                            self.assertIn(lang.text("CGYRO Transfer Explorer"), app.root.title())
                            self.assertEqual(app.option_tabs.tab(0, "text"), lang.text("Trace settings"))
                            self.assertEqual(app.option_tabs.tab(1, "text"), lang.text("Figure style"))
                            views = app.table.master.master
                            self.assertEqual(views.tab(0, "text"), lang.text("Stacked figure"))
                            self.assertEqual(views.tab(1, "text"), lang.text("Trace coordinates and values"))
                self.assertEqual(json.loads(lang.path.read_text())["language"], "en")

    def test_chinese_switch_and_saved_preference_do_not_depend_on_font_detection(self):
        from cgyro_ui_language import get_language, UiLanguage
        for kind in ("main", "trace"):
            with self.subTest(kind=kind), self.app(kind) as app:
                # Reproduce a missed CJK font even on machines with Noto installed.
                app.root._cgyro_theme["cjk"] = False
                lang = get_language(app.root)
                lang.variable.set("zh")
                lang.apply()
                app.root.update_idletasks()
                self.assertTrue(lang.chinese)
                if kind == "main":
                    self.assertEqual(app.plot_button.cget("text"), "绘制")
                    self.assertEqual(app.plot_type_combo.get(), "频率")
                else:
                    self.assertEqual(app.option_tabs.tab(0, "text"), "追踪设置")
                    self.assertIn("CGYRO 转移分析", app.root.title())
                self.assertEqual(json.loads(lang.path.read_text())["language"], "zh")
                with patch.dict(os.environ):
                    os.environ.pop("CGYRO_UI_LANGUAGE", None)
                    restored = UiLanguage(app.root)
                    self.assertTrue(restored.chinese)
                    self.assertEqual(restored.text("File"), "文件")
                lang.variable.set("en")
                lang.apply()
                self.assertFalse(lang.chinese)

    def test_language_switch_preserves_rendered_trace_and_selected_result(self):
        from dataclasses import replace
        from cgyro_fullt_trace_core import TraceSettings, trace_case
        from cgyro_ui_language import get_language
        from test_fullt_trace import make_case
        with self.app("trace") as app:
            lang = get_language(app.root)
            case = make_case(lang.path.parent / "synthetic", compressed=False)
            result = trace_case(case, TraceSettings(ky=.1, steps=2, t_start=100., t_end=300.))
            app.results = [result, replace(result, case="second case")]
            app._refresh_result_choices()
            app.result_combo.current(1)
            app.redraw()
            canvas = app.canvas
            figure = canvas.figure
            with patch("cgyro_fullt_trace_ui.build_figure", side_effect=AssertionError("language must not redraw")), \
                    patch("cgyro_fullt_trace_ui.trace_case", side_effect=AssertionError("language must not reread data")):
                for code in ("zh", "en"):
                    lang.variable.set(code)
                    lang.refresh()
                    app.root.update_idletasks()
                    self.assertEqual(app.result_combo.current(), 1)
                    self.assertIn("second case", app.result_var.get())
                    self.assertIn("个面板" if lang.chinese else "panels", app.result_var.get())
                    self.assertIs(app.canvas, canvas)
                    self.assertIs(app.canvas.figure, figure)
                    row = app.table.item(app.table.get_children()[0], "values")
                    self.assertEqual(row[1:6], ("0", "0.1", "0.1", "0.1", "20"))
                    self.assertEqual(row[-1], "否" if lang.chinese else "No")
                    status = app._result_status(result, figure._fullt_trace_scale)
                    self.assertIn("倍率" if lang.chinese else "multipliers", status)

    def test_translated_combo_choices_update_model_and_dynamic_choices(self):
        from cgyro_ui_language import get_language
        with self.app() as app:
            lang = get_language(app.root)
            lang.variable.set("zh")
            lang.apply()
            app.plot_type_combo.current(2)
            self.assertEqual(app.plot_type_var.get(), "Flux")
            app.species_combo["values"] = ["Main Ion (Z=1, M=2)", "Electron"]
            app.species_combo.current(1)
            self.assertEqual(app.species_var.get(), "Electron")
            app.plot_type_var.set("Growth Rate")
            self.assertEqual(app.plot_type_combo.get(), lang.text("Growth Rate"))

    def test_formula_notes_stay_lazy_until_opened(self):
        with self.app() as app:
            app.plot_type_var.set("Fluctuation 1D")
            app.fluc_xaxis_var.set("v.s theta")
            app.update_options()
            app.root.update_idletasks()
            self.assertFalse(app.fluc_formula_frame.expanded)
            self.assertFalse(hasattr(app.fluc_formula_canvas, "_formula_render_key"))
            app.fluc_formula_frame.toggle()
            app.root.update_idletasks()
            self.assertTrue(hasattr(app.fluc_formula_canvas, "_formula_render_key"))

    def test_time_mode_uses_active_controls_and_preserves_manual_input(self):
        import numpy as np
        with self.app() as app:
            app.t_start_var.set("300")
            app.t_end_var.set("600")
            time = np.arange(1001.)
            app.time_mode_var.set("Last percent")
            app.time_percent_var.set("20")
            self.assertEqual(app._get_time_indices(time)[1:], (800., 1000.))
            self.assertTrue(app.t_start_entry.instate(["disabled"]))
            app.time_mode_var.set("Full range")
            self.assertEqual(app._get_time_indices(time)[1:], (0., 1000.))
            app.time_mode_var.set("Manual Start/End")
            self.assertEqual(app._get_time_indices(time)[1:], (300., 600.))
            self.assertFalse(app.t_start_entry.instate(["disabled"]))

    def test_transfer_forms_fit_and_language_keeps_scroll_and_inputs(self):
        from cgyro_ui_language import get_language
        for scaling in (1., 2.):
            with self.app("trace", scaling=scaling) as app:
                app.root.geometry("920x640")
                app.root.update_idletasks()
                app.main_pane.sashpos(0, 320)
                for tab in (0, 1):
                    app.option_tabs.select(tab)
                    for code in ("en", "zh"):
                        lang = get_language(app.root)
                        lang.variable.set(code)
                        lang.refresh()
                        app.root.update_idletasks()
                        for form in app._forms:
                            self.assert_form_fits(form)
                self.assertFalse(app.export_button.instate(["!disabled"]))


if __name__ == "__main__":
    unittest.main()
