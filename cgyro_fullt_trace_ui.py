"""Independent Tk interface for FULLT tracing; workers never call Tk."""

from dataclasses import asdict
from pathlib import Path
import json
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

try:
    from cgyro_fullt_trace_core import TraceCancelled, TraceSettings, trace_case
    from cgyro_fullt_trace_plot import FigureSettings, build_figure, export_bundle, suggested_bundle_name
    from cgyro_io import write_json_atomic
    from cgyro_ui_language import get_language, LocalizedCombobox
    from cgyro_ui_common import configure_theme, ui_text, FormLayout, CollapsibleSection, bind_local_scroll, show_interface_guide
except ImportError:
    from .cgyro_fullt_trace_core import TraceCancelled, TraceSettings, trace_case
    from .cgyro_fullt_trace_plot import FigureSettings, build_figure, export_bundle, suggested_bundle_name
    from .cgyro_io import write_json_atomic
    from .cgyro_ui_language import get_language, LocalizedCombobox
    from .cgyro_ui_common import configure_theme, ui_text, FormLayout, CollapsibleSection, bind_local_scroll, show_interface_guide


DIRECTIONS = {"Maximum positive T": "positive", "Minimum negative T": "negative"}
DIAGNOSTICS = {"FULLT_ASYM": "asym", "FULLT": "fullt"}
CYCLES = {"Stop at repeat": "stop", "Allow repeats": "allow", "Skip visited modes": "skip"}
GLYPHS = {"Smooth (classic FULLT)": "smooth", "Diamonds (grid samples)": "diamonds", "Cells": "cells"}
PANEL_SCALING = {"Auto (1, 2, 5, 10...)": "auto", "None (original ratios)": "none", "Custom multipliers": "custom"}


class FullTTraceApp:
    def __init__(self, root, cases=(), trace=None, figure=None, output=None, formats=None):
        self.root = root
        configure_theme(root)
        version = (Path(__file__).resolve().parent / "VERSION").read_text().strip()
        root.title("CGYRO Transfer Explorer — v" + version)
        root.geometry("1380x900")
        root.minsize(920, 640)
        self.cases, self.results, self.vars, self.busy_widgets = [], [], {}, []
        self._setting_options = False
        self.events = queue.Queue()
        self.cancel = threading.Event()
        self.busy = False
        self.canvas = self.toolbar = None
        self.trace_settings = trace or TraceSettings()
        self.figure_settings = figure or FigureSettings()
        self.formats = tuple(formats or ("pdf", "svg", "png"))
        self._layout(output)
        self._set_options(self.trace_settings, self.figure_settings)
        for case in cases:
            self._add_path(case)
        self._sync_actions()
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.after(100, self._poll)

    def _layout(self, output):
        tr = lambda zh, en: ui_text(self.root, zh, en)
        menu = tk.Menu(self.root)
        file_menu = tk.Menu(menu, tearoff=False)
        file_menu.add_command(label=tr("打开会话…", "Load session…"), command=self.load_session, accelerator="Ctrl+O")
        file_menu.add_command(label=tr("保存会话…", "Save session…"), command=self.save_session, accelerator="Ctrl+S")
        file_menu.add_separator()
        file_menu.add_command(label=tr("退出", "Exit"), command=self.close)
        menu.add_cascade(label=tr("文件", "File"), menu=file_menu)
        help_menu = tk.Menu(menu, tearoff=False)
        help_menu.add_command(label=tr("追踪说明…", "Transfer guide…"), command=self.show_help)
        help_menu.add_command(label="Interface Guide...", command=lambda: show_interface_guide(self.root))
        menu.add_cascade(label=tr("帮助", "Help"), menu=help_menu)
        self.root.configure(menu=menu)
        self.root.bind("<Control-o>", lambda _e: self.load_session())
        self.root.bind("<Control-s>", lambda _e: self.save_session())
        self.status = tk.StringVar(value=tr("就绪 · 各面板共用色标", "Ready · one shared color scale"))
        ttk.Label(self.root, textvariable=self.status, style="Status.TLabel", padding=(12, 7)).pack(side="bottom", fill="x")
        body = ttk.Frame(self.root, padding=(10, 8, 10, 0))
        body.pack(fill="both", expand=True)
        self.main_pane = ttk.PanedWindow(body, orient="horizontal")
        self.main_pane.pack(fill="both", expand=True)
        self.left_panel = left = ttk.Frame(self.main_pane, width=410)
        self.main_pane.add(left, weight=0)
        left.columnconfigure(0, weight=1)
        left.rowconfigure(2, weight=1)
        ttk.Label(left, text="CGYRO Transfer Explorer", style="SidebarTitle.TLabel").grid(row=0, column=0, sticky="w", pady=(2, 10))
        case_frame = ttk.LabelFrame(left, text=tr("1  选择算例", "1  Select cases"), style="Card.TLabelframe")
        case_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10), padx=(0, 8))
        case_frame.columnconfigure(0, weight=1)
        self.case_list = tk.Listbox(case_frame, height=4, width=12, selectmode="extended",
                                    exportselection=False, relief="flat", activestyle="none")
        self.case_list.grid(row=0, column=0, columnspan=2, sticky="ew")
        ybar = ttk.Scrollbar(case_frame, command=self.case_list.yview)
        ybar.grid(row=0, column=2, sticky="ns")
        xbar = ttk.Scrollbar(case_frame, orient="horizontal", command=self.case_list.xview)
        xbar.grid(row=1, column=0, columnspan=2, sticky="ew")
        self.case_list.configure(yscrollcommand=ybar.set, xscrollcommand=xbar.set)
        self._button(case_frame, tr("添加算例…", "Add case…"), self.add_case).grid(row=2, column=0, sticky="ew", pady=(6, 0), padx=(0, 5))
        self._button(case_frame, tr("移除", "Remove"), self.remove_cases).grid(row=2, column=1, sticky="ew", pady=(6, 0))
        ttk.Label(case_frame, text=tr("未选择时追踪全部算例", "No selection uses all cases"), style="Hint.TLabel").grid(row=3, column=0, columnspan=2, sticky="w", pady=(4, 0))

        self.option_tabs = tabs = ttk.Notebook(left)
        tabs.grid(row=2, column=0, sticky="nsew", padx=(0, 8))
        self.option_canvases = []
        self._forms, self.field_widgets = [], {}
        trace_tab = self._scroll_tab(tabs, tr("2  追踪设置", "2  Trace settings"))
        figure_tab = self._scroll_tab(tabs, tr("3  出图样式", "3  Figure style"))

        def group(parent, title, collapsed=False):
            if collapsed:
                section = CollapsibleSection(parent, title)
                section.pack(fill="x", pady=(3, 7))
                frame = section.body
            else:
                frame = ttk.LabelFrame(parent, text=title, style="Card.TLabelframe")
                frame.pack(fill="x", pady=(0, 10))
            self._forms.append(FormLayout(frame))
            return frame

        mode = group(trace_tab, tr("起点与追踪", "Initial mode and trace"))
        for row, (key, label, choices) in enumerate([
            ("diagnostic", tr("诊断", "Diagnostic"), DIAGNOSTICS),
            ("kx", tr("kx · 起点", "Initial kx"), None),
            ("ky", tr("ky · 固定", "Fixed ky"), None),
            ("steps", tr("面板数", "Panels"), None),
            ("direction", tr("极值方向", "Direction"), DIRECTIONS)]):
            self._entry(mode, row, key, label, choices)
        time_group = group(trace_tab, tr("时间平均", "Time average"))
        self._entry(time_group, 0, "t_start", "Time start")
        self._entry(time_group, 1, "t_end", "Time end")
        ttk.Label(time_group, text=tr("留空：末段 50% / 最后时刻", "Blank: last 50% / last sample"), style="Hint.TLabel").grid(row=2, column=0, columnspan=2, sticky="w")
        advanced = group(trace_tab, tr("重复模式与计算说明", "Repeated modes and calculation notes"), collapsed=True)
        self._entry(advanced, 0, "cycles", tr("重复模式", "Repeated mode"), CYCLES)
        ttk.Label(advanced, text=tr("每一步固定 ky，只改变 kx。极值在 ky′ = ky 的行上选取，各步使用相同的 Re(T) 时间平均区间。", "Fixed ky at every step; only kx changes. Extrema use the ky′ = ky row and the same mean Re(T) time window."),
                  style="Hint.TLabel", wraplength=265).grid(row=1, column=0, columnspan=2, sticky="w", pady=5)

        size = group(figure_tab, tr("尺寸与文字", "Size and typography"))
        for row, (key, label) in enumerate([("width_mm", tr("图宽 · mm", "Width · mm")),
                    ("panel_height_mm", tr("面板高 · mm", "Panel height · mm")), ("font_pt", tr("字号 · pt", "Font · pt"))]):
            self._entry(size, row, key, label)
        colors = group(figure_tab, tr("渲染与共用色标", "Rendering and shared color scale"))
        for row, (key, label, choices) in enumerate([
            ("glyph", tr("渲染", "Rendering"), GLYPHS),
            ("color_max", tr("色标上限", "Color max"), None),
            ("panel_scaling", tr("面板倍率", "Panel scaling"), PANEL_SCALING),
            ("panel_multipliers", tr("自定义倍率", "Multipliers"), None)]):
            self._entry(colors, row, key, label, choices)
        self.vars["panel_multipliers"].trace_add("write", self._multipliers_changed)
        self.vars["normalize"] = tk.BooleanVar(value=False)
        ttk.Checkbutton(colors, text=tr("除以全局 max |T|", "Divide by global max |T|"), variable=self.vars["normalize"]).grid(row=4, column=0, columnspan=2, sticky="w")
        ttk.Label(colors, text=tr("色标留空时由原始数据确定；倍率不改变色标范围。", "Blank color max uses raw maps; multipliers keep that range."),
                  style="Hint.TLabel", wraplength=265).grid(row=5, column=0, columnspan=2, sticky="w", pady=4)
        limits = group(figure_tab, tr("坐标范围", "Axis limits"))
        self._entry(limits, 0, "xlim", "kx′ min, max")
        self._entry(limits, 1, "ylim", "ky′ min, max")
        details = group(figure_tab, tr("标题与导出精度", "Titles and export resolution"), collapsed=True)
        self._entry(details, 0, "title", tr("总标题", "Context title"))
        self._entry(details, 1, "colorbar_label", tr("色标标题", "Colorbar label"))
        self._entry(details, 2, "dpi", "DPI")
        self.vars["hide_title"] = tk.BooleanVar(value=False)
        ttk.Checkbutton(details, text=tr("隐藏算例与时间标题", "Hide case/time title"), variable=self.vars["hide_title"]).grid(row=3, column=0, columnspan=2, sticky="w")
        for form in self._forms:
            form.adopt_grid()
        for frame, canvas in self.option_canvases:
            bind_local_scroll(frame, canvas)

        actions = ttk.Frame(left, padding=(0, 8, 8, 0))
        actions.grid(row=3, column=0, sticky="ew")
        actions.columnconfigure(0, weight=1)
        self.run_button = self._button(actions, tr("4  开始追踪", "4  Trace cases"), self.run_trace)
        self.run_button.configure(style="Accent.TButton")
        self.run_button.grid(row=0, column=0, sticky="ew")
        self.cancel_button = ttk.Button(actions, text=tr("取消", "Cancel"), command=self.cancel.set, state="disabled")
        self.cancel_button.grid(row=0, column=1, padx=(5, 0))
        self.apply_button = self._button(actions, tr("应用出图样式", "Apply figure style"), self.redraw)
        self.apply_button.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(5, 0))

        self.right_panel = right = ttk.Frame(self.main_pane, padding=(10, 0, 0, 0))
        self.main_pane.add(right, weight=1)
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)
        results_header = ttk.Frame(right, padding=(0, 0, 0, 8))
        results_header.grid(row=0, column=0, sticky="ew")
        results_header.columnconfigure(1, weight=1)
        ttk.Label(results_header, text=tr("5  结果", "5  Results"), style="SectionLabel.TLabel").grid(row=0, column=0, padx=(0, 12))
        self.result_var = tk.StringVar()
        self.result_combo = LocalizedCombobox(results_header, textvariable=self.result_var, width=12, state="readonly")
        self.result_combo.grid(row=0, column=1, sticky="ew")
        self.result_combo.bind("<<ComboboxSelected>>", lambda _e: self.redraw())
        views = ttk.Notebook(right)
        views.grid(row=1, column=0, sticky="nsew")
        preview, table_frame = ttk.Frame(views), ttk.Frame(views)
        views.add(preview, text=tr("多面板预览", "Stacked figure"))
        views.add(table_frame, text=tr("追踪坐标与数值", "Trace coordinates and values"))
        preview.columnconfigure(0, weight=1)
        preview.rowconfigure(1, weight=1)
        self.toolbar_frame = ttk.Frame(preview)
        self.toolbar_frame.grid(row=0, column=0, columnspan=2, sticky="ew")
        self.scroll_canvas = tk.Canvas(preview, background="#e8ecf1", highlightthickness=0)
        self.scroll_canvas.grid(row=1, column=0, sticky="nsew")
        bar = ttk.Scrollbar(preview, command=self.scroll_canvas.yview)
        bar.grid(row=1, column=1, sticky="ns")
        xbar = ttk.Scrollbar(preview, orient="horizontal", command=self.scroll_canvas.xview)
        xbar.grid(row=2, column=0, sticky="ew")
        self.scroll_canvas.configure(yscrollcommand=bar.set, xscrollcommand=xbar.set)
        self.figure_frame = ttk.Frame(self.scroll_canvas)
        self.scroll_canvas.create_window(0, 0, window=self.figure_frame, anchor="nw")
        self.figure_frame.bind("<Configure>", lambda _e: self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox("all")))
        ttk.Label(self.figure_frame, text=tr("添加算例，设置起点，然后开始追踪。\n预览保留出图比例，可滚动查看全部面板。", "Add cases, set the initial mode, then Trace.\nThe preview keeps the publication aspect ratio."),
                  padding=24, style="Muted.TLabel", wraplength=340).pack()
        columns = ("step", "kx", "ky", "kx_prime", "ky_prime", "T", "next_kx", "next_ky", "mirror")
        captions = ("Step", "kx", "ky", "Peak kx′", "Peak ky′", "Re(T)", "Next kx", "Next ky", "Conjugate")
        self.table = ttk.Treeview(table_frame, columns=columns, show="headings")
        for key, title in zip(columns, captions):
            self.table.heading(key, text=title)
            self.table.column(key, width=100, minwidth=75, anchor="center", stretch=False)
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)
        self.table.grid(row=0, column=0, sticky="nsew")
        table_xbar = ttk.Scrollbar(table_frame, orient="horizontal", command=self.table.xview)
        table_xbar.grid(row=1, column=0, sticky="ew")
        table_ybar = ttk.Scrollbar(table_frame, command=self.table.yview)
        table_ybar.grid(row=0, column=1, sticky="ns")
        self.table.configure(xscrollcommand=table_xbar.set, yscrollcommand=table_ybar.set)

        export = ttk.Frame(right, padding=(0, 8, 0, 0))
        export.grid(row=2, column=0, sticky="ew")
        export.columnconfigure(0, weight=1)
        self.output_var = tk.StringVar(value=str(output or (Path.cwd() / "FULLT_exports")))
        self.tiff_var = tk.BooleanVar(value="tiff" in self.formats)
        export_settings = CollapsibleSection(export, tr("导出目录与格式", "Export folder and formats"))
        export_settings.grid(row=0, column=0, sticky="ew", padx=(0, 12))
        export_settings.body.columnconfigure(0, weight=1)
        ttk.Entry(export_settings.body, textvariable=self.output_var, width=12).grid(row=0, column=0, sticky="ew")
        ttk.Button(export_settings.body, text="…", width=3, command=self.choose_output).grid(row=0, column=1, padx=(5, 0))
        ttk.Checkbutton(export_settings.body, text="PDF + SVG + PNG; + TIFF", variable=self.tiff_var).grid(row=1, column=0, columnspan=2, sticky="w")
        self.export_button = self._button(export, tr("导出图片与数据…", "Export figures + data…"), self.export)
        self.export_button.grid(row=0, column=1, sticky="ne")
        self.main_pane.bind("<Configure>", self._fit_panes, add="+")
        self.main_pane.bind("<ButtonRelease-1>", self._fit_panes, add="+")
        self.root.after_idle(lambda: self.main_pane.sashpos(0, 410))
        get_language(self.root).install_menu(menu, self.status, self._language_changed)

    def _language_changed(self):
        """Refresh result labels in place without rebuilding the figure or reading data."""
        if self.results:
            self._refresh_result_choices()
            self._fill_result_table(self.results[max(0, self.result_combo.current())])

    def _refresh_result_choices(self):
        index = max(0, self.result_combo.current())
        template = ui_text(self.root, "{i}. {case} — {n} 个面板", "{i}. {case} — {n} panels")
        self.result_combo.configure(values=[template.format(i=i + 1, case=r.case, n=len(r.panels))
                                           for i, r in enumerate(self.results)])
        if self.results:
            self.result_combo.current(min(index, len(self.results) - 1))

    def _fill_result_table(self, result):
        self.table.delete(*self.table.get_children())
        for i, panel in enumerate(result.panels, 1):
            peak = panel.peak
            values = [i, panel.mode.kx, panel.mode.ky] + (
                [peak.kx, peak.ky, peak.value, peak.next_mode.kx, peak.next_mode.ky,
                 ui_text(self.root, "是", "Yes") if peak.conjugated else ui_text(self.root, "否", "No")]
                if peak else [""] * 6)
            self.table.insert("", "end", values=[f"{v:.6g}" if isinstance(v, float) else v for v in values])

    def _result_status(self, result, scale):
        factors = ", ".join(f"{x:g}" for x in scale["panel_multipliers"])
        summary = ui_text(self.root, "{n} 个面板 | 倍率：{factors} | 色标范围：±{limit:.5g}",
                          "{n} panels | multipliers: {factors} | color range: ±{limit:.5g}")
        reason = result.stop_reason
        if result.stop_code == "completed":
            reason = ui_text(self.root, "已完成请求的全部面板。", "Completed all requested panels.")
        elif result.stop_code == "cycle":
            reason = ui_text(self.root, "遇到重复模式而停止；可选择允许重复或跳过已访问模式。",
                             "Stopped at a repeated mode; choose Allow repeats or Skip visited modes to continue.")
        elif result.stop_code == "no_extremum":
            reason = ui_text(self.root, "固定 ky 行上没有符合条件的有限极值。",
                             "No eligible finite extremum remains on the fixed ky row.")
        return summary.format(n=len(result.panels), factors=factors, limit=scale["limit"]) + " — " + reason

    def _fit_panes(self, _event=None):
        width = self.main_pane.winfo_width()
        if width > 700:
            current = self.main_pane.sashpos(0)
            wanted = min(max(320, current), width - 380)
            if current != wanted:
                self.main_pane.sashpos(0, wanted)

    def _button(self, parent, text, callback):
        button = ttk.Button(parent, text=text, command=callback)
        self.busy_widgets.append(button)
        return button

    def _sync_actions(self):
        self.run_button.configure(state="normal" if self.cases and not self.busy else "disabled")
        state = "normal" if self.results and not self.busy else "disabled"
        self.apply_button.configure(state=state)
        self.export_button.configure(state=state)

    def _scroll_tab(self, tabs, title):
        shell = ttk.Frame(tabs)
        tabs.add(shell, text=title)
        shell.columnconfigure(0, weight=1)
        shell.rowconfigure(0, weight=1)
        canvas = tk.Canvas(shell, width=330, height=100, highlightthickness=0, background=configure_theme(self.root)["background"])
        canvas.grid(row=0, column=0, sticky="nsew")
        bar = ttk.Scrollbar(shell, command=canvas.yview)
        bar.grid(row=0, column=1, sticky="ns")
        canvas.configure(yscrollcommand=bar.set)
        frame = ttk.Frame(canvas, padding=10)
        item = canvas.create_window(0, 0, window=frame, anchor="nw")
        def resize_content(_event):
            canvas.configure(scrollregion=canvas.bbox("all"))
        frame.bind("<Configure>", resize_content)
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(item, width=e.width))
        self.option_canvases.append((frame, canvas))
        return frame

    def _entry(self, parent, row, key, label, choices=None):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=5, padx=(0, 8))
        self.vars[key] = tk.StringVar()
        widget = (LocalizedCombobox(parent, textvariable=self.vars[key], values=list(choices), state="readonly", width=20)
                  if choices else ttk.Entry(parent, textvariable=self.vars[key], width=21))
        widget.grid(row=row, column=1, sticky="ew", pady=5)
        self.field_widgets[key] = widget

    def _set_options(self, trace, figure):
        values = {**asdict(trace), **asdict(figure)}
        mappings = {"direction": DIRECTIONS, "diagnostic": DIAGNOSTICS, "cycles": CYCLES,
                    "glyph": GLYPHS, "panel_scaling": PANEL_SCALING}
        self._setting_options = True
        try:
            for key, value in values.items():
                if key in mappings:
                    value = next(label for label, code in mappings[key].items() if code == value)
                if key in ("xlim", "ylim", "panel_multipliers") and value is not None:
                    value = ", ".join(map(str, value))
                self.vars[key].set("" if value is None else value)
            self.vars["hide_title"].set(figure.title == "")
        finally:
            self._setting_options = False

    def _multipliers_changed(self, *_args):
        """Typing or pasting manual factors must not leave them ignored in Auto."""
        if not self._setting_options and self.vars["panel_multipliers"].get().strip():
            self.vars["panel_scaling"].set("Custom multipliers")
            self.status.set(ui_text(self.root, "已选择自定义倍率，点击“应用出图样式”更新预览", "Manual multipliers selected — click Apply figure style to update."))

    def _options(self):
        def number(key, optional=False):
            value = self.vars[key].get().strip()
            return None if optional and not value else float(value)
        def limits(key):
            text = self.vars[key].get().replace(",", " ").strip()
            return tuple(map(float, text.split())) if text else None
        trace = TraceSettings(kx=number("kx"), ky=number("ky"), steps=int(self.vars["steps"].get()),
                              direction=DIRECTIONS[self.vars["direction"].get()], diagnostic=DIAGNOSTICS[self.vars["diagnostic"].get()],
                              t_start=number("t_start", True), t_end=number("t_end", True), cycles=CYCLES[self.vars["cycles"].get()])
        figure = FigureSettings(width_mm=number("width_mm"), panel_height_mm=number("panel_height_mm"), font_pt=number("font_pt"),
                                color_max=number("color_max", True), normalize=self.vars["normalize"].get(), glyph=GLYPHS[self.vars["glyph"].get()],
                                xlim=limits("xlim"), ylim=limits("ylim"), title="" if self.vars["hide_title"].get() else self.vars["title"].get().strip() or None,
                                colorbar_label=self.vars["colorbar_label"].get(), dpi=int(self.vars["dpi"].get()),
                                panel_scaling=PANEL_SCALING[self.vars["panel_scaling"].get()],
                                panel_multipliers=limits("panel_multipliers"))
        trace.validate()
        figure.validate()
        return trace, figure

    def _add_path(self, path):
        path = str(Path(path).expanduser().absolute())
        if path not in self.cases:
            self.cases.append(path)
            self.case_list.insert("end", Path(path).name)
        self.case_list.selection_set(self.cases.index(path))
        self._sync_actions()

    def add_case(self):
        path = filedialog.askdirectory(parent=self.root, title="Select a case containing FULLT files")
        if path:
            self._add_path(path)

    def remove_cases(self):
        for index in reversed(self.case_list.curselection()):
            self.case_list.delete(index)
            del self.cases[index]
        self._sync_actions()

    def choose_output(self):
        path = filedialog.askdirectory(parent=self.root, title="Choose export parent folder")
        if path:
            self.output_var.set(path)

    def _start_worker(self, task, cancellable=True):
        if self.busy:
            return
        self.busy = True
        self.cancel.clear()
        for button in self.busy_widgets:
            button.configure(state="disabled")
        self.result_combo.configure(state="disabled")
        self.cancel_button.configure(state="normal" if cancellable else "disabled")
        def work():
            try:
                self.events.put(("done", task()))
            except Exception as exc:
                self.events.put(("error", exc))
        threading.Thread(target=work, daemon=True).start()

    def run_trace(self):
        try:
            trace, figure = self._options()
            indices = self.case_list.curselection() or range(len(self.cases))
            cases = [self.cases[i] for i in indices]
            if not cases:
                raise ValueError("Add at least one case first.")
        except (ValueError, KeyError) as exc:
            messagebox.showerror("Trace options", str(exc), parent=self.root)
            return
        self.trace_settings, self.figure_settings = trace, figure
        self.status.set(ui_text(self.root, "正在读取选中算例的 FULLT 数据…", "Reading selected FULLT diagnostics…"))
        def task():
            results, errors = [], []
            for case in cases:
                try:
                    results.append(trace_case(case, trace, progress=lambda text: self.events.put(("progress", text)), cancel=self.cancel))
                except TraceCancelled:
                    raise
                except (ValueError, OSError, RuntimeError) as exc:
                    errors.append(f"{Path(case).name}: {exc}")
            return "trace", results, errors
        self._start_worker(task)

    def redraw(self):
        if self.busy or not self.results:
            return
        try:
            _, figure = self._options()
            index = max(0, self.result_combo.current())
            result = self.results[index]
            rendered = build_figure(result, figure)
            for child in self.figure_frame.winfo_children():
                child.destroy()
            for child in self.toolbar_frame.winfo_children():
                child.destroy()
            self.canvas = FigureCanvasTkAgg(rendered, master=self.figure_frame)
            self.canvas.get_tk_widget().pack()
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbar_frame, pack_toolbar=False)
            self.toolbar.pack(fill="x")
            self.canvas.draw()
            self.figure_settings = figure
            self._fill_result_table(result)
            scale = rendered._fullt_trace_scale
            self.status.set(self._result_status(result, scale))
        except (ValueError, KeyError, OSError) as exc:
            messagebox.showerror("Figure options", str(exc), parent=self.root)

    def export(self):
        if not self.results:
            messagebox.showinfo("Export", "Trace at least one case first.", parent=self.root)
            return
        try:
            _, figure = self._options()
        except (ValueError, KeyError) as exc:
            messagebox.showerror("Figure options", str(exc), parent=self.root)
            return
        parent = Path(self.output_var.get()).expanduser()
        formats = ("pdf", "svg", "png", "tiff") if self.tiff_var.get() else ("pdf", "svg", "png")
        results = list(self.results)
        self.status.set(ui_text(self.root, "正在导出图片与追踪数据…", "Exporting publication figures and trace data…"))
        def task():
            paths = [export_bundle(r, figure, parent / suggested_bundle_name(r), formats) for r in results]
            return "export", paths, []
        self._start_worker(task, cancellable=False)

    def _poll(self):
        try:
            while True:
                kind, payload = self.events.get_nowait()
                if kind == "progress":
                    self.status.set(payload)
                    continue
                self.busy = False
                for button in self.busy_widgets:
                    button.configure(state="normal")
                self.result_combo.configure(state="readonly")
                self.cancel_button.configure(state="disabled")
                if kind == "error":
                    self.status.set(str(payload))
                    self._sync_actions()
                    if not isinstance(payload, TraceCancelled):
                        messagebox.showerror("CGYRO Transfer Explorer", str(payload), parent=self.root)
                    continue
                action, items, errors = payload
                if action == "trace" and items:
                    self.results = items
                    self._refresh_result_choices()
                    self.result_combo.current(0)
                    self.redraw()
                elif action == "export":
                    self.status.set(ui_text(self.root, "已导出 {n} 组图片与数据到 {path}",
                                            "Exported {n} bundle(s) to {path}").format(n=len(items), path=items[0].parent))
                if errors:
                    self.status.set(ui_text(self.root, "部分算例追踪失败；如果全部失败，将保留原有结果。",
                                            "Some cases could not be traced; previous results are retained if none succeeded."))
                    messagebox.showwarning("Case diagnostics", "\n".join(errors), parent=self.root)
                self._sync_actions()
        except queue.Empty:
            pass
        self.root.after(100, self._poll)

    def save_session(self):
        try:
            trace, figure = self._options()
            path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".json", filetypes=[("JSON", "*.json")])
            if path:
                write_json_atomic(path, {"format": "fullt-trace-1", "cases": self.cases, "trace": asdict(trace), "figure": asdict(figure), "output": self.output_var.get()})
                self.status.set(ui_text(self.root, "已保存会话：", "Session saved: ") + path)
        except (ValueError, KeyError, OSError) as exc:
            messagebox.showerror("Save session", str(exc), parent=self.root)

    def load_session(self):
        if self.busy:
            return
        path = filedialog.askopenfilename(parent=self.root, filetypes=[("JSON", "*.json")])
        if not path:
            return
        try:
            config = json.loads(Path(path).read_text(encoding="utf-8-sig"))
            if not isinstance(config, dict) or config.get("format") != "fullt-trace-1":
                raise ValueError("This is not a FULLT Trace session.")
            trace, figure = TraceSettings(**config["trace"]), FigureSettings(**config["figure"])
            trace.validate()
            figure.validate()
            cases = config["cases"]
            if not isinstance(cases, list) or not all(isinstance(p, str) for p in cases):
                raise ValueError("Invalid session case list.")
            self._set_options(trace, figure)
            self.case_list.delete(0, "end")
            self.cases.clear()
            for case in cases:
                self._add_path(case)
            self.output_var.set(config.get("output", self.output_var.get()))
            self.status.set(ui_text(self.root, "会话已加载；重新追踪以更新结果", "Session loaded. Trace again to refresh the figure."))
        except (ValueError, TypeError, KeyError, OSError) as exc:
            messagebox.showerror("Load session", str(exc), parent=self.root)

    def show_help(self):
        window = tk.Toplevel(self.root)
        window.title("CGYRO Transfer Explorer — Guide")
        window.geometry("850x650")
        text = tk.Text(window, wrap="word", padx=16, pady=12)
        scroll = ttk.Scrollbar(window, command=text.yview)
        scroll.pack(side="right", fill="y")
        text.configure(yscrollcommand=scroll.set)
        text.pack(fill="both", expand=True)
        path = Path(__file__).resolve().parent / "FULLT_TRACE_GUIDE.md"
        text.insert("1.0", path.read_text(encoding="utf-8") if path.is_file() else "See FULLT_TRACE_GUIDE.md.")
        text.configure(state="disabled")

    def close(self):
        self.cancel.set()
        self.root.destroy()


def launch(cases=(), trace=None, figure=None, output=None, formats=None):
    root = tk.Tk()
    FullTTraceApp(root, cases, trace, figure, output, formats)
    root.mainloop()
