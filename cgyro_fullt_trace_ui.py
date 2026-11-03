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
except ImportError:
    from .cgyro_fullt_trace_core import TraceCancelled, TraceSettings, trace_case
    from .cgyro_fullt_trace_plot import FigureSettings, build_figure, export_bundle, suggested_bundle_name
    from .cgyro_io import write_json_atomic


DIRECTIONS = {"Maximum positive T": "positive", "Minimum negative T": "negative"}
DIAGNOSTICS = {"FULLT_ASYM": "asym", "FULLT": "fullt"}
CYCLES = {"Stop at repeat": "stop", "Allow repeats": "allow", "Skip visited modes": "skip"}
GLYPHS = {"Smooth (classic FULLT)": "smooth", "Diamonds (grid samples)": "diamonds", "Cells": "cells"}
PANEL_SCALING = {"Auto (1, 2, 5, 10...)": "auto", "None (original ratios)": "none", "Custom multipliers": "custom"}


class FullTTraceApp:
    def __init__(self, root, cases=(), trace=None, figure=None, output=None, formats=None):
        self.root = root
        version = (Path(__file__).resolve().parent / "VERSION").read_text().strip()
        root.title("FULLT Trace Studio — v" + version)
        root.geometry("1380x900")
        root.minsize(1040, 720)
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
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.after(100, self._poll)

    def _layout(self, output):
        menu = tk.Menu(self.root)
        file_menu = tk.Menu(menu, tearoff=False)
        file_menu.add_command(label="Load session…", command=self.load_session)
        file_menu.add_command(label="Save session…", command=self.save_session)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.close)
        menu.add_cascade(label="File", menu=file_menu)
        menu.add_command(label="Help", command=self.show_help)
        self.root.configure(menu=menu)
        body = ttk.Frame(self.root, padding=10)
        body.pack(fill="both", expand=True)
        body.columnconfigure(1, weight=1)
        body.rowconfigure(0, weight=1)
        left = ttk.Frame(body, width=350)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 12))
        left.columnconfigure(0, weight=1)
        left.rowconfigure(2, weight=1)
        ttk.Label(left, text="FULLT Trace Studio", font=("TkDefaultFont", 15, "bold")).grid(sticky="w", pady=(0, 9))
        case_frame = ttk.LabelFrame(left, text="Cases", padding=8)
        case_frame.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        case_frame.columnconfigure(0, weight=1)
        self.case_list = tk.Listbox(case_frame, height=4, selectmode="extended", exportselection=False)
        self.case_list.grid(row=0, column=0, columnspan=2, sticky="ew")
        ybar = ttk.Scrollbar(case_frame, command=self.case_list.yview)
        ybar.grid(row=0, column=2, sticky="ns")
        self.case_list.configure(yscrollcommand=ybar.set)
        self._button(case_frame, "Add case…", self.add_case).grid(row=1, column=0, sticky="ew", pady=(6, 0))
        self._button(case_frame, "Remove", self.remove_cases).grid(row=1, column=1, sticky="ew", padx=(5, 0), pady=(6, 0))
        tabs = ttk.Notebook(left)
        self.option_tabs = tabs
        self.option_canvases = []
        tabs.grid(row=2, column=0, sticky="nsew")
        trace_tab = self._scroll_tab(tabs, "Trace")
        figure_tab = self._scroll_tab(tabs, "Publication")
        for tab in (trace_tab, figure_tab):
            tab.columnconfigure(1, weight=1)
        trace_fields = [("diagnostic", "Diagnostic", DIAGNOSTICS), ("kx", "Initial kx", None),
                        ("ky", "Fixed ky", None), ("steps", "Panels / steps", None),
                        ("direction", "Follow", DIRECTIONS), ("cycles", "On repeated mode", CYCLES),
                        ("t_start", "Time start", None), ("t_end", "Time end", None)]
        for row, (key, label, choices) in enumerate(trace_fields):
            self._entry(trace_tab, row, key, label, choices)
        ttk.Label(trace_tab, text="Blank time start: last 50% of the chosen span.\nBlank time end: last saved sample.\n\nFixed ky at every step; only kx changes.\nExtrema use the row ky' = fixed ky\nand the same mean Re(T) time window.",
                  wraplength=315, foreground="#586575").grid(row=8, column=0, columnspan=2, sticky="w", pady=12)
        figure_fields = [("width_mm", "Width (mm)", None), ("panel_height_mm", "Panel height (mm)", None),
                         ("font_pt", "Font (pt)", None), ("color_max", "Shared color max", None),
                         ("glyph", "Map rendering", GLYPHS),
                         ("panel_scaling", "Panel scaling", PANEL_SCALING),
                         ("panel_multipliers", "Multipliers: 1, 1, 10", None),
                         ("xlim", "x limits: min, max", None),
                         ("ylim", "y limits: min, max", None), ("title", "Context title", None),
                         ("colorbar_label", "Colorbar label", None), ("dpi", "Export DPI", None)]
        for row, (key, label, choices) in enumerate(figure_fields):
            self._entry(figure_tab, row, key, label, choices)
        self.vars["panel_multipliers"].trace_add("write", self._multipliers_changed)
        self.vars["normalize"] = tk.BooleanVar(value=False)
        ttk.Checkbutton(figure_tab, text="Divide every panel by ONE global max |T|", variable=self.vars["normalize"]).grid(
            row=len(figure_fields), column=0, columnspan=2, sticky="w", pady=6)
        self.vars["hide_title"] = tk.BooleanVar(value=False)
        ttk.Checkbutton(figure_tab, text="Hide case/time title", variable=self.vars["hide_title"]).grid(row=len(figure_fields) + 1, column=0, columnspan=2, sticky="w")
        ttk.Label(figure_tab, text="Entering multipliers selects Custom automatically. Click Apply to update.\nBlank color max: fixed range from the unmultiplied maps.\nChanging multipliers keeps this color range; values beyond it saturate.\nBlank limits: full grid; tracing keeps the same ky.",
                  wraplength=315, foreground="#586575").grid(row=len(figure_fields) + 2, column=0, columnspan=2, sticky="w", pady=7)
        for frame, canvas in self.option_canvases:
            self._bind_option_scroll(frame, canvas)
        actions = ttk.Frame(left)
        actions.grid(row=3, column=0, sticky="ew", pady=(8, 0))
        actions.columnconfigure(0, weight=1)
        self._button(actions, "Trace selected cases", self.run_trace).grid(row=0, column=0, sticky="ew")
        self.cancel_button = ttk.Button(actions, text="Cancel", command=self.cancel.set, state="disabled")
        self.cancel_button.grid(row=0, column=1, padx=(5, 0))
        self._button(actions, "Apply publication style to preview", self.redraw).grid(row=1, column=0, columnspan=2, sticky="ew", pady=5)
        export = ttk.LabelFrame(left, text="Export all traced results", padding=8)
        export.grid(row=4, column=0, sticky="ew", pady=(3, 0))
        export.columnconfigure(0, weight=1)
        self.output_var = tk.StringVar(value=str(output or (Path.cwd() / "FULLT_exports")))
        ttk.Entry(export, textvariable=self.output_var, width=27).grid(row=0, column=0, sticky="ew")
        ttk.Button(export, text="…", width=3, command=self.choose_output).grid(row=0, column=1, padx=(4, 0))
        self.tiff_var = tk.BooleanVar(value="tiff" in self.formats)
        ttk.Checkbutton(export, text="PDF + SVG + PNG; also TIFF", variable=self.tiff_var).grid(row=1, column=0, columnspan=2, sticky="w", pady=4)
        self._button(export, "Export figure + trace data…", self.export).grid(row=2, column=0, columnspan=2, sticky="ew")
        right = ttk.Frame(body)
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)
        self.result_var = tk.StringVar()
        self.result_combo = ttk.Combobox(right, textvariable=self.result_var, state="readonly")
        self.result_combo.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        self.result_combo.bind("<<ComboboxSelected>>", lambda _e: self.redraw())
        views = ttk.Notebook(right)
        views.grid(row=1, column=0, sticky="nsew")
        preview = ttk.Frame(views)
        table_frame = ttk.Frame(views)
        views.add(preview, text="Stacked figure")
        views.add(table_frame, text="Trace coordinates")
        preview.columnconfigure(0, weight=1)
        preview.rowconfigure(1, weight=1)
        self.toolbar_frame = ttk.Frame(preview)
        self.toolbar_frame.grid(row=0, column=0, columnspan=2, sticky="ew")
        self.scroll_canvas = tk.Canvas(preview, background="#e8ecf1", highlightthickness=0)
        self.scroll_canvas.grid(row=1, column=0, sticky="nsew")
        bar = ttk.Scrollbar(preview, orient="vertical", command=self.scroll_canvas.yview)
        bar.grid(row=1, column=1, sticky="ns")
        xbar = ttk.Scrollbar(preview, orient="horizontal", command=self.scroll_canvas.xview)
        xbar.grid(row=2, column=0, sticky="ew")
        self.scroll_canvas.configure(yscrollcommand=bar.set, xscrollcommand=xbar.set)
        self.figure_frame = ttk.Frame(self.scroll_canvas)
        self.scroll_canvas.create_window(0, 0, window=self.figure_frame, anchor="nw")
        self.figure_frame.bind("<Configure>", lambda _e: self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox("all")))
        ttk.Label(self.figure_frame, text="Load a case, choose the initial mode, then Trace.\nThe preview keeps the publication aspect ratio; scroll to see further panels.", padding=25).pack()
        columns = ("step", "kx", "ky", "kx_prime", "ky_prime", "T", "next_kx", "next_ky", "mirror")
        self.table = ttk.Treeview(table_frame, columns=columns, show="headings")
        for key in columns:
            self.table.heading(key, text=key)
            self.table.column(key, width=85, anchor="center")
        self.table.pack(fill="both", expand=True)
        table_bar = ttk.Scrollbar(table_frame, orient="horizontal", command=self.table.xview)
        table_bar.pack(fill="x")
        self.table.configure(xscrollcommand=table_bar.set)
        self.status = tk.StringVar(value="Ready — all panels use one shared color scale.")
        ttk.Label(self.root, textvariable=self.status, anchor="w", padding=(12, 6)).pack(fill="x")

    def _button(self, parent, text, callback):
        button = ttk.Button(parent, text=text, command=callback)
        self.busy_widgets.append(button)
        return button

    def _scroll_tab(self, tabs, title):
        shell = ttk.Frame(tabs)
        tabs.add(shell, text=title)
        shell.columnconfigure(0, weight=1)
        shell.rowconfigure(0, weight=1)
        canvas = tk.Canvas(shell, width=360, height=200, highlightthickness=0)
        canvas.grid(row=0, column=0, sticky="nsew")
        bar = ttk.Scrollbar(shell, command=canvas.yview)
        bar.grid(row=0, column=1, sticky="ns")
        canvas.configure(yscrollcommand=bar.set)
        frame = ttk.Frame(canvas, padding=10)
        item = canvas.create_window(0, 0, window=frame, anchor="nw")
        def resize_content(_event):
            canvas.configure(scrollregion=canvas.bbox("all"))
            # Respect DPI/font changes: grow the sidebar instead of clipping
            # checkbox labels or combobox selections inside a narrow viewport.
            width = max(360, frame.winfo_reqwidth())
            if int(float(canvas.cget("width"))) != width:
                canvas.configure(width=width)
        frame.bind("<Configure>", resize_content)
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(item, width=e.width))
        self.option_canvases.append((frame, canvas))
        return frame

    def _bind_option_scroll(self, widget, canvas):
        """Local bindings keep case-list wheel events out of the settings scroll."""
        def scroll(event):
            if canvas.yview() != (0., 1.):
                upward = getattr(event, "num", None) == 4 or getattr(event, "delta", 0) > 0
                canvas.yview_scroll(-2 if upward else 2, "units")
            return "break"
        for event in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            widget.bind(event, scroll)
        for child in widget.winfo_children():
            self._bind_option_scroll(child, canvas)

    def _entry(self, parent, row, key, label, choices=None):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=5, padx=(0, 8))
        self.vars[key] = tk.StringVar()
        widget = (ttk.Combobox(parent, textvariable=self.vars[key], values=list(choices), state="readonly", width=20)
                  if choices else ttk.Entry(parent, textvariable=self.vars[key], width=21))
        widget.grid(row=row, column=1, sticky="ew", pady=5)

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
            self.status.set("Manual multipliers selected — click Apply publication style to preview.")

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

    def add_case(self):
        path = filedialog.askdirectory(parent=self.root, title="Select a case containing FULLT files")
        if path:
            self._add_path(path)

    def remove_cases(self):
        for index in reversed(self.case_list.curselection()):
            self.case_list.delete(index)
            del self.cases[index]

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
        self.status.set("Reading selected FULLT diagnostics…")
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
            self.table.delete(*self.table.get_children())
            for i, panel in enumerate(result.panels, 1):
                peak = panel.peak
                values = [i, panel.mode.kx, panel.mode.ky] + ([peak.kx, peak.ky, peak.value, peak.next_mode.kx, peak.next_mode.ky, peak.conjugated]
                                                           if peak else [""] * 6)
                self.table.insert("", "end", values=[f"{v:.6g}" if isinstance(v, float) else v for v in values])
            scale = rendered._fullt_trace_scale
            factors = ", ".join(f"{x:g}" for x in scale["panel_multipliers"])
            self.status.set(f"{len(result.panels)} panels | multipliers: {factors} | color range: +/-{scale['limit']:.5g} — {result.stop_reason}")
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
        self.status.set("Exporting publication figures and trace data…")
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
                    if not isinstance(payload, TraceCancelled):
                        messagebox.showerror("FULLT Trace", str(payload), parent=self.root)
                    continue
                action, items, errors = payload
                if action == "trace" and items:
                    self.results = items
                    self.result_combo.configure(values=[f"{i + 1}. {r.case} — {len(r.panels)} panels" for i, r in enumerate(items)])
                    self.result_combo.current(0)
                    self.redraw()
                elif action == "export":
                    self.status.set(f"Exported {len(items)} bundle(s) to {items[0].parent}")
                    messagebox.showinfo("Export complete", "\n".join(map(str, items)), parent=self.root)
                if errors:
                    self.status.set("Some cases could not be traced; previous results are retained if none succeeded.")
                    messagebox.showwarning("Case diagnostics", "\n".join(errors), parent=self.root)
        except queue.Empty:
            pass
        self.root.after(100, self._poll)

    def save_session(self):
        try:
            trace, figure = self._options()
            path = filedialog.asksaveasfilename(parent=self.root, defaultextension=".json", filetypes=[("JSON", "*.json")])
            if path:
                write_json_atomic(path, {"format": "fullt-trace-1", "cases": self.cases, "trace": asdict(trace), "figure": asdict(figure), "output": self.output_var.get()})
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
            self.status.set("Session loaded. Trace again to refresh the figure.")
        except (ValueError, TypeError, KeyError, OSError) as exc:
            messagebox.showerror("Load session", str(exc), parent=self.root)

    def show_help(self):
        window = tk.Toplevel(self.root)
        window.title("FULLT Trace Studio — Guide")
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
