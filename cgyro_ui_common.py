"""Shared, offline-only Tk presentation helpers. No simulation data lives here."""

import json
import os
from pathlib import Path
import tkinter as tk
from tkinter import font as tkfont, messagebox, simpledialog, ttk

try:
    from cgyro_io import write_json_atomic
    from cgyro_ui_language import get_language
    from cgyro_fonts import GUI_FONT_SIZE, gui_font_family, configure_plot_fonts
except ImportError:
    from .cgyro_io import write_json_atomic
    from .cgyro_ui_language import get_language
    from .cgyro_fonts import GUI_FONT_SIZE, gui_font_family, configure_plot_fonts


def configure_theme(root):
    """Resolve installed fonts once per interpreter; never download at startup."""
    if hasattr(root, "_cgyro_theme"):
        return root._cgyro_theme
    family, cjk = gui_font_family(root)
    configure_plot_fonts()
    fonts = {}
    for name, size, weight in (("body", GUI_FONT_SIZE, "normal"), ("small", 10, "normal"),
                               ("section", GUI_FONT_SIZE, "bold"), ("title", 13, "bold")):
        fonts[name] = tkfont.Font(root=root, family=family, size=size, weight=weight)
    for name in ("TkDefaultFont", "TkTextFont", "TkFixedFont", "TkMenuFont", "TkHeadingFont", "TkTooltipFont",
                 "TkCaptionFont", "TkSmallCaptionFont", "TkIconFont"):
        tkfont.nametofont(name, root=root).configure(family=family, size=GUI_FONT_SIZE)
    style = ttk.Style(root)
    if "clam" in style.theme_names():
        style.theme_use("clam")
    bg, fg, muted, border, accent = "#f3f5f8", "#243246", "#536176", "#cbd3df", "#2364ad"
    root.configure(background=bg)
    root.option_add("*Font", "TkDefaultFont")
    root.option_add("*Listbox.background", "#ffffff")
    root.option_add("*Listbox.foreground", fg)
    root.option_add("*Listbox.selectBackground", accent)
    root.option_add("*Listbox.selectForeground", "#ffffff")
    root.option_add("*Listbox.highlightBackground", border)
    root.option_add("*Text.background", "#ffffff")
    root.option_add("*Text.foreground", fg)
    style.configure(".", font=fonts["body"], background=bg, foreground=fg,
                    bordercolor=border, lightcolor=bg, darkcolor=border)
    style.configure("TEntry", padding=(5, 4), fieldbackground="white")
    style.configure("TCombobox", padding=(5, 4), arrowsize=14)
    style.map("TCombobox", fieldbackground=[("readonly", "white"), ("disabled", "#e8ecf1")],
              foreground=[("disabled", muted), ("readonly", fg)])
    style.configure("TButton", padding=(9, 5))
    style.configure("Compact.TButton", padding=(7, 4))
    style.configure("Accent.TButton", font=fonts["section"], background=accent,
                    foreground="white", padding=(10, 8))
    style.map("Accent.TButton", background=[("disabled", "#d7dfe9"), ("active", "#174d87")],
              foreground=[("disabled", "#536176"), ("!disabled", "white")])
    for widget in ("TButton", "TCheckbutton", "TEntry", "TRadiobutton"):
        style.map(widget, foreground=[("disabled", muted)])
    style.configure("TCheckbutton", padding=(0, 3))
    style.configure("TNotebook", tabmargins=(0, 2, 0, 0))
    style.configure("TNotebook.Tab", padding=(12, 7))
    style.map("TNotebook.Tab", background=[("selected", "#ffffff")],
              foreground=[("selected", accent)])
    style.configure("Treeview", rowheight=fonts["body"].metrics("linespace") + 10,
                    fieldbackground="white", background="white")
    style.configure("Treeview.Heading", font=fonts["section"], padding=(7, 5))
    for name in ("AppTitle.TLabel", "SidebarTitle.TLabel"):
        style.configure(name, font=fonts["title"])
    style.configure("SectionLabel.TLabel", font=fonts["section"])
    for name in ("Muted.TLabel", "Hint.TLabel", "Status.TLabel"):
        style.configure(name, font=fonts["small"], foreground=muted)
    style.configure("Placeholder.TEntry", foreground="#68788f")
    for name in ("Card.TLabelframe", "Inner.TLabelframe"):
        style.configure(name, padding=10, borderwidth=1, relief="solid")
        style.configure(name + ".Label", font=fonts["section"])
    style.configure("Disclosure.TButton", anchor="w", padding=(2, 5), relief="flat",
                    foreground=muted, font=fonts["small"])
    root._cgyro_theme = {"style": style, "fonts": fonts, "cjk": cjk, "background": bg}
    return root._cgyro_theme


def ui_text(root, chinese, english):
    """Use readable English on minimal Linux desktops without a CJK font."""
    configure_theme(root)
    return get_language(root).pair(chinese, english)


def wrapped_text(text, font, width):
    lines = []
    for paragraph in text.split("\n"):
        line = ""
        # Break Latin text at words, CJK text at characters.
        words = paragraph.split(" ") if " " in paragraph else list(paragraph)
        separator = " " if " " in paragraph else ""
        for word in words:
            candidate = line + (separator if line else "") + word
            if line and font.measure(candidate) > width:
                lines.append(line)
                line = word
            else:
                line = candidate
        lines.append(line)
    return "\n".join(lines)


class FormLayout:
    """Lay out complete fields together; reflow without replacing widgets/state."""
    def __init__(self, parent):
        self.parent, self.rows, self._pending, self._key = parent, [], None, None
        parent._cgyro_form = self
        parent.bind("<Configure>", self.schedule, add="+")

    def set_rows(self, rows):
        self.rows = [tuple(row) if isinstance(row, (list, tuple)) else (row,) for row in rows]
        self._key = None
        self.schedule()

    def adopt_grid(self):
        """Adapt legacy forms, preserving row order and label/control ownership."""
        widgets = sorted(self.parent.grid_slaves(), key=lambda w: (int(w.grid_info()["row"]), int(w.grid_info()["column"])))
        rows, i = [], 0
        while i < len(widgets):
            widget = widgets[i]
            if (widget.winfo_class() == "TLabel" and i + 1 < len(widgets)
                    and widgets[i + 1].winfo_class() in ("TEntry", "TCombobox", "TSpinbox")):
                rows.append((widget, widgets[i + 1]))
                i += 2
            else:
                rows.append(widget)
                i += 1
        self.set_rows(rows)

    def schedule(self, _event=None):
        if self._pending is None:
            self._pending = self.parent._root().after_idle(self._layout)

    def _layout(self):
        self._pending = None
        if not self.parent.winfo_exists():
            return
        width = max(80, self.parent.winfo_width() - 22)
        font = tkfont.nametofont("TkDefaultFont", root=self.parent)
        key = (width, tuple(tuple(str(w) for w in row) for row in self.rows),
               tuple(font.actual().items()))
        if key == self._key:
            return
        self._key = key
        pairs = [row for row in self.rows if len(row) > 1]
        label_width = max([font.measure(str(row[0].cget("text"))) for row in pairs] or [0])
        stacked = width < min(label_width, font.measure("M") * 14) + font.measure("0") * 17 + 28
        for col in range(max(4, self.parent.grid_size()[0])):
            self.parent.columnconfigure(col, weight=0, minsize=0, uniform="")
        self.parent.columnconfigure(0 if stacked else 1, weight=1)
        if not stacked:
            self.parent.columnconfigure(0, minsize=min(label_width + 10, int(width * .44)))
        y = 0
        for row in self.rows:
            if len(row) == 1:
                widget = row[0]
                if widget.winfo_class() in ("TCheckbutton", "TRadiobutton"):
                    if not hasattr(widget, "_form_text"):
                        widget._form_text = str(widget.cget("text"))
                    style_font = ttk.Style(widget).lookup(widget.cget("style") or "TCheckbutton", "font")
                    check_font = tkfont.Font(root=widget, font=style_font) if style_font else font
                    widget.configure(text=wrapped_text(widget._form_text, check_font, width - 30))
                elif widget.winfo_class() == "TLabel":
                    widget.configure(wraplength=width)
                widget.grid(row=y, column=0, columnspan=2, rowspan=1, sticky="ew", padx=0, pady=3)
                y += 1
            else:
                label, control = row[:2]
                label.configure(wraplength=width if stacked else max(70, int(width * .42)))
                label.grid(row=y, column=0, columnspan=2 if stacked else 1, rowspan=1,
                           sticky="w", padx=(0, 8), pady=(5, 1) if stacked else 3)
                control.grid(row=y + 1 if stacked else y, column=0 if stacked else 1,
                             columnspan=2 if stacked else 1, rowspan=1, sticky="ew", padx=0, pady=3)
                y += 2 if stacked else 1


class CollapsibleSection(ttk.Frame):
    def __init__(self, parent, title, *, expanded=False):
        super().__init__(parent)
        self.title, self.expanded, self.on_open = title, expanded, None
        self.columnconfigure(0, weight=1)
        self.button = ttk.Button(self, style="Disclosure.TButton", command=self.toggle)
        self.button.grid(row=0, column=0, sticky="ew")
        self.body = ttk.Frame(self, padding=(0, 4, 0, 2))
        self._sync()

    def _sync(self):
        self.button.configure(text=("−  " if self.expanded else "+  ") + self.title)
        if self.expanded:
            self.body.grid(row=1, column=0, sticky="ew")
            if self.on_open:
                self.on_open()
        else:
            self.body.grid_remove()

    def toggle(self):
        self.expanded = not self.expanded
        self._sync()


def bind_local_scroll(widget, canvas):
    """Consume wheel events locally, including comboboxes; never bind_all."""
    def scroll(event):
        if canvas.yview() != (0., 1.):
            up = getattr(event, "num", None) == 4 or getattr(event, "delta", 0) > 0
            canvas.yview_scroll(-2 if up else 2, "units")
        return "break"
    for sequence in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
        widget.bind(sequence, scroll)
    for child in widget.winfo_children():
        # Formula viewports own their scrolling, rather than moving the form.
        target = child if isinstance(child, tk.Canvas) and hasattr(child, "_formula_window_id") else canvas
        bind_local_scroll(child, target)


def show_interface_guide(root):
    """Open the current UI guide without launching a browser."""
    window = tk.Toplevel(root)
    window.title(ui_text(root, "界面与语言说明", "Interface and language guide"))
    window.geometry("820x620")
    window.rowconfigure(0, weight=1)
    window.columnconfigure(0, weight=1)
    view = tk.Text(window, wrap="word", padx=18, pady=14, relief="flat")
    view.grid(row=0, column=0, sticky="nsew")
    scroll = ttk.Scrollbar(window, command=view.yview)
    scroll.grid(row=0, column=1, sticky="ns")
    view.configure(yscrollcommand=scroll.set)
    try:
        filename = "GUI_GUIDE.md" if get_language(root).chinese else "GUI_GUIDE_EN.md"
        content = (Path(__file__).resolve().parent / filename).read_text(encoding="utf-8")
    except OSError as exc:
        content = str(exc)
    view.insert("1.0", content)
    view.configure(state="disabled")


class PresetStore:
    """One data-only JSON file per app; malformed files are never overwritten."""
    def __init__(self, name, path=None):
        base = Path(os.environ.get("CGYRO_PRESET_DIR", str(Path.home() / ".cgyro_comparison_tool")))
        self.path = Path(path) if path else base / (name + "_presets.json")

    def load(self):
        if not self.path.exists():
            return {}
        data = json.loads(self.path.read_text(encoding="utf-8-sig"))
        if not isinstance(data, dict) or data.get("format") != "cgyro-presets-1" or not isinstance(data.get("presets"), dict):
            raise ValueError("Invalid preset file: " + str(self.path))
        if not all(isinstance(k, str) and isinstance(v, dict) for k, v in data["presets"].items()):
            raise ValueError("Invalid preset entries: " + str(self.path))
        return data["presets"]

    def save(self, name, state):
        name = name.strip()
        if not name:
            raise ValueError("Enter a preset name.")
        presets = self.load()
        presets[name] = state
        write_json_atomic(self.path, {"format": "cgyro-presets-1", "presets": presets})


class PresetBar(ttk.Frame):
    def __init__(self, parent, *, root, name, capture, apply, status):
        super().__init__(parent)
        self.root, self.capture, self.apply, self.status = root, capture, apply, status
        self.store = PresetStore(name)
        self.columnconfigure(0, weight=1)
        self.choice = tk.StringVar()
        self.combo = ttk.Combobox(self, textvariable=self.choice, width=12, state="readonly")
        self.combo.grid(row=0, column=0, sticky="ew", padx=(0, 5))
        self.combo.bind("<<ComboboxSelected>>", self._apply)
        self.save_button = ttk.Button(self, text=ui_text(root, "保存预设", "Save preset"), command=self._save)
        self.save_button.grid(row=0, column=1)
        self.refresh()

    def refresh(self):
        try:
            self.presets = self.store.load()
        except (ValueError, OSError) as exc:
            self.presets = {}
            self.status.set(str(exc))
        self.combo.configure(values=list(self.presets))
        if not self.choice.get():
            self.choice.set(ui_text(self.root, "选择预设…", "Choose preset…"))

    def retranslate(self):
        if self.choice.get() not in self.presets:
            self.choice.set(ui_text(self.root, "选择预设…", "Choose preset…"))

    def _apply(self, _event=None):
        try:
            self.apply(self.presets[self.choice.get()])
            self.status.set(ui_text(self.root, "已应用预设，可继续修改", "Preset applied; settings remain editable."))
        except (ValueError, TypeError, KeyError, tk.TclError) as exc:
            messagebox.showerror("Preset", str(exc), parent=self.root)

    def _save(self):
        try:
            state = self.capture()
            name = simpledialog.askstring(ui_text(self.root, "保存预设", "Save preset"),
                                          ui_text(self.root, "预设名称（同名覆盖）", "Name (same name replaces preset)"), parent=self.root)
            if name and name.strip():
                self.store.save(name, state)
                self.refresh()
                self.choice.set(name.strip())
                self.status.set(ui_text(self.root, "已保存预设：", "Preset saved: ") + name.strip())
        except (ValueError, TypeError, KeyError, OSError, tk.TclError) as exc:
            messagebox.showerror("Preset", str(exc), parent=self.root)
