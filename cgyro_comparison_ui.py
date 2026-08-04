"""
UI and case-management mixin for CGYRO comparison GUI.
"""

import os
import re
import json
import sys
import subprocess
import tempfile
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np

try:
    from cgyro_comparison_bootstrap import (
        cgyrodata,
        cgyrodata_plot,
        DEFAULT_APP_TITLE,
        DEFAULT_WINDOW_GEOMETRY,
        DEFAULT_CASE_PICKER_ROOT,
        DEFAULT_LINEAR_GAMMA_FILE,
        default_share_dir,
    )
except ImportError:
    from .cgyro_comparison_bootstrap import (
        cgyrodata,
        cgyrodata_plot,
        DEFAULT_APP_TITLE,
        DEFAULT_WINDOW_GEOMETRY,
        DEFAULT_CASE_PICKER_ROOT,
        DEFAULT_LINEAR_GAMMA_FILE,
        default_share_dir,
    )

try:
    from cgyro_update import APP_VERSION, REPOSITORY_URL, UpdateCheckError, check_for_updates
except ImportError:
    from .cgyro_update import APP_VERSION, REPOSITORY_URL, UpdateCheckError, check_for_updates


class CgyroUiMixin:
    # Keep UI option literals centralized.  Plot dispatch below depends on
    # exact strings, so collecting them here makes later menu changes less
    # likely to desynchronize controls and backend handlers.
    _FLUX_TYPE_OPTIONS = ("Energy", "Particle")
    _FLUX_XAXIS_OPTIONS = ("v.s ky", "v.s kx (estimated)", "v.s Time", "v.s ky_time", "v.s 2D")
    _FLUC_FIELD_OPTIONS = ("Phi", "Apar", "Bpar")
    _FLUC_XAXIS_OPTIONS = ("v.s ky", "v.s kx", "v.s Time", "v.s theta", "fft")
    _FLUC_MOMENT_OPTIONS = ("Phi", "Density", "Energy", "Temperature", "Apar", "Bpar")
    _FLUC2D_VIEW_OPTIONS = ("vs xy", "vs xt", "vs kxky")
    # User-facing Fluctuation-2D views are translated into exact internal
    # plot-type tokens.  The exact token matters because "vs kxky" also
    # contains "vs kx"; exact dispatch must catch it before generic 1D logic.
    _FLUC2D_VIEW_PLOT_TYPES = {
        "vs xy": ("Fluctuation 2D", "Fluctuation 2D (vs xy)"),
        "vs xt": ("Fluctuation 2D vs xt", "Fluctuation 2D (vs xt)"),
        "vs kxky": ("Fluctuation 2D vs kxky", "Fluctuation 2D (vs kxky)"),
    }
    _FFT_MODE_OPTIONS = ("Nonlinear", "Linear")
    _FFT_VIEW_OPTIONS = ("Omega vs ky", "Omega vs kx")
    _FFT_SPECTRUM_OPTIONS = ("Amplitude", "Power")
    _ZF_XAXIS_OPTIONS = ("vs Time", "vs kx", "phi vs kx(theta=0)", "vs gamma_lin")
    _ENERGY_BALANCE_MODE_OPTIONS = (
        "Entropy balance",
        "ZF energy balance",
        "Effective growth rate",
        "Single plot",
        "FULLT transfer map",
        "2D scan",
    )
    _ENERGY_BALANCE_SPEC_OPTIONS = ("Total (-1)", "Main ion (0)", "Electron (1)")
    _ENERGY_BALANCE_SINGLE_QUANTITY_OPTIONS = ("T", "N", "T-N", "Dr", "Dtheta", "Dc", "DZ", "entropy")
    _ENERGY_BALANCE_SINGLE_XAXIS_OPTIONS = ("vs Time", "vs ky", "vs kx", "vs kxky")
    _ENERGY_BALANCE_SINGLE_NORM_OPTIONS = ("None", "|min(T)|", "|max(T)|")
    _ENERGY_BALANCE_SINGLE_ENTROPY_NORM_OPTIONS = ("None", "min entropy", "max entropy", "ky=0 entropy")
    _ENERGY_BALANCE_TRANSFER_XAXIS_OPTIONS = ("vs kxky", "vs kx", "trace pxqx")
    _ENERGY_BALANCE_TRANSFER_QUANTITY_OPTIONS = ("Re", "Im", "Abs")
    _OTHERS_PLOT_OPTIONS = ("Error", "rcorr_phi", "POD_parity")
    _OTHERS_FIELD_OPTIONS = ("Phi", "Apar", "Bpar")
    _OTHERS_POD_FIELD_OPTIONS = ("Apar", "Phi")
    # Workspace files persist Tk variables, not widget objects.  Cases are
    # stored separately as paths so a workspace remains portable between GUI
    # sessions and does not attempt to pickle pygacode data objects.
    _WORKSPACE_STATE_VARS = [
        "lock_case_selection_var",
        "fluc_average_var",
        "time_mode_var",
        "time_percent_var",
        "time_duration_var",
        "plot_type_var",
        "t_start_var",
        "t_end_var",
        "log_x_var",
        "log_y_var",
        "axis_kx_min_var",
        "axis_kx_max_var",
        "axis_ky_min_var",
        "axis_ky_max_var",
        "norm_ky_var",
        "flux_type_var",
        "flux_xaxis_var",
        "flux_scan_xparam_var",
        "flux_decomp_var",
        "flux_2d_errorbar_var",
        "flux_norm_real_ion_var",
        "fluc_field_var",
        "fluc_xaxis_var",
        "fluc_norm_max_var",
        "fluc_advanced_var",
        "fluc_theta_kx_var",
        "fluc_theta_ky_var",
        "species_var",
        "plot_all_species_var",
        "moment_var",
        "fluc2d_view_var",
        "fluc2d_x_elec_var",
        "fluc2d_log_z_var",
        "linear_gamma_file_var",
        "fft_mode_var",
        "fft_view_var",
        "fft_spectrum_var",
        "fft_overlay_freq_var",
        "zf_xaxis_var",
        "zf_gamma_lin_ky_var",
        "energy_balance_mode_var",
        "energy_balance_n_var",
        "energy_balance_spec_var",
        "energy_balance_single_quantity_var",
        "energy_balance_single_xaxis_var",
        "energy_balance_single_norm_var",
        "energy_balance_single_norm_entropy_var",
        "energy_balance_transfer_quantity_var",
        "energy_balance_transfer_xaxis_var",
        "energy_balance_transfer_ky_var",
        "energy_balance_transfer_kx_var",
        "energy_balance_transfer_asym_var",
        "energy_balance_transfer_norm_max_var",
        "others_plot_var",
        "others_rcorr_field_var",
        "others_rcorr_theta_var",
        "others_pod_field_var",
        "others_pod_kx_var",
        "others_pod_ky_var",
    ]

    def _create_formula_panel(self, master, figsize=(4.0, 2.4), dpi=100):
        """Create a scrollable matplotlib formula panel with consistent behavior."""
        frame = ttk.Frame(master)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)

        fig = plt.Figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111)
        ax.axis("off")
        fig.patch.set_facecolor("white")

        viewport = tk.Canvas(frame, highlightthickness=0, borderwidth=0, bg="white")
        viewport.grid(row=0, column=0, sticky=tk.NSEW, padx=(0, 2), pady=(0, 2))

        canvas = FigureCanvasTkAgg(fig, master=viewport)
        widget = canvas.get_tk_widget()
        window_id = viewport.create_window((0, 0), window=widget, anchor=tk.NW)
        viewport._formula_inner_widget = widget
        viewport._formula_window_id = window_id

        vscroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=viewport.yview)
        vscroll.grid(row=0, column=1, sticky=tk.NS)
        hscroll = ttk.Scrollbar(frame, orient=tk.HORIZONTAL, command=viewport.xview)
        hscroll.grid(row=1, column=0, sticky=tk.EW, padx=(0, 2))
        viewport.configure(xscrollcommand=hscroll.set, yscrollcommand=vscroll.set)

        return frame, fig, ax, canvas, viewport, vscroll, hscroll

    def __init__(self, root):
        """Initialize UI state, widgets, and menus."""
        self.root = root
        self.root.title(DEFAULT_APP_TITLE)
        self.root.geometry(DEFAULT_WINDOW_GEOMETRY)
        self.fluc_average_var = tk.StringVar(value="Root Mean Square")
        self.time_mode_var = tk.StringVar(value="Manual Start/End")
        self.time_percent_var = tk.StringVar(value="50")
        self.time_duration_var = tk.StringVar(value="")
        self.axis_kx_min_var = tk.StringVar()
        self.axis_kx_max_var = tk.StringVar()
        self.axis_ky_min_var = tk.StringVar()
        self.axis_ky_max_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Ready")
        self.case_summary_var = tk.StringVar(value="0 cases loaded")
        self.lock_case_selection_var = tk.BooleanVar(value=False)

        self.cases = {}  # Dictionary to store loaded cases: {name: cgyrodata_object}
        self._locked_case_names = []
        self._restoring_case_selection = False
        # Optional per-case selectors used by Fluctuation 1D Advanced mode.
        # Each blank value means average over that spectral axis.
        self._fluc_advanced_case_values = {}
        self._fluc_advanced_seeded = False
        self.ani = None # Animation object
        self.is_paused = False
        self.current_frame = 0
        self.total_frames = 0
        self._drag_start_index = None
        self._manual_pager_active = False
        self._manual_pager_label = "Page"
        self._update_check_in_progress = False
        self._update_check_menu = None
        self._auto_workspace_path = os.environ.get("CGYRO_AUTO_WORKSPACE", "").strip()
        
        self._create_layout()
        self._create_menu()
        if self._auto_workspace_path:
            self.root.after(100, self._restore_auto_workspace)
        # Do not hold up startup or show an error for an offline workstation.
        self.root.after(1500, self._check_for_updates_silently)

    def _configure_ui_styles(self):
        """Apply a small, theme-friendly visual system to the main window."""
        style = ttk.Style(self.root)
        self._ui_style = style
        try:
            style.configure("App.TFrame", padding=0)
            style.configure("AppTitle.TLabel", font=("Segoe UI", 13, "bold"))
            style.configure("SidebarTitle.TLabel", font=("Segoe UI", 12, "bold"))
            style.configure("SectionLabel.TLabel", font=("Segoe UI", 9, "bold"))
            style.configure("Muted.TLabel", foreground="#667085", font=("Segoe UI", 9))
            style.configure("Hint.TLabel", foreground="#7a8493", font=("Segoe UI", 8))
            style.configure("Status.TLabel", foreground="#667085", font=("Segoe UI", 9))
            style.configure("Card.TLabelframe", padding=8)
            style.configure("Card.TLabelframe.Label", font=("Segoe UI", 10, "bold"))
            style.configure("Inner.TLabelframe", padding=6)
            style.configure("Inner.TLabelframe.Label", font=("Segoe UI", 9, "bold"))
            style.configure("Compact.TButton", padding=(7, 4))
            style.configure("Accent.TButton", font=("Segoe UI", 10, "bold"), padding=(10, 7))
        except tk.TclError:
            # Some platform themes expose fewer style options.  The default
            # ttk theme is still fully usable when a custom option is rejected.
            pass

    def _set_initial_pane_position(self):
        """Set a comfortable initial control-panel width after geometry settles."""
        try:
            pane_width = int(self.main_pane.winfo_width())
            if pane_width <= 0:
                self.root.after(50, self._set_initial_pane_position)
                return
            left_min = int(getattr(self, "_left_pane_min_width", 340))
            right_min = int(getattr(self, "_right_pane_min_width", 500))
            max_position = max(left_min, pane_width - right_min)
            desired_position = min(440, max(left_min, int(pane_width * 0.32)))
            position = min(desired_position, max_position)
            self.main_pane.sashpos(0, position)
            self._enforce_pane_minimums()
        except (AttributeError, tk.TclError, ValueError):
            pass

    def _enforce_pane_minimums(self, _event=None):
        """Keep both panes usable on Tk versions without a ``minsize`` option."""
        try:
            pane_width = int(self.main_pane.winfo_width())
            if pane_width <= 0:
                return
            left_min = int(getattr(self, "_left_pane_min_width", 340))
            right_min = int(getattr(self, "_right_pane_min_width", 500))
            current_position = int(self.main_pane.sashpos(0))
            max_position = max(left_min, pane_width - right_min)
            position = min(max(current_position, left_min), max_position)
            if position != current_position:
                self.main_pane.sashpos(0, position)
        except (AttributeError, tk.TclError, ValueError):
            pass

    def _refresh_case_summary(self, status_prefix="Ready"):
        """Refresh the compact case count and status-bar summary."""
        try:
            count = int(self.case_listbox.size())
            selected = len(self.case_listbox.curselection())
        except (AttributeError, tk.TclError):
            count = len(getattr(self, "cases", {}))
            selected = 0

        case_word = "case" if count == 1 else "cases"
        summary = f"{count} {case_word} loaded"
        if selected:
            summary += f" | {selected} selected"
        self.case_summary_var.set(summary)

        plot_var = getattr(self, "plot_type_var", None)
        plot_name = str(plot_var.get()).strip() if plot_var is not None else ""
        status = f"{status_prefix} | {summary}"
        if plot_name:
            status += f" | {plot_name}"
        self.status_var.set(status)

    def _create_layout(self):
        """Create top-level panels, control widgets, and matplotlib canvas."""
        self._configure_ui_styles()

        # Keep a small, always-visible status strip separate from the scrollable
        # control panel so long option lists never hide important feedback.
        status_frame = ttk.Frame(self.root, padding=(10, 4))
        status_frame.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Separator(status_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, side=tk.TOP)
        ttk.Label(
            status_frame,
            textvariable=self.status_var,
            style="Status.TLabel",
        ).pack(anchor=tk.W, pady=(3, 0))

        # A PanedWindow makes the control area resizable instead of locking it
        # to a fixed 420 px column on every screen size.
        main_frame = ttk.Frame(self.root, style="App.TFrame", padding=(8, 8, 8, 0))
        main_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self._left_pane_min_width = 340
        self._right_pane_min_width = 500
        self.main_pane = ttk.PanedWindow(main_frame, orient=tk.HORIZONTAL)
        self.main_pane.pack(fill=tk.BOTH, expand=True)
        self.main_pane.bind("<Configure>", self._enforce_pane_minimums, add="+")
        self.main_pane.bind("<B1-Motion>", self._enforce_pane_minimums, add="+")
        self.main_pane.bind("<ButtonRelease-1>", self._enforce_pane_minimums, add="+")

        # Left sidebar uses three independent zones: case selection stays fixed
        # at the top, only plot options scroll in the middle, and the primary
        # actions remain available at the bottom.  This avoids losing context
        # when a plot mode exposes a long set of dynamic controls.
        left_container = ttk.Frame(self.main_pane, width=400)
        self.main_pane.add(left_container, weight=0)
        left_container.pack_propagate(False)
        left_container.columnconfigure(0, weight=1)
        left_container.rowconfigure(2, weight=1)
        self.left_container = left_container

        # Fixed sidebar header.
        header = ttk.Frame(left_container, padding=(12, 10, 12, 6), style="App.TFrame")
        header.grid(row=0, column=0, sticky=tk.EW)
        self.sidebar_header = header
        header.columnconfigure(0, weight=1)
        title_row = ttk.Frame(header, style="App.TFrame")
        title_row.grid(row=0, column=0, sticky=tk.EW)
        title_row.columnconfigure(0, weight=1)
        ttk.Label(
            title_row,
            text=DEFAULT_APP_TITLE,
            style="SidebarTitle.TLabel",
        ).grid(row=0, column=0, sticky=tk.W)
        ttk.Label(
            title_row,
            text=f"v{APP_VERSION}",
            style="Muted.TLabel",
        ).grid(row=0, column=1, sticky=tk.E)
        ttk.Label(
            header,
            textvariable=self.case_summary_var,
            style="Muted.TLabel",
        ).grid(row=1, column=0, sticky=tk.W, pady=(2, 0))

        # Fixed case selection card.
        case_section = ttk.LabelFrame(
            left_container,
            text="Cases",
            padding=8,
            style="Card.TLabelframe",
        )
        case_section.grid(row=1, column=0, sticky=tk.EW, padx=12, pady=(0, 8))
        self.case_section = case_section

        case_list_frame = ttk.Frame(case_section)
        case_list_frame.pack(fill=tk.X)
        case_list_frame.columnconfigure(0, weight=1)
        case_list_frame.rowconfigure(0, weight=1)
        self.case_listbox = tk.Listbox(
            case_list_frame,
            selectmode=tk.EXTENDED,
            height=6,
            activestyle="none",
            relief=tk.FLAT,
            borderwidth=1,
            highlightthickness=1,
            highlightcolor="#7aa2d6",
            selectbackground="#2f6fad",
            selectforeground="white",
            exportselection=not bool(self.lock_case_selection_var.get()),
        )
        self.case_listbox.grid(row=0, column=0, sticky=tk.NSEW)
        self.case_list_scrollbar = ttk.Scrollbar(
            case_list_frame,
            orient=tk.VERTICAL,
            command=self.case_listbox.yview,
        )
        self.case_list_scrollbar.grid(row=0, column=1, sticky=tk.NS, padx=(4, 0))
        self.case_listbox.configure(yscrollcommand=self.case_list_scrollbar.set)

        # Enable drag and drop reordering.
        self.case_listbox.bind('<Button-1>', self._on_drag_start)
        self.case_listbox.bind('<B1-Motion>', self._on_drag_motion)
        self.case_listbox.bind('<<ListboxSelect>>', self._on_case_listbox_select)
        # A bind_all wheel handler is used for the outer options panel. Handle
        # the case list first so it scrolls independently.
        self.case_listbox.bind(
            '<MouseWheel>', self._on_case_listbox_mousewheel, add='+'
        )
        self.case_listbox.bind(
            '<Button-4>', self._on_case_listbox_mousewheel, add='+'
        )
        self.case_listbox.bind(
            '<Button-5>', self._on_case_listbox_mousewheel, add='+'
        )

        selection_row = ttk.Frame(case_section)
        selection_row.pack(fill=tk.X, pady=(5, 5))
        self.lock_case_selection_check = ttk.Checkbutton(
            selection_row,
            text="Lock selection",
            variable=self.lock_case_selection_var,
            command=self._on_lock_case_selection_toggle,
        )
        self.lock_case_selection_check.pack(side=tk.LEFT)
        ttk.Label(
            selection_row,
            text="No selection = all cases",
            style="Hint.TLabel",
        ).pack(side=tk.RIGHT, padx=(8, 0))

        case_buttons = ttk.Frame(case_section)
        case_buttons.pack(fill=tk.X)
        for column in range(3):
            case_buttons.columnconfigure(column, weight=1, uniform="case_actions")
        ttk.Button(
            case_buttons,
            text="Add Case",
            command=self.add_case_single,
            style="Compact.TButton",
        ).grid(row=0, column=0, padx=(0, 2), pady=(0, 3), sticky=tk.EW)
        ttk.Button(
            case_buttons,
            text="Add Multiple",
            command=self.add_case_multiple,
            style="Compact.TButton",
        ).grid(row=0, column=1, padx=2, pady=(0, 3), sticky=tk.EW)
        ttk.Button(
            case_buttons,
            text="Add Group",
            command=self.add_group,
            style="Compact.TButton",
        ).grid(row=0, column=2, padx=(2, 0), pady=(0, 3), sticky=tk.EW)
        ttk.Button(
            case_buttons,
            text="Remove",
            command=self.remove_case,
            style="Compact.TButton",
        ).grid(row=1, column=0, padx=(0, 2), sticky=tk.EW)
        ttk.Button(
            case_buttons,
            text="Reload",
            command=self.reload_cases,
            style="Compact.TButton",
        ).grid(row=1, column=1, padx=2, sticky=tk.EW)
        ttk.Button(
            case_buttons,
            text="Remove All",
            command=self.remove_all_cases,
            style="Compact.TButton",
        ).grid(row=1, column=2, padx=(2, 0), sticky=tk.EW)

        # Scrollable plot-options region.  Cases and actions deliberately live
        # outside this canvas, so only the potentially long dynamic form moves.
        scroll_host = ttk.Frame(left_container)
        scroll_host.grid(row=2, column=0, sticky=tk.NSEW)
        self.options_scroll_host = scroll_host
        scroll_host.columnconfigure(0, weight=1)
        scroll_host.rowconfigure(0, weight=1)
        self.left_scrollbar = ttk.Scrollbar(scroll_host, orient=tk.VERTICAL)
        self.left_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.left_canvas = tk.Canvas(
            scroll_host,
            height=80,
            highlightthickness=0,
            borderwidth=0,
        )
        self.left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.left_canvas.configure(yscrollcommand=self.left_scrollbar.set)
        self.left_scrollbar.configure(command=self.left_canvas.yview)
        try:
            canvas_bg = self._ui_style.lookup("TFrame", "background")
            if canvas_bg:
                self.left_canvas.configure(background=canvas_bg)
        except (AttributeError, tk.TclError):
            pass

        left_panel = ttk.Frame(
            self.left_canvas,
            padding=(12, 0, 8, 10),
            style="App.TFrame",
        )
        self.left_panel = left_panel
        self._left_panel_window = self.left_canvas.create_window((0, 0), window=left_panel, anchor=tk.NW)
        left_panel.bind("<Configure>", self._on_left_panel_configure)
        self.left_canvas.bind("<Configure>", self._on_left_canvas_configure)
        self.left_canvas.bind("<Enter>", self._on_left_panel_enter)
        self.left_canvas.bind("<Leave>", self._on_left_panel_leave)
        left_panel.bind("<Enter>", self._on_left_panel_enter)
        left_panel.bind("<Leave>", self._on_left_panel_leave)

        self.root.after_idle(self._set_initial_pane_position)

        # Plot controls are the only content in the scrollable middle region.
        plot_section = ttk.LabelFrame(left_panel, text="Plot setup", padding=8, style="Card.TLabelframe")
        plot_section.pack(fill=tk.X)
        ttk.Label(plot_section, text="Plot type:", style="SectionLabel.TLabel").pack(anchor=tk.W)
        
        self.plot_type_var = tk.StringVar(value="Frequency")
        plot_types = [
            "Frequency",
            "Growth Rate",
            "Flux",
            "Fluctuation 1D",
            "Fluctuation 2D",
            "Energy balance",
            "Zonal ExB Shearing Rate",
            "Others",
        ]
        self.plot_type_combo = ttk.Combobox(plot_section, textvariable=self.plot_type_var, values=plot_types, state="readonly")
        self.plot_type_combo.pack(fill=tk.X, pady=(3, 0))
        self.plot_type_combo.bind("<<ComboboxSelected>>", self.update_options)

        # Dynamic Options Frame
        self.options_frame = ttk.LabelFrame(plot_section, text="Options", padding=7, style="Inner.TLabelframe")
        self.options_frame.pack(fill=tk.X, pady=(8, 0))
        self._init_options()
        self._refresh_case_summary()

        # Fixed action footer.  Plot and navigation remain reachable even when
        # the current option form is much taller than the window.
        action_footer = ttk.Frame(left_container, padding=(12, 6, 12, 10))
        action_footer.grid(row=3, column=0, sticky=tk.EW)
        self.action_footer = action_footer
        ttk.Separator(action_footer, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=(0, 6))
        action_section = ttk.LabelFrame(
            action_footer,
            text="Actions",
            padding=8,
            style="Card.TLabelframe",
        )
        action_section.pack(fill=tk.X)
        self.action_section = action_section
        action_section.columnconfigure(0, weight=1)
        action_section.columnconfigure(1, weight=1)
        action_section.columnconfigure(2, weight=1)
        self.plot_button = ttk.Button(
            action_section,
            text="Plot",
            command=self.plot_comparison,
            style="Accent.TButton",
        )
        self.plot_button.grid(row=0, column=0, columnspan=3, sticky=tk.EW, pady=(0, 5))
        ttk.Button(
            action_section,
            text="Case info",
            command=self.plot_case_info,
            style="Compact.TButton",
        ).grid(row=1, column=0, padx=(0, 2), sticky=tk.EW)
        ttk.Button(
            action_section,
            text="Diff input",
            command=self.plot_input_diff,
            style="Compact.TButton",
        ).grid(row=1, column=1, padx=2, sticky=tk.EW)
        ttk.Button(
            action_section,
            text="Clear plot",
            command=self.clear_plot,
            style="Compact.TButton",
        ).grid(row=1, column=2, padx=(2, 0), sticky=tk.EW)

        ttk.Separator(action_section, orient=tk.HORIZONTAL).grid(
            row=2, column=0, columnspan=3, sticky=tk.EW, pady=(7, 5)
        )
        navigation_header = ttk.Frame(action_section)
        navigation_header.grid(row=3, column=0, columnspan=3, sticky=tk.EW, pady=(0, 4))
        navigation_header.columnconfigure(0, weight=1)
        ttk.Label(
            navigation_header,
            text="Animation / pages",
            style="SectionLabel.TLabel",
        ).grid(row=0, column=0, sticky=tk.W)
        ttk.Label(
            navigation_header,
            text="enabled when available",
            style="Hint.TLabel",
        ).grid(row=0, column=1, sticky=tk.E)

        self.anim_controls_frame = ttk.Frame(action_section)
        self.anim_controls_frame.grid(row=4, column=0, columnspan=3, sticky=tk.EW)
        for column in range(3):
            self.anim_controls_frame.columnconfigure(column, weight=1, uniform="animation")

        self.btn_prev = ttk.Button(self.anim_controls_frame, text="< Prev", command=self.prev_frame, state="disabled")
        self.btn_prev.grid(row=0, column=0, padx=(0, 2), sticky=tk.EW)
        
        self.btn_pause = ttk.Button(self.anim_controls_frame, text="Pause", command=self.toggle_pause, state="disabled")
        self.btn_pause.grid(row=0, column=1, padx=2, sticky=tk.EW)
        
        self.btn_next = ttk.Button(self.anim_controls_frame, text="Next >", command=self.next_frame, state="disabled")
        self.btn_next.grid(row=0, column=2, padx=(2, 0), sticky=tk.EW)

        # Right panel: Plot Area
        right_panel = ttk.Frame(self.main_pane, padding=10)
        self.main_pane.add(right_panel, weight=1)
        self.right_panel = right_panel
        right_panel.rowconfigure(0, weight=1)
        right_panel.rowconfigure(1, weight=0, minsize=38)
        right_panel.columnconfigure(0, weight=1)

        self.fig = plt.Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right_panel)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.NSEW)

        toolbar_frame = ttk.Frame(right_panel)
        toolbar_frame.grid(row=1, column=0, sticky=tk.EW)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side=tk.LEFT, fill=tk.X, expand=True)

    def _bind_menu_shortcut(self, sequence, command, widget=None):
        """Bind one menu shortcut and stop duplicate widget processing."""
        target = widget if widget is not None else self.root

        def invoke(_event=None):
            command()
            return "break"

        target.bind(sequence, invoke, add="+")

    def _add_menu_command(
        self,
        menu,
        label,
        command,
        accelerator=None,
        shortcut=None,
        shortcut_widget=None,
    ):
        """Add a consistently labelled menu command and optional shortcut."""
        options = {"label": label, "command": command}
        if accelerator:
            options["accelerator"] = accelerator
        menu.add_command(**options)
        if shortcut:
            self._bind_menu_shortcut(shortcut, command, widget=shortcut_widget)

    def _create_menu(self):
        """Create a compact menu bar grouped by user workflow."""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        self.menubar = menubar

        # File contains workspace-level operations only.
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu, underline=0)
        self._add_menu_command(
            file_menu,
            "Open Workspace...",
            self.load_workspace,
            accelerator="Ctrl+O",
            shortcut="<Control-o>",
        )
        self._add_menu_command(
            file_menu,
            "Save Workspace...",
            self.save_workspace,
            accelerator="Ctrl+S",
            shortcut="<Control-s>",
        )
        file_menu.add_separator()
        self._add_menu_command(
            file_menu,
            "Exit",
            self.root.quit,
            accelerator="Ctrl+Q",
            shortcut="<Control-q>",
        )
        self.file_menu = file_menu

        # Case management is kept separate from workspace file operations.
        cases_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Cases", menu=cases_menu, underline=0)
        self._add_menu_command(
            cases_menu,
            "Add Case...",
            self.add_case_single,
            accelerator="Ctrl+N",
            shortcut="<Control-n>",
        )
        self._add_menu_command(
            cases_menu,
            "Add Multiple Cases...",
            self.add_case_multiple,
            accelerator="Ctrl+Shift+N",
            shortcut="<Control-Shift-N>",
        )
        self._add_menu_command(cases_menu, "Add Case Group...", self.add_group)
        cases_menu.add_separator()
        self._add_menu_command(
            cases_menu,
            "Remove Selected",
            self.remove_case,
            accelerator="Delete",
            shortcut="<Delete>",
            shortcut_widget=self.case_listbox,
        )
        self._add_menu_command(cases_menu, "Remove All...", self.remove_all_cases)
        self._add_menu_command(
            cases_menu,
            "Reload Cases",
            self.reload_cases,
            accelerator="F5",
            shortcut="<F5>",
        )
        self.cases_menu = cases_menu

        data_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Data", menu=data_menu, underline=0)
        self._add_menu_command(
            data_menu,
            "Convert Binary to Readable...",
            self.transfer_bin_to_readable,
        )
        self._add_menu_command(
            data_menu,
            "Export Current Plot Data...",
            self.save_current_plot_data,
            accelerator="Ctrl+Shift+E",
            shortcut="<Control-Shift-E>",
        )
        self.data_menu = data_menu

        plot_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Plot", menu=plot_menu, underline=0)
        self._add_menu_command(
            plot_menu,
            "Plot Now",
            self.plot_comparison,
            accelerator="Ctrl+Enter",
            shortcut="<Control-Return>",
        )
        self._add_menu_command(plot_menu, "Clear Plot", self.clear_plot)
        plot_menu.add_separator()

        time_menu = tk.Menu(plot_menu, tearoff=0)
        plot_menu.add_cascade(label="Time Window", menu=time_menu)
        time_menu.add_radiobutton(
            label="Manual Start / End",
            variable=self.time_mode_var,
            value="Manual Start/End",
        )
        time_menu.add_radiobutton(
            label="Last Percentage...",
            variable=self.time_mode_var,
            value="Last percent",
            command=self._set_time_last_percent,
        )
        time_menu.add_radiobutton(
            label="Last Duration...",
            variable=self.time_mode_var,
            value="Last duration",
            command=self._set_time_last_duration,
        )
        time_menu.add_radiobutton(
            label="Full Range",
            variable=self.time_mode_var,
            value="Full range",
            command=self._clear_simple_time_entries,
        )
        time_menu.add_separator()
        time_menu.add_command(label="Reset Time Window", command=self.clear_time_range)

        average_menu = tk.Menu(plot_menu, tearoff=0)
        plot_menu.add_cascade(label="Averaging", menu=average_menu)
        average_menu.add_radiobutton(
            label="Mean Absolute",
            variable=self.fluc_average_var,
            value="Mean Absolute",
            command=self.update_options,
        )
        average_menu.add_radiobutton(
            label="Root Mean Square",
            variable=self.fluc_average_var,
            value="Root Mean Square",
            command=self.update_options,
        )

        axis_menu = tk.Menu(plot_menu, tearoff=0)
        plot_menu.add_cascade(label="Axis", menu=axis_menu)
        axis_menu.add_command(
            label="Set Axis Limits...",
            command=self.open_axis_limits_dialog,
        )
        axis_menu.add_command(label="Clear Axis Limits", command=self.clear_axis_limits)
        self.plot_menu = plot_menu
        self.time_menu = time_menu
        self.average_menu = average_menu
        self.axis_menu = axis_menu

        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu, underline=0)
        help_menu.add_command(label="Check for Updates...", command=self.check_for_updates)
        help_menu.add_separator()
        help_menu.add_command(label="About", command=self.show_about)
        self.help_menu = help_menu
        self._update_check_menu = help_menu

    def _check_for_updates_silently(self):
        """Run the startup check without interrupting offline users."""
        self._start_update_check(silent=True)

    def check_for_updates(self):
        """Start a user-requested update check in the background."""
        self._start_update_check(silent=False)

    def _start_update_check(self, silent):
        if self._update_check_in_progress:
            if not silent:
                messagebox.showinfo(
                    "Check for Updates",
                    "An update check is already in progress.",
                    parent=self.root,
                )
            return

        self._update_check_in_progress = True
        if self._update_check_menu is not None:
            self._update_check_menu.entryconfig("Check for Updates...", state="disabled")

        worker = threading.Thread(
            target=self._update_check_worker,
            args=(silent,),
            name="cgyro-update-check",
            daemon=True,
        )
        worker.start()

    def _update_check_worker(self, silent):
        try:
            result = check_for_updates(APP_VERSION)
        except UpdateCheckError as exc:
            self._post_update_check_result(silent, None, str(exc))
        except Exception as exc:  # Keep an unexpected network/parser error from killing the worker.
            self._post_update_check_result(silent, None, str(exc))
        else:
            self._post_update_check_result(silent, result, None)

    def _post_update_check_result(self, silent, result, error):
        """Hand the worker result back to Tk unless the window was closed."""
        try:
            self.root.after(0, self._finish_update_check, silent, result, error)
        except tk.TclError:
            pass

    def _finish_update_check(self, silent, result, error):
        self._update_check_in_progress = False
        if self._update_check_menu is not None:
            self._update_check_menu.entryconfig("Check for Updates...", state="normal")

        if error:
            if not silent:
                messagebox.showwarning("Check for Updates", error, parent=self.root)
            return

        if result.update_available:
            release_label = f"\n{result.release_name}" if result.release_name else ""
            update_now = messagebox.askyesno(
                "Update Available",
                f"A newer version of CGYRO Comparison Tool is available.\n\n"
                f"Current version: {result.current_version}\n"
                f"Latest version: {result.latest_version}{release_label}\n\n"
                "Download the update with git pull now?\n"
                "After downloading, the application will restart and restore the current workspace.",
                parent=self.root,
            )
            if update_now:
                self._pull_update_from_git(result.latest_version)
        elif not silent:
            messagebox.showinfo(
                "Check for Updates",
                f"You are using the latest version ({result.current_version}).",
                parent=self.root,
            )

    def _pull_update_from_git(self, latest_version):
        """Fast-forward the current Git checkout after explicit user approval."""
        app_dir = os.path.dirname(os.path.abspath(__file__))
        try:
            repo_check = subprocess.run(
                ["git", "rev-parse", "--show-toplevel"],
                cwd=app_dir,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                timeout=10,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            self._show_manual_git_update(app_dir, f"Could not locate Git: {exc}")
            return

        if repo_check.returncode != 0:
            self._show_manual_git_update(
                app_dir,
                "This application is not running from a Git checkout.",
            )
            return

        repo_root = repo_check.stdout.strip() or app_dir
        branch_check = subprocess.run(
            ["git", "symbolic-ref", "HEAD"],
            cwd=repo_root,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=10,
        )
        branch_ref = branch_check.stdout.strip()
        branch_prefix = "refs/heads/"
        branch = (
            branch_ref[len(branch_prefix):]
            if branch_ref.startswith(branch_prefix)
            else ""
        )
        if branch != "main":
            self._show_manual_git_update(
                repo_root,
                f"The current Git branch is '{branch or 'detached HEAD'}', not 'main'.",
            )
            return

        status_check = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=repo_root,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=10,
        )
        if status_check.returncode != 0 or status_check.stdout.strip():
            self._show_manual_git_update(
                repo_root,
                "Local Git changes were found, so the update was not applied automatically.",
            )
            return

        try:
            auto_workspace_path = self._save_auto_workspace()
        except Exception as exc:
            self._show_manual_git_update(
                repo_root,
                f"Could not save the current workspace before updating: {exc}",
            )
            return

        try:
            pull_result = subprocess.run(
                ["git", "pull", "--ff-only", "origin", "main"],
                cwd=repo_root,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                timeout=60,
            )
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            self._remove_auto_workspace(auto_workspace_path)
            self._show_manual_git_update(repo_root, f"The Git update failed: {exc}")
            return

        output = pull_result.stdout.strip()
        if pull_result.returncode == 0:
            try:
                self._restart_after_update(auto_workspace_path, latest_version)
            except (OSError, ValueError, subprocess.SubprocessError) as exc:
                self._remove_auto_workspace(auto_workspace_path)
                messagebox.showerror(
                    "Update Restart Failed",
                    f"Version {latest_version} was downloaded, but the application could not restart automatically.\n\n"
                    f"{exc}\n\nPlease restart the application manually.\n\n{output}",
                    parent=self.root,
                )
        else:
            self._remove_auto_workspace(auto_workspace_path)
            self._show_manual_git_update(
                repo_root,
                f"The Git update failed:\n{output or 'unknown error'}",
            )

    def _save_auto_workspace(self):
        """Save a restart workspace in the system temporary directory."""
        path = os.path.join(
            tempfile.gettempdir(),
            f"cgyro_comparison_autorestore_{os.getpid()}.json",
        )
        try:
            self._write_workspace_file(path, self._build_workspace_payload())
        except Exception:
            self._remove_auto_workspace(path)
            raise
        return path

    @staticmethod
    def _remove_auto_workspace(path):
        """Remove a temporary restart workspace when it is no longer needed."""
        if not path:
            return
        try:
            os.remove(path)
        except (FileNotFoundError, OSError):
            pass

    def _restart_after_update(self, workspace_path, latest_version):
        """Start the updated checkout and close this pre-update process."""
        app_dir = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.abspath(sys.argv[0])
        if not os.path.isfile(script_path):
            script_path = os.path.join(app_dir, "cgyro_comparison.py")
        if not os.path.isfile(script_path):
            raise OSError(f"Could not locate the application entry point: {script_path}")

        env = os.environ.copy()
        env["CGYRO_AUTO_WORKSPACE"] = str(workspace_path)
        env["CGYRO_AUTO_UPDATE_VERSION"] = str(latest_version)
        command = [sys.executable, script_path, *sys.argv[1:]]
        popen_kwargs = {
            "cwd": app_dir,
            "env": env,
            "stdin": subprocess.DEVNULL,
        }
        if os.name == "nt":
            popen_kwargs["creationflags"] = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
        else:
            popen_kwargs["start_new_session"] = True

        subprocess.Popen(command, **popen_kwargs)
        self.status_var.set(f"Updated to {latest_version}; restarting and restoring workspace...")
        self.root.after(250, self.root.destroy)

    def _show_manual_git_update(self, repo_dir, reason):
        """Show a copyable command when an automatic fast-forward is unsafe."""
        command = f'cd "{repo_dir}"\ngit pull --ff-only origin main'
        messagebox.showwarning(
            "Update Not Applied",
            f"{reason}\n\nRun this command manually after checking the repository:\n\n{command}",
            parent=self.root,
        )

    def show_about(self):
        """Show the local version and project page."""
        messagebox.showinfo(
            "About CGYRO Comparison Tool",
            f"CGYRO Comparison Tool\nVersion {APP_VERSION}\n\n{REPOSITORY_URL}",
            parent=self.root,
        )

    def _init_options(self):
        """Initialize all plot-option widgets and bind dynamic callbacks."""
        self.options_frame.columnconfigure(0, weight=1)
        self.options_frame.columnconfigure(1, weight=1)
        self.options_frame.columnconfigure(2, weight=0)
        self.options_frame.columnconfigure(3, weight=0)

        # --- Persistent Options (Time Range) ---
        ttk.Label(self.options_frame, text="Time Start:").grid(row=0, column=0, sticky=tk.W)
        self.t_start_var = tk.StringVar()
        ttk.Entry(self.options_frame, textvariable=self.t_start_var, width=9).grid(row=0, column=1, sticky=tk.W)

        ttk.Label(self.options_frame, text="Time End:").grid(row=1, column=0, sticky=tk.W)
        self.t_end_var = tk.StringVar()
        ttk.Entry(self.options_frame, textvariable=self.t_end_var, width=9).grid(row=1, column=1, sticky=tk.W)

        # --- Global Log Scale Options ---
        self.log_x_var = tk.BooleanVar(value=False)
        self.log_x_check = ttk.Checkbutton(self.options_frame, text="Log X", variable=self.log_x_var)
        self.log_x_check.grid(row=0, column=2, sticky=tk.W, padx=(6, 0))

        self.log_y_var = tk.BooleanVar(value=False)
        self.log_y_check = ttk.Checkbutton(self.options_frame, text="Log Y", variable=self.log_y_var)
        self.log_y_check.grid(row=1, column=2, sticky=tk.W, padx=(6, 0))

        ttk.Button(
            self.options_frame,
            text="Clear",
            width=6,
            command=self.clear_time_range,
        ).grid(row=0, column=3, rowspan=2, padx=(6, 0), sticky=tk.EW)

        # Separator
        ttk.Separator(self.options_frame, orient=tk.HORIZONTAL).grid(row=2, column=0, columnspan=4, sticky="ew", pady=5)

        # --- Dynamic Options ---
        # 1. Divided by ky (Frequency, Growth rate)
        self.norm_ky_var = tk.BooleanVar(value=False)
        self.norm_ky_check = ttk.Checkbutton(self.options_frame, text="Divided by ky", variable=self.norm_ky_var)
        
        # 2. Flux Options
        self.flux_type_var = tk.StringVar(value="Energy")
        self.flux_type_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.flux_type_var,
            values=list(self._FLUX_TYPE_OPTIONS),
            state="readonly",
            width=15,
        )
        self.flux_xaxis_var = tk.StringVar(value="v.s ky")
        self.flux_xaxis_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.flux_xaxis_var,
            values=list(self._FLUX_XAXIS_OPTIONS),
            state="readonly",
            width=18
        )
        self.flux_scan_xparam_label = ttk.Label(self.options_frame, text="2D X Param:")
        self.flux_scan_xparam_var = tk.StringVar(value="")
        self.flux_scan_xparam_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.flux_scan_xparam_var,
            values=[],
            state="readonly",
            width=18,
        )
        
        self.flux_decomp_var = tk.BooleanVar(value=False)
        self.flux_decomp_check = ttk.Checkbutton(self.options_frame, text="Decompose by Field", variable=self.flux_decomp_var)
        self.flux_2d_errorbar_var = tk.BooleanVar(value=False)
        self.flux_2d_errorbar_check = ttk.Checkbutton(
            self.options_frame,
            text="Error bars",
            variable=self.flux_2d_errorbar_var,
        )
        self.flux_norm_real_ion_var = tk.BooleanVar(value=False)
        self.flux_norm_real_ion_check = ttk.Checkbutton(
            self.options_frame,
            text="Normalized by real ion (max-density ion)",
            variable=self.flux_norm_real_ion_var,
        )
        (
            self.flux_formula_frame,
            self.flux_formula_fig,
            self.flux_formula_ax,
            self.flux_formula_canvas,
            self.flux_formula_widget,
            self.flux_formula_vscroll,
            self.flux_formula_hscroll,
        ) = self._create_formula_panel(self.options_frame, figsize=(4.0, 2.4))

        # 3. Fluctuation 1D Options
        self.fluc_field_var = tk.StringVar(value="Phi")
        self.fluc_field_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.fluc_field_var,
            values=list(self._FLUC_FIELD_OPTIONS),
            state="readonly",
            width=10,
        )
        
        self.fluc_xaxis_var = tk.StringVar(value="v.s ky")
        self.fluc_xaxis_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.fluc_xaxis_var,
            values=list(self._FLUC_XAXIS_OPTIONS),
            state="readonly",
            width=15,
        )
        self.fluc_norm_max_var = tk.BooleanVar(value=False)
        self.fluc_norm_max_check = ttk.Checkbutton(
            self.options_frame,
            text="Normalize by max value",
            variable=self.fluc_norm_max_var,
        )
        self.fluc_theta_kx_label = ttk.Label(self.options_frame, text="kx (blank=avg):")
        self.fluc_theta_kx_var = tk.StringVar(value="")
        self.fluc_theta_kx_entry = ttk.Entry(
            self.options_frame,
            textvariable=self.fluc_theta_kx_var,
            width=12,
        )
        self.fluc_theta_ky_label = ttk.Label(self.options_frame, text="ky (blank=avg):")
        self.fluc_theta_ky_var = tk.StringVar(value="")
        self.fluc_theta_ky_entry = ttk.Entry(
            self.options_frame,
            textvariable=self.fluc_theta_ky_var,
            width=12,
        )
        self.fluc_advanced_var = tk.BooleanVar(value=False)
        self.fluc_advanced_check = ttk.Checkbutton(
            self.options_frame,
            text="Advanced: per-case kx/ky",
            variable=self.fluc_advanced_var,
            command=self._on_fluc_advanced_toggle,
        )
        self.fluc_advanced_button = ttk.Button(
            self.options_frame,
            text="Edit per-case kx/ky...",
            command=self._open_advanced_fluc_selector,
        )
        (
            self.fluc_formula_frame,
            self.fluc_formula_fig,
            self.fluc_formula_ax,
            self.fluc_formula_canvas,
            self.fluc_formula_widget,
            self.fluc_formula_vscroll,
            self.fluc_formula_hscroll,
        ) = self._create_formula_panel(self.options_frame, figsize=(4.0, 2.2))

        # 4. Species Selection (Flux, Fluctuation 2D)
        self.species_label = ttk.Label(self.options_frame, text="Species:")
        self.species_var = tk.StringVar()
        self.species_combo = ttk.Combobox(self.options_frame, textvariable=self.species_var, width=20, state="readonly")
        
        self.plot_all_species_var = tk.BooleanVar(value=False)
        self.plot_all_species_check = ttk.Checkbutton(self.options_frame, text="Plot All Species (First Case Only)", variable=self.plot_all_species_var)

        # 5. Fluctuation 2D options
        self.moment_label = ttk.Label(self.options_frame, text="Moment:")
        self.moment_var = tk.StringVar(value="Phi")
        self.moment_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.moment_var,
            values=list(self._FLUC_MOMENT_OPTIONS),
            state="readonly",
        )

        self.fluc2d_view_label = ttk.Label(self.options_frame, text="View:")
        self.fluc2d_view_var = tk.StringVar(value="vs xy")
        self.fluc2d_view_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.fluc2d_view_var,
            values=list(self._FLUC2D_VIEW_OPTIONS),
            state="readonly",
            width=10
        )
        self.fluc2d_x_elec_var = tk.BooleanVar(value=False)
        self.fluc2d_x_elec_check = ttk.Checkbutton(
            self.options_frame,
            text=r"Spatial axes normalize to electron scale (x,y)/\rho_e",
            variable=self.fluc2d_x_elec_var,
        )
        self.fluc2d_log_z_var = tk.BooleanVar(value=False)
        self.fluc2d_log_z_check = ttk.Checkbutton(
            self.options_frame,
            text="Log z",
            variable=self.fluc2d_log_z_var,
        )

        # 6. FFT Options (reused from before, but managed dynamically)
        self.fft_options_frame = ttk.LabelFrame(self.options_frame, text="FFT Settings", padding=5)
        # Shared linear spectrum file path (used by both FFT overlay and Zonal-vs-Gamma_Linear)
        self.linear_gamma_file_var = tk.StringVar(value=DEFAULT_LINEAR_GAMMA_FILE)
        # Analysis Mode
        ttk.Label(self.fft_options_frame, text="Mode:").grid(row=0, column=0, sticky=tk.W)
        self.fft_mode_var = tk.StringVar(value="Nonlinear")
        self.fft_mode_combo = ttk.Combobox(
            self.fft_options_frame,
            textvariable=self.fft_mode_var,
            values=list(self._FFT_MODE_OPTIONS),
            state="readonly",
            width=15,
        )
        self.fft_mode_combo.grid(row=0, column=1, sticky=tk.W)
        # View Mode
        ttk.Label(self.fft_options_frame, text="View:").grid(row=1, column=0, sticky=tk.W)
        self.fft_view_var = tk.StringVar(value="Omega vs ky")
        self.fft_view_combo = ttk.Combobox(
            self.fft_options_frame,
            textvariable=self.fft_view_var,
            values=list(self._FFT_VIEW_OPTIONS),
            state="readonly",
            width=15,
        )
        self.fft_view_combo.grid(row=1, column=1, sticky=tk.W)
        # Spectrum type
        ttk.Label(self.fft_options_frame, text="Spectrum:").grid(row=2, column=0, sticky=tk.W)
        self.fft_spectrum_var = tk.StringVar(value="Amplitude")
        self.fft_spectrum_combo = ttk.Combobox(
            self.fft_options_frame,
            textvariable=self.fft_spectrum_var,
            values=list(self._FFT_SPECTRUM_OPTIONS),
            state="readonly",
            width=15,
        )
        self.fft_spectrum_combo.grid(row=2, column=1, sticky=tk.W)
        # Overlay Frequency Checkbox (from exported omega/gamma-vs-ky file)
        self.fft_overlay_freq_var = tk.BooleanVar(value=False)
        self.fft_overlay_freq_check = ttk.Checkbutton(
            self.fft_options_frame,
            text="Overlay Linear Freq",
            variable=self.fft_overlay_freq_var,
        )
        self.fft_overlay_freq_check.grid(row=3, column=0, columnspan=2, sticky=tk.W)
        ttk.Label(self.fft_options_frame, text="Linear File:").grid(row=4, column=0, sticky=tk.W)
        self.fft_linear_file_entry = ttk.Entry(
            self.fft_options_frame,
            textvariable=self.linear_gamma_file_var,
            width=20,
        )
        self.fft_linear_file_entry.grid(row=4, column=1, sticky=tk.W)
        self.fft_linear_file_browse = ttk.Button(
            self.fft_options_frame,
            text="Browse",
            command=self._browse_linear_gamma_file,
            width=8,
        )
        self.fft_linear_file_browse.grid(row=5, column=1, sticky=tk.W, pady=(2, 0))
        (
            self.fft_formula_frame,
            self.fft_formula_fig,
            self.fft_formula_ax,
            self.fft_formula_canvas,
            self.fft_formula_widget,
            self.fft_formula_vscroll,
            self.fft_formula_hscroll,
        ) = self._create_formula_panel(self.options_frame, figsize=(4.0, 2.2))

        # 7. Zonal ExB shearing options
        self.zf_xaxis_label = ttk.Label(self.options_frame, text="Mode:")
        self.zf_xaxis_var = tk.StringVar(value="vs Time")
        self.zf_xaxis_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.zf_xaxis_var,
            values=list(self._ZF_XAXIS_OPTIONS),
            state="readonly",
            width=15
        )
        self.linear_gamma_file_label = ttk.Label(self.options_frame, text="Linear File:")
        self.linear_gamma_file_entry = ttk.Entry(
            self.options_frame,
            textvariable=self.linear_gamma_file_var,
            width=30
        )
        self.linear_gamma_file_browse = ttk.Button(
            self.options_frame,
            text="Browse",
            command=self._browse_linear_gamma_file
        )
        self.zf_gamma_lin_ky_label = ttk.Label(self.options_frame, text="ky for ratio:")
        self.zf_gamma_lin_ky_var = tk.StringVar(value="")
        self.zf_gamma_lin_ky_entry = ttk.Entry(
            self.options_frame,
            textvariable=self.zf_gamma_lin_ky_var,
            width=12,
        )
        (
            self.zf_formula_frame,
            self.zf_formula_fig,
            self.zf_formula_ax,
            self.zf_formula_canvas,
            self.zf_formula_widget,
            self.zf_formula_vscroll,
            self.zf_formula_hscroll,
        ) = self._create_formula_panel(self.options_frame, figsize=(4.0, 2.2))

        # 8. Energy-balance options (triad_v2-based)
        self.energy_balance_options_frame = ttk.Frame(self.options_frame)
        self.energy_balance_options_frame.columnconfigure(0, weight=0)
        self.energy_balance_options_frame.columnconfigure(1, weight=1)
        self.energy_balance_options_frame.columnconfigure(2, weight=0)
        self.energy_balance_options_frame.columnconfigure(3, weight=1)

        self.energy_balance_mode_label = ttk.Label(self.energy_balance_options_frame, text="Mode:")
        self.energy_balance_mode_var = tk.StringVar(value="Entropy balance")
        self.energy_balance_mode_combo = ttk.Combobox(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_mode_var,
            values=list(self._ENERGY_BALANCE_MODE_OPTIONS),
            state="readonly",
            width=24,
        )
        self.energy_balance_n_label = ttk.Label(self.energy_balance_options_frame, text="n index:")
        self.energy_balance_n_var = tk.StringVar(value="0")
        self.energy_balance_n_entry = ttk.Entry(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_n_var,
            width=10,
        )
        self.energy_balance_spec_label = ttk.Label(self.energy_balance_options_frame, text="Species:")
        self.energy_balance_spec_var = tk.StringVar(value="Total (-1)")
        self.energy_balance_spec_combo = ttk.Combobox(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_spec_var,
            values=list(self._ENERGY_BALANCE_SPEC_OPTIONS),
            state="readonly",
            width=16,
        )
        self.energy_balance_single_quantity_label = ttk.Label(self.energy_balance_options_frame, text="Quantity:")
        self.energy_balance_single_quantity_var = tk.StringVar(value="T-N")
        self.energy_balance_single_quantity_combo = ttk.Combobox(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_single_quantity_var,
            values=list(self._ENERGY_BALANCE_SINGLE_QUANTITY_OPTIONS),
            state="readonly",
            width=12,
        )
        self.energy_balance_single_xaxis_label = ttk.Label(self.energy_balance_options_frame, text="X-axis:")
        self.energy_balance_single_xaxis_var = tk.StringVar(value="vs Time")
        self.energy_balance_single_xaxis_combo = ttk.Combobox(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_single_xaxis_var,
            values=list(self._ENERGY_BALANCE_SINGLE_XAXIS_OPTIONS),
            state="readonly",
            width=12,
        )
        self.energy_balance_single_norm_label = ttk.Label(self.energy_balance_options_frame, text="Normalize:")
        self.energy_balance_single_norm_var = tk.StringVar(value="None")
        self.energy_balance_single_norm_combo = ttk.Combobox(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_single_norm_var,
            values=list(self._ENERGY_BALANCE_SINGLE_NORM_OPTIONS),
            state="readonly",
            width=12,
        )
        # Kept only so older workspace JSON files that saved this BooleanVar
        # still restore cleanly.  New UI state uses energy_balance_single_norm_var.
        self.energy_balance_single_norm_entropy_var = tk.BooleanVar(value=False)
        self.energy_balance_transfer_quantity_label = ttk.Label(self.energy_balance_options_frame, text="Quantity:")
        self.energy_balance_transfer_quantity_var = tk.StringVar(value="Re")
        self.energy_balance_transfer_quantity_combo = ttk.Combobox(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_transfer_quantity_var,
            values=list(self._ENERGY_BALANCE_TRANSFER_QUANTITY_OPTIONS),
            state="readonly",
            width=10,
        )
        self.energy_balance_transfer_xaxis_label = ttk.Label(self.energy_balance_options_frame, text="X-axis:")
        self.energy_balance_transfer_xaxis_var = tk.StringVar(value="vs kxky")
        self.energy_balance_transfer_xaxis_combo = ttk.Combobox(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_transfer_xaxis_var,
            values=list(self._ENERGY_BALANCE_TRANSFER_XAXIS_OPTIONS),
            state="readonly",
            width=10,
        )
        self.energy_balance_transfer_ky_label = ttk.Label(self.energy_balance_options_frame, text="fixed source ky:")
        self.energy_balance_transfer_ky_var = tk.StringVar(value="0")
        self.energy_balance_transfer_ky_entry = ttk.Entry(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_transfer_ky_var,
            width=10,
        )
        self.energy_balance_transfer_kx_label = ttk.Label(self.energy_balance_options_frame, text="fixed source kx:")
        self.energy_balance_transfer_kx_var = tk.StringVar(value="0")
        self.energy_balance_transfer_kx_entry = ttk.Entry(
            self.energy_balance_options_frame,
            textvariable=self.energy_balance_transfer_kx_var,
            width=10,
        )
        self.energy_balance_transfer_asym_var = tk.BooleanVar(value=False)
        self.energy_balance_transfer_asym_check = ttk.Checkbutton(
            self.energy_balance_options_frame,
            text="Asym",
            variable=self.energy_balance_transfer_asym_var,
            command=self.update_options,
        )
        self.energy_balance_transfer_norm_max_var = tk.BooleanVar(value=False)
        self.energy_balance_transfer_norm_max_check = ttk.Checkbutton(
            self.energy_balance_options_frame,
            text="Normalized by max T",
            variable=self.energy_balance_transfer_norm_max_var,
            command=self.update_options,
        )
        (
            self.energy_balance_formula_frame,
            self.energy_balance_formula_fig,
            self.energy_balance_formula_ax,
            self.energy_balance_formula_canvas,
            self.energy_balance_formula_widget,
            self.energy_balance_formula_vscroll,
            self.energy_balance_formula_hscroll,
        ) = self._create_formula_panel(self.options_frame, figsize=(4.4, 4.0))

        # 9. Others options
        self.others_plot_label = ttk.Label(self.options_frame, text="Others:")
        self.others_plot_var = tk.StringVar(value="Error")
        self.others_plot_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.others_plot_var,
            values=list(self._OTHERS_PLOT_OPTIONS),
            state="readonly",
            width=15,
        )
        self.others_rcorr_field_label = ttk.Label(self.options_frame, text="Field:")
        self.others_rcorr_field_var = tk.StringVar(value="Phi")
        self.others_rcorr_field_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.others_rcorr_field_var,
            values=list(self._OTHERS_FIELD_OPTIONS),
            state="readonly",
            width=10,
        )
        self.others_rcorr_theta_label = ttk.Label(self.options_frame, text="Theta Idx:")
        self.others_rcorr_theta_var = tk.StringVar(value="-1")
        self.others_rcorr_theta_entry = ttk.Entry(
            self.options_frame,
            textvariable=self.others_rcorr_theta_var,
            width=10,
        )
        self.others_pod_field_label = ttk.Label(self.options_frame, text="Field:")
        self.others_pod_field_var = tk.StringVar(value="Apar")
        self.others_pod_field_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.others_pod_field_var,
            values=list(self._OTHERS_POD_FIELD_OPTIONS),
            state="readonly",
            width=8,
        )
        self.others_pod_kx_label = ttk.Label(self.options_frame, text="kx:")
        self.others_pod_kx_var = tk.StringVar(value="0")
        self.others_pod_kx_entry = ttk.Entry(
            self.options_frame,
            textvariable=self.others_pod_kx_var,
            width=12,
        )
        self.others_pod_ky_label = ttk.Label(self.options_frame, text="ky:")
        self.others_pod_ky_var = tk.StringVar(value="0")
        self.others_pod_ky_entry = ttk.Entry(
            self.options_frame,
            textvariable=self.others_pod_ky_var,
            width=12,
        )
        (
            self.pod_formula_frame,
            self.pod_formula_fig,
            self.pod_formula_ax,
            self.pod_formula_canvas,
            self.pod_formula_widget,
            self.pod_formula_vscroll,
            self.pod_formula_hscroll,
        ) = self._create_formula_panel(self.options_frame, figsize=(4.0, 2.4))

        # Bind once; avoid repeated callback registration inside update_options.
        self.flux_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.flux_type_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.fluc_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.fluc2d_view_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.moment_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.zf_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.energy_balance_mode_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.energy_balance_single_quantity_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.energy_balance_single_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.energy_balance_single_norm_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.energy_balance_transfer_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.fft_mode_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.fft_view_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.fft_spectrum_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.others_plot_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.others_rcorr_field_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.others_pod_field_combo.bind("<<ComboboxSelected>>", self.update_options)

        # Initialize layout
        self.update_options()

    def clear_time_range(self):
        """Clear GUI time-range entries and revert to automatic time-window logic."""
        self.time_mode_var.set("Manual Start/End")
        self.time_percent_var.set("50")
        self.time_duration_var.set("")
        self.t_start_var.set("")
        self.t_end_var.set("")

    def clear_axis_limits(self):
        """Clear manual plot-axis limit entries."""
        self.axis_kx_min_var.set("")
        self.axis_kx_max_var.set("")
        self.axis_ky_min_var.set("")
        self.axis_ky_max_var.set("")
        try:
            self._apply_manual_axis_limits()
            self.canvas.draw_idle()
        except Exception:
            pass

    def open_axis_limits_dialog(self):
        """Open a top-menu Axis dialog for persistent manual plot limits."""
        if hasattr(self, "_axis_dialog") and self._axis_dialog.winfo_exists():
            self._axis_dialog.lift()
            self._axis_dialog.focus_set()
            return

        dialog = tk.Toplevel(self.root)
        self._axis_dialog = dialog
        dialog.title("Axis")
        dialog.transient(self.root)
        dialog.resizable(False, False)

        frame = ttk.Frame(dialog, padding=10)
        frame.grid(row=0, column=0, sticky=tk.NSEW)
        for col in (1, 2):
            frame.columnconfigure(col, weight=1)

        ttk.Label(frame, text="Min").grid(row=0, column=1, sticky=tk.W, padx=(0, 6))
        ttk.Label(frame, text="Max").grid(row=0, column=2, sticky=tk.W)

        ttk.Label(frame, text="x lim:").grid(row=1, column=0, sticky=tk.W, padx=(0, 8), pady=(4, 0))
        kx_min_entry = ttk.Entry(frame, textvariable=self.axis_kx_min_var, width=12)
        kx_min_entry.grid(row=1, column=1, sticky=tk.EW, padx=(0, 6), pady=(4, 0))
        ttk.Entry(frame, textvariable=self.axis_kx_max_var, width=12).grid(
            row=1, column=2, sticky=tk.EW, pady=(4, 0)
        )

        ttk.Label(frame, text="y lim:").grid(row=2, column=0, sticky=tk.W, padx=(0, 8), pady=(4, 0))
        ttk.Entry(frame, textvariable=self.axis_ky_min_var, width=12).grid(
            row=2, column=1, sticky=tk.EW, padx=(0, 6), pady=(4, 0)
        )
        ttk.Entry(frame, textvariable=self.axis_ky_max_var, width=12).grid(
            row=2, column=2, sticky=tk.EW, pady=(4, 0)
        )

        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=3, column=0, columnspan=3, sticky=tk.E, pady=(10, 0))

        def apply_limits():
            try:
                self._manual_axis_limits()
            except Exception as e:
                messagebox.showerror("Axis", f"Invalid axis limit: {e}", parent=dialog)
                return
            try:
                self._apply_manual_axis_limits()
                self.canvas.draw_idle()
            except Exception:
                pass

        ttk.Button(btn_frame, text="Apply", command=apply_limits).pack(side=tk.LEFT, padx=(0, 6))
        ttk.Button(btn_frame, text="Clear", command=self.clear_axis_limits).pack(side=tk.LEFT, padx=(0, 6))
        ttk.Button(btn_frame, text="Close", command=dialog.destroy).pack(side=tk.LEFT)

        dialog.bind("<Return>", lambda _event: apply_limits())
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        kx_min_entry.focus_set()

    def _clear_simple_time_entries(self):
        """Clear left-panel time entries so the top Time menu controls averaging."""
        self.t_start_var.set("")
        self.t_end_var.set("")

    def _set_time_last_percent(self):
        """Prompt for the final time percentage used by the top Time menu."""
        try:
            current = float(self.time_percent_var.get())
        except Exception:
            current = 50.0
        value = simpledialog.askfloat(
            "Time average",
            "Average over last percent of available time:",
            parent=self.root,
            initialvalue=current,
            minvalue=0.0,
            maxvalue=100.0,
        )
        if value is None:
            return
        self.time_percent_var.set(f"{value:g}")
        self.time_mode_var.set("Last percent")
        self._clear_simple_time_entries()

    def _set_time_last_duration(self):
        """Prompt for the final physical-time duration used by the top Time menu."""
        try:
            current = float(self.time_duration_var.get())
        except Exception:
            current = 0.0
        value = simpledialog.askfloat(
            "Time average",
            "Average over last time duration:",
            parent=self.root,
            initialvalue=current if current > 0.0 else 1.0,
            minvalue=0.0,
        )
        if value is None:
            return
        self.time_duration_var.set(f"{value:g}")
        self.time_mode_var.set("Last duration")
        self._clear_simple_time_entries()

    def _browse_linear_gamma_file(self):
        """Open file chooser and update the linear omega/gamma reference file path."""
        initial_dir = default_share_dir()
        try:
            current = self.linear_gamma_file_var.get().strip()
            if current:
                if os.path.isabs(current):
                    initial_dir = os.path.dirname(current) or initial_dir
                elif os.path.isdir(current):
                    initial_dir = current
        except Exception:
            pass

        path = filedialog.askopenfilename(
            title="Select omega/gamma vs ky file",
            initialdir=initial_dir,
            filetypes=[("Text files", "*.txt *.dat *.csv"), ("All files", "*.*")]
        )
        if path:
            self.linear_gamma_file_var.set(path)

    def _fluc_advanced_enabled(self):
        """Return whether Fluctuation 1D uses per-case kx/ky selections."""
        var = getattr(self, "fluc_advanced_var", None)
        try:
            return bool(var.get()) if var is not None else False
        except Exception:
            return False

    def _ensure_fluc_advanced_case_values(self):
        """Keep the Advanced kx/ky table synchronized with loaded case names."""
        values = getattr(self, "_fluc_advanced_case_values", None)
        if not isinstance(values, dict):
            values = {}
            self._fluc_advanced_case_values = values

        case_names = []
        try:
            case_names = [
                str(self.case_listbox.get(index))
                for index in range(self.case_listbox.size())
            ]
        except (AttributeError, tk.TclError):
            pass

        for case_name in case_names:
            current = values.get(case_name)
            if not isinstance(current, dict):
                current = {}
            values[case_name] = {
                "kx": str(current.get("kx", "")),
                "ky": str(current.get("ky", "")),
            }

        for case_name in list(values):
            if case_name not in case_names:
                values.pop(case_name, None)
        return values

    def _seed_fluc_advanced_case_values(self):
        """Initialize new per-case entries from the shared selector values."""
        values = self._ensure_fluc_advanced_case_values()
        try:
            default_kx = str(self.fluc_theta_kx_var.get()).strip()
        except Exception:
            default_kx = ""
        try:
            default_ky = str(self.fluc_theta_ky_var.get()).strip()
        except Exception:
            default_ky = ""

        for entry in values.values():
            if not str(entry.get("kx", "")).strip():
                entry["kx"] = default_kx
            if not str(entry.get("ky", "")).strip():
                entry["ky"] = default_ky

    def _on_fluc_advanced_toggle(self):
        """Switch between shared and per-case Fluctuation 1D selectors."""
        if self._fluc_advanced_enabled() and not getattr(self, "_fluc_advanced_seeded", False):
            self._seed_fluc_advanced_case_values()
            self._fluc_advanced_seeded = True
        self.update_options()

    def _get_fluc_case_axis_text(self, case_label, axis_name):
        """Return the effective kx/ky text for one case."""
        if not self._fluc_advanced_enabled():
            var_name = "fluc_theta_kx_var" if axis_name == "kx" else "fluc_theta_ky_var"
            var = getattr(self, var_name, None)
            try:
                return str(var.get()).strip() if var is not None else ""
            except Exception:
                return ""

        values = self._ensure_fluc_advanced_case_values()
        entry = values.get(str(case_label), {})
        if not isinstance(entry, dict):
            return ""
        return str(entry.get(axis_name, "")).strip()

    def _open_advanced_fluc_selector(self):
        """Open the per-case kx/ky editor used by Fluctuation 1D Advanced mode."""
        try:
            case_names = [
                str(self.case_listbox.get(index))
                for index in range(self.case_listbox.size())
            ]
        except (AttributeError, tk.TclError):
            case_names = []
        if not case_names:
            messagebox.showwarning(
                "Advanced kx/ky",
                "Load at least one case before editing per-case selections.",
                parent=self.root,
            )
            return

        values = self._ensure_fluc_advanced_case_values()
        dialog = tk.Toplevel(self.root)
        dialog.title("Advanced kx/ky selection")
        dialog.transient(self.root)
        dialog.geometry("680x460")
        dialog.minsize(520, 300)

        content = ttk.Frame(dialog, padding=10)
        content.pack(fill=tk.BOTH, expand=True)
        content.columnconfigure(0, weight=1)
        content.rowconfigure(1, weight=1)

        ttk.Label(
            content,
            text=(
                "Choose a value for each case. Leave a cell blank to average "
                "over that axis; use idx:n or val:x when needed."
            ),
            wraplength=640,
            justify=tk.LEFT,
        ).grid(row=0, column=0, sticky=tk.W, pady=(0, 8))

        table_frame = ttk.Frame(content)
        table_frame.grid(row=1, column=0, sticky=tk.NSEW)
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(0, weight=1)

        canvas = tk.Canvas(table_frame, highlightthickness=0, borderwidth=0)
        scrollbar = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=canvas.yview)
        inner = ttk.Frame(canvas)
        window_id = canvas.create_window((0, 0), window=inner, anchor=tk.NW)
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.grid(row=0, column=0, sticky=tk.NSEW)
        scrollbar.grid(row=0, column=1, sticky=tk.NS)

        inner.columnconfigure(0, weight=1)
        inner.bind(
            "<Configure>",
            lambda _event: canvas.configure(scrollregion=canvas.bbox("all")),
        )
        canvas.bind(
            "<Configure>",
            lambda event: canvas.itemconfigure(window_id, width=event.width),
        )

        ttk.Label(inner, text="Case").grid(row=0, column=0, sticky=tk.W, padx=(2, 12), pady=(0, 4))
        ttk.Label(inner, text="kx (blank=avg)").grid(row=0, column=1, sticky=tk.W, padx=4, pady=(0, 4))
        ttk.Label(inner, text="ky (blank=avg)").grid(row=0, column=2, sticky=tk.W, padx=4, pady=(0, 4))

        editors = []
        for row, case_name in enumerate(case_names, start=1):
            entry = values.get(case_name, {})
            kx_var = tk.StringVar(value=str(entry.get("kx", "")))
            ky_var = tk.StringVar(value=str(entry.get("ky", "")))
            ttk.Label(inner, text=case_name).grid(
                row=row, column=0, sticky=tk.W, padx=(2, 12), pady=3
            )
            ttk.Entry(inner, textvariable=kx_var, width=18).grid(
                row=row, column=1, sticky=tk.EW, padx=4, pady=3
            )
            ttk.Entry(inner, textvariable=ky_var, width=18).grid(
                row=row, column=2, sticky=tk.EW, padx=4, pady=3
            )
            editors.append((case_name, kx_var, ky_var))
        inner.columnconfigure(1, weight=1)
        inner.columnconfigure(2, weight=1)

        buttons = ttk.Frame(content)
        buttons.grid(row=2, column=0, sticky=tk.EW, pady=(10, 0))

        def on_apply(event=None):
            for case_name, kx_var, ky_var in editors:
                values[case_name] = {
                    "kx": kx_var.get().strip(),
                    "ky": ky_var.get().strip(),
                }
            dialog.destroy()
            self.update_options()

        def on_cancel(event=None):
            dialog.destroy()

        ttk.Button(buttons, text="Cancel", command=on_cancel).pack(side=tk.RIGHT)
        ttk.Button(buttons, text="Apply", command=on_apply).pack(
            side=tk.RIGHT, padx=(0, 8)
        )
        dialog.bind("<Return>", on_apply)
        dialog.bind("<Escape>", on_cancel)
        dialog.grab_set()
        dialog.focus_set()
        dialog.wait_window()

    def _hide_dynamic_options(self):
        """Hide all dynamic option widgets before re-laying out current mode controls."""
        widgets = [
            self.norm_ky_check,
            self.flux_type_combo, self.flux_xaxis_combo, self.flux_decomp_check,
            self.flux_2d_errorbar_check, self.flux_norm_real_ion_check,
            self.flux_scan_xparam_label, self.flux_scan_xparam_combo,
            self.flux_formula_frame,
            self.fluc_field_combo, self.fluc_xaxis_combo,
            self.fluc_norm_max_check,
            self.fluc_advanced_check, self.fluc_advanced_button,
            self.fluc_theta_kx_label, self.fluc_theta_kx_entry,
            self.fluc_theta_ky_label, self.fluc_theta_ky_entry,
            self.fluc_formula_frame,
            self.species_label, self.species_combo, self.plot_all_species_check,
            self.fluc2d_view_label, self.fluc2d_view_combo,
            self.fluc2d_x_elec_check, self.fluc2d_log_z_check,
            self.moment_label, self.moment_combo,
            self.fft_options_frame,
            self.fft_formula_frame,
            self.zf_xaxis_label, self.zf_xaxis_combo,
            self.linear_gamma_file_label, self.linear_gamma_file_entry, self.linear_gamma_file_browse,
            self.zf_gamma_lin_ky_label, self.zf_gamma_lin_ky_entry,
            self.zf_formula_frame,
            self.energy_balance_options_frame,
            self.energy_balance_mode_label, self.energy_balance_mode_combo,
            self.energy_balance_n_label, self.energy_balance_n_entry,
            self.energy_balance_spec_label, self.energy_balance_spec_combo,
            self.energy_balance_single_quantity_label, self.energy_balance_single_quantity_combo,
            self.energy_balance_single_xaxis_label, self.energy_balance_single_xaxis_combo,
            self.energy_balance_single_norm_label, self.energy_balance_single_norm_combo,
            self.energy_balance_transfer_quantity_label, self.energy_balance_transfer_quantity_combo,
            self.energy_balance_transfer_xaxis_label, self.energy_balance_transfer_xaxis_combo,
            self.energy_balance_transfer_ky_label, self.energy_balance_transfer_ky_entry,
            self.energy_balance_transfer_kx_label, self.energy_balance_transfer_kx_entry,
            self.energy_balance_transfer_asym_check,
            self.energy_balance_transfer_norm_max_check,
            self.energy_balance_formula_frame,
            self.others_plot_label, self.others_plot_combo,
            self.others_rcorr_field_label, self.others_rcorr_field_combo,
            self.others_rcorr_theta_label, self.others_rcorr_theta_entry,
            self.others_pod_field_label, self.others_pod_field_combo,
            self.others_pod_ky_label, self.others_pod_ky_entry,
            self.others_pod_kx_label, self.others_pod_kx_entry,
            self.pod_formula_frame,
        ]
        for w in widgets:
            w.grid_remove()

    @staticmethod
    def _safe_mathtext_line(line):
        """Normalize formula-panel strings for older Matplotlib mathtext parsers."""
        text = str(line)
        text = re.sub(r"\\ge(?![A-Za-z])", r"\\geq", text)
        text = re.sub(r"\\le(?![A-Za-z])", r"\\leq", text)
        return text

    @staticmethod
    def _formula_visual_length(line):
        """Rough display-length estimate used to shrink long formula notes."""
        text = str(line)
        text = re.sub(r"\$+", "", text)
        text = re.sub(r"\\mathrm\{([^{}]*)\}", r"\1", text)
        text = re.sub(r"\\bf\{([^{}]*)\}", r"\1", text)
        text = re.sub(r"\\[A-Za-z]+", "X", text)
        text = re.sub(r"[{}_^]", "", text)
        return max(1, len(text))

    def _formula_font_size(self, line, base_fs, width_chars=62, min_fs=5.8):
        """Scale formula-panel font size so long lines are not clipped."""
        visual_len = self._formula_visual_length(line)
        if visual_len <= width_chars:
            return base_fs
        return max(min_fs, min(base_fs, base_fs * float(width_chars) / float(visual_len)))

    def _draw_formula_panel(
        self,
        fig,
        ax,
        canvas,
        lines,
        *,
        widget=None,
        base_fig_w=4.4,
        min_fig_h=2.0,
        title_fs=8.8,
        line_fs=7.8,
        title_dy=0.115,
        line_dy=0.105,
        frac_fs=7.2,
        frac_dy=0.125,
        width_chars=62,
        min_fs=5.8,
    ):
        """Render formula notes in a wide, unclipped axes with sane scroll extents."""
        def line_style(index, text):
            if index == 0:
                return title_fs, title_dy
            if (
                r"\dfrac" in text
                or r"\frac" in text
                or r"\sqrt" in text
                or r"\sum" in text
                or r"\left" in text
                or r"\langle" in text
            ):
                return frac_fs, frac_dy
            return line_fs, line_dy

        # Work in pixels first. Mathtext rows with sqrt/sum/fractions are much
        # taller than ordinary text; estimating in axes-fraction units causes
        # overlap when the Tk widget is narrow.
        dpi = float(fig.dpi)
        top_px = 24.0
        bottom_px = 24.0
        row_specs = []
        total_px = top_px + bottom_px
        max_visual_len = 1
        for idx, text in enumerate(lines):
            fs_base, dy_factor = line_style(idx, text)
            fs = self._formula_font_size(text, fs_base, width_chars=width_chars, min_fs=min_fs)
            visual_len = self._formula_visual_length(text)
            max_visual_len = max(max_visual_len, visual_len)
            fs_px = float(fs) * dpi / 72.0
            if idx == 0:
                row_px = max(24.0, fs_px * 2.4)
            else:
                row_px = max(20.0, fs_px * 2.0)
            if any(token in text for token in (r"\sqrt", r"\dfrac", r"\frac")):
                row_px = max(row_px, fs_px * 6.2)
            if r"\sum" in text:
                row_px = max(row_px, fs_px * 5.4)
            if any(token in text for token in (r"\left", r"\langle")):
                row_px = max(row_px, fs_px * 3.1)
            row_specs.append((text, fs, row_px))
            total_px += row_px

        fig_h = max(min_fig_h, total_px / dpi)
        fig_w = max(base_fig_w, min(9.0, 0.9 + 0.047 * float(max_visual_len)))
        try:
            fig.set_size_inches(fig_w, fig_h, forward=True)
            ax.set_position([0.05, 0.04, 0.92, 0.93])
        except Exception:
            pass

        ax.clear()
        ax.axis("off")
        current_px = top_px
        fig_h_px_for_text = max(1.0, fig_h * dpi)
        for line, fs, row_px in row_specs:
            y = 1.0 - current_px / fig_h_px_for_text
            ax.text(
                0.01,
                y,
                self._safe_mathtext_line(line),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=fs,
                clip_on=False,
            )
            current_px += row_px

        try:
            pad_px = 36
            fig_w_px = int(fig.get_figwidth() * fig.dpi)
            fig_h_px = int(fig.get_figheight() * fig.dpi)
            target = widget if widget is not None else canvas.get_tk_widget()
            inner = getattr(target, "_formula_inner_widget", None)
            window_id = getattr(target, "_formula_window_id", None)
            if inner is not None:
                inner.configure(width=fig_w_px, height=fig_h_px)
            if window_id is not None and hasattr(target, "itemconfigure"):
                target.itemconfigure(window_id, width=fig_w_px, height=fig_h_px)
            if hasattr(target, "configure"):
                target.configure(scrollregion=(0, 0, fig_w_px + pad_px, fig_h_px + pad_px))
            if hasattr(target, "xview_moveto"):
                target.xview_moveto(0.0)
            if hasattr(target, "yview_moveto"):
                target.yview_moveto(0.0)
        except Exception:
            pass

        canvas.draw_idle()

    def _render_flux_kx_formula_math(self):
        """Render math notes for the estimated `Flux vs kx` mode."""
        flux_type = self.flux_type_var.get().strip().lower()
        if flux_type == "particle":
            moment_note = r"Particle mode: $m_{\mathrm{ES}}=n,\;m_A=v$ (flutter proxy), $m_B=n$."
        else:
            moment_note = r"Energy mode: $m_{\mathrm{ES}}=e,\;m_A=e,\;m_B=e$ (estimated EM proxies)."
        norm_note = (
            r"Optional: real-ion normalization uses ion with max $n_i$ ($Z>0$) "
            r"as GyroBohm reference."
        )

        lines = [
            r"Flux vs $k_x$ (estimated, follows code path)",
            r"Input fields loaded as $[k_x,\theta,k_y,t]$ (species fields first summed over selected species).",
            r"$\Gamma_{\mathrm{ES}}(k_x,t)=\left\langle \sum_{k_y}"
            r"\Re\!\left[m_{\mathrm{ES}}^*(i k_y)\!\left(-\phi\right)\right]\right\rangle_{\theta}$",
            r"$\Gamma_{\mathrm{EM},A}(k_x,t)=\left\langle \sum_{k_y}"
            r"\Re\!\left[m_A^*(i k_y)\!\left(A_{\parallel}\right)\right]\right\rangle_{\theta}$",
            r"$\Gamma_{\mathrm{EM},B}(k_x,t)=\left\langle \sum_{k_y}"
            r"\Re\!\left[m_B^*(i k_y)B_{\parallel}\right]\right\rangle_{\theta}$",
            r"$\Gamma_{\mathrm{tot}}(k_x,t)=\Gamma_{\mathrm{ES}}+\Gamma_{\mathrm{EM},A}+\Gamma_{\mathrm{EM},B}$",
            r"Time average in code: $\bar{\Gamma}(k_x)=\langle\Gamma_{\mathrm{tot}}(k_x,t)\rangle_{t\in[t_0,t_1]}$;",
            r"if no valid time window, uses last time slice $t=t_{\mathrm{end}}$.",
            r"Final display keeps non-negative branch only: $k_x \geq 0$.",
            moment_note,
            norm_note,
        ]

        self._draw_formula_panel(
            self.flux_formula_fig,
            self.flux_formula_ax,
            self.flux_formula_canvas,
            lines,
            widget=self.flux_formula_widget,
            base_fig_w=3.8,
            min_fig_h=2.4,
            title_fs=8.8,
            line_fs=7.6,
            line_dy=0.08,
            width_chars=56,
        )

    def _render_fluctuation_1d_formula_math(self):
        """Render formula notes for Fluctuation 1D averaging definition."""
        mode = str(self.fluc_average_var.get()).strip().lower()
        xaxis = str(self.fluc_xaxis_var.get()).strip().lower()
        if xaxis == "v.s theta":
            if mode == "mean absolute":
                lines = [
                    r"Phi vs theta: Mean Absolute (code definition)",
                    r"Field slice: $F(k_x,\theta,k_y,t)=\phi/\rho_s$.",
                    r"Blank kx or ky selection averages over that spectral axis.",
                    r"Use $idx:n$ for an index or $val:x$ for the nearest physical kx/ky value.",
                    r"Optional max normalization divides the final profile by its finite maximum.",
                    r"$A(\theta)=\left\langle\mathrm{mean}_{k_x,k_y}|F(k_x,\theta,k_y,t)|\right\rangle_t$.",
                    r"The time window follows the shared Time controls.",
                ]
            else:
                lines = [
                    r"Phi vs theta: Root Mean Square (code definition)",
                    r"Field slice: $F(k_x,\theta,k_y,t)=\phi/\rho_s$.",
                    r"Blank kx or ky selection averages over that spectral axis.",
                    r"Use $idx:n$ for an index or $val:x$ for the nearest physical kx/ky value.",
                    r"Optional max normalization divides the final profile by its finite maximum.",
                    r"$A(\theta)=\sqrt{\left\langle\mathrm{mean}_{k_x,k_y}|F(k_x,\theta,k_y,t)|^2\right\rangle_t}$.",
                    r"The time window follows the shared Time controls.",
                ]
        elif mode == "mean absolute":
            lines = [
                r"Fluctuation 1D: Mean Absolute (code definition)",
                r"Field slice used in code: $F(k_x,k_y,t)$ from midplane $\theta$ and radial index $[1:]$.",
                r"Normalization in code: $F\leftarrow F/\rho_s$ (uses case `rho`; if invalid, fallback $1$).",
                r"Blank fixed-axis selection averages; enter $idx:n$ or $val:x$ to select a mode.",
                r"Optional max normalization divides the final profile by its finite maximum.",
                r"vs $k_y$: $A(k_y)=\left\langle\mathrm{mean}_{k_x}|F(k_x,k_y,t)|\right\rangle_{t\in[t_0,t_1]}$; fixed kx selects one mode.",
                r"vs $k_x$: $A(k_x)=\left\langle\mathrm{mean}_{k_y}|F(k_x,k_y,t)|\right\rangle_{t\in[t_0,t_1]}$; fixed ky selects one mode.",
                r"vs Time: $n=0$ and $n>0$ channels are split, each sums over $k_x$ (and $k_y\neq0$ for $n>0$).",
                r"If selected window is empty, code uses last time slice.",
            ]
        else:
            lines = [
                r"Fluctuation 1D: Root Mean Square (code definition)",
                r"Field slice used in code: $F(k_x,k_y,t)$ from midplane $\theta$ and radial index $[1:]$.",
                r"Normalization in code: $F\leftarrow F/\rho_s$ (uses case `rho`; if invalid, fallback $1$).",
                r"Blank fixed-axis selection averages; enter $idx:n$ or $val:x$ to select a mode.",
                r"Optional max normalization divides the final profile by its finite maximum.",
                r"vs $k_y$: $A(k_y)=\sqrt{\left\langle\mathrm{mean}_{k_x}|F(k_x,k_y,t)|^2\right\rangle_{t\in[t_0,t_1]}}$; fixed kx selects one mode.",
                r"vs $k_x$: $A(k_x)=\sqrt{\left\langle\mathrm{mean}_{k_y}|F(k_x,k_y,t)|^2\right\rangle_{t\in[t_0,t_1]}}$; fixed ky selects one mode.",
                r"vs Time: $A_{n=0}(t)=\sqrt{\sum_{k_x}|F(k_x,k_y=0,t)|^2}$,",
                r"$A_{n>0}(t)=\sqrt{\sum_{k_x}\sum_{k_y\neq0}|F(k_x,k_y,t)|^2}$.",
                r"If selected window is empty, code uses last time slice.",
            ]

        self._draw_formula_panel(
            self.fluc_formula_fig,
            self.fluc_formula_ax,
            self.fluc_formula_canvas,
            lines,
            widget=self.fluc_formula_widget,
            base_fig_w=3.8,
            min_fig_h=1.8,
            title_fs=8.8,
            line_fs=7.3,
            line_dy=0.11,
            frac_dy=0.18,
            width_chars=56,
        )

    def _render_fft_formula_math(self):
        """Render math notes for FFT amplitude/power interpretation."""
        spectrum_mode = self.fft_spectrum_var.get().strip().lower()
        if spectrum_mode == "power":
            lines = [
                r"FFT Spectrum (Power)",
                r"$F(\omega)=\mathcal{F}_t[f(t)]$",
                r"$P(\omega)=|F(\omega)|^2$",
                r"Linear coherent: $P_{\mathrm{coh}}=|\sum_k F_k|^2$",
                r"Nonlinear incoherent: $P_{\mathrm{inc}}=\sum_k |F_k|^2$",
                r"Power mode in this GUI: no extra normalization.",
            ]
        else:
            lines = [
                r"FFT Spectrum (Amplitude)",
                r"$F(\omega)=\mathcal{F}_t[f(t)]$",
                r"$A(\omega)=|F(\omega)|$",
                r"Linear coherent: $A_{\mathrm{coh}}=|\sum_k F_k|$",
                r"Nonlinear incoherent: $A_{\mathrm{inc}}=\sum_k |F_k|$",
                r"Displayed as normalized amplitude:",
                r"$\widetilde{A}(\omega)=A(\omega)/\max_{\omega}A(\omega)$",
            ]

        self._draw_formula_panel(
            self.fft_formula_fig,
            self.fft_formula_ax,
            self.fft_formula_canvas,
            lines,
            widget=self.fft_formula_widget,
            base_fig_w=3.8,
            min_fig_h=2.0,
            title_fs=8.8,
            line_fs=8.0,
            line_dy=0.11,
            width_chars=56,
        )

    def _render_zf_formula_math(self, mode_text):
        """Render math notes for selected zonal ExB shearing sub-mode."""
        mode = (mode_text or "").strip().lower()
        if mode == "vs kx":
            lines = [
                r"Zonal ExB Shearing Rate (vs $k_x$, exact code path)",
                r"Source in code: $k_y$ index $0$,",
                r"$\theta$ index $i_\theta=4N_\theta/8$, radial slice $[1:]$.",
                r"Code normalization: $\phi_{abs}(k_x,t)=|\phi(k_x,t)|/\rho_s$ (uses case `rho`).",
                r"$\phi_{rms}(k_x)=\sqrt{\left\langle\phi_{abs}(k_x,t)^2\right\rangle_{t\in[t_0,t_1]}}$.",
                r"$\omega_{ZF}(k_x)=2k_x^2\phi_{rms}(k_x)$.",
                r"If window empty, uses last time slice; display keeps $k_x \geq 0$ branch.",
            ]
        elif mode == "phi vs kx(theta=0)":
            lines = [
                r"$\phi_{ZF}$ vs $k_x$ (theta label, exact code path)",
                r"Code uses $k_y$ index $0$, $\theta$ index $i_\theta=4N_\theta/8$, radial slice $[1:]$.",
                r"$\phi_{rms}(k_x)=\sqrt{\left\langle |\phi(k_x,t)|^2\right\rangle_{t\in[t_0,t_1]}}/\rho_s$.",
                r"If window empty, uses last time slice.",
                r"Plotted on full $k_x$ range (negative and positive).",
            ]
        elif mode == "vs gamma_lin":
            lines = [
                r"vs $\gamma_{lin}$ (code compares on shared axis with $k_x=k_y$)",
                r"$\phi(k_x,t)$ source in code: $k_y$ index $0$,",
                r"$\theta$ index $i_\theta=4N_\theta/8$, radial slice $[1:]$.",
                r"Code normalization: $\phi_{abs}(k_x,t)=|\phi(k_x,t)|/\rho$.",
                r"If time window is valid:",
                r"$\phi_{rms}(k_x)=\sqrt{\left\langle \phi_{abs}(k_x,t)^2\right\rangle_{t\in[t_0,t_1]}}$.",
                r"Else: $\phi_{rms}(k_x)=\phi_{abs}(k_x,t_{last})$.",
                r"$\omega_{ZF}(k_x)=2k_x^2\,\phi_{rms}(k_x)$.",
                r"$V_{ZF}^{mean}=0.5\sqrt{\sum_{k_x}|k_x\,\phi_{rms}(k_x)|^2}$.",
                r"Plotted third curve: $k_y V_{ZF}^{mean}$ on the same $k_x=k_y$ grid.",
                r"$\gamma_{lin}(k_y)$ is read from file and plotted directly.",
                r"finite values only, then $k_y>0$ is kept for the positive branch comparison.",
                r"If $k_y^\star$ is set, the ratios use interpolation of $\gamma_{lin}$ at $k_y^\star$.",
            ]
        else:
            lines = [
                r"Zonal ExB Shearing Rate (vs Time, exact code path)",
                r"Source in code: $k_y$ index $0$,",
                r"$\theta$ index $i_\theta=4N_\theta/8$, radial slice $[1:]$.",
                r"Normalization in code:",
                r"$\phi_{ZF}(k_x,t)=|\phi(k_x,t)|/\rho_s$ (uses case `rho`).",
                r"$\omega_{ZF}(t)=\sum_{k_x}k_x^2\phi_{ZF}(k_x,t)$.",
            ]

        self._draw_formula_panel(
            self.zf_formula_fig,
            self.zf_formula_ax,
            self.zf_formula_canvas,
            lines,
            widget=self.zf_formula_widget,
            base_fig_w=4.4,
            min_fig_h=3.8 if mode == "vs gamma_lin" else 2.2,
            title_fs=8.5,
            line_fs=7.4,
            line_dy=0.075 if mode == "vs gamma_lin" else 0.105,
            frac_dy=0.095 if mode == "vs gamma_lin" else 0.125,
            width_chars=60,
        )

    def _render_pod_formula_math(self):
        """Render math notes for POD parity decomposition."""
        field = str(self.others_pod_field_var.get()).strip()
        field_tex = r"\phi" if field.lower() == "phi" else r"A_{\parallel}"
        lines = [
            r"POD-Parity Decomposition",
            rf"$F(p,\theta,t)={field_tex}(k_x(p),k_y,\theta,t)$",
            r"$\widetilde{F}(p,\theta,t)=(-1)^p\,F(p,\theta,t)$  (alt-phase unwrap)",
            r"$X_{(p,\theta),t}=\widetilde{F}-\langle \widetilde{F}\rangle_t$",
            r"$X=U\Sigma V^\dagger,\;\;a_n(z)\leftarrow U_n$",
            r"$P_{\mathrm{even}}=\dfrac{\int |A(z)+A(-z)|^2\,dz}{4\int |A(z)|^2\,dz}$",
            r"$P_{\mathrm{odd}}=\dfrac{\int |A(z)-A(-z)|^2\,dz}{4\int |A(z)|^2\,dz}$",
            r"Label rule: $P_{\mathrm{even}}>0.7$ tearing, $P_{\mathrm{odd}}>0.7$ ballooning, else mixed.",
            r"Requirement: continuous parallel data only ($\theta>1$).",
        ]
        self._draw_formula_panel(
            self.pod_formula_fig,
            self.pod_formula_ax,
            self.pod_formula_canvas,
            lines,
            widget=self.pod_formula_widget,
            base_fig_w=3.8,
            min_fig_h=2.2,
            title_fs=8.8,
            line_fs=7.5,
            line_dy=0.095,
            frac_dy=0.13,
            width_chars=56,
        )

    def _energy_balance_mode_key(self):
        """Return a stable internal key for the Energy-balance mode combobox."""
        mode = str(self.energy_balance_mode_var.get()).strip().lower()
        if "zf" in mode:
            return "zf"
        if (
            "effective growth" in mode
            or "gamma_eff" in mode
            or r"\gamma_{eff}" in mode
        ):
            return "gamma_eff"
        if "single" in mode:
            return "single"
        if "fullt" in mode:
            return "fullt"
        if "2d" in mode:
            return "2d"
        return "entropy"

    def _render_energy_balance_formula_math(self):
        """Render math notes for Energy-balance sub-modes (triad_v2 style)."""
        mode_key = self._energy_balance_mode_key()
        is_gamma_eff_mode = (mode_key == "gamma_eff")
        is_energy_2d_mode = (mode_key == "2d")
        is_fullt_mode = (mode_key == "fullt")
        if mode_key == "zf":
            lines = [
                r"$\bf{ZF\ energy\ balance\ (code\ equations)}$",
                r"$f[s,r,c,n,t]=f_{\Re}+i\,f_{\Im}$",
                r"$N_0(t)=\sum_r \Re\{f[r,1,n,t]\}$",
                r"$N_{\mathrm{total}}(t)=\sum_a\sum_r \Re\{f[a,r,1,n,t]\}$",
                r"$\mathrm{idx}_3(t)\equiv\sum_r \Re\{f[r,2,n,t]\}=\sum_r \frac{d\!\left(T_a\,\delta S_a\right)}{dt}$",
                r"$\mathrm{idx}_4(t)=\sum_r \Re\{f[r,3,n,t]\}=\frac{dW_{em}}{dt}$",
                r"$W_{es}(t)=\sum_a \frac{n_a z_a^2}{2T_a}\sum_{k_x}\!\left(|\phi(k_x,k_y\!=\!0,t)|^2\,[1-\Gamma_{0,a}(k_x)]\right)$",
                r"$\Gamma_{0,a}(k_x)=I_0(b_a)e^{-b_a},\ \ b_a=(|k_x|\rho_a)^2$",
                r"$\frac{dW_{es}}{dt}=\frac{d}{dt}W_{es}(t)$",
                r"$L_{Z,\mathrm{total}}(t)=\frac{dW_{es}}{dt}-N_{\mathrm{total}}(t)$",
                r"$D_r=\sum_r\Re\{f[r,5,n,t]\},\ D_\theta=\sum_r\Re\{f[r,6,n,t]\},\ D_c=\sum_r\Re\{f[r,7,n,t]\}$",
                r"$D'_Z(t)=0$",
                r"$\mathrm{Plotted:}\ \mathcal{N}^{NZ\rightarrow Z},\ D'_Z,\ L_{Z,\mathrm{total}},\ L_{T_e}\,\frac{dW_{es}}{dt}$",
            ]
        elif is_gamma_eff_mode:
            lines = [
                r"$\bf{Triad\ v3:\ \gamma_{eff}\ spectra}$",
                r"$\mathrm{avg}_t\equiv\langle\cdot\rangle_{[t_0,t_1]},\ \ \sum_a:\mathrm{sum\ over\ species}$",
                r"$\mathrm{idx}\leftrightarrow\mathrm{physics}:\ 1\!\to\!T_{a,k_\perp},\ 2\!\to\![k_y\!=\!0]\mathcal{N}_{a,k_\perp}\ \mathrm{or}\ [k_y\!\ne\!0]T^{*}_{a,k_\perp},\ 5\!\to\!\delta S_{a,k_\perp}$",
                r"$\mathcal{E}(k_x,k_y)=\left\langle \sum_a \delta S_{a,k_\perp}\right\rangle_t$",
                r"$A_{NZ}(k_x,k_y)=\left\langle \sum_a [k_y\neq 0]\,T^{*}_{a,k_\perp}\right\rangle_t,\ \ A_{NZ}(k_x,k_y\!=\!0)=0$",
                r"$A_{Z}(k_x,k_y)=\left\langle \sum_a T_{a,k_\perp}\right\rangle_t - A_{NZ}(k_x,k_y)$",
                r"$\gamma_{eff}^{NZ}(k_x,k_y)= -\dfrac{A_{NZ}}{2\,\mathcal{E}},\ \ \gamma_{eff}^{Z}(k_x,k_y)= -\dfrac{A_{Z}}{2\,\mathcal{E}}$",
                r"$\gamma_{eff}^{NL}=\gamma_{eff}^{NZ}+\gamma_{eff}^{Z}$",
                r"$\mathrm{If\ }\gamma_{lin}\mathrm{\ file\ provided:}\ \gamma_{eff}^{stable}=\gamma_{lin}-\gamma_{eff}^{NL}$",
                r"$\mathrm{Plotted\ (no\ linear\ file):}\ \gamma_{eff}^{NZ},\ \gamma_{eff}^{Z},\ \gamma_{eff}^{NL}$",
                r"$\mathrm{Plotted\ (with\ linear\ file):}\ \gamma_{lin},\ \gamma_{eff}^{NZ},\ \gamma_{eff}^{Z},\ \gamma_{eff}^{stable}$",
            ]
        elif is_energy_2d_mode:
            lines = [
                r"$\bf{Energy\ balance\ vs\ 2D\ scan}$",
                r"$\mathrm{For\ each\ case,\ average\ over\ selected}\ [t_0,t_1]$",
                r"$T_{a}^{NZ\rightarrow Z}(t)=\sum_{k_x}\Re\{f[a,k_x,\mathrm{idx1},n,t]\}$",
                r"$N_{a}^{NZ\rightarrow Z}(t)=\sum_{k_x}\Re\{f[a,k_x,\mathrm{idx2},n,t]\}$",
                r"$\delta S_a(k_x,n,t)\equiv \Re\{f[a,k_x,\mathrm{idx5},n,t]\}$",
                r"$S(t)=\sum_{a}\sum_{k_x}\sum_{n\neq 0}\delta S_a(k_x,n,t)$",
                r"$\left(T_{a}^{NZ\rightarrow Z}/S\right)=\langle T_{a}^{NZ\rightarrow Z}\rangle_t/\langle S\rangle_t$",
                r"$\left(N_{a}^{NZ\rightarrow Z}/S\right)=\langle N_{a}^{NZ\rightarrow Z}\rangle_t/\langle S\rangle_t$",
                r"$\mathrm{Plotted:}\ \mathcal{N}_D/S,\ \mathcal{T}_D/S,\ \mathcal{N}_e/S,\ \mathcal{T}_e/S$",
            ]
        elif is_fullt_mode:
            use_asym = False
            try:
                use_asym = bool(self.energy_balance_transfer_asym_var.get())
            except Exception:
                use_asym = False
            try:
                view_txt = str(self.energy_balance_transfer_xaxis_var.get()).strip()
            except Exception:
                view_txt = "vs kxky"
            transfer_name = r"\mathrm{FULLT}_{\mathrm{asym}}" if use_asym else r"\mathrm{FULLT}"
            heading = r"$\bf{FULLT\ asym\ transfer\ map}$" if use_asym else r"$\bf{FULLT\ transfer\ map}$"
            map_name = r"\mathrm{FULLT}_{\mathrm{asym}}" if use_asym else r"\mathrm{FULLT}"
            sum_note = (
                r"$\mathrm{Asym\ diagnostic:}\ C(\psi' H''-\psi'' H')H_k^*,\ "
                r"\mathrm{not\ summed\ against\ original\ triad.}$"
                if use_asym
                else r"$\sum_{k^{\prime}}T_k^\Phi(k^{\prime})\simeq\mathrm{triad\ ch1}(k_x=0,k_{y,\mathrm{sel}}^{\mathrm{target}})$"
            )
            lines = [
                heading,
                rf"${transfer_name}[\Re/\Im,\ k_y^{{\mathrm{{target}}}},$",
                r"$\quad k_x^{\mathrm{source}},\ k_y^{\mathrm{source}},\ channel,\ t]$",
                r"$\mathrm{fixed:}\ k=(0,k_y^{\mathrm{target}})$",
                r"$\mathrm{plot:}\ k^{\prime}=(k_x^{\mathrm{source}},k_y^{\mathrm{source}})$",
                r"$\mathrm{Plotted\ quantity:}\ \Re\{T_k^\Phi(k^{\prime})\}\quad(\mathrm{fixed})$",
                rf"$\mathrm{{View:}}\ {view_txt}$",
                rf"$\mathrm{{Map:}}\ \langle{map_name}(k_{{y,\mathrm{{sel}}}}^{{\mathrm{{target}}}},$",
                r"$\quad k_x^{\prime},k_y^{\prime})\rangle_t=\langle T_k^\Phi(k^{\prime})\rangle_t$",
                r"$\mathrm{vs}\ k_x:\ \mathrm{take\ the}\ k_y^{\prime}=k_y^{\mathrm{target}}\ \mathrm{slice\ through\ the\ marker.}$",
                r"$\mathrm{trace}\ p_xq_x:\ \mathrm{fix}\ k_y^{source}=k_y^{target},\ \mathrm{plot}\ T(p_x,q_x).$",
                r"$\mathrm{zonal\ chain:}\ q_x=p_x\pm k_{zf};\ \ q_x\ \mathrm{becomes\ the\ next}\ p_x.$",
                r"$\mathrm{Optional\ display:}\ T>0\ \mathrm{divided\ by}\ \max(T),\quad T<0\ \mathrm{divided\ by}\ |\min(T)|$",
                sum_note,
                r"$T_k^\Phi>0:\ k^{\prime}\rightarrow k;\ \ \ T_k^\Phi<0:\ k\rightarrow k^{\prime}.$",
                r"$\mathrm{marker:}\ k=(0,k_y^{\mathrm{target}})$",
            ]
        elif mode_key == "single":
            lines = [
                r"$\bf{Energy\ balance\ single\ plot}$",
                r"$f[s,r,c,n,t]=f_{\Re}+i\,f_{\Im}$",
                r"$T_0(t)=\sum_r \Re\{f[r,0,n,t]\},\ \ N_0(t)=\sum_r \Re\{f[r,1,n,t]\}$",
                r"$\delta S_a(k_x,k_y,t)\equiv \Re\{f[a,k_x,\mathrm{idx5},k_y,t]\}$",
                r"$\mathrm{entropy}_{a}(k_y)=\log\!\left(\sum_{k_x}\left\langle \delta S_a(k_x,k_y,t)\right\rangle_t\right)$",
                r"$\mathrm{(for\ vs\ ky,\ entropy\ uses\ idx5,\ not\ idx3)}$",
                r"$\mathrm{Quantity} \in \{T,\ N,\ T\!-\!N,\ D_Z,\ entropy\}$",
                r"$\mathrm{X\!-\!axis} \in \{t,\ k_y,\ k_x,\ (k_x,k_y)\}$",
                r"$\mathrm{vs}\ t,\ k_x:\ \mathrm{use\ ky\ scan\ to\ match\ nearest\ stored}\ k_y$",
                r"$\mathrm{vs}\ k_x:\ \mathrm{plot}\ \langle T(k_x,k_y,t)\rangle_t$",
                r"$\mathrm{Normalize:}\ \{T,N,T\!-\!N\}/|\min(T)|\ \mathrm{or}\ /|\max(T)|;$",
                r"$\mathrm{entropy:}\ S_g/|\min(S_g)|,\ S_g/|\max(S_g)|,\ \mathrm{or}\ S_g/S_g(k_y=0).$",
            ]
        else:
            lines = [
                r"$\bf{Gyro\!-\!center\ entropy\ balance\ (code\ equations)}$",
                r"$f[s,r,c,n,t]=f_{\Re}+i\,f_{\Im}$",
                r"$T_0(t)=\sum_r \Re\{f[r,0,n,t]\},\ \ N_0(t)=\sum_r \Re\{f[r,1,n,t]\}$",
                r"$N_{\mathrm{total}}(t)=\sum_a\sum_r \Re\{f[a,r,1,n,t]\}$",
                r"$\mathrm{idx}_{3,\mathrm{total}}(t)\equiv\sum_a\sum_r \Re\{f[a,r,2,n,t]\}=\sum_a\sum_r \frac{d\!\left(T_a\,\delta S_a\right)}{dt}$",
                r"$\mathrm{idx}_{4}(t)=\sum_r \Re\{f[r,3,n,t]\}=\frac{dW_{em}}{dt}$",
                r"$D_r=\sum_r\Re\{f[r,5,n,t]\},\ D_\theta=\sum_r\Re\{f[r,6,n,t]\},\ D_c=\sum_r\Re\{f[r,7,n,t]\}$",
                r"$D_Z(t)=D_r+D_\theta+D_c$",
                r"$W_{es}(t)=\sum_a \frac{n_a z_a^2}{2T_a}\sum_{k_x}\!\left(|\phi(k_x,k_y\!=\!0,t)|^2\,[1-\Gamma_{0,a}(k_x)]\right)$",
                r"$\Gamma_{0,a}(k_x)=I_0(b_a)e^{-b_a},\ \ b_a=(|k_x|\rho_a)^2$",
                r"$\frac{dW_{es}}{dt}=\frac{d}{dt}W_{es}(t),\ \ \frac{dS_g}{dt}=\mathrm{idx}_{3,\mathrm{total}}-\frac{dW_{es}}{dt}$",
                r"$L_{Z,\mathrm{total}}(t)=\frac{dW_{es}}{dt}-N_{\mathrm{total}}(t)$",
                r"$\mathrm{Plotted:}\ (\mathcal{T}-\mathcal{N})^{NZ\rightarrow Z},\ D_Z,\ -L_{Z,\mathrm{total}},\ \frac{dS_g}{dt}$",
            ]
        self._draw_formula_panel(
            self.energy_balance_formula_fig,
            self.energy_balance_formula_ax,
            self.energy_balance_formula_canvas,
            lines,
            widget=self.energy_balance_formula_widget,
            base_fig_w=4.4,
            min_fig_h=4.0,
            title_fs=8.8,
            line_fs=7.4 if is_fullt_mode else 7.8,
            line_dy=0.074 if is_fullt_mode else 0.105,
            frac_fs=7.0 if is_fullt_mode else 7.2,
            frac_dy=0.084 if is_fullt_mode else 0.125,
            width_chars=58,
            min_fs=5.4 if is_fullt_mode else 5.8,
        )

    @staticmethod
    def _pod_theta_resolution_error_text():
        return (
            "Error: Insufficient theta resolution. "
            "Parity decomposition requires continuous parallel field data (theta > 1)."
        )

    def _show_pod_theta_resolution_warning(self):
        msg = self._pod_theta_resolution_error_text()
        try:
            messagebox.showwarning("Error", msg)
        except Exception:
            pass
        print(msg)

    @staticmethod
    def _case_theta_resolution_is_insufficient(data):
        """Return True when known theta resolution is insufficient for parity POD."""
        theta_candidates = []
        for attr in ("theta_plot", "n_theta"):
            v = getattr(data, attr, None)
            try:
                iv = int(v)
                if iv > 0:
                    theta_candidates.append(iv)
            except Exception:
                pass
        if theta_candidates:
            return min(theta_candidates) <= 1

        # Optional fallback from already-loaded field arrays.
        for attr in ("kxky_apar", "kxky_phi"):
            raw = getattr(data, attr, None)
            if raw is None:
                continue
            arr = np.asarray(raw)
            if arr.ndim == 5 and arr.shape[0] == 2:
                n_th = int(arr.shape[2])  # [2,nr,theta,ky,t]
                return n_th <= 1
            if arr.ndim == 4:
                n_th = int(arr.shape[1])  # [nr,theta,ky,t]
                return n_th <= 1
        return False

    def update_options(self, event=None):
        """Refresh dynamic option layout according to currently selected plot type."""
        plot_type_changed = (
            event is not None
            and getattr(event, "widget", None) is getattr(self, "plot_type_combo", None)
        )
        self._refresh_case_summary()
        self._hide_dynamic_options()
        plot_type = self.plot_type_var.get()
        row = 3 # Start below the shared Time Start/End and log-scale controls.

        if plot_type in ["Frequency", "Growth Rate"]:
            self.norm_ky_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
            
        elif plot_type == "Flux":
            flux_mode = self.flux_xaxis_var.get().strip().lower()
            self.flux_type_combo.grid(row=row, column=0, columnspan=2, sticky=tk.W + tk.E)
            row += 1
            self.flux_xaxis_combo.grid(row=row, column=0, columnspan=2, sticky=tk.W + tk.E)
            row += 1
            self.species_label.grid(row=row, column=0, sticky=tk.W)
            self.species_combo.grid(row=row, column=1, sticky=tk.W + tk.E)
            row += 1
            if flux_mode == "v.s 2d":
                self._refresh_flux_2d_param_choices()
                self.flux_scan_xparam_label.grid(row=row, column=0, sticky=tk.W)
                self.flux_scan_xparam_combo.grid(row=row, column=1, sticky=tk.W + tk.E)
                row += 1
                self.flux_2d_errorbar_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
                row += 1
            else:
                self.flux_decomp_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
                row += 1
            self.flux_norm_real_ion_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
            row += 1
            if flux_mode not in ["v.s 2d", "v.s ky_time"]:
                self.plot_all_species_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
                row += 1
            if flux_mode == "v.s kx (estimated)":
                self._render_flux_kx_formula_math()
                self.flux_formula_frame.grid(row=row, column=0, columnspan=4, sticky=tk.W + tk.E, pady=(4, 0))

        elif plot_type == "Fluctuation 1D":
            self.fluc_field_combo.grid(row=row, column=0, sticky=tk.W)
            self.fluc_xaxis_combo.grid(row=row, column=1, sticky=tk.W)
            row += 1

            fluc_xaxis = self.fluc_xaxis_var.get()
            self.fluc_advanced_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
            row += 1

            advanced = self._fluc_advanced_enabled()
            if advanced and fluc_xaxis in ("v.s ky", "v.s kx", "v.s theta"):
                self._ensure_fluc_advanced_case_values()
                self.fluc_advanced_button.grid(
                    row=row, column=0, columnspan=2, sticky=tk.W
                )
                row += 1
            else:
                if fluc_xaxis == "v.s ky":
                    self.fluc_theta_kx_label.configure(text="fixed kx (blank=avg):")
                    self.fluc_theta_kx_label.grid(row=row, column=0, sticky=tk.W)
                    self.fluc_theta_kx_entry.grid(row=row, column=1, sticky=tk.W)
                    row += 1
                elif fluc_xaxis == "v.s kx":
                    self.fluc_theta_ky_label.configure(text="fixed ky (blank=avg):")
                    self.fluc_theta_ky_label.grid(row=row, column=0, sticky=tk.W)
                    self.fluc_theta_ky_entry.grid(row=row, column=1, sticky=tk.W)
                    row += 1
                elif fluc_xaxis == "v.s theta":
                    self.fluc_theta_kx_label.configure(text="kx (blank=avg):")
                    self.fluc_theta_ky_label.configure(text="ky (blank=avg):")
                    self.fluc_theta_kx_label.grid(row=row, column=0, sticky=tk.W)
                    self.fluc_theta_kx_entry.grid(row=row, column=1, sticky=tk.W)
                    self.fluc_theta_ky_label.grid(row=row, column=2, sticky=tk.W, padx=(8, 0))
                    self.fluc_theta_ky_entry.grid(row=row, column=3, sticky=tk.W)
                    row += 1

            if fluc_xaxis in ("v.s ky", "v.s kx", "v.s theta"):
                self.fluc_norm_max_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
                row += 1
            
            # Check if FFT is selected in the sub-option
            if self.fluc_xaxis_var.get() == "fft":
                 self.fft_options_frame.grid(row=row, column=0, columnspan=2, sticky=tk.W+tk.E, pady=5)
                 row += 1
                 self._render_fft_formula_math()
                 self.fft_formula_frame.grid(row=row, column=0, columnspan=4, sticky=tk.W + tk.E, pady=(4, 0))
            else:
                 self._render_fluctuation_1d_formula_math()
                 self.fluc_formula_frame.grid(row=row, column=0, columnspan=4, sticky=tk.W + tk.E, pady=(4, 0))

        elif plot_type == "Fluctuation 2D":
            self.fluc2d_view_label.grid(row=row, column=0, sticky=tk.W)
            self.fluc2d_view_combo.grid(row=row, column=1, sticky=tk.W)
            row += 1
            self.moment_label.grid(row=row, column=0, sticky=tk.W)
            self.moment_combo.grid(row=row, column=1, sticky=tk.W)
            row += 1
            # Species selection needed for Density/Energy/Temperature moments
            if self.moment_var.get() in ["Density", "Energy", "Temperature"]:
                 self.species_label.grid(row=row, column=0, sticky=tk.W)
                 self.species_combo.grid(row=row, column=1, sticky=tk.W)
                 row += 1
            if self._is_fluc2d_kxky_view(self.fluc2d_view_var.get()):
                self.fluc2d_log_z_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
            else:
                self.fluc2d_x_elec_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)

        elif plot_type == "Zonal ExB Shearing Rate":
            self.zf_xaxis_label.grid(row=row, column=0, sticky=tk.W)
            self.zf_xaxis_combo.grid(row=row, column=1, sticky=tk.W)
            row += 1
            zf_mode = self.zf_xaxis_var.get().strip().lower()
            if zf_mode == "vs gamma_lin":
                self.linear_gamma_file_label.grid(row=row, column=0, sticky=tk.W)
                self.linear_gamma_file_entry.grid(row=row, column=1, columnspan=2, sticky=tk.W + tk.E)
                self.linear_gamma_file_browse.grid(row=row, column=3, sticky=tk.W)
                row += 1
                self.zf_gamma_lin_ky_label.grid(row=row, column=0, sticky=tk.W)
                self.zf_gamma_lin_ky_entry.grid(row=row, column=1, sticky=tk.W)
                row += 1
            self._render_zf_formula_math(zf_mode)
            self.zf_formula_frame.grid(row=row, column=0, columnspan=4, sticky=tk.W + tk.E, pady=(4, 0))

        elif plot_type == "Energy balance":
            self.energy_balance_options_frame.grid(
                row=row, column=0, columnspan=4, sticky=tk.W + tk.E
            )
            erow = 0
            self.energy_balance_mode_label.grid(
                row=erow, column=0, sticky=tk.W, padx=(0, 6), pady=(0, 4)
            )
            self.energy_balance_mode_combo.grid(
                row=erow, column=1, columnspan=3, sticky=tk.W + tk.E, pady=(0, 4)
            )
            erow += 1
            row += 1
            mode_key = self._energy_balance_mode_key()
            is_gamma_eff_mode = (mode_key == "gamma_eff")
            is_energy_2d_mode = (mode_key == "2d")
            is_single_plot_mode = (mode_key == "single")
            is_fullt_mode = (mode_key == "fullt")
            xaxis_single = str(self.energy_balance_single_xaxis_var.get()).strip().lower()
            if is_gamma_eff_mode:
                n_label_txt = "ky:"
            elif is_single_plot_mode:
                n_label_txt = "ky scan:"
            else:
                n_label_txt = "n index:"
            self.energy_balance_n_label.configure(text=n_label_txt)
            if not is_fullt_mode:
                self.energy_balance_n_label.grid(
                    row=erow, column=0, sticky=tk.W, padx=(0, 6), pady=2
                )
                self.energy_balance_n_entry.grid(
                    row=erow, column=1, sticky=tk.W + tk.E, pady=2
                )
                self.energy_balance_spec_label.grid(
                    row=erow, column=2, sticky=tk.W, padx=(10, 6), pady=2
                )
                self.energy_balance_spec_combo.grid(
                    row=erow, column=3, sticky=tk.W + tk.E, pady=2
                )
                erow += 1
            if is_single_plot_mode:
                qty_single = str(self.energy_balance_single_quantity_var.get()).strip().lower()
                self.energy_balance_single_quantity_label.grid(
                    row=erow, column=0, sticky=tk.W, padx=(0, 6), pady=2
                )
                self.energy_balance_single_quantity_combo.grid(
                    row=erow, column=1, sticky=tk.W + tk.E, pady=2
                )
                self.energy_balance_single_xaxis_label.grid(
                    row=erow, column=2, sticky=tk.W, padx=(10, 6), pady=2
                )
                self.energy_balance_single_xaxis_combo.grid(
                    row=erow, column=3, sticky=tk.W + tk.E, pady=2
                )
                erow += 1
                if qty_single == "entropy":
                    self.energy_balance_single_norm_combo.configure(
                        values=list(self._ENERGY_BALANCE_SINGLE_ENTROPY_NORM_OPTIONS)
                    )
                    current_norm = str(self.energy_balance_single_norm_var.get()).strip().lower()
                    if current_norm not in ["none", "min entropy", "max entropy", "ky=0 entropy"]:
                        self.energy_balance_single_norm_var.set("None")
                    norm_state = "readonly"
                else:
                    self.energy_balance_single_norm_combo.configure(
                        values=list(self._ENERGY_BALANCE_SINGLE_NORM_OPTIONS)
                    )
                    current_norm = str(self.energy_balance_single_norm_var.get()).strip()
                    if current_norm.lower() in ["min entropy", "max entropy", "ky=0 entropy"]:
                        self.energy_balance_single_norm_var.set("None")
                    norm_state = tk.DISABLED if qty_single in ["dr", "dtheta", "dc", "dz"] else "readonly"
                self.energy_balance_single_norm_label.grid(
                    row=erow, column=0, sticky=tk.W, padx=(0, 6), pady=2
                )
                self.energy_balance_single_norm_combo.configure(state=norm_state)
                self.energy_balance_single_norm_combo.grid(
                    row=erow, column=1, sticky=tk.W + tk.E, pady=2
                )
                erow += 1
            if is_fullt_mode:
                self.energy_balance_transfer_xaxis_label.grid(
                    row=erow, column=0, sticky=tk.W, padx=(0, 6), pady=2
                )
                self.energy_balance_transfer_xaxis_combo.grid(
                    row=erow, column=1, sticky=tk.W + tk.E, pady=2
                )
                self.energy_balance_transfer_ky_label.grid(
                    row=erow, column=2, sticky=tk.W, padx=(10, 6), pady=2
                )
                self.energy_balance_transfer_ky_entry.grid(
                    row=erow, column=3, sticky=tk.W + tk.E, pady=2
                )
                erow += 1
                transfer_view = str(self.energy_balance_transfer_xaxis_var.get()).strip().lower()
                if transfer_view == "trace pxqx":
                    self.energy_balance_transfer_kx_label.configure(text="fixed zonal kx:")
                else:
                    self.energy_balance_transfer_kx_label.configure(text="fixed source kx:")
                self.energy_balance_transfer_kx_label.grid(
                    row=erow, column=2, sticky=tk.W, padx=(10, 6), pady=2
                )
                self.energy_balance_transfer_kx_entry.grid(
                    row=erow, column=3, sticky=tk.W + tk.E, pady=2
                )
                self.energy_balance_transfer_asym_check.grid(
                    row=erow, column=0, sticky=tk.W, pady=2
                )
                self.energy_balance_transfer_norm_max_check.grid(
                    row=erow, column=1, sticky=tk.W, padx=(8, 0), pady=2
                )
                erow += 1
            if is_gamma_eff_mode:
                self.linear_gamma_file_label.grid(row=row, column=0, sticky=tk.W)
                self.linear_gamma_file_entry.grid(row=row, column=1, columnspan=2, sticky=tk.W + tk.E)
                self.linear_gamma_file_browse.grid(row=row, column=3, sticky=tk.W)
                row += 1
            if is_energy_2d_mode:
                self._refresh_flux_2d_param_choices()
                self.flux_scan_xparam_label.grid(row=row, column=0, sticky=tk.W)
                self.flux_scan_xparam_combo.grid(row=row, column=1, columnspan=3, sticky=tk.W + tk.E)
                row += 1
            self._render_energy_balance_formula_math()
            self.energy_balance_formula_frame.grid(
                row=row, column=0, columnspan=4, sticky=tk.W + tk.E, pady=(8, 0)
            )

        elif plot_type == "Others":
            self.others_plot_label.grid(row=row, column=0, sticky=tk.W)
            self.others_plot_combo.grid(row=row, column=1, sticky=tk.W)
            other_mode = self.others_plot_var.get().strip().lower()
            if other_mode == "rcorr_phi":
                row += 1
                self.others_rcorr_field_label.grid(row=row, column=0, sticky=tk.W)
                self.others_rcorr_field_combo.grid(row=row, column=1, sticky=tk.W)
                row += 1
                self.others_rcorr_theta_label.grid(row=row, column=0, sticky=tk.W)
                self.others_rcorr_theta_entry.grid(row=row, column=1, sticky=tk.W)
            elif other_mode == "pod_parity":
                self.others_pod_field_label.grid(row=row, column=2, sticky=tk.W, padx=(10, 0))
                self.others_pod_field_combo.grid(row=row, column=3, sticky=tk.W)
                row += 1
                self.others_pod_kx_label.grid(row=row, column=0, sticky=tk.W)
                self.others_pod_kx_entry.grid(row=row, column=1, sticky=tk.W)
                row += 1
                self.others_pod_ky_label.grid(row=row, column=0, sticky=tk.W)
                self.others_pod_ky_entry.grid(row=row, column=1, sticky=tk.W)
                row += 1
                self._render_pod_formula_math()
                self.pod_formula_frame.grid(row=row, column=0, columnspan=4, sticky=tk.W + tk.E, pady=(4, 0))

        if plot_type_changed:
            # A newly selected mode should start at its first option instead of
            # retaining a deep scroll offset from the previous, longer form.
            self.root.after_idle(self._scroll_plot_options_to_top)

    def _get_case_species(self, data):
        """Extract species (Z, Mass) tuples from a case object."""
        try:
            z = None
            m = None
            
            if hasattr(data, 'z') and hasattr(data, 'mass'):
                z = data.z
                m = data.mass
            
            if z is not None and m is not None:
                # Combine into list of tuples, rounded to avoid float issues
                return [(round(float(zi), 2), round(float(mi), 2)) for zi, mi in zip(z, m)]
            
            return []
        except Exception as e:
            print(f"Error getting species: {e}")
            return []

    def _collect_flux_scan_case_scalars(self, selected_case_names=None):
        """
        Return per-case scalar dictionaries parsed from each case `input.cgyro`.

        Output format: list of `(case_name, scalars_dict)`.
        """
        if selected_case_names is None:
            selected_case_names = self._get_selected_case_names()

        out = []
        for case_name in selected_case_names:
            data = self.cases.get(case_name, None)
            if data is None:
                continue
            try:
                case_dir = self._resolve_case_dir(data)
            except Exception:
                case_dir = None
            scalars = self._read_input_cgyro_scalars(case_dir) if case_dir else {}
            if hasattr(scalars, 'keys') and len(scalars) > 0:
                out.append((case_name, scalars))
        return out

    def _get_flux_2d_varying_params(self, selected_case_names=None):
        """
        Return sorted scalar keys that vary across selected cases.
        """
        case_scalars = self._collect_flux_scan_case_scalars(selected_case_names)
        if len(case_scalars) <= 0:
            return []

        common_keys = None
        for _name, scalars in case_scalars:
            keys = set(str(k) for k in scalars.keys())
            common_keys = keys if common_keys is None else (common_keys & keys)
        if not common_keys:
            return []

        varying = []
        for key in sorted(common_keys):
            if not self._is_flux_2d_physical_param_key(key):
                continue
            values = []
            ok = True
            for _name, scalars in case_scalars:
                try:
                    values.append(float(scalars.get(key)))
                except Exception:
                    ok = False
                    break
            if not ok:
                continue
            if len(values) > 1 and (max(values) - min(values) > 1.0e-12):
                varying.append(str(key))
        return varying

    @staticmethod
    def _is_flux_2d_physical_param_key(key):
        """
        True when `key` is treated as a physical input parameter for Flux-vs-2D scans.

        This filter excludes numerical/solver controls (grid size, time-step, flags, etc.)
        and keeps plasma/geometry/profile/species physics knobs.
        """
        k = str(key).strip().upper()
        if len(k) == 0:
            return False

        # Species-indexed physics parameters, e.g. Z_1, MASS_2, DENS_3, TEMP_1 ...
        if re.match(
            r'^(Z|MASS|DENS|TEMP|DLNNDR|DLNTDR|SDLNNDR|SDLNTDR|NU)_\d+$',
            k
        ):
            return True

        physical_exact = {
            # Core physics scalars
            'KY', 'BETAE_UNIT', 'BETAE', 'BETA_STAR', 'LAMBDA_STAR',
            'GAMMA_E', 'GAMMA_P', 'MACH', 'Z_EFF', 'NU_EE',
            # Geometry / equilibrium
            'RMIN', 'RMAJ', 'Q', 'S', 'SHIFT',
            'KAPPA', 'S_KAPPA', 'DELTA', 'S_DELTA', 'ZETA', 'S_ZETA',
            'ZMAG', 'DZMAG',
            # Shape coefficients
            'SHAPE_SIN3', 'SHAPE_S_SIN3', 'SHAPE_SIN4', 'SHAPE_S_SIN4',
            'SHAPE_SIN5', 'SHAPE_S_SIN5', 'SHAPE_SIN6', 'SHAPE_S_SIN6',
            'SHAPE_COS0', 'SHAPE_S_COS0', 'SHAPE_COS1', 'SHAPE_S_COS1',
            'SHAPE_COS2', 'SHAPE_S_COS2', 'SHAPE_COS3', 'SHAPE_S_COS3',
            'SHAPE_COS4', 'SHAPE_S_COS4', 'SHAPE_COS5', 'SHAPE_S_COS5',
            'SHAPE_COS6', 'SHAPE_S_COS6',
            # Model-physics choices often used in scans
            'N_SPECIES', 'TEMP_AE', 'DENS_AE', 'MASS_AE',
            'DLNTDR_AE', 'DLNNDR_AE',
        }
        if k in physical_exact:
            return True

        # Common alternate species-parameter naming styles in some input files.
        if re.match(r'^(AS_|TAUS_|MACHS_|ZEFFS_)\d+$', k):
            return True

        # Everything else defaults to non-physical (numerical/control) for this GUI mode.
        return False

    def _refresh_flux_2d_param_choices(self):
        """Refresh Flux-vs-2D x-parameter combo from currently selected cases."""
        try:
            varying = self._get_flux_2d_varying_params()
        except Exception:
            varying = []

        if len(varying) <= 0:
            self.flux_scan_xparam_combo['values'] = []
            self.flux_scan_xparam_var.set("")
            return

        self.flux_scan_xparam_combo['values'] = varying
        current = str(self.flux_scan_xparam_var.get()).strip()
        if current not in varying:
            self.flux_scan_xparam_var.set(varying[0])

    def _get_species_name(self, z, m):
        """Return a descriptive name for the species based on Z and M."""
        # Note: Mass is normalized to Deuterium (m=2.014 u), so we multiply by 2 to get approx AMU
        real_m = m * 2.0
        
        if abs(z - (-1)) < 0.1: return "Electron"
        if abs(z - 1) < 0.1:
            if abs(real_m - 1) < 0.2: return "Hydrogen" # m=0.5
            if abs(real_m - 2) < 0.2: return "Deuterium" # m=1.0
            if abs(real_m - 3) < 0.2: return "Tritium" # m=1.5
        if abs(z - 2) < 0.1:
            if abs(real_m - 3) < 0.2: return "Helium-3" # m=1.5
            if abs(real_m - 4) < 0.2: return "Helium-4" # m=2.0
        if abs(z - 3) < 0.1: return "Lithium"
        if abs(z - 4) < 0.1: return "Beryllium"
        if abs(z - 5) < 0.1: return "Boron"
        if abs(z - 6) < 0.1: return "Carbon"
        if abs(z - 7) < 0.1: return "Nitrogen"
        if abs(z - 8) < 0.1: return "Oxygen"
        if abs(z - 10) < 0.1: return "Neon"
        if abs(z - 18) < 0.1: return "Argon"
        if abs(z - 28) < 0.1: return "Nickel"
        if abs(z - 74) < 0.1: return "Tungsten"
        
        return f"Ion (Z={z}, M={real_m:.1f})"

    def _get_species_short_name(self, z, m):
        """Return a compact species token for plot labels."""
        name = self._get_species_name(z, m)
        compact = {
            "Electron": "e",
            "Hydrogen": "H",
            "Deuterium": "D",
            "Tritium": "T",
            "Helium-3": "He3",
            "Helium-4": "He4",
            "Lithium": "Li",
            "Beryllium": "Be",
            "Boron": "B",
            "Carbon": "C",
            "Nitrogen": "N",
            "Oxygen": "O",
            "Neon": "Ne",
            "Argon": "Ar",
            "Nickel": "Ni",
            "Tungsten": "W",
        }
        if name in compact:
            return compact[name]
        if abs(z - round(z)) < 1.0e-8:
            return f"Z{int(round(z))}"
        return f"Z{z:.3g}"

    def _get_case_species_densities(self, data, n_species=None):
        """Return normalized DENS values for the case species."""
        specs = self._get_case_species(data)
        if n_species is None:
            n_species = len(specs)
        dens = np.full(int(n_species), np.nan, dtype=float)

        try:
            dens_attr = np.asarray(getattr(data, 'dens', []), dtype=float).reshape(-1)
            ncopy = min(dens.size, dens_attr.size)
            if ncopy > 0:
                dens[:ncopy] = dens_attr[:ncopy]
        except Exception:
            pass

        if np.any(~np.isfinite(dens)):
            scalars = {}
            try:
                case_dir = self._resolve_case_dir(data)
                scalars = self._read_input_cgyro_scalars(case_dir)
            except Exception:
                scalars = {}
            for i in range(dens.size):
                if np.isfinite(dens[i]):
                    continue
                try:
                    val = scalars.get(f"DENS_{i+1}", None)
                    if val is not None:
                        dens[i] = float(val)
                except Exception:
                    pass

        return dens

    def _format_species_join_label(self, data, indices):
        """Format selected species as compact tokens joined by '+'."""
        specs = self._get_case_species(data)
        labels = []
        for idx in indices:
            ii = int(idx)
            if 0 <= ii < len(specs):
                z, m = specs[ii]
                labels.append(self._get_species_short_name(z, m))
            else:
                labels.append(f"s{ii+1}")
        return "+".join(labels)

    def _update_species_list(self):
        """Find common species across all loaded cases and update the dropdown."""
        if not self.cases:
            self.species_combo['values'] = []
            self.species_var.set("")
            return

        common_species = None
        has_main_ion = False
        has_any_ion = False
        
        for data in self.cases.values():
            specs = set(self._get_case_species(data))
            if common_species is None:
                common_species = specs
            else:
                common_species &= specs
            
            # Main ion is defined per case by normalized DENS > 0.3 among ions.
            if self._get_main_ion_indices(data):
                has_main_ion = True
            for z, m in specs:
                if z > 0:
                    has_any_ion = True
        
        if common_species:
            # Sort by Z (usually ions first, then electrons if Z<0)
            sorted_species = sorted(list(common_species), key=lambda x: (x[0], x[1]), reverse=True)
            
            values = []
            
            # Add Main Ion option if applicable
            if has_main_ion:
                values.append("Main Ion (DENS>0.3)")
            if has_any_ion:
                values.append("All Ions")

            for z, m in sorted_species:
                name = self._get_species_name(z, m)
                values.append(f"{name} (Z={z}, M={m})")
            
            self.species_combo['values'] = values
            
            # Preserve selection if possible, else select first
            current = self.species_var.get()
            if current not in values and values:
                self.species_combo.current(0)
        else:
            self.species_combo['values'] = []
            self.species_var.set("No common species")
        try:
            self._refresh_flux_2d_param_choices()
        except Exception:
            pass

    def _get_main_ion_indices(self, data):
        """Return ion indices whose normalized DENS is greater than 0.3."""
        indices = []
        specs = self._get_case_species(data)
        dens = self._get_case_species_densities(data, n_species=len(specs))
        for i, (z, m) in enumerate(specs):
            if i >= dens.size:
                continue
            if z > 0 and np.isfinite(dens[i]) and dens[i] > 0.3:
                indices.append(i)
        return indices

    def _get_all_ion_indices(self, data):
        """Return indices of all ion species (Z > 0)."""
        indices = []
        specs = self._get_case_species(data)
        for i, (z, _m) in enumerate(specs):
            if z > 0:
                indices.append(i)
        return indices

    def _build_effective_plot_type(self):
        """Translate GUI selections into internal plot-type token and display label."""
        # `plot_type` is a stable backend command string; `display_plot_type`
        # is allowed to be more descriptive.  Keeping both lets old plotting
        # handlers remain unchanged while the UI text stays readable.
        plot_type_selection = self.plot_type_var.get()
        plot_type = plot_type_selection
        display_plot_type = plot_type_selection

        if plot_type_selection == "Flux":
            flux_type = self.flux_type_var.get()
            flux_xaxis = self.flux_xaxis_var.get()
            xaxis_str = flux_xaxis.replace("v.s", "vs")
            plot_type = f"{flux_type} Flux {xaxis_str}"
            display_plot_type = plot_type
            if self.flux_decomp_var.get() and flux_xaxis.strip().lower() != "v.s 2d":
                plot_type += " (Decomp)"
        elif plot_type_selection == "Fluctuation 1D":
            field = self.fluc_field_var.get()
            xaxis = self.fluc_xaxis_var.get()
            if xaxis == "fft":
                plot_type = f"{field} FFT"
            else:
                xaxis_str = xaxis.replace("v.s", "vs")
                plot_type = f"{field} {xaxis_str}"
                if xaxis in ("v.s ky", "v.s kx", "v.s theta"):
                    if self._fluc_advanced_enabled():
                        if xaxis == "v.s ky":
                            display_plot_type = f"{field} vs ky (per-case kx)"
                        elif xaxis == "v.s kx":
                            display_plot_type = f"{field} vs kx (per-case ky)"
                        else:
                            display_plot_type = f"{field} vs theta (per-case kx/ky)"
                    else:
                        kx_text = self.fluc_theta_kx_var.get().strip() or "average"
                        ky_text = self.fluc_theta_ky_var.get().strip() or "average"
                        if xaxis == "v.s ky":
                            display_plot_type = f"{field} vs ky (kx={kx_text})"
                        elif xaxis == "v.s kx":
                            display_plot_type = f"{field} vs kx (ky={ky_text})"
                        else:
                            display_plot_type = f"{field} vs theta (kx={kx_text}, ky={ky_text})"
        elif plot_type_selection == "Fluctuation 2D":
            view = self.fluc2d_view_var.get().strip().lower()
            # Use the exact mapping above instead of composing strings here;
            # this avoids accidental overlap with generic "vs kx" 1D routing.
            plot_type, display_plot_type = self._FLUC2D_VIEW_PLOT_TYPES.get(
                view,
                self._FLUC2D_VIEW_PLOT_TYPES["vs xy"],
            )
        elif plot_type_selection == "Zonal ExB Shearing Rate":
            zf_xaxis = self.zf_xaxis_var.get().strip().lower()
            if zf_xaxis == "vs kx":
                plot_type = "ZF ExB Shearing Spectrum"
                display_plot_type = "Zonal ExB Shearing Rate (vs kx)"
            elif zf_xaxis == "phi vs kx(theta=0)":
                plot_type = "ZF Phi Spectrum (theta0)"
                display_plot_type = "Zonal ExB Shearing Rate (phi vs kx, theta=0)"
            elif zf_xaxis == "vs gamma_lin":
                plot_type = "ZF ExB Fig4 (kx=ky)"
                display_plot_type = "Zonal ExB Shearing Rate (vs gamma_lin)"
            else:
                plot_type = "ZF ExB Shearing Rate"
                display_plot_type = "Zonal ExB Shearing Rate (vs Time)"
        elif plot_type_selection == "Energy balance":
            mode_key = self._energy_balance_mode_key()
            if mode_key == "zf":
                plot_type = "Energy Balance ZF"
                display_plot_type = "Energy balance: ZF energy balance"
            elif mode_key == "gamma_eff":
                plot_type = "Energy Balance Gamma Eff"
                display_plot_type = r"Energy balance: vs $\gamma_{eff}^Z$ and $\gamma_{eff}^{NZ}$"
            elif mode_key == "single":
                qty = str(self.energy_balance_single_quantity_var.get()).strip()
                xax = str(self.energy_balance_single_xaxis_var.get()).strip()
                xax_plot = xax.replace("v.s", "vs")
                plot_type = f"Energy Balance Single {qty} {xax_plot}"
                display_plot_type = f"Energy balance: Single plot ({qty}, {xax})"
            elif mode_key == "fullt":
                ky = str(self.energy_balance_transfer_ky_var.get()).strip()
                kx = str(self.energy_balance_transfer_kx_var.get()).strip()
                xax = str(self.energy_balance_transfer_xaxis_var.get()).strip()
                use_asym = False
                try:
                    use_asym = bool(self.energy_balance_transfer_asym_var.get())
                except Exception:
                    use_asym = False
                transfer_name = "FULLT asym" if use_asym else "FULLT"
                plot_type = "Energy Balance FULLT"
                if xax.strip().lower() == "trace pxqx":
                    display_plot_type = f"Energy balance: {transfer_name} trace ({xax}, Re, source ky={ky}, zonal kx={kx})"
                else:
                    display_plot_type = f"Energy balance: {transfer_name} transfer ({xax}, Re, source ky={ky}, source kx={kx})"
            elif mode_key == "2d":
                plot_type = "Energy Balance vs 2D"
                display_plot_type = "Energy balance: vs 2D"
            else:
                plot_type = "Energy Balance Entropy"
                display_plot_type = "Energy balance: gyro-center entropy balance"
        elif plot_type_selection == "Others":
            other = self.others_plot_var.get().strip().lower()
            if other == "rcorr_phi":
                plot_type = "Radial Correlation (rcorr_phi)"
                display_plot_type = "Others: rcorr_phi"
            elif other == "pod_parity":
                plot_type = "POD Parity"
                display_plot_type = "Others: POD parity"
            else:
                plot_type = "Integration Error"
                display_plot_type = "Others: Error"

        return plot_type_selection, plot_type, display_plot_type

    def _case_selection_lock_enabled(self):
        """Return whether explicit case selection should survive focus changes."""
        var = getattr(self, "lock_case_selection_var", None)
        try:
            return bool(var.get()) if var is not None else False
        except Exception:
            return False

    def _case_listbox_selected_names(self):
        """Return only the selection currently reported by the Tk listbox."""
        try:
            return [
                str(self.case_listbox.get(index))
                for index in self.case_listbox.curselection()
            ]
        except (AttributeError, tk.TclError):
            return []

    def _valid_locked_case_names(self):
        """Return cached locked cases that are still present, in listbox order."""
        locked = {str(name) for name in getattr(self, "_locked_case_names", [])}
        if not locked:
            return []
        try:
            return [
                str(self.case_listbox.get(index))
                for index in range(self.case_listbox.size())
                if str(self.case_listbox.get(index)) in locked
            ]
        except (AttributeError, tk.TclError):
            return []

    def _restore_locked_case_selection(self):
        """Restore the cached locked selection in the visible case list."""
        if not self._case_selection_lock_enabled():
            return
        names = self._valid_locked_case_names()
        if not names:
            return

        selected = set(self._case_listbox_selected_names())
        if selected == set(names):
            return

        self._restoring_case_selection = True
        try:
            self.case_listbox.selection_clear(0, tk.END)
            names_set = set(names)
            for index in range(self.case_listbox.size()):
                if str(self.case_listbox.get(index)) in names_set:
                    self.case_listbox.selection_set(index)
        except (AttributeError, tk.TclError):
            pass
        finally:
            self._restoring_case_selection = False

    def _apply_case_selection_lock_state(self, capture=True):
        """Apply the lock option to Tk and optionally capture current cases."""
        enabled = self._case_selection_lock_enabled()
        try:
            self.case_listbox.configure(exportselection=not enabled)
        except (AttributeError, tk.TclError):
            pass

        if enabled:
            selected = self._case_listbox_selected_names()
            if capture and selected:
                self._locked_case_names = selected
            elif self._valid_locked_case_names():
                self._restore_locked_case_selection()
        elif capture:
            self._locked_case_names = []
        return enabled

    def _on_lock_case_selection_toggle(self):
        """Enable or disable persistent explicit case selection."""
        enabled = self._apply_case_selection_lock_state(capture=True)
        status = "Selection locked" if enabled else "Selection unlocked"
        self._refresh_case_summary(status)

    def _on_case_listbox_select(self, event=None):
        """Track an explicit case selection while the lock option is enabled."""
        if not getattr(self, "_restoring_case_selection", False):
            selected = self._case_listbox_selected_names()
            if self._case_selection_lock_enabled():
                if selected:
                    self._locked_case_names = selected
                elif self._valid_locked_case_names():
                    self._restore_locked_case_selection()
        self.update_options(event)

    def _get_selected_case_names(self):
        """Return explicit/locked cases, or all cases when no choice exists."""
        selected = self._case_listbox_selected_names()
        if selected:
            if self._case_selection_lock_enabled():
                self._locked_case_names = list(selected)
            return selected

        if self._case_selection_lock_enabled():
            locked = self._valid_locked_case_names()
            if locked:
                return locked

        return [self.case_listbox.get(i) for i in range(self.case_listbox.size())]

    def _get_workspace_case_entries(self):
        """Return loaded case entries in current listbox order."""
        entries = []
        for i in range(self.case_listbox.size()):
            case_name = self.case_listbox.get(i)
            data = self.cases.get(case_name, None)
            case_dir = None
            if data is not None:
                # Different pygacode wrappers expose the case directory through
                # different attributes; resolve through the shared helper first,
                # then fall back to raw attributes for older wrapper objects.
                try:
                    case_dir = self._resolve_case_dir(data)
                except Exception:
                    case_dir = None
                if not case_dir:
                    try:
                        case_dir = getattr(data, 'dir', None) or getattr(data, 'path', None)
                    except Exception:
                        case_dir = None
            entries.append({
                "name": str(case_name),
                "path": str(case_dir) if case_dir else "",
            })
        return entries

    def _get_workspace_state(self):
        """Collect Tk variable values that define the current plotting state."""
        state = {}
        for attr in self._WORKSPACE_STATE_VARS:
            var = getattr(self, attr, None)
            if var is None or not hasattr(var, "get"):
                continue
            try:
                state[attr] = var.get()
            except Exception:
                pass
        state["fluc_advanced_case_values"] = {
            str(case_name): {
                "kx": str(values.get("kx", "")),
                "ky": str(values.get("ky", "")),
            }
            for case_name, values in getattr(self, "_fluc_advanced_case_values", {}).items()
            if isinstance(values, dict)
        }
        return state

    def _restore_workspace_state(self, state):
        """Restore saved Tk variable values when corresponding controls exist."""
        if not isinstance(state, dict):
            return
        for attr in self._WORKSPACE_STATE_VARS:
            if attr not in state:
                continue
            var = getattr(self, attr, None)
            if var is None or not hasattr(var, "set"):
                continue
            try:
                # Ignore variables removed by newer versions of the tool; this
                # keeps old workspace JSON files loadable after UI evolution.
                var.set(state[attr])
            except Exception:
                pass
        self._apply_case_selection_lock_state(capture=False)
        raw_advanced = state.get("fluc_advanced_case_values", {})
        if isinstance(raw_advanced, dict):
            self._fluc_advanced_case_values = {
                str(case_name): {
                    "kx": str(values.get("kx", "")),
                    "ky": str(values.get("ky", "")),
                }
                for case_name, values in raw_advanced.items()
                if isinstance(values, dict)
            }
        self._fluc_advanced_seeded = bool(state.get("fluc_advanced_var", False))
        self._ensure_fluc_advanced_case_values()
        if (
            "energy_balance_single_norm_var" not in state
            and bool(state.get("energy_balance_single_norm_entropy_var", False))
            and hasattr(self, "energy_balance_single_norm_var")
        ):
            try:
                self.energy_balance_single_norm_var.set("|min(T)|")
            except Exception:
                pass

    def _default_workspace_dir(self):
        """Return preferred workspace directory, shared with data export when available."""
        try:
            if hasattr(self, "_default_data_export_dir"):
                return self._default_data_export_dir()
        except Exception:
            pass
        return default_share_dir()

    def _has_plot_content(self):
        """Return whether the main axes currently contain rendered content."""
        try:
            return any(
                bool(getattr(self.ax, attr, []))
                for attr in ("lines", "collections", "images", "patches", "texts", "artists")
            )
        except Exception:
            return False

    def _build_workspace_payload(self):
        """Build a serializable snapshot of the current GUI workspace."""
        selected_names = self._case_listbox_selected_names()
        if not selected_names and self._case_selection_lock_enabled():
            selected_names = self._valid_locked_case_names()
        workspace = {
            "version": 1,
            "tool": "cgyro_comparison_tool",
            "cases": self._get_workspace_case_entries(),
            "selected_cases": selected_names,
            "state": self._get_workspace_state(),
            "replot": self._has_plot_content(),
        }
        try:
            workspace["window_geometry"] = self.root.geometry()
        except (AttributeError, tk.TclError):
            pass
        try:
            workspace["pane_position"] = int(self.main_pane.sashpos(0))
        except (AttributeError, tk.TclError, TypeError, ValueError):
            pass
        return workspace

    @staticmethod
    def _write_workspace_file(out_path, workspace):
        """Write a workspace JSON file, creating its parent directory if needed."""
        parent_dir = os.path.dirname(os.path.abspath(out_path))
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(workspace, f, indent=2)

    def _restore_auto_workspace(self):
        """Restore the workspace passed by the update-restarted process."""
        path = self._auto_workspace_path
        self._auto_workspace_path = ""
        if not path:
            return

        try:
            result = self._load_workspace_file(path, show_message=False)
        except Exception as exc:
            print(f"Failed to restore the workspace after update: {exc}")
            result = None
        finally:
            self._remove_auto_workspace(path)

        if result is None:
            self.status_var.set("Updated successfully, but the previous workspace could not be restored.")
            return

        loaded_count = result["loaded_count"]
        failed = result["failed"]
        self._refresh_case_summary("Restored")
        if failed:
            self.status_var.set(
                f"Restored {loaded_count} case(s); skipped {len(failed)} unavailable case(s)."
            )
        else:
            self.status_var.set(f"Restored {loaded_count} case(s) after update.")

        if result["workspace"].get("replot") and loaded_count > 0:
            self.root.after_idle(self._restore_auto_plot)

    def _restore_auto_plot(self):
        """Re-render the last workspace after its cases and options are restored."""
        try:
            self.plot_comparison()
        except Exception as exc:
            print(f"Failed to restore the previous plot after update: {exc}")

    def save_workspace(self):
        """Save current case list, selection, and plotting options to JSON."""
        if not hasattr(self, "cases") or len(self.cases) == 0:
            messagebox.showwarning("Save workspace", "No loaded cases to save.")
            return

        out_path = filedialog.asksaveasfilename(
            title="Save CGYRO comparison workspace",
            defaultextension=".json",
            filetypes=[
                ("Workspace files", "*.json"),
                ("All files", "*.*"),
            ],
            initialdir=self._default_workspace_dir(),
        )
        if not out_path:
            return

        try:
            self._write_workspace_file(out_path, self._build_workspace_payload())
            messagebox.showinfo("Save workspace", f"Workspace saved:\n{out_path}")
        except Exception as e:
            messagebox.showerror("Save workspace", f"Failed to save workspace:\n{e}")

    def _load_workspace_file(self, in_path, show_message=True):
        """Load and apply a workspace file, optionally without user dialogs."""
        try:
            with open(in_path, "r", encoding="utf-8") as f:
                workspace = json.load(f)
        except Exception as e:
            if show_message:
                messagebox.showerror("Load workspace", f"Failed to read workspace:\n{e}")
            return None

        if not isinstance(workspace, dict):
            if show_message:
                messagebox.showerror("Load workspace", "Invalid workspace: expected a JSON object.")
            return None

        cases = workspace.get("cases", [])
        if not isinstance(cases, list):
            if show_message:
                messagebox.showerror("Load workspace", "Invalid workspace: missing case list.")
            return None

        self.cases.clear()
        self._locked_case_names = []
        self._fluc_advanced_case_values.clear()
        self._fluc_advanced_seeded = False
        self.case_listbox.delete(0, tk.END)

        loaded_count = 0
        failed = []
        for entry in cases:
            if isinstance(entry, dict):
                case_dir = entry.get("path", "")
            else:
                # Backward-compatible path-only entries are accepted for
                # hand-written or older workspace files.
                case_dir = str(entry)
            case_dir = str(case_dir).strip()
            if not case_dir:
                failed.append("(empty path)")
                continue
            if self._load_case_from_dir(case_dir, silent=True):
                loaded_count += 1
            else:
                failed.append(case_dir)

        self._update_species_list()
        self._restore_workspace_state(workspace.get("state", {}))
        try:
            self.update_options()
        except Exception:
            pass

        selected = set(str(x) for x in workspace.get("selected_cases", []))
        if selected:
            self.case_listbox.selection_clear(0, tk.END)
            for i in range(self.case_listbox.size()):
                if self.case_listbox.get(i) in selected:
                    self.case_listbox.selection_set(i)
        if self._case_selection_lock_enabled():
            self._locked_case_names = self._case_listbox_selected_names()
            self._restore_locked_case_selection()

        geometry = workspace.get("window_geometry")
        if geometry:
            try:
                self.root.geometry(str(geometry))
            except (tk.TclError, TypeError, ValueError):
                pass
        pane_position = workspace.get("pane_position")
        if pane_position is not None:
            try:
                pane_position = int(pane_position)
                self.root.after_idle(lambda: self.main_pane.sashpos(pane_position))
            except (tk.TclError, TypeError, ValueError):
                pass

        result = {
            "workspace": workspace,
            "loaded_count": loaded_count,
            "failed": failed,
        }
        if show_message:
            msg = f"Loaded {loaded_count} case(s) from workspace."
            if failed:
                msg += "\n\nSkipped:\n- " + "\n- ".join(failed[:20])
                if len(failed) > 20:
                    msg += f"\n- ... ({len(failed) - 20} more)"
            messagebox.showinfo("Load workspace", msg)
        return result

    def load_workspace(self):
        """Load case list, selection, and plotting options from a JSON workspace."""
        in_path = filedialog.askopenfilename(
            title="Load CGYRO comparison workspace",
            filetypes=[
                ("Workspace files", "*.json"),
                ("All files", "*.*"),
            ],
            initialdir=self._default_workspace_dir(),
        )
        if not in_path:
            return
        self._load_workspace_file(in_path, show_message=True)

    @staticmethod
    def _is_fluc2d_kxky_view(view):
        return str(view).strip().lower() == "vs kxky"

    def _is_contour_like_plot(self, plot_type):
        """True when plot type is contour/multi-panel style (single-case rendering)."""
        is_fullt_kxky = False
        if plot_type == "Energy Balance FULLT":
            try:
                fullt_xaxis = str(self.energy_balance_transfer_xaxis_var.get()).strip().lower()
            except Exception:
                fullt_xaxis = "vs kxky"
            is_fullt_kxky = (fullt_xaxis != "vs kx")

        return (
            ("FFT" in plot_type)
            or plot_type.startswith("Fluctuation 2D")
            or plot_type == "POD Parity"
            or is_fullt_kxky
            or (
                plot_type.startswith("Energy Balance Single ")
                and "vs kxky" in plot_type.lower()
            )
            or ("vs kxky" in plot_type.lower())
            or ("vs ky_time" in plot_type)
        )

    def _is_standard_line_plot(self, plot_type):
        """True when plot type is standard line rendering."""
        return not self._is_contour_like_plot(plot_type)

    def _resolve_species_indices(
        self,
        data,
        n_species,
        case_label,
        species_override_index=None,
        fallback_first=True,
        main_ion_policy="all",
        single_species_only=False,
    ):
        """
        Resolve species selection from GUI selection or explicit override.
        Returns (indices, species_label). Empty indices means failure.
        """
        indices = []
        spec_label = ""
        selected_spec_str = ""

        if species_override_index is not None:
            indices = [int(species_override_index)]
            try:
                specs = self._get_case_species(data)
                if 0 <= int(species_override_index) < len(specs):
                    z, m = specs[int(species_override_index)]
                    spec_label = self._get_species_name(z, m)
                else:
                    spec_label = f"s{int(species_override_index)+1}"
            except Exception:
                spec_label = f"s{int(species_override_index)+1}"
        else:
            try:
                selected_spec_str = self.species_var.get().strip()
            except Exception:
                selected_spec_str = ""

            if selected_spec_str.startswith("Main Ion"):
                indices = self._get_main_ion_indices(data)
                if not indices:
                    print(f"Warning: No Main Ion (Z>0 and DENS>0.3) found in {case_label}")
                    return [], "Main Ion"
                if main_ion_policy == "first" and len(indices) > 1:
                    indices = [indices[0]]
                    print("Note: Main Ion sum is reduced to first main-ion species for this plot.")
                spec_label = self._format_species_join_label(data, indices)
            elif selected_spec_str == "All Ions":
                indices = self._get_all_ion_indices(data)
                if not indices:
                    print(f"Warning: No ions (Z>0) found in {case_label}")
                    return [], "All Ions"
                spec_label = "All Ions"
            else:
                match = re.search(r"Z=([-\d\.]+), M=([-\d\.]+)", selected_spec_str)
                if match:
                    target_z = float(match.group(1))
                    target_m = float(match.group(2))
                    current_specs = self._get_case_species(data)
                    for idx, (z, m) in enumerate(current_specs):
                        if abs(z - target_z) < 0.01 and abs(m - target_m) < 0.01:
                            indices = [idx]
                            spec_label = self._get_species_name(z, m)
                            break
                if not indices and fallback_first:
                    indices = [0]
                    spec_label = "s1"

        indices = [i for i in indices if 0 <= int(i) < int(n_species)]
        if not indices:
            print(f"No valid species index for {case_label}")
            return [], spec_label

        if single_species_only and len(indices) > 1:
            print(
                f"Note: multi-species selection is reduced to first species index "
                f"{indices[0]} for {case_label}."
            )
            indices = [indices[0]]
            if selected_spec_str.startswith("Main Ion") or selected_spec_str == "All Ions":
                spec_label = self._format_species_join_label(data, indices)

        return indices, spec_label

    # ------------------------------------------------------------------
    # Case Discovery / Loading
    # ------------------------------------------------------------------

    @staticmethod
    def _get_case_picker_initial_dir():
        """Return preferred initial directory for case selection dialogs."""
        if os.path.exists(DEFAULT_CASE_PICKER_ROOT):
            return DEFAULT_CASE_PICKER_ROOT
        return default_share_dir()

    @staticmethod
    def _looks_like_case_dir(dir_path):
        """Heuristically decide whether a directory contains CGYRO case markers."""
        markers = ("input.cgyro", "out.cgyro.freq", "input.gacode")
        return any(os.path.exists(os.path.join(dir_path, m)) for m in markers)

    @staticmethod
    def _instantiate_cgyrodata_light(dir_path):
        """Instantiate cgyrodata without preloading heavy binary diagnostics."""
        try:
            return cgyrodata(dir_path, fast=True)
        except TypeError:
            return cgyrodata(dir_path)

    def _load_case_from_dir(self, dir_path, silent=False):
        """Load one CGYRO case directory and register it in case list."""
        dir_path = dir_path.replace('\\', '/')
        if not dir_path.endswith('/'):
            dir_path = dir_path + '/'
            
        case_name = os.path.basename(os.path.dirname(dir_path))
        
        # Ensure unique name
        if case_name in self.cases:
            i = 1
            while f"{case_name}_{i}" in self.cases:
                i += 1
            case_name = f"{case_name}_{i}"

        try:
            try:
                data = self._instantiate_cgyrodata_light(dir_path)
            except Exception as e:
                data = cgyrodata_plot(dir_path)
            
            self.cases[case_name] = data
            self.case_listbox.insert(tk.END, case_name)
            self._ensure_fluc_advanced_case_values()
            self._refresh_case_summary("Loaded")
            if not silent:
                print(f"Loaded case: {case_name} from {dir_path}")
            return True
        except Exception as e:
            if not silent:
                print(f"Failed to load case from {dir_path}: {e}")
            return False

    def _load_case_dirs_batch(self, dir_paths):
        """Load multiple directories and return `(loaded_count, failed_dirs)`."""
        loaded_count = 0
        failed_dirs = []
        for d in dir_paths:
            if self._load_case_from_dir(d, silent=True):
                loaded_count += 1
            else:
                failed_dirs.append(d)
        return loaded_count, failed_dirs

    def _select_multiple_case_dirs(self, initial_dir):
        """Select multiple case directories from one parent folder."""
        parent_dir = filedialog.askdirectory(
            title="Select Parent Directory",
            initialdir=initial_dir,
        )
        if not parent_dir:
            return []

        try:
            subdirs_all = [
                os.path.join(parent_dir, name)
                for name in sorted(os.listdir(parent_dir))
                if os.path.isdir(os.path.join(parent_dir, name))
            ]
        except Exception as e:
            messagebox.showerror("Error", f"Failed to list subfolders:\n{e}")
            return []

        # Prefer showing only case-like folders; fallback to all subfolders.
        case_like_subdirs = [d for d in subdirs_all if self._looks_like_case_dir(d)]
        subdirs = case_like_subdirs if case_like_subdirs else subdirs_all

        if not subdirs:
            messagebox.showwarning("Warning", "No subfolders found in selected parent directory.")
            return []

        selected_dirs = []
        dialog = tk.Toplevel(self.root)
        dialog.title("Select Case Folders (Ctrl/Shift multi-select)")
        dialog.transient(self.root)
        dialog.grab_set()
        dialog.geometry("700x500")

        ttk.Label(dialog, text=f"Parent: {parent_dir}").pack(anchor=tk.W, padx=10, pady=(10, 4))
        ttk.Label(dialog, text="Tip: use Ctrl/Shift to select multiple folders.").pack(
            anchor=tk.W, padx=10, pady=(0, 6)
        )

        list_frame = ttk.Frame(dialog)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        yscroll = ttk.Scrollbar(list_frame, orient=tk.VERTICAL)
        xscroll = ttk.Scrollbar(list_frame, orient=tk.HORIZONTAL)
        listbox = tk.Listbox(
            list_frame,
            selectmode=tk.EXTENDED,
            yscrollcommand=yscroll.set,
            xscrollcommand=xscroll.set,
        )
        yscroll.config(command=listbox.yview)
        xscroll.config(command=listbox.xview)
        listbox.grid(row=0, column=0, sticky=tk.NSEW)
        yscroll.grid(row=0, column=1, sticky=tk.NS)
        xscroll.grid(row=1, column=0, sticky=tk.EW)
        list_frame.rowconfigure(0, weight=1)
        list_frame.columnconfigure(0, weight=1)

        for d in subdirs:
            listbox.insert(tk.END, os.path.basename(d))

        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

        def on_select_all():
            """Select all listed directories."""
            listbox.selection_set(0, tk.END)

        def on_ok(event=None):
            """Accept selection and close the chooser."""
            nonlocal selected_dirs
            indices = listbox.curselection()
            selected_dirs = [subdirs[i] for i in indices]
            dialog.destroy()

        def on_cancel(event=None):
            """Cancel selection and close the chooser."""
            nonlocal selected_dirs
            selected_dirs = []
            dialog.destroy()

        ttk.Button(btn_frame, text="Select All", command=on_select_all).pack(side=tk.LEFT)
        ttk.Button(btn_frame, text="Cancel", command=on_cancel).pack(side=tk.RIGHT)
        ttk.Button(btn_frame, text="OK", command=on_ok).pack(side=tk.RIGHT, padx=(0, 8))

        dialog.bind("<Return>", on_ok)
        dialog.bind("<Escape>", on_cancel)
        listbox.focus_set()
        dialog.wait_window()

        return selected_dirs

    def add_case_single(self):
        """Open one-folder chooser and add a single CGYRO case directory."""
        initial_dir = self._get_case_picker_initial_dir()
        dir_path = filedialog.askdirectory(
            title="Select Case Directory (Single)",
            initialdir=initial_dir,
        )
        if not dir_path:
            return

        if self._load_case_from_dir(dir_path, silent=False):
            self._update_species_list()
            messagebox.showinfo("Load Complete", "Loaded 1 case.")
        else:
            messagebox.showerror("Error", "Failed to load selected case directory.")

    def add_case_multiple(self):
        """Open multi-folder chooser and add one or more CGYRO case directories."""
        initial_dir = self._get_case_picker_initial_dir()
        dir_paths = self._select_multiple_case_dirs(initial_dir)
        if not dir_paths:
            return

        loaded_count, failed_dirs = self._load_case_dirs_batch(dir_paths)

        if loaded_count > 0:
            self._update_species_list()
            msg = f"Loaded {loaded_count} case(s)."
            if failed_dirs:
                msg += f"\nFailed: {len(failed_dirs)}"
            messagebox.showinfo("Load Complete", msg)
        else:
            messagebox.showerror("Error", "Failed to load selected case directories.")

    def add_case(self):
        """Backward-compatible alias: keep old call path mapped to multiple-case loading."""
        self.add_case_multiple()

    def add_group(self):
        """Load all valid CGYRO-like subfolders under one parent directory."""
        initial_dir = self._get_case_picker_initial_dir()
        parent_dir = filedialog.askdirectory(title="Select Parent Directory (Load All Subfolders)", initialdir=initial_dir)
        
        if parent_dir:
            candidate_dirs = []
            for item in os.listdir(parent_dir):
                sub_path = os.path.join(parent_dir, item)
                if os.path.isdir(sub_path) and self._looks_like_case_dir(sub_path):
                    candidate_dirs.append(sub_path)

            loaded_count, _failed_dirs = self._load_case_dirs_batch(candidate_dirs)
            
            if loaded_count > 0:
                self._update_species_list()
                messagebox.showinfo("Success", f"Loaded {loaded_count} cases from {parent_dir}")
            else:
                messagebox.showwarning("Warning", "No valid cases found in selected directory.")

    def remove_case(self):
        """Remove currently selected cases from list and internal cache."""
        selected_indices = self.case_listbox.curselection()
        self._restoring_case_selection = True
        try:
            for index in reversed(selected_indices):
                case_name = self.case_listbox.get(index)
                del self.cases[case_name]
                self._fluc_advanced_case_values.pop(str(case_name), None)
                self.case_listbox.delete(index)
        finally:
            self._restoring_case_selection = False
        self._locked_case_names = [
            name for name in self._locked_case_names if name in self.cases
        ]
        self._restore_locked_case_selection()
        self._update_species_list()
        self._refresh_case_summary("Removed")
        if hasattr(self, "_clear_current_plot_data"):
            self._clear_current_plot_data()

    def remove_all_cases(self):
        """Clear all loaded cases after user confirmation."""
        if messagebox.askyesno("Confirm", "Are you sure you want to remove all loaded cases?"):
            self.cases.clear()
            self._locked_case_names = []
            self._fluc_advanced_case_values.clear()
            self.case_listbox.delete(0, tk.END)
            self._update_species_list()
            self._refresh_case_summary("Removed all")
            if hasattr(self, "_clear_current_plot_data"):
                self._clear_current_plot_data()

    def reload_cases(self):
        """Reload all currently registered cases from their original directories."""
        loaded_count = 0
        failed_cases = []
        for case_name, data in list(self.cases.items()):
            # Find the directory
            dir_path = getattr(data, 'dir', None)
            if not dir_path:
                dir_path = getattr(data, 'path', None)
            
            if dir_path:
                try:
                    try:
                        new_data = self._instantiate_cgyrodata_light(dir_path)
                    except Exception as e:
                        new_data = cgyrodata_plot(dir_path)
                    self.cases[case_name] = new_data
                    loaded_count += 1
                except Exception as e:
                    print(f"Failed to reload {case_name} from {dir_path}: {e}")
                    failed_cases.append(case_name)
            else:
                print(f"Could not find directory path for {case_name}")
                failed_cases.append(case_name)
                
        if loaded_count > 0:
            self._update_species_list()
            if hasattr(self, "_clear_current_plot_data"):
                self._clear_current_plot_data()
            message = f"Successfully reloaded {loaded_count} cases."
            if failed_cases:
                message += f"\nFailed to reload: {', '.join(failed_cases)}"
            messagebox.showinfo("Reload Complete", message)
        elif failed_cases:
            messagebox.showwarning("Reload Failed", f"Failed to reload: {', '.join(failed_cases)}")

    def _on_drag_start(self, event):
        """Record drag-start index for case-list reordering."""
        if self.case_listbox.size() > 0:
            self._drag_start_index = self.case_listbox.nearest(event.y)

    def _on_drag_motion(self, event):
        """Reorder case-list entries while dragging with left mouse button."""
        if self.case_listbox.size() == 0:
            return "break"
            
        new_index = self.case_listbox.nearest(event.y)
        if new_index != self._drag_start_index:
            if 0 <= new_index < self.case_listbox.size():
                text = self.case_listbox.get(self._drag_start_index)
                is_selected = self.case_listbox.selection_includes(self._drag_start_index)
                
                self.case_listbox.delete(self._drag_start_index)
                self.case_listbox.insert(new_index, text)
                
                if is_selected:
                    self.case_listbox.selection_set(new_index)
                
                self._drag_start_index = new_index
        return "break"

    def _on_case_listbox_mousewheel(self, event):
        """Scroll only the case list and stop the outer-panel wheel binding."""
        try:
            delta = int(getattr(event, "delta", 0))
        except Exception:
            delta = 0

        if delta != 0:
            step = -1 * int(delta / 120) if abs(delta) >= 120 else (-1 if delta > 0 else 1)
            self.case_listbox.yview_scroll(step, "units")
            return "break"

        num = int(getattr(event, "num", 0))
        if num == 4:
            self.case_listbox.yview_scroll(-1, "units")
        elif num == 5:
            self.case_listbox.yview_scroll(1, "units")
        return "break"

    def _stop_animation(self):
        """Stop active matplotlib animation and reset animation control state."""
        if self.ani:
            self.ani.event_source.stop()
            self.ani = None
        scroll_cid = getattr(self, "_case_info_scroll_cid", None)
        if scroll_cid is not None:
            try:
                self.canvas.mpl_disconnect(scroll_cid)
            except Exception:
                pass
        self._case_info_scroll_cid = None
        self._case_info_scroll_active = False
        self._case_info_lines = []
        self._case_info_line_colors = []
        self._case_info_case_name = ""
        self._case_info_file_path = ""
        self._manual_pager_active = False
        self._manual_pager_label = "Page"
        self.is_paused = False
        self.anim_update_func = None
        self.current_frame = 0
        self.total_frames = 0
        if hasattr(self, 'btn_pause'):
            self.btn_pause.config(text="Pause", state="disabled")
            self.btn_prev.config(state="disabled")
            self.btn_next.config(state="disabled")

    def _on_left_panel_configure(self, _event=None):
        """Update left scroll-region when control panel content changes size."""
        try:
            self.left_canvas.configure(scrollregion=self.left_canvas.bbox("all"))
        except Exception:
            pass

    def _scroll_plot_options_to_top(self):
        """Show the beginning of the dynamic options form after mode changes."""
        try:
            self.left_canvas.yview_moveto(0.0)
        except (AttributeError, tk.TclError):
            pass

    def _on_left_canvas_configure(self, event):
        """Keep embedded left control panel width synced to canvas width."""
        try:
            self.left_canvas.itemconfigure(self._left_panel_window, width=event.width)
            self.left_canvas.configure(scrollregion=self.left_canvas.bbox("all"))
        except Exception:
            pass

    def _on_left_panel_mousewheel(self, event):
        """Mouse-wheel scrolling for left control panel."""
        try:
            delta = int(getattr(event, "delta", 0))
        except Exception:
            delta = 0

        if delta != 0:
            step = -1 * int(delta / 120) if abs(delta) >= 120 else (-1 if delta > 0 else 1)
            self.left_canvas.yview_scroll(step, "units")
            return "break"

        num = int(getattr(event, "num", 0))
        if num == 4:
            self.left_canvas.yview_scroll(-1, "units")
            return "break"
        if num == 5:
            self.left_canvas.yview_scroll(1, "units")
            return "break"
        return None

    def _on_left_panel_enter(self, _event=None):
        """Activate global wheel binding while pointer is over left control area."""
        self.left_canvas.bind_all("<MouseWheel>", self._on_left_panel_mousewheel)
        self.left_canvas.bind_all("<Button-4>", self._on_left_panel_mousewheel)
        self.left_canvas.bind_all("<Button-5>", self._on_left_panel_mousewheel)

    def _on_left_panel_leave(self, _event=None):
        """Release global wheel binding after pointer leaves left control area."""
        try:
            self.left_canvas.unbind_all("<MouseWheel>")
            self.left_canvas.unbind_all("<Button-4>")
            self.left_canvas.unbind_all("<Button-5>")
        except Exception:
            pass

    def _update_manual_pager_status(self):
        """Update pager status text for manual text-page mode."""
        if not getattr(self, "_manual_pager_active", False):
            return
        total = int(getattr(self, "total_frames", 0))
        current = int(getattr(self, "current_frame", 0))
        if total <= 0:
            self.btn_pause.config(text=self._manual_pager_label, state="disabled")
            self.btn_prev.config(state="disabled")
            self.btn_next.config(state="disabled")
            return
        self.btn_pause.config(text=f"{self._manual_pager_label} {current + 1}/{total}", state="disabled")
        if total > 1:
            self.btn_prev.config(state="normal")
            self.btn_next.config(state="normal")
        else:
            self.btn_prev.config(state="disabled")
            self.btn_next.config(state="disabled")

    def _enable_manual_pager(self, total_frames, update_func, start_frame=0, label="Page"):
        """
        Enable Prev/Next stepping for non-animation paged content.
        """
        self._manual_pager_active = True
        self._manual_pager_label = str(label)
        self.ani = None
        self.is_paused = True
        self.anim_update_func = update_func
        self.total_frames = max(0, int(total_frames))
        if self.total_frames > 0:
            self.current_frame = int(start_frame) % self.total_frames
        else:
            self.current_frame = 0
        self._update_manual_pager_status()

    def toggle_pause(self):
        """Toggle animation pause/play state."""
        if getattr(self, "_manual_pager_active", False):
            return
        if self.ani:
            if self.is_paused:
                self.ani.event_source.start()
                self.btn_pause.config(text="Pause")
                self.is_paused = False
            else:
                self.ani.event_source.stop()
                self.btn_pause.config(text="Play")
                self.is_paused = True

    def prev_frame(self):
        """Step to previous animation frame when paused."""
        if getattr(self, "_manual_pager_active", False) and hasattr(self, 'anim_update_func') and self.anim_update_func:
            if int(getattr(self, "total_frames", 0)) <= 0:
                return
            self.current_frame = max(0, self.current_frame - 1)
            self.anim_update_func(self.current_frame)
            self._update_manual_pager_status()
            self.canvas.draw()
            return
        if self.ani and self.is_paused and hasattr(self, 'anim_update_func'):
            self.current_frame = (self.current_frame - 1) % self.total_frames
            self.anim_update_func(self.current_frame)
            self.canvas.draw()

    def next_frame(self):
        """Step to next animation frame when paused."""
        if getattr(self, "_manual_pager_active", False) and hasattr(self, 'anim_update_func') and self.anim_update_func:
            if int(getattr(self, "total_frames", 0)) <= 0:
                return
            self.current_frame = min(self.total_frames - 1, self.current_frame + 1)
            self.anim_update_func(self.current_frame)
            self._update_manual_pager_status()
            self.canvas.draw()
            return
        if self.ani and self.is_paused and hasattr(self, 'anim_update_func'):
            self.current_frame = (self.current_frame + 1) % self.total_frames
            self.anim_update_func(self.current_frame)
            self.canvas.draw()

    def clear_plot(self):
        """Clear current figure and reset to one empty axes."""
        self._reset_plot_area()
        self.canvas.draw()

    # ------------------------------------------------------------------
    # Shared Math / Plot Utilities
    # ------------------------------------------------------------------


