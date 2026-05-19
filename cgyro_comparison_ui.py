"""
UI and case-management mixin for CGYRO comparison GUI.
"""

import os
import re
import getpass
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

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
    )
except ImportError:
    from .cgyro_comparison_bootstrap import (
        cgyrodata,
        cgyrodata_plot,
        DEFAULT_APP_TITLE,
        DEFAULT_WINDOW_GEOMETRY,
        DEFAULT_CASE_PICKER_ROOT,
        DEFAULT_LINEAR_GAMMA_FILE,
    )


class CgyroUiMixin:
    def __init__(self, root):
        """Initialize UI state, widgets, and menus."""
        self.root = root
        self.root.title(DEFAULT_APP_TITLE)
        self.root.geometry(DEFAULT_WINDOW_GEOMETRY)

        self.cases = {}  # Dictionary to store loaded cases: {name: cgyrodata_object}
        self.ani = None # Animation object
        self.is_paused = False
        self.current_frame = 0
        self.total_frames = 0
        self._drag_start_index = None
        self._manual_pager_active = False
        self._manual_pager_label = "Page"
        
        self._create_layout()
        self._create_menu()

    def _create_layout(self):
        """Create top-level panels, control widgets, and matplotlib canvas."""
        # Main container
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Left panel (scrollable): Controls and Case List
        left_container = ttk.Frame(main_frame, width=420)
        left_container.pack(side=tk.LEFT, fill=tk.Y)
        left_container.pack_propagate(False)

        self.left_scrollbar = ttk.Scrollbar(left_container, orient=tk.VERTICAL)
        self.left_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.left_canvas = tk.Canvas(left_container, highlightthickness=0, borderwidth=0)
        self.left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.left_canvas.configure(yscrollcommand=self.left_scrollbar.set)
        self.left_scrollbar.configure(command=self.left_canvas.yview)

        left_panel = ttk.Frame(self.left_canvas, padding=10)
        self._left_panel_window = self.left_canvas.create_window((0, 0), window=left_panel, anchor=tk.NW)
        left_panel.bind("<Configure>", self._on_left_panel_configure)
        self.left_canvas.bind("<Configure>", self._on_left_canvas_configure)
        self.left_canvas.bind("<Enter>", self._on_left_panel_enter)
        self.left_canvas.bind("<Leave>", self._on_left_panel_leave)
        left_panel.bind("<Enter>", self._on_left_panel_enter)
        left_panel.bind("<Leave>", self._on_left_panel_leave)
        
        # Case List
        ttk.Label(left_panel, text="Loaded Cases:").pack(anchor=tk.W)
        self.case_listbox = tk.Listbox(left_panel, selectmode=tk.EXTENDED, height=10)
        self.case_listbox.pack(fill=tk.X, pady=5)
        
        # Enable drag and drop reordering
        self.case_listbox.bind('<Button-1>', self._on_drag_start)
        self.case_listbox.bind('<B1-Motion>', self._on_drag_motion)
        self.case_listbox.bind('<<ListboxSelect>>', self.update_options)
        
        btn_frame_load = ttk.Frame(left_panel)
        btn_frame_load.pack(fill=tk.X, pady=5)
        btn_frame_load.columnconfigure(0, weight=1)
        btn_frame_load.columnconfigure(1, weight=1)
        ttk.Button(
            btn_frame_load, text="Add Case (Single)", command=self.add_case_single
        ).grid(row=0, column=0, padx=2, pady=2, sticky=tk.EW)
        ttk.Button(
            btn_frame_load, text="Add Case (Multiple)", command=self.add_case_multiple
        ).grid(row=0, column=1, padx=2, pady=2, sticky=tk.EW)
        ttk.Button(btn_frame_load, text="Add Group", command=self.add_group).grid(
            row=1, column=0, columnspan=2, padx=2, pady=2, sticky=tk.EW
        )

        btn_frame_manage = ttk.Frame(left_panel)
        btn_frame_manage.pack(fill=tk.X, pady=2)
        btn_frame_manage.columnconfigure(0, weight=1)
        btn_frame_manage.columnconfigure(1, weight=1)
        ttk.Button(btn_frame_manage, text="Remove", command=self.remove_case).grid(
            row=0, column=0, padx=2, pady=2, sticky=tk.EW
        )
        ttk.Button(btn_frame_manage, text="Remove All", command=self.remove_all_cases).grid(
            row=0, column=1, padx=2, pady=2, sticky=tk.EW
        )
        ttk.Button(btn_frame_manage, text="Reload", command=self.reload_cases).grid(
            row=1, column=0, columnspan=2, padx=2, pady=2, sticky=tk.EW
        )

        # Plot Controls
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        ttk.Label(left_panel, text="Plot Type:").pack(anchor=tk.W)
        
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
        self.plot_type_combo = ttk.Combobox(left_panel, textvariable=self.plot_type_var, values=plot_types, state="readonly")
        self.plot_type_combo.pack(fill=tk.X, pady=5)
        self.plot_type_combo.bind("<<ComboboxSelected>>", self.update_options)

        # Dynamic Options Frame
        self.options_frame = ttk.LabelFrame(left_panel, text="Options", padding=5)
        self.options_frame.pack(fill=tk.X, pady=10)
        self._init_options()

        # Action Buttons
        ttk.Button(left_panel, text="Plot", command=self.plot_comparison).pack(fill=tk.X, pady=10)
        ttk.Button(left_panel, text="Check out.cgyro.info", command=self.plot_case_info).pack(fill=tk.X, pady=(0, 6))
        ttk.Button(left_panel, text="Diff input.cgyro", command=self.plot_input_diff).pack(fill=tk.X, pady=(0, 6))
        ttk.Button(left_panel, text="Clear Plot", command=self.clear_plot).pack(fill=tk.X)

        # Animation Controls
        self.anim_controls_frame = ttk.Frame(left_panel)
        self.anim_controls_frame.pack(fill=tk.X, pady=10)
        
        self.btn_prev = ttk.Button(self.anim_controls_frame, text="< Prev", command=self.prev_frame, state="disabled")
        self.btn_prev.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        self.btn_pause = ttk.Button(self.anim_controls_frame, text="Pause", command=self.toggle_pause, state="disabled")
        self.btn_pause.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        self.btn_next = ttk.Button(self.anim_controls_frame, text="Next >", command=self.next_frame, state="disabled")
        self.btn_next.pack(side=tk.LEFT, fill=tk.X, expand=True)

        # Right panel: Plot Area
        right_panel = ttk.Frame(main_frame, padding=10)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
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

    def _create_menu(self):
        """Create application menu bar and file actions."""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Add Case (Single)", command=self.add_case_single)
        file_menu.add_command(label="Add Case (Multiple)", command=self.add_case_multiple)
        file_menu.add_command(label="Add Group", command=self.add_group)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)

        data_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Data", menu=data_menu)
        data_menu.add_command(
            label="transfer bin to readable",
            command=self.transfer_bin_to_readable,
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
        self.flux_type_combo = ttk.Combobox(self.options_frame, textvariable=self.flux_type_var, values=["Energy", "Particle"], state="readonly", width=15)
        self.flux_xaxis_var = tk.StringVar(value="v.s ky")
        self.flux_xaxis_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.flux_xaxis_var,
            values=["v.s ky", "v.s kx (estimated)", "v.s Time", "v.s ky_time", "v.s 2D"],
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
        self.flux_norm_real_ion_var = tk.BooleanVar(value=False)
        self.flux_norm_real_ion_check = ttk.Checkbutton(
            self.options_frame,
            text="Normalized by real ion (max-density ion)",
            variable=self.flux_norm_real_ion_var,
        )
        self.flux_formula_frame = ttk.Frame(self.options_frame)
        self.flux_formula_fig = plt.Figure(figsize=(3.2, 2.4), dpi=100)
        self.flux_formula_ax = self.flux_formula_fig.add_subplot(111)
        self.flux_formula_ax.axis("off")
        self.flux_formula_fig.patch.set_facecolor("white")
        self.flux_formula_canvas = FigureCanvasTkAgg(self.flux_formula_fig, master=self.flux_formula_frame)
        self.flux_formula_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # 3. Fluctuation 1D Options
        self.fluc_field_var = tk.StringVar(value="Phi")
        self.fluc_field_combo = ttk.Combobox(self.options_frame, textvariable=self.fluc_field_var, values=["Phi", "Apar", "Bpar"], state="readonly", width=10)
        
        self.fluc_xaxis_var = tk.StringVar(value="v.s ky")
        self.fluc_xaxis_combo = ttk.Combobox(self.options_frame, textvariable=self.fluc_xaxis_var, values=["v.s ky", "v.s kx", "v.s Time", "fft"], state="readonly", width=15)

        # 4. Species Selection (Flux, Fluctuation 2D)
        self.species_label = ttk.Label(self.options_frame, text="Species:")
        self.species_var = tk.StringVar()
        self.species_combo = ttk.Combobox(self.options_frame, textvariable=self.species_var, width=20, state="readonly")
        
        self.plot_all_species_var = tk.BooleanVar(value=False)
        self.plot_all_species_check = ttk.Checkbutton(self.options_frame, text="Plot All Species (First Case Only)", variable=self.plot_all_species_var)

        # 5. Fluctuation 2D options
        self.moment_label = ttk.Label(self.options_frame, text="Moment:")
        self.moment_var = tk.StringVar(value="Phi")
        self.moment_combo = ttk.Combobox(self.options_frame, textvariable=self.moment_var,
                                         values=["Phi", "Density", "Energy", "Temperature", "Apar", "Bpar"], state="readonly")

        self.fluc2d_view_label = ttk.Label(self.options_frame, text="View:")
        self.fluc2d_view_var = tk.StringVar(value="vs xy")
        self.fluc2d_view_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.fluc2d_view_var,
            values=["vs xy", "vs xt"],
            state="readonly",
            width=10
        )
        self.fluc2d_x_elec_var = tk.BooleanVar(value=False)
        self.fluc2d_x_elec_check = ttk.Checkbutton(
            self.options_frame,
            text=r"Spatial axes normalize to electron scale (x,y)/\rho_e",
            variable=self.fluc2d_x_elec_var,
        )

        # 6. FFT Options (reused from before, but managed dynamically)
        self.fft_options_frame = ttk.LabelFrame(self.options_frame, text="FFT Settings", padding=5)
        # Shared linear spectrum file path (used by both FFT overlay and Zonal-vs-Gamma_Linear)
        self.linear_gamma_file_var = tk.StringVar(value=DEFAULT_LINEAR_GAMMA_FILE)
        # Analysis Mode
        ttk.Label(self.fft_options_frame, text="Mode:").grid(row=0, column=0, sticky=tk.W)
        self.fft_mode_var = tk.StringVar(value="Nonlinear")
        self.fft_mode_combo = ttk.Combobox(self.fft_options_frame, textvariable=self.fft_mode_var,
                                         values=["Nonlinear", "Linear"], state="readonly", width=15)
        self.fft_mode_combo.grid(row=0, column=1, sticky=tk.W)
        # View Mode
        ttk.Label(self.fft_options_frame, text="View:").grid(row=1, column=0, sticky=tk.W)
        self.fft_view_var = tk.StringVar(value="Omega vs ky")
        self.fft_view_combo = ttk.Combobox(self.fft_options_frame, textvariable=self.fft_view_var,
                                         values=["Omega vs ky", "Omega vs kx"], state="readonly", width=15)
        self.fft_view_combo.grid(row=1, column=1, sticky=tk.W)
        # Spectrum type
        ttk.Label(self.fft_options_frame, text="Spectrum:").grid(row=2, column=0, sticky=tk.W)
        self.fft_spectrum_var = tk.StringVar(value="Amplitude")
        self.fft_spectrum_combo = ttk.Combobox(
            self.fft_options_frame,
            textvariable=self.fft_spectrum_var,
            values=["Amplitude", "Power"],
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
        self.fft_formula_frame = ttk.Frame(self.options_frame)
        self.fft_formula_fig = plt.Figure(figsize=(3.2, 2.0), dpi=100)
        self.fft_formula_ax = self.fft_formula_fig.add_subplot(111)
        self.fft_formula_ax.axis("off")
        self.fft_formula_fig.patch.set_facecolor("white")
        self.fft_formula_canvas = FigureCanvasTkAgg(self.fft_formula_fig, master=self.fft_formula_frame)
        self.fft_formula_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # 7. Zonal ExB shearing options
        self.zf_xaxis_label = ttk.Label(self.options_frame, text="Mode:")
        self.zf_xaxis_var = tk.StringVar(value="vs Time")
        self.zf_xaxis_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.zf_xaxis_var,
            values=["vs Time", "vs kx", "vs gamma_lin"],
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
        self.zf_formula_frame = ttk.Frame(self.options_frame)
        self.zf_formula_frame.columnconfigure(0, weight=1)
        self.zf_formula_frame.rowconfigure(0, weight=1)
        self.zf_formula_fig = plt.Figure(figsize=(3.2, 1.8), dpi=100)
        self.zf_formula_ax = self.zf_formula_fig.add_subplot(111)
        self.zf_formula_ax.axis("off")
        self.zf_formula_fig.patch.set_facecolor("white")
        self.zf_formula_canvas = FigureCanvasTkAgg(self.zf_formula_fig, master=self.zf_formula_frame)
        self.zf_formula_widget = self.zf_formula_canvas.get_tk_widget()
        self.zf_formula_widget.grid(row=0, column=0, sticky=tk.NSEW, padx=(0, 2), pady=(0, 2))
        self.zf_formula_vscroll = ttk.Scrollbar(
            self.zf_formula_frame,
            orient=tk.VERTICAL,
            command=self.zf_formula_widget.yview,
        )
        self.zf_formula_vscroll.grid(row=0, column=1, sticky=tk.NS)
        self.zf_formula_hscroll = ttk.Scrollbar(
            self.zf_formula_frame,
            orient=tk.HORIZONTAL,
            command=self.zf_formula_widget.xview,
        )
        self.zf_formula_hscroll.grid(row=1, column=0, sticky=tk.EW, padx=(0, 2))
        self.zf_formula_widget.configure(
            xscrollcommand=self.zf_formula_hscroll.set,
            yscrollcommand=self.zf_formula_vscroll.set,
        )

        # 8. Energy-balance options (triad_v2-based)
        self.energy_balance_mode_label = ttk.Label(self.options_frame, text="Mode:")
        self.energy_balance_mode_var = tk.StringVar(value="gyro-center entropy balance")
        self.energy_balance_mode_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.energy_balance_mode_var,
            values=[
                "gyro-center entropy balance",
                "ZF energy balance",
                r"vs $\gamma_{eff}^Z$ and $\gamma_{eff}^{NZ}$",
                "Single plot",
                "v.s 2D",
            ],
            state="readonly",
            width=28,
        )
        self.energy_balance_n_label = ttk.Label(self.options_frame, text="n index:")
        self.energy_balance_n_var = tk.StringVar(value="0")
        self.energy_balance_n_entry = ttk.Entry(
            self.options_frame,
            textvariable=self.energy_balance_n_var,
            width=10,
        )
        self.energy_balance_spec_label = ttk.Label(self.options_frame, text="Species:")
        self.energy_balance_spec_var = tk.StringVar(value="Total (-1)")
        self.energy_balance_spec_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.energy_balance_spec_var,
            values=["Total (-1)", "Main ion (0)", "Electron (1)"],
            state="readonly",
            width=16,
        )
        self.energy_balance_single_quantity_label = ttk.Label(self.options_frame, text="Quantity:")
        self.energy_balance_single_quantity_var = tk.StringVar(value="T-N")
        self.energy_balance_single_quantity_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.energy_balance_single_quantity_var,
            values=["T", "N", "T-N", "entropy"],
            state="readonly",
            width=12,
        )
        self.energy_balance_single_xaxis_label = ttk.Label(self.options_frame, text="X-axis:")
        self.energy_balance_single_xaxis_var = tk.StringVar(value="vs Time")
        self.energy_balance_single_xaxis_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.energy_balance_single_xaxis_var,
            values=["vs Time", "vs ky"],
            state="readonly",
            width=12,
        )
        self.energy_balance_formula_frame = ttk.Frame(self.options_frame)
        self.energy_balance_formula_frame.columnconfigure(0, weight=1)
        self.energy_balance_formula_frame.rowconfigure(0, weight=1)
        self.energy_balance_formula_fig = plt.Figure(figsize=(4.4, 4.0), dpi=100)
        self.energy_balance_formula_ax = self.energy_balance_formula_fig.add_subplot(111)
        self.energy_balance_formula_ax.axis("off")
        self.energy_balance_formula_fig.patch.set_facecolor("white")
        self.energy_balance_formula_canvas = FigureCanvasTkAgg(
            self.energy_balance_formula_fig,
            master=self.energy_balance_formula_frame,
        )
        self.energy_balance_formula_widget = self.energy_balance_formula_canvas.get_tk_widget()
        self.energy_balance_formula_widget.grid(
            row=0, column=0, sticky=tk.NSEW, padx=(0, 2), pady=(0, 2)
        )
        self.energy_balance_formula_vscroll = ttk.Scrollbar(
            self.energy_balance_formula_frame,
            orient=tk.VERTICAL,
            command=self.energy_balance_formula_widget.yview,
        )
        self.energy_balance_formula_vscroll.grid(row=0, column=1, sticky=tk.NS)
        self.energy_balance_formula_hscroll = ttk.Scrollbar(
            self.energy_balance_formula_frame,
            orient=tk.HORIZONTAL,
            command=self.energy_balance_formula_widget.xview,
        )
        self.energy_balance_formula_hscroll.grid(row=1, column=0, sticky=tk.EW, padx=(0, 2))
        self.energy_balance_formula_widget.configure(
            xscrollcommand=self.energy_balance_formula_hscroll.set,
            yscrollcommand=self.energy_balance_formula_vscroll.set,
        )

        # 9. Others options
        self.others_plot_label = ttk.Label(self.options_frame, text="Others:")
        self.others_plot_var = tk.StringVar(value="Error")
        self.others_plot_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.others_plot_var,
            values=["Error", "rcorr_phi", "POD_parity"],
            state="readonly",
            width=15,
        )
        self.others_rcorr_field_label = ttk.Label(self.options_frame, text="Field:")
        self.others_rcorr_field_var = tk.StringVar(value="Phi")
        self.others_rcorr_field_combo = ttk.Combobox(
            self.options_frame,
            textvariable=self.others_rcorr_field_var,
            values=["Phi", "Apar", "Bpar"],
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
            values=["Apar", "Phi"],
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
        self.pod_formula_frame = ttk.Frame(self.options_frame)
        self.pod_formula_fig = plt.Figure(figsize=(3.2, 2.2), dpi=100)
        self.pod_formula_ax = self.pod_formula_fig.add_subplot(111)
        self.pod_formula_ax.axis("off")
        self.pod_formula_fig.patch.set_facecolor("white")
        self.pod_formula_canvas = FigureCanvasTkAgg(self.pod_formula_fig, master=self.pod_formula_frame)
        self.pod_formula_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Bind once; avoid repeated callback registration inside update_options.
        self.flux_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.flux_type_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.fluc_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.moment_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.zf_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.energy_balance_mode_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.energy_balance_single_quantity_combo.bind("<<ComboboxSelected>>", self.update_options)
        self.energy_balance_single_xaxis_combo.bind("<<ComboboxSelected>>", self.update_options)
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
        self.t_start_var.set("")
        self.t_end_var.set("")

    def _browse_linear_gamma_file(self):
        """Open file chooser and update the linear omega/gamma reference file path."""
        user_name = (
            os.environ.get("USER", "").strip()
            or os.environ.get("USERNAME", "").strip()
            or getpass.getuser().strip()
        )
        user_share_dir = os.path.join("/data/share", user_name) if user_name else "/data/share"
        if os.path.isdir(user_share_dir):
            initial_dir = user_share_dir
        elif os.path.isdir("/data/share"):
            initial_dir = "/data/share"
        else:
            initial_dir = os.getcwd()
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

    def _hide_dynamic_options(self):
        """Hide all dynamic option widgets before re-laying out current mode controls."""
        widgets = [
            self.norm_ky_check,
            self.flux_type_combo, self.flux_xaxis_combo, self.flux_decomp_check,
            self.flux_norm_real_ion_check,
            self.flux_scan_xparam_label, self.flux_scan_xparam_combo,
            self.flux_formula_frame,
            self.fluc_field_combo, self.fluc_xaxis_combo,
            self.species_label, self.species_combo, self.plot_all_species_check,
            self.fluc2d_view_label, self.fluc2d_view_combo,
            self.fluc2d_x_elec_check,
            self.moment_label, self.moment_combo,
            self.fft_options_frame,
            self.fft_formula_frame,
            self.zf_xaxis_label, self.zf_xaxis_combo,
            self.linear_gamma_file_label, self.linear_gamma_file_entry, self.linear_gamma_file_browse,
            self.zf_gamma_lin_ky_label, self.zf_gamma_lin_ky_entry,
            self.zf_formula_frame,
            self.energy_balance_mode_label, self.energy_balance_mode_combo,
            self.energy_balance_n_label, self.energy_balance_n_entry,
            self.energy_balance_spec_label, self.energy_balance_spec_combo,
            self.energy_balance_single_quantity_label, self.energy_balance_single_quantity_combo,
            self.energy_balance_single_xaxis_label, self.energy_balance_single_xaxis_combo,
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
            r"Flux vs $k_x$ (estimated)",
            r"$\Gamma_{\mathrm{ES}}(k_x,t)=\left\langle \sum_{k_y}"
            r"\Re\!\left[m_{\mathrm{ES}}^*(i k_y)\!\left(-\phi\right)\right]\right\rangle_{\theta}$",
            r"$\Gamma_{\mathrm{EM},A}(k_x,t)=\left\langle \sum_{k_y}"
            r"\Re\!\left[m_A^*(i k_y)\!\left(A_{\parallel}\right)\right]\right\rangle_{\theta}$",
            r"$\Gamma_{\mathrm{EM},B}(k_x,t)=\left\langle \sum_{k_y}"
            r"\Re\!\left[m_B^*(i k_y)B_{\parallel}\right]\right\rangle_{\theta}$",
            r"$\Gamma_{\mathrm{tot}}=\Gamma_{\mathrm{ES}}+\Gamma_{\mathrm{EM},A}+\Gamma_{\mathrm{EM},B}$",
            r"$\bar{\Gamma}(k_x)=\left\langle \Gamma_{\mathrm{tot}}(k_x,t)"
            r"\right\rangle_{t\in[t_1,t_2]},\;\;k_x\geq 0$",
            moment_note,
            norm_note,
        ]

        ax = self.flux_formula_ax
        ax.clear()
        ax.axis("off")
        y = 0.98
        dy = 0.14
        for i, line in enumerate(lines):
            fs = 9 if i == 0 else 8.5
            ax.text(0.01, y, line, transform=ax.transAxes, ha="left", va="top", fontsize=fs)
            y -= dy
        self.flux_formula_canvas.draw_idle()

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

        ax = self.fft_formula_ax
        ax.clear()
        ax.axis("off")
        y = 0.98
        dy = 0.16
        for i, line in enumerate(lines):
            fs = 9 if i == 0 else 8.5
            ax.text(0.01, y, line, transform=ax.transAxes, ha="left", va="top", fontsize=fs)
            y -= dy
        self.fft_formula_canvas.draw_idle()

    def _render_zf_formula_math(self, mode_text):
        """Render math notes for selected zonal ExB shearing sub-mode."""
        mode = (mode_text or "").strip().lower()
        if mode == "vs kx":
            lines = [
                r"Zonal ExB Shearing Rate (vs $k_x$)",
                r"$\phi_{ZF}(k_x,t)=\langle \phi(k_x,k_y\!\approx\!0,\theta,t)\rangle_{\theta}$",
                r"$\Omega(k_x,t)=-k_x^2\,\phi_{ZF}(k_x,t)$",
                r"$\bar{\omega}_{ZF}(k_x)=\left\langle |\Omega(k_x,t)|\right\rangle_{t\in[t_1,t_2]}$",
                r"Plotted on non-negative branch: $k_x\geq 0$.",
            ]
        elif mode == "vs gamma_lin":
            lines = [
                r"vs $\gamma_{lin}$ (set $k_x=k_y$)",
                r"Curves: $\gamma_{lin}(k_y)$, $\langle \omega_{ZF}(k_x)\rangle$, $k_x\langle V_{ZF}\rangle$",
                r"$\phi_{ZF}\!\rightarrow\!\Omega(k_x,t)=-k_x^2\phi_{ZF}(k_x,t)$",
                r"$\omega_{ZF}(k_x,t)=\left|(k_x\rho_D)^2\,\delta\bar{\phi}_{k_x,k_\theta=0}(t)\right|$",
                r"$V_{ZF}(t)=0.5\sqrt{\sum_{k_x}\left|k_x\rho_D\bar{\phi}_{k_x,k_\theta=0}(t)\right|^2}$",
                r"Optional input: $k_y^\star$ for ratios $\omega_{ZF}(k_y^\star)/\gamma_{lin}(k_y^\star)$",
                r"and $k_xV_{ZF}(k_y^\star)/\gamma_{lin}(k_y^\star)$.",
                r"Time-average in selected window; compare on shared $k\rho_s$ axis ($k_x=k_y$).",
            ]
        else:
            lines = [
                r"Zonal ExB Shearing Rate (vs Time)",
                r"$\phi_{ZF}(k_x,t)=\langle \phi(k_x,k_y\!\approx\!0,\theta,t)\rangle_{\theta}$",
                r"$\Omega(k_x,t)=-k_x^2\,\phi_{ZF}(k_x,t)$",
                r"$\omega_{ZF}(t)=\left(\sum_{k_x}|\Omega(k_x,t)|^2\right)^{1/2}$",
            ]

        ax = self.zf_formula_ax
        ax.clear()
        ax.axis("off")

        # Keep long FIG4 notes readable in narrow option panes.
        fig_h = 1.8
        y = 0.96
        dy = 0.18
        title_fs = 9
        line_fs = 8.5
        if mode == "vs gamma_lin":
            # Keep FIG4 notes compact to prevent bottom clipping in small windows.
            fig_h = 2.25
            dy = 0.13
            title_fs = 8.4
            line_fs = 7.7
        try:
            self.zf_formula_fig.set_size_inches(3.2, fig_h, forward=True)
            self.zf_formula_fig.subplots_adjust(left=0.02, right=0.97, top=0.97, bottom=0.08)
        except Exception:
            pass

        # Keep the last line visible regardless of line count / math-text size.
        if len(lines) > 1:
            dy_auto = min(dy, 0.80 / float(len(lines) - 1))
        else:
            dy_auto = dy

        for i, line in enumerate(lines):
            fs = title_fs if i == 0 else line_fs
            ax.text(0.01, y, line, transform=ax.transAxes, ha="left", va="top", fontsize=fs)
            y -= dy_auto
        try:
            self.zf_formula_widget.configure(scrollregion=self.zf_formula_widget.bbox("all"))
            self.zf_formula_widget.xview_moveto(0.0)
            self.zf_formula_widget.yview_moveto(0.0)
        except Exception:
            pass
        self.zf_formula_canvas.draw_idle()

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
        ax = self.pod_formula_ax
        ax.clear()
        ax.axis("off")
        y = 0.97
        for i, line in enumerate(lines):
            if i == 0:
                fs = 9.0
                dy_line = 0.12
            elif r"\dfrac" in line or r"\frac" in line:
                # Fraction lines are taller; reserve extra vertical spacing.
                fs = 7.4
                dy_line = 0.18
            else:
                fs = 8.0
                dy_line = 0.11
            ax.text(0.01, y, line, transform=ax.transAxes, ha="left", va="top", fontsize=fs)
            y -= dy_line
        self.pod_formula_canvas.draw_idle()

    def _render_energy_balance_formula_math(self):
        """Render math notes for Energy-balance sub-modes (triad_v2 style)."""
        mode = str(self.energy_balance_mode_var.get()).strip().lower()
        is_gamma_eff_mode = (r"\gamma_{eff}" in mode) or ("gamma_eff" in mode)
        is_energy_2d_mode = (mode == "v.s 2d")
        if mode == "zf energy balance":
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
                r"$\mathrm{idx3}(a,k_x,n,t)=\dfrac{d\!\left(T_a\,\delta S_{a,k_\perp}\right)}{dt}$",
                r"$\dfrac{dS}{dt}(t)=\sum_{a}\sum_{k_x}\sum_{n\neq 0}\Re\{f[a,k_x,\mathrm{idx3},n,t]\}$",
                r"$S(t)=\int^t \dfrac{dS}{dt}(t')\,dt'\ \propto\ \sum_{a}\sum_{k_x}\sum_{n\neq 0}T_a\delta S_{a,k_\perp}$",
                r"$\left(T_{a}^{NZ\rightarrow Z}/S\right)=\langle T_{a}^{NZ\rightarrow Z}\rangle_t/\langle S\rangle_t$",
                r"$\left(N_{a}^{NZ\rightarrow Z}/S\right)=\langle N_{a}^{NZ\rightarrow Z}\rangle_t/\langle S\rangle_t$",
                r"$\mathrm{Plotted:}\ \mathcal{N}_D/S,\ \mathcal{T}_D/S,\ \mathcal{N}_e/S,\ \mathcal{T}_e/S$",
            ]
        elif mode == "single plot":
            lines = [
                r"$\bf{Energy\ balance\ single\ plot}$",
                r"$f[s,r,c,n,t]=f_{\Re}+i\,f_{\Im}$",
                r"$T_0(t)=\sum_r \Re\{f[r,0,n,t]\},\ \ N_0(t)=\sum_r \Re\{f[r,1,n,t]\}$",
                r"$\delta S_a(k_x,k_y,t)\equiv \Re\{f[a,k_x,\mathrm{idx5},k_y,t]\}$",
                r"$\mathrm{entropy}_{a}(k_y)=\log\!\left(\sum_{k_x}\left\langle \delta S_a(k_x,k_y,t)\right\rangle_t\right)$",
                r"$\mathrm{(for\ vs\ ky,\ entropy\ uses\ idx5,\ not\ idx3)}$",
                r"$\mathrm{Quantity} \in \{T,\ N,\ T\!-\!N,\ entropy\}$",
                r"$\mathrm{X\!-\!axis} \in \{t,\ k_y\}$;  k_y\ \mathrm{mode\ uses\ time\ average\ over\ selected\ window.}$",
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
        ax = self.energy_balance_formula_ax
        ax.clear()
        ax.axis("off")
        y = 0.98
        for i, line in enumerate(lines):
            if i == 0:
                fs = 8.8
                dy_line = 0.115
            elif r"\frac" in line or r"\sum" in line or r"\left(" in line:
                fs = 7.2
                dy_line = 0.125
            else:
                fs = 7.8
                dy_line = 0.105
            ax.text(0.01, y, line, transform=ax.transAxes, ha="left", va="top", fontsize=fs)
            y -= dy_line
        try:
            fig_w_px = int(self.energy_balance_formula_fig.get_figwidth() * self.energy_balance_formula_fig.dpi)
            fig_h_px = int(self.energy_balance_formula_fig.get_figheight() * self.energy_balance_formula_fig.dpi)
            self.energy_balance_formula_widget.configure(scrollregion=(0, 0, fig_w_px, fig_h_px))
            self.energy_balance_formula_widget.xview_moveto(0.0)
            self.energy_balance_formula_widget.yview_moveto(0.0)
        except Exception:
            pass
        self.energy_balance_formula_canvas.draw_idle()

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
        self._hide_dynamic_options()
        plot_type = self.plot_type_var.get()
        row = 3 # Start gridding from row 3

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
            
            # Check if FFT is selected in the sub-option
            if self.fluc_xaxis_var.get() == "fft":
                 self.fft_options_frame.grid(row=row, column=0, columnspan=2, sticky=tk.W+tk.E, pady=5)
                 row += 1
                 self._render_fft_formula_math()
                 self.fft_formula_frame.grid(row=row, column=0, columnspan=4, sticky=tk.W + tk.E, pady=(4, 0))

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
            self.energy_balance_mode_label.grid(row=row, column=0, sticky=tk.W)
            self.energy_balance_mode_combo.grid(row=row, column=1, columnspan=2, sticky=tk.W + tk.E)
            row += 1
            mode = self.energy_balance_mode_var.get().strip().lower()
            is_gamma_eff_mode = (r"\gamma_{eff}" in mode) or ("gamma_eff" in mode)
            is_energy_2d_mode = (mode == "v.s 2d")
            is_single_plot_mode = (mode == "single plot")
            xaxis_single = str(self.energy_balance_single_xaxis_var.get()).strip().lower()
            if is_gamma_eff_mode:
                n_label_txt = "ky:"
            elif is_single_plot_mode and xaxis_single == "vs ky":
                n_label_txt = "n index:"
            else:
                n_label_txt = "n index:"
            self.energy_balance_n_label.configure(text=n_label_txt)
            self.energy_balance_n_label.grid(row=row, column=0, sticky=tk.W)
            self.energy_balance_n_entry.grid(row=row, column=1, sticky=tk.W)
            self.energy_balance_spec_label.grid(row=row, column=2, sticky=tk.W, padx=(10, 0))
            self.energy_balance_spec_combo.grid(row=row, column=3, sticky=tk.W)
            row += 1
            if is_single_plot_mode:
                self.energy_balance_single_quantity_label.grid(row=row, column=0, sticky=tk.W)
                self.energy_balance_single_quantity_combo.grid(row=row, column=1, sticky=tk.W)
                self.energy_balance_single_xaxis_label.grid(row=row, column=2, sticky=tk.W, padx=(10, 0))
                self.energy_balance_single_xaxis_combo.grid(row=row, column=3, sticky=tk.W)
                row += 1
            if is_gamma_eff_mode:
                self.linear_gamma_file_label.grid(row=row, column=0, sticky=tk.W)
                self.linear_gamma_file_entry.grid(row=row, column=1, columnspan=2, sticky=tk.W + tk.E)
                self.linear_gamma_file_browse.grid(row=row, column=3, sticky=tk.W)
                row += 1
            if is_energy_2d_mode:
                self._refresh_flux_2d_param_choices()
                self.flux_scan_xparam_label.grid(row=row, column=0, sticky=tk.W)
                self.flux_scan_xparam_combo.grid(row=row, column=1, sticky=tk.W + tk.E)
                row += 1
            self._render_energy_balance_formula_math()
            self.energy_balance_formula_frame.grid(
                row=row, column=0, columnspan=4, sticky=tk.W + tk.E, pady=(4, 0)
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
            
            # Check for Main Ion (D or T)
            d_or_t = False
            for z, m in specs:
                name = self._get_species_name(z, m)
                if name in ["Deuterium", "Tritium"]:
                    d_or_t = True
                if z > 0:
                    has_any_ion = True
            if d_or_t:
                has_main_ion = True
        
        if common_species:
            # Sort by Z (usually ions first, then electrons if Z<0)
            sorted_species = sorted(list(common_species), key=lambda x: (x[0], x[1]), reverse=True)
            
            values = []
            
            # Add Main Ion option if applicable
            if has_main_ion:
                values.append("Main Ion (D+T)")
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
        """Return indices of Deuterium and Tritium species."""
        indices = []
        specs = self._get_case_species(data)
        for i, (z, m) in enumerate(specs):
            name = self._get_species_name(z, m)
            if name in ["Deuterium", "Tritium"]:
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
        elif plot_type_selection == "Fluctuation 2D":
            view = self.fluc2d_view_var.get().strip().lower()
            if view == "vs xt":
                plot_type = "Fluctuation 2D vs xt"
                display_plot_type = "Fluctuation 2D (vs xt)"
            else:
                plot_type = "Fluctuation 2D"
                display_plot_type = "Fluctuation 2D (vs xy)"
        elif plot_type_selection == "Zonal ExB Shearing Rate":
            zf_xaxis = self.zf_xaxis_var.get().strip().lower()
            if zf_xaxis == "vs kx":
                plot_type = "ZF ExB Shearing Spectrum"
                display_plot_type = "Zonal ExB Shearing Rate (vs kx)"
            elif zf_xaxis == "vs gamma_lin":
                plot_type = "ZF ExB Fig4 (kx=ky)"
                display_plot_type = "Zonal ExB Shearing Rate (vs gamma_lin)"
            else:
                plot_type = "ZF ExB Shearing Rate"
                display_plot_type = "Zonal ExB Shearing Rate (vs Time)"
        elif plot_type_selection == "Energy balance":
            mode = self.energy_balance_mode_var.get().strip().lower()
            if mode == "zf energy balance":
                plot_type = "Energy Balance ZF"
                display_plot_type = "Energy balance: ZF energy balance"
            elif (r"\gamma_{eff}" in mode) or ("gamma_eff" in mode):
                plot_type = "Energy Balance Gamma Eff"
                display_plot_type = r"Energy balance: vs $\gamma_{eff}^Z$ and $\gamma_{eff}^{NZ}$"
            elif mode == "single plot":
                qty = str(self.energy_balance_single_quantity_var.get()).strip()
                xax = str(self.energy_balance_single_xaxis_var.get()).strip()
                xax_plot = xax.replace("v.s", "vs")
                plot_type = f"Energy Balance Single {qty} {xax_plot}"
                display_plot_type = f"Energy balance: Single plot ({qty}, {xax})"
            elif mode == "v.s 2d":
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

    def _get_selected_case_names(self):
        """Return selected case names, or all loaded cases when nothing is selected."""
        selected_indices = self.case_listbox.curselection()
        if not selected_indices:
            return [self.case_listbox.get(i) for i in range(self.case_listbox.size())]
        return [self.case_listbox.get(i) for i in selected_indices]

    @staticmethod
    def _is_contour_like_plot(plot_type):
        """True when plot type is contour/multi-panel style (single-case rendering)."""
        return (
            ("FFT" in plot_type)
            or plot_type.startswith("Fluctuation 2D")
            or plot_type == "POD Parity"
            or ("vs ky_time" in plot_type)
        )

    @classmethod
    def _is_standard_line_plot(cls, plot_type):
        """True when plot type is standard line rendering."""
        return not cls._is_contour_like_plot(plot_type)

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
            selected_spec_str = ""
            try:
                selected_spec_str = self.species_var.get().strip()
            except Exception:
                selected_spec_str = ""

            if selected_spec_str == "Main Ion (D+T)":
                indices = self._get_main_ion_indices(data)
                if not indices:
                    print(f"Warning: No Main Ion (D/T) found in {case_label}")
                    return [], "Main Ion"
                if main_ion_policy == "first" and len(indices) > 1:
                    indices = [indices[0]]
                    print("Note: Main Ion sum is reduced to first main-ion species for this plot.")
                spec_label = "Main Ion"
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

        return indices, spec_label

    # ------------------------------------------------------------------
    # Case Discovery / Loading
    # ------------------------------------------------------------------

    @staticmethod
    def _get_case_picker_initial_dir():
        """Return preferred initial directory for case selection dialogs."""
        if os.path.exists(DEFAULT_CASE_PICKER_ROOT):
            return DEFAULT_CASE_PICKER_ROOT
        return os.getcwd()

    @staticmethod
    def _looks_like_case_dir(dir_path):
        """Heuristically decide whether a directory contains CGYRO case markers."""
        markers = ("input.cgyro", "out.cgyro.freq", "input.gacode")
        return any(os.path.exists(os.path.join(dir_path, m)) for m in markers)

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
                data = cgyrodata(dir_path)
            except Exception as e:
                data = cgyrodata_plot(dir_path)
            
            self.cases[case_name] = data
            self.case_listbox.insert(tk.END, case_name)
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
        for index in reversed(selected_indices):
            case_name = self.case_listbox.get(index)
            del self.cases[case_name]
            self.case_listbox.delete(index)
        self._update_species_list()

    def remove_all_cases(self):
        """Clear all loaded cases after user confirmation."""
        if messagebox.askyesno("Confirm", "Are you sure you want to remove all loaded cases?"):
            self.cases.clear()
            self.case_listbox.delete(0, tk.END)
            self._update_species_list()

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
                        new_data = cgyrodata(dir_path)
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


