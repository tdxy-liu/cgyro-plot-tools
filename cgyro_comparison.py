
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import sys
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
from matplotlib import ticker
import matplotlib.animation as animation
import traceback

# Ensure pygacode can be imported
sys.path.append(os.path.join(os.path.dirname(__file__), '../gacode/f2py'))
try:
    from pygacode.cgyro.data import cgyrodata
    from pygacode.cgyro.data_plot import cgyrodata_plot
except ImportError:
    print("Error: Could not import cgyrodata. Please ensure pygacode is available.")
    # Mock for development if imports fail completely (remove in production)
    class cgyrodata:
        def __init__(self, path):
            self.path = path
            self.ky = np.linspace(0,1,5)
            self.freq = np.zeros((2,5,1))
            self.t = np.linspace(0,10,10)
            self.ky_flux = np.zeros((1, 2, 2, 5, 10)) # [species, moment, field, ky, t]
            self.z = np.array([1.0])
            self.mass = np.array([2.0])
            self.n_radial = 1
            self.theta_plot = 1
        def getflux(self):
            pass
        def getbigfield(self):
            pass
    class cgyrodata_plot:
        pass

class CGYRO_Comparison:
    def __init__(self, root):
        self.root = root
        self.root.title("CGYRO Comparison Tool")
        self.root.geometry("1200x800")

        self.cases = {}  # Dictionary to store loaded cases: {name: cgyrodata_object}
        self.ani = None # Animation object
        self.is_paused = False
        self.current_frame = 0
        self.total_frames = 0
        self._drag_start_index = None
        
        self._create_layout()
        self._create_menu()

    def _create_layout(self):
        # Main container
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Left panel: Controls and Case List
        left_panel = ttk.Frame(main_frame, width=300, padding=10)
        left_panel.pack(side=tk.LEFT, fill=tk.Y)
        
        # Case List
        ttk.Label(left_panel, text="Loaded Cases:").pack(anchor=tk.W)
        self.case_listbox = tk.Listbox(left_panel, selectmode=tk.EXTENDED, height=10)
        self.case_listbox.pack(fill=tk.X, pady=5)
        
        # Enable drag and drop reordering
        self.case_listbox.bind('<Button-1>', self._on_drag_start)
        self.case_listbox.bind('<B1-Motion>', self._on_drag_motion)
        
        btn_frame = ttk.Frame(left_panel)
        btn_frame.pack(fill=tk.X, pady=5)
        ttk.Button(btn_frame, text="Add Case", command=self.add_case).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="Add Group", command=self.add_group).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="Remove", command=self.remove_case).pack(side=tk.LEFT, padx=2)
        ttk.Button(btn_frame, text="Remove All", command=self.remove_all_cases).pack(side=tk.LEFT, padx=2)

        # Plot Controls
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        ttk.Label(left_panel, text="Plot Type:").pack(anchor=tk.W)
        
        self.plot_type_var = tk.StringVar(value="Frequency")
        plot_types = ["Frequency", "Growth Rate", "Flux", "Phi", "Fluctuation 2D"]
        self.plot_type_combo = ttk.Combobox(left_panel, textvariable=self.plot_type_var, values=plot_types, state="readonly")
        self.plot_type_combo.pack(fill=tk.X, pady=5)
        self.plot_type_combo.bind("<<ComboboxSelected>>", self.update_options)

        # Dynamic Options Frame
        self.options_frame = ttk.LabelFrame(left_panel, text="Options", padding=5)
        self.options_frame.pack(fill=tk.X, pady=10)
        self._init_options()

        # Action Buttons
        ttk.Button(left_panel, text="Plot", command=self.plot_comparison).pack(fill=tk.X, pady=10)
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

        self.fig = plt.Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right_panel)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(self.canvas, right_panel)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _create_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Add Case Directory", command=self.add_case)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)

    def _init_options(self):
        # --- Persistent Options (Time Range) ---
        ttk.Label(self.options_frame, text="Time Start:").grid(row=0, column=0, sticky=tk.W)
        self.t_start_var = tk.StringVar()
        ttk.Entry(self.options_frame, textvariable=self.t_start_var, width=10).grid(row=0, column=1, sticky=tk.W)

        ttk.Label(self.options_frame, text="Time End:").grid(row=1, column=0, sticky=tk.W)
        self.t_end_var = tk.StringVar()
        ttk.Entry(self.options_frame, textvariable=self.t_end_var, width=10).grid(row=1, column=1, sticky=tk.W)

        # Separator
        ttk.Separator(self.options_frame, orient=tk.HORIZONTAL).grid(row=2, column=0, columnspan=2, sticky="ew", pady=5)

        # --- Dynamic Options ---
        self.dynamic_widgets = []

        # 1. Divided by ky (Frequency, Growth rate)
        self.norm_ky_var = tk.BooleanVar(value=False)
        self.norm_ky_check = ttk.Checkbutton(self.options_frame, text="Divided by ky", variable=self.norm_ky_var)
        
        # 2. Flux Options
        self.flux_type_var = tk.StringVar(value="Energy")
        self.flux_type_combo = ttk.Combobox(self.options_frame, textvariable=self.flux_type_var, values=["Energy", "Particle"], state="readonly", width=15)
        self.flux_xaxis_var = tk.StringVar(value="v.s ky")
        self.flux_xaxis_combo = ttk.Combobox(self.options_frame, textvariable=self.flux_xaxis_var, values=["v.s ky", "v.s Time"], state="readonly", width=15)
        
        # 3. Phi Options
        self.phi_xaxis_var = tk.StringVar(value="v.s ky")
        self.phi_xaxis_combo = ttk.Combobox(self.options_frame, textvariable=self.phi_xaxis_var, values=["v.s ky", "v.s Time", "fft"], state="readonly", width=15)

        # 4. Species Selection (Flux, Fluctuation 2D)
        self.species_label = ttk.Label(self.options_frame, text="Species:")
        self.species_var = tk.StringVar()
        self.species_combo = ttk.Combobox(self.options_frame, textvariable=self.species_var, width=20, state="readonly")
        
        self.plot_all_species_var = tk.BooleanVar(value=False)
        self.plot_all_species_check = ttk.Checkbutton(self.options_frame, text="Plot All Species (First Case Only)", variable=self.plot_all_species_var)

        # 5. Fluctuation 2D Moment
        self.moment_label = ttk.Label(self.options_frame, text="Moment:")
        self.moment_var = tk.StringVar(value="Phi")
        self.moment_combo = ttk.Combobox(self.options_frame, textvariable=self.moment_var,
                                         values=["Phi", "Density", "Energy", "Temperature", "Apar", "Bpar"], state="readonly")

        # 6. Phi FFT Options (reused from before, but managed dynamically)
        self.fft_options_frame = ttk.LabelFrame(self.options_frame, text="Phi FFT Settings", padding=5)
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
        # Overlay Frequency Checkbox
        self.fft_overlay_freq_var = tk.BooleanVar(value=False)
        self.fft_overlay_freq_check = ttk.Checkbutton(self.fft_options_frame, text="Overlay Linear Freq", variable=self.fft_overlay_freq_var)
        self.fft_overlay_freq_check.grid(row=2, column=0, columnspan=2, sticky=tk.W)

        # Initialize layout
        self.update_options()

    def _hide_dynamic_options(self):
        # Helper to grid_remove all dynamic widgets
        widgets = [
            self.norm_ky_check,
            self.flux_type_combo, self.flux_xaxis_combo,
            self.phi_xaxis_combo,
            self.species_label, self.species_combo, self.plot_all_species_check,
            self.moment_label, self.moment_combo,
            self.fft_options_frame
        ]
        for w in widgets:
            w.grid_remove()

    def update_options(self, event=None):
        self._hide_dynamic_options()
        plot_type = self.plot_type_var.get()
        row = 3 # Start gridding from row 3

        if plot_type in ["Frequency", "Growth Rate"]:
            self.norm_ky_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)
            
        elif plot_type == "Flux":
            self.flux_type_combo.grid(row=row, column=0, sticky=tk.W)
            self.flux_xaxis_combo.grid(row=row, column=1, sticky=tk.W)
            row += 1
            self.species_label.grid(row=row, column=0, sticky=tk.W)
            self.species_combo.grid(row=row, column=1, sticky=tk.W)
            row += 1
            self.plot_all_species_check.grid(row=row, column=0, columnspan=2, sticky=tk.W)

        elif plot_type == "Phi":
            self.phi_xaxis_combo.grid(row=row, column=0, columnspan=2, sticky=tk.W)
            row += 1
            
            # Check if FFT is selected in the sub-option
            if self.phi_xaxis_var.get() == "fft":
                 self.fft_options_frame.grid(row=row, column=0, columnspan=2, sticky=tk.W+tk.E, pady=5)
            
            # Bind change event for phi_xaxis_combo to update FFT frame visibility immediately
            self.phi_xaxis_combo.bind("<<ComboboxSelected>>", lambda e: self.update_options())

        elif plot_type == "Fluctuation 2D":
            self.moment_label.grid(row=row, column=0, sticky=tk.W)
            self.moment_combo.grid(row=row, column=1, sticky=tk.W)
            row += 1
            # Species selection needed for Density/Energy moments
            if self.moment_var.get() in ["Density", "Energy", "Temperature"]:
                 self.species_label.grid(row=row, column=0, sticky=tk.W)
                 self.species_combo.grid(row=row, column=1, sticky=tk.W)
            
            self.moment_combo.bind("<<ComboboxSelected>>", lambda e: self.update_options())

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
            if d_or_t:
                has_main_ion = True
        
        if common_species:
            # Sort by Z (usually ions first, then electrons if Z<0)
            sorted_species = sorted(list(common_species), key=lambda x: (x[0], x[1]), reverse=True)
            
            values = []
            
            # Add Main Ion option if applicable
            if has_main_ion:
                values.append("Main Ion (D+T)")

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

    def _get_main_ion_indices(self, data):
        """Return indices of Deuterium and Tritium species."""
        indices = []
        specs = self._get_case_species(data)
        for i, (z, m) in enumerate(specs):
            name = self._get_species_name(z, m)
            if name in ["Deuterium", "Tritium"]:
                indices.append(i)
        return indices

    def _load_case_from_dir(self, dir_path, silent=False):
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

    def add_case(self):
        initial_dir = "/data/share/jinyue_liu" if os.path.exists("/data/share/jinyue_liu") else os.getcwd()
        dir_path = filedialog.askdirectory(title="Select CGYRO Case Directory", initialdir=initial_dir)
        
        if dir_path:
            if self._load_case_from_dir(dir_path):
                self._update_species_list()
            else:
                messagebox.showerror("Error", f"Failed to load case from {dir_path}")

    def add_group(self):
        initial_dir = "/data/share/jinyue_liu" if os.path.exists("/data/share/jinyue_liu") else os.getcwd()
        parent_dir = filedialog.askdirectory(title="Select Parent Directory (Load All Subfolders)", initialdir=initial_dir)
        
        if parent_dir:
            loaded_count = 0
            for item in os.listdir(parent_dir):
                sub_path = os.path.join(parent_dir, item)
                if os.path.isdir(sub_path):
                    # Check if it looks like a case (optional, but good for filtering)
                    # For now, just try to load it. cgyrodata might fail if files missing, which is caught.
                    # Or we can check for input.cgyro / out.cgyro.freq
                    if os.path.exists(os.path.join(sub_path, 'input.cgyro')) or \
                       os.path.exists(os.path.join(sub_path, 'out.cgyro.freq')) or \
                       os.path.exists(os.path.join(sub_path, 'input.gacode')):
                        if self._load_case_from_dir(sub_path, silent=True):
                            loaded_count += 1
            
            if loaded_count > 0:
                self._update_species_list()
                messagebox.showinfo("Success", f"Loaded {loaded_count} cases from {parent_dir}")
            else:
                messagebox.showwarning("Warning", "No valid cases found in selected directory.")

    def remove_case(self):
        selected_indices = self.case_listbox.curselection()
        for index in reversed(selected_indices):
            case_name = self.case_listbox.get(index)
            del self.cases[case_name]
            self.case_listbox.delete(index)
        self._update_species_list()

    def remove_all_cases(self):
        if messagebox.askyesno("Confirm", "Are you sure you want to remove all loaded cases?"):
            self.cases.clear()
            self.case_listbox.delete(0, tk.END)
            self._update_species_list()

    def _on_drag_start(self, event):
        if self.case_listbox.size() > 0:
            self._drag_start_index = self.case_listbox.nearest(event.y)

    def _on_drag_motion(self, event):
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
        if self.ani:
            self.ani.event_source.stop()
            self.ani = None
        self.is_paused = False
        self.anim_update_func = None
        if hasattr(self, 'btn_pause'):
            self.btn_pause.config(text="Pause", state="disabled")
            self.btn_prev.config(state="disabled")
            self.btn_next.config(state="disabled")

    def toggle_pause(self):
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
        if self.ani and self.is_paused and hasattr(self, 'anim_update_func'):
            self.current_frame = (self.current_frame - 1) % self.total_frames
            self.anim_update_func(self.current_frame)
            self.canvas.draw()

    def next_frame(self):
        if self.ani and self.is_paused and hasattr(self, 'anim_update_func'):
            self.current_frame = (self.current_frame + 1) % self.total_frames
            self.anim_update_func(self.current_frame)
            self.canvas.draw()

    def clear_plot(self):
        self._stop_animation()
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.canvas.draw()

    def _get_time_indices(self, t):
        try:
            t_start_str = self.t_start_var.get()
            t_end_str = self.t_end_var.get()
            
            if len(t) == 0:
                return [], 0, 0

            t_min = t[0]
            t_max = t[-1]
            
            if t_start_str: t_start = float(t_start_str)
            else: t_start = t_max * 0.5
            
            if t_end_str: t_end = float(t_end_str)
            else: t_end = t_max

            # If specified range is completely outside data, fallback to last 50%
            # Also fallback if user specified end time extends beyond simulation time
            if t_start >= t_max or t_end <= t_min or t_end > t_max:
                t_start = t_max * 0.5
                t_end = t_max

            indices = np.where((t >= t_start) & (t <= t_end))[0]
            return indices, t_start, t_end
        except Exception as e:
            print(f"Error parsing time range: {e}")
            return [], 0, 0

    @staticmethod
    def _maptoreal_fft(nr, nn, nx, ny, c):
        # Storage for numpy inverse real transform (irfft2)
        d = np.zeros([nx, nn], dtype=complex)

        for i in range(nr):
            p = i - nr // 2
            # k is the "standard FFT index"
            if -p < 0:
                k = -p + nx
            else:
                k = -p
            # Use identity f(p,-n) = f(-p,n)*
            d[k, 0:nn] = np.conj(c[i, 0:nn])

        # 2D inverse real Hermitian transform
        # NOTE: using inverse FFT with convention exp(ipx+iny), so need n -> -n
        # NOTE: need factor of 0.5 to match half-sum method of slow maptoreal()
        f = np.fft.irfft2(d, s=[nx, ny], norm='forward') * 0.5
        return f

    def _plot_1d(self, x, y, label, plot_type):
        if len(x) > 0 and len(y) > 0:
            if len(x) == len(y):
                if "vs Time" in plot_type:
                    self.ax.plot(x, y, label=label) # No marker for time traces
                else:
                    self.ax.plot(x, y, marker='o', label=label)
            else:
                print(f"Dimension mismatch for {label}: x={len(x)}, y={len(y)}")

    def plot_comparison(self):
        # Stop any existing animation
        self._stop_animation()
        
        # Clear figure to remove colorbars if they exist
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        
        plot_type_selection = self.plot_type_var.get()
        
        # Synthesize the effective plot_type string expected by the plotting methods
        plot_type = plot_type_selection
        
        if plot_type_selection == "Flux":
            flux_type = self.flux_type_var.get() # Energy or Particle
            flux_xaxis = self.flux_xaxis_var.get() # v.s ky or v.s time
            
            # Construct string like "Energy Flux vs ky" or "Particle Flux vs Time"
            # Note: The dropdown values are "v.s ky" but previous code expected "vs ky"
            xaxis_str = flux_xaxis.replace("v.s", "vs")
            plot_type = f"{flux_type} Flux {xaxis_str}"
            
        elif plot_type_selection == "Phi":
            phi_xaxis = self.phi_xaxis_var.get()
            if phi_xaxis == "fft":
                plot_type = "Phi FFT"
            else:
                xaxis_str = phi_xaxis.replace("v.s", "vs")
                plot_type = f"Phi {xaxis_str}"

        selected_indices = self.case_listbox.curselection()
        
        if not selected_indices:
            # If nothing selected, plot all (in listbox order)
            selected_cases = [self.case_listbox.get(i) for i in range(self.case_listbox.size())]
        else:
            selected_cases = [self.case_listbox.get(i) for i in selected_indices]

        # Special handling for contour plots (Phi FFT, Fluctuation 2D)
        if plot_type in ["Phi FFT", "Fluctuation 2D"]:
            if len(selected_cases) > 1:
                messagebox.showwarning("Warning", f"{plot_type} only supports one case at a time. Plotting the first selected case.")
                selected_cases = list(selected_cases)[:1]

        # Handle "Plot All Species" logic
        if self.plot_all_species_var.get() and "Flux" in plot_type and len(selected_cases) > 0:
            # Only plot the first selected case, but iterate over all its species
            # selected_cases might be a dict_keys object if from self.cases.keys()
            # Convert to list to be safe
            case_list = list(selected_cases)
            case_name = case_list[0]
            data = self.cases[case_name]
            
            # Ensure flux data is loaded to get species list
            if not hasattr(data, 'ky_flux'):
                try:
                    data.getflux()
                except:
                    pass
            
            specs = self._get_case_species(data)
            if not specs:
                print(f"No species found for {case_name}")
                return

            for i, (z, m) in enumerate(specs):
                # Create a temporary override for species selection
                try:
                    self._plot_single_case(data, case_name, plot_type, species_override_index=i)
                except Exception as e:
                    print(f"Error plotting species {i} for {case_name}: {e}")
                    traceback.print_exc()

        else:
            # Standard plotting loop
            for case_name in selected_cases:
                data = self.cases[case_name]
                try:
                    self._plot_single_case(data, case_name, plot_type)
                except Exception as e:
                    print(f"Error plotting {case_name}: {e}")
                    traceback.print_exc()

        if plot_type != "Phi FFT":
            handles, labels = self.ax.get_legend_handles_labels()
            if handles:
                self.ax.legend()
            self.ax.grid(True)
            
        title = f"Comparison: {plot_type}"
        if self.norm_ky_var.get() and plot_type_selection in ["Frequency", "Growth Rate"]:
             title = f"{title}/ky"
        self.ax.set_title(title)
        self.canvas.draw()

    def _plot_frequency_growth(self, data, label, plot_type, t_indices, t_start, t_end):
        x = []
        y = []
        if plot_type == "Frequency":
            if hasattr(data, 'freq'):
                x = data.ky
                n_ky = len(x)
                nt = len(data.t) if hasattr(data, 't') else 0
                
                # Check if freq has time dimension and handle possible shape variations
                if data.freq.ndim == 3:
                    shape = data.freq.shape
                    # Expected: (2, n_ky, nt) or (2, nt, n_ky)
                    
                    # Case 1: (2, n_ky, nt) - Standard
                    if shape[1] == n_ky and shape[2] == nt:
                        if len(t_indices) > 0 and t_indices.max() < shape[2]:
                            y = np.mean(data.freq[0, :, t_indices], axis=0)
                            label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
                        else:
                            y = data.freq[0, :, -1]
                    
                    # Case 2: (2, nt, n_ky) - Swapped
                    elif shape[1] == nt and shape[2] == n_ky:
                            if len(t_indices) > 0 and t_indices.max() < shape[1]:
                                y = np.mean(data.freq[0, t_indices, :], axis=1)
                                label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
                            else:
                                y = data.freq[0, -1, :]
                    
                    # Fallback: Try to match dimensions
                    else:
                        if shape[1] == n_ky:
                            # Assume (2, n_ky, nt)
                            if len(t_indices) > 0 and t_indices.max() < shape[2]:
                                y = np.mean(data.freq[0, :, t_indices], axis=0)
                                label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
                            else:
                                y = data.freq[0, :, -1]
                        elif shape[2] == n_ky:
                            # Assume (2, nt, n_ky)
                            if len(t_indices) > 0 and t_indices.max() < shape[1]:
                                y = np.mean(data.freq[0, t_indices, :], axis=1)
                                label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
                            else:
                                y = data.freq[0, -1, :]
                        else:
                            # Last resort: assume standard
                            y = data.freq[0, :, -1]
                else:
                    # 2D array [2, ky] (no time dimension)
                    y = data.freq[0, :]

        elif plot_type == "Growth Rate":
            if hasattr(data, 'freq'):
                x = data.ky
                n_ky = len(x)
                nt = len(data.t) if hasattr(data, 't') else 0

                if data.freq.ndim == 3:
                    shape = data.freq.shape
                    
                    # Case 1: (2, n_ky, nt) - Standard
                    if shape[1] == n_ky and shape[2] == nt:
                        if len(t_indices) > 0 and t_indices.max() < shape[2]:
                            y = np.mean(data.freq[1, :, t_indices], axis=0)
                            label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
                        else:
                            y = data.freq[1, :, -1]
                    
                    # Case 2: (2, nt, n_ky) - Swapped
                    elif shape[1] == nt and shape[2] == n_ky:
                            if len(t_indices) > 0 and t_indices.max() < shape[1]:
                                y = np.mean(data.freq[1, t_indices, :], axis=1)
                                label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
                            else:
                                y = data.freq[1, -1, :]
                            
                    # Fallback
                    else:
                        if shape[1] == n_ky:
                            # Assume (2, n_ky, nt)
                            if len(t_indices) > 0 and t_indices.max() < shape[2]:
                                y = np.mean(data.freq[1, :, t_indices], axis=0)
                                label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
                            else:
                                y = data.freq[1, :, -1]
                        elif shape[2] == n_ky:
                            # Assume (2, nt, n_ky)
                            if len(t_indices) > 0 and t_indices.max() < shape[1]:
                                y = np.mean(data.freq[1, t_indices, :], axis=1)
                                label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
                            else:
                                y = data.freq[1, -1, :]
                        else:
                            y = data.freq[1, :, -1]
                else:
                    y = data.freq[1, :]
        
        # Normalize by ky if requested
        if self.norm_ky_var.get() and len(x) > 0 and len(y) > 0:
            # Avoid division by zero
            with np.errstate(divide='ignore', invalid='ignore'):
                y = np.divide(y, x)
        
        self._plot_1d(x, y, label, plot_type)

    def _plot_fluctuation_2d(self, data, label, plot_type, t_indices, t_start, t_end):
        if not hasattr(data, 'kxky_phi'):
            try:
                print(f"Loading big field data for {label}...")
                data.getbigfield()
            except Exception as e:
                print(f"Could not load big field data for {label}: {e}")
                return

        moment = self.moment_var.get()
        target_indices = []
        
        if moment in ["Density", "Energy"]:
            selected_spec_str = self.species_var.get()
            if selected_spec_str == "Main Ion (D+T)":
                 # Fluctuation 2D does not support Main Ion sum, default to first main ion or first species
                 target_indices = self._get_main_ion_indices(data)
                 if target_indices:
                      target_indices = [target_indices[0]] # Just use the first one (e.g. D)
                      print(f"Note: Fluctuation 2D does not support 'Main Ion' sum. Plotting first main ion species.")
                 else:
                      print(f"Warning: No Main Ion (D/T) found in {label}")
                      return
            else:
                match = re.search(r"Z=([-\d\.]+), M=([-\d\.]+)", selected_spec_str)
                if match:
                    target_z = float(match.group(1))
                    target_m = float(match.group(2))
                    current_specs = self._get_case_species(data)
                    for idx, (z, m) in enumerate(current_specs):
                            if abs(z - target_z) < 0.01 and abs(m - target_m) < 0.01:
                                target_indices = [idx]
                                break
                else:
                     target_indices = [0]
        else:
             # For field moments (Phi, Apar, Bpar), species doesn't matter, just use 0 as placeholder if needed
             target_indices = [0]

        if not target_indices:
             target_indices = [0]
        
        # Select data based on moment
        c_field = None
        if moment == "Phi":
            if hasattr(data, 'kxky_phi'):
                    c_field = data.kxky_phi
        elif moment == "Apar":
            if hasattr(data, 'kxky_apar'):
                    c_field = data.kxky_apar
        elif moment == "Bpar":
            if hasattr(data, 'kxky_bpar'):
                    c_field = data.kxky_bpar
        elif moment == "Density":
            if hasattr(data, 'kxky_n'):
                    c_field = data.kxky_n
        elif moment == "Energy":
            if hasattr(data, 'kxky_e'):
                    c_field = data.kxky_e
        elif moment == "Temperature":
            if hasattr(data, 'kxky_n') and hasattr(data, 'kxky_e'):
                c_field = (data.kxky_n, data.kxky_e)
            else:
                print(f"Data for Temperature (n and e) not available in {label}")
                return

        if c_field is None and moment != "Temperature":
            print(f"Data for {moment} not available in {label}")
            return

        # Theta index: midplane
        itheta = data.theta_plot * 4 // 8
        if itheta >= data.theta_plot: itheta = 0

        # Define calculation function
        def get_frame_data(ti):
            if moment == "Temperature":
                # Only support single species for Temperature for now
                species_idx = target_indices[0]
                
                # c_field is tuple (kxky_n, kxky_e)
                n_data = c_field[0]
                e_data = c_field[1]
                
                # Reconstruct complex arrays
                # data shape is [2, nr, theta, species, nn, nt]
                # we want [2, nr, nn]
                ndim = data.kxky_n.ndim
                shape = data.kxky_n.shape
                # print(f"DEBUG: kxky_phi shape: {shape}")

                if ndim == 6:
                    # Standard: [2, nr, theta, species, nn, nt]
                    n_real = n_data[0, :, itheta, species_idx, :, ti]
                    n_imag = n_data[1, :, itheta, species_idx, :, ti]
                    e_real = e_data[0, :, itheta, species_idx, :, ti]
                    e_imag = e_data[1, :, itheta, species_idx, :, ti]
                
                    cn = n_real + 1j * n_imag
                    ce = e_real + 1j * e_imag
        
                elif ndim == 5:
                    # Standard: [nr, theta, species, nn, nt]
                    cn = n_data[:, itheta, species_idx, :, ti]
                    ce = e_data[:, itheta, species_idx, :, ti]
                else:
                    print(f"Error: Unsupported kxky_phi dimensions: {ndim} {shape}")
                    return

                t0 = data.temp[species_idx]
                n0 = data.dens[species_idx]
                c = ((2.0/3.0) * ce - t0 * cn) / n0
                
            elif moment in ["Density", "Energy"]:
                    # [2, nr, theta, species, nn, nt]
                    ndim = data.kxky_n.ndim
                    shape = data.kxky_n.shape
                    
                    c = 0
                    
                    for s_idx in target_indices:
                        if ndim == 6:
                        # Standard: [2, nr, theta, nn, nt]
                            c_real = c_field[0, :, itheta, s_idx, :, ti]
                            c_imag = c_field[1, :, itheta, s_idx, :, ti]
                            c += c_real + 1j * c_imag
            
                        elif ndim == 5:
                        # Standard: [nr, theta, species, nn, nt]
                            c += c_field[:, itheta, s_idx, :, ti]
                        else:
                            print(f"Error: Unsupported kxky_phi dimensions: {ndim} {shape}")
                            return

            else:
                    # Fields are species independent
                    ndim = data.kxky_phi.ndim
                    shape = data.kxky_phi.shape

                    if ndim == 5:
                        c_real = c_field[0, :, itheta, :, ti]
                        c_imag = c_field[1, :, itheta, :, ti]
                        c = c_real + 1j * c_imag
                    elif ndim == 4:
                    # Standard: [nr, theta, nn, nt]
                        c = c_field[:, itheta, :, ti]
                    else:
                        print(f"Error: Unsupported kxky_phi dimensions: {ndim} {shape}")
                        return
            
            # FFT params
            nr = data.n_radial
            nn = data.n_n
            nx = nr + 1
            ny = 2 * nn - 1
            
            f = self._maptoreal_fft(nr, nn, nx, ny, c)
            f = 2 * f # Correct for half-sum

            # Grid
            x = np.arange(nx) * 2 * np.pi / nx
            y = np.arange(ny) * 2 * np.pi / ny
            
            xmax = data.length
            if len(data.ky) > 1 and abs(data.ky[1]) > 1e-9:
                    ymax = (2 * np.pi) / np.abs(data.ky[1])
            else:
                    ymax = 100.0 # Fallback
            
            xp = x / (2 * np.pi) * xmax
            yp = y / (2 * np.pi) * ymax
            
            # Periodic extensions
            xp = np.append(xp, xmax)
            yp = np.append(yp, ymax)
            fp = np.zeros([nx + 1, ny + 1])
            fp[0:nx, 0:ny] = f[:, :]
            fp[-1, :] = fp[0, :]
            fp[:, -1] = fp[:, 0]
            
            return xp, yp, fp, xmax, ymax

        # Determine if animation is needed
        # Only animate if user explicitly provided a time range (start or end)
        user_specified_time = self.t_start_var.get().strip() or self.t_end_var.get().strip()
        should_animate = user_specified_time and len(t_indices) > 1

        if should_animate:
            # Animation Logic
            self.total_frames = len(t_indices)
            self.current_frame = 0

            def update(frame):
                # Handle loop or bounds if called manually
                frame = frame % self.total_frames
                self.current_frame = frame
                
                ti = t_indices[frame]
                xp, yp, fp, xmax, ymax = get_frame_data(ti)
                
                self.ax.clear()
                levels = 50
                self.ax.contourf(xp, yp, np.transpose(fp), levels, cmap='RdBu_r')
                # Note: Colorbar in animation might be tricky if range changes.
                # Ideally set vmin/vmax fixed based on global min/max, but here we let it scale.
                
                self.ax.set_xlabel(r'$x/\rho_s$')
                self.ax.set_ylabel(r'$y/\rho_s$')
                self.ax.set_title(f'{moment} Fluctuation: {label} (t={data.t[ti]:.2f})')
                self.ax.set_aspect('equal')
                self.ax.set_xlim([0, xmax])
                self.ax.set_ylim([0, ymax])

            self.anim_update_func = update
            self.ani = animation.FuncAnimation(self.fig, update, frames=self.total_frames, interval=100, repeat=True)
            
            # Enable controls
            self.btn_pause.config(state="normal", text="Pause")
            self.btn_prev.config(state="normal")
            self.btn_next.config(state="normal")
            
            self.canvas.draw()
            
        else:
            # Static Plot (Last time point or end of range)
            ti = -1
            if len(t_indices) > 0:
                ti = t_indices[-1]

            xp, yp, fp, xmax, ymax = get_frame_data(ti)
            
            levels = 50
            cs = self.ax.contourf(xp, yp, np.transpose(fp), levels, cmap='RdBu_r')
            self.fig.colorbar(cs, ax=self.ax)
            
            self.ax.set_xlabel(r'$x/\rho_s$')
            self.ax.set_ylabel(r'$y/\rho_s$')
            self.ax.set_title(f'{moment} Fluctuation: {label} (t={data.t[ti]:.2f})')
            self.ax.set_aspect('equal')
            self.ax.set_xlim([0, xmax])
            self.ax.set_ylim([0, ymax])

    def _plot_phi_fft(self, data, label, t_indices):
        # Implementation of phi_f_fft logic
        # Need to compute FFT of phi over time
        
        # Check if we have big field data
        if not hasattr(data, 'kxky_phi'):
            try:
                print(f"Loading big field data for {label}...")
                data.getbigfield()
            except Exception as e:
                print(f"Could not load big field data for {label}: {e}")
                return

        if not hasattr(data, 'kxky_phi'):
            print(f"No kxky_phi data available for {label}")
            return

        # kxky_phi: [2, n_radial, theta_plot, n_n, nt]
        # We need phi(ky, t) -> FFT -> phi(ky, omega)
        
        ndim = data.kxky_phi.ndim
        shape = data.kxky_phi.shape
        # print(f"DEBUG: kxky_phi shape: {shape}")

        phi_t_all = None

        if ndim == 5:
            # Standard: [2, n_radial, theta_plot, n_n, nt]
            i_theta = data.theta_plot * 4 // 8 # Outboard midplane approx
            phi_t_all = data.kxky_phi[0, :, i_theta, :, :] + 1j * data.kxky_phi[1, :, i_theta, :, :]
        
        elif ndim == 4:
            # Standard: [n_radial, theta_plot, n_n, nt]
            i_theta = data.theta_plot * 4 // 8 # Outboard midplane approx
            phi_t_all = data.kxky_phi[:, i_theta, :, :]

        else:
            print(f"Error: Unsupported kxky_phi dimensions: {ndim} {shape}")
            return
        # phi_t_all shape: [n_radial, n_n, nt]
        
        # Slice time if needed
        if len(t_indices) > 0:
            phi_t_all = phi_t_all[..., t_indices]
            t_array = data.t[t_indices]
        else:
            t_array = data.t
        
        nt_slice = len(t_array)
        if nt_slice < 2:
            print("Not enough time points for FFT")
            return

        # Perform FFT along time axis
        # omega = 2*pi*f
        dt = t_array[1] - t_array[0]
        # Parse out.cgyro.info for ion direction to ensure correctness as per user request
        # Default to -1.0 (Standard CGYRO e^-iwt vs FFT e^-i2pift)
        freq_mult = -1.0
        
        try:
            info_path = os.path.join(data.dir, 'out.cgyro.info')
            if os.path.exists(info_path):
                with open(info_path, 'r') as f:
                    info_content = f.read()
                    if "Ion direction: omega > 0" in info_content:
                        # Ions are positive freq. FFT of e^-i|w|t is at -|w|.
                        # Mult by -1 gives +|w|. Correct.
                        freq_mult = -1.0
                    elif "Ion direction: omega < 0" in info_content:
                        # Ions are negative freq. FFT of e^i|w|t is at +|w|.
                        # Mult by -1 gives -|w|. Correct.
                        freq_mult = -1.0
        except Exception as e:
            print(f"Warning checking info file: {e}")

        # Standard FFT freq first
        omega = np.fft.fftfreq(nt_slice, d=dt) * 2 * np.pi
        phi_omega_all = np.fft.fft(phi_t_all, axis=-1)
        
        # Shift zero frequency to center
        omega_shifted = np.fft.fftshift(omega)
        phi_omega_all_shifted = np.fft.fftshift(phi_omega_all, axes=-1)

        # Apply sign correction
        omega_shifted *= freq_mult
        
        # Flip to ensure ascending axis for contourf
        # Since freq_mult is -1.0, the axis is now descending. Flip it.
        if freq_mult < 0:
            omega_shifted = np.flip(omega_shifted)
            phi_omega_all_shifted = np.flip(phi_omega_all_shifted, axis=-1)
        
        # 1. kx=0 component
        i_r_0 = data.n_radial // 2
        phi_kx0 = phi_omega_all_shifted[i_r_0, :, :]
        mag_kx0 = np.abs(phi_kx0)

        # 3. ky=0 component (omega vs kx)
        # phi_omega_all_shifted: [n_radial, n_ky, n_omega]
        # We need [n_radial, n_omega] for ky=0
        # ky=0 is usually at index 0, but let's verify
        ky_idx_0 = 0
        if abs(data.ky[0]) > 1e-6:
                # Search for 0
                idx = np.where(np.abs(data.ky) < 1e-6)[0]
                if len(idx) > 0:
                    ky_idx_0 = idx[0]
                else:
                    print("Warning: ky=0 not found. Using first ky index.")
                    ky_idx_0 = 0
        
        phi_ky0 = phi_omega_all_shifted[:, ky_idx_0, :] # [n_radial, n_omega]
        mag_ky0 = np.abs(phi_ky0)
        
        # 2. Calculate quantities based on View Mode
        view_mode = self.fft_view_var.get()
        analysis_mode = self.fft_mode_var.get()
        
        # Pre-calculate kx array
        # data.kx usually contains only positive modes [0, 1, ..., n_radial-1] * dkx
        # But for FFT shifting we need the full negative to positive range.
        # Reconstruct full kx array based on data.kx[1] (dkx) and n_radial
        
        n_radial = data.n_radial
        
        # Check if data.kx is sufficient
        if len(data.kx) > 1:
            dkx = data.kx[1] - data.kx[0] # Assuming uniform spacing and starts at 0
            if dkx == 0 and len(data.kx) > 2: dkx = data.kx[2] - data.kx[1]
        else:
            # Fallback if kx is too short (e.g. n_radial=1)
            dkx = 2 * np.pi / data.length if data.length > 0 else 1.0
        
        # Standard FFT frequency ordering: 0, 1, ..., n/2-1, -n/2, ..., -1
        kx_full = np.fft.fftfreq(n_radial, d=1.0) * dkx * n_radial # fftfreq returns f/fs, we want k = 2*pi*m/L
        # Wait, fftfreq(n, d) returns f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)
        # We want kx = 2*pi * f
        # So kx = 2*pi * fftfreq(n_radial, d=1.0) ? No, physical length L matters.
        # kx = 2*pi * m / L.
        # dkx = 2*pi/L.
        # fftfreq(n, d=L/n) returns [0, 1/L, ...]. Multiply by 2*pi -> [0, 2pi/L, ...]
        
        # Let's use dkx directly.
        # kx_indices = np.fft.fftfreq(n_radial, d=1.0) * n_radial # [0, 1, ..., -n/2, ...]
        # kx_full = kx_indices * dkx
        
        # But data.kx might already be [0, dkx, 2dkx...].
        # If we assume standard FFT layout for the data being plotted (phi_shifted_radial is shifted later):
        
        # Actually, simpler approach:
        # We want the axis for the shifted plot.
        # The shifted plot has n_radial points.
        # The grid should center at 0.
        
        kx_full = np.fft.fftfreq(n_radial, d=1.0) * n_radial * dkx
        kx_shifted = np.fft.fftshift(kx_full)

        if view_mode == "Omega vs ky":
            # --- Logic for Omega vs ky ---
            # 1. kx=0 component
            # 2. Sum over all kx
            
            mag_plot1 = mag_kx0
            
            if analysis_mode == "Linear":
                # Coherent sum
                complex_sum = np.sum(phi_omega_all_shifted, axis=0) # Sum over radial (kx)
                mag_plot2 = np.abs(complex_sum)
                title2_suffix = "(Sum all kx/x=0)"
                
                # Per-ky normalization
                for i in range(mag_plot1.shape[0]):
                    if not np.all(np.isfinite(mag_plot1[i, :])): mag_plot1[i, :] = 0.0
                    else:
                        m = np.max(mag_plot1[i, :])
                        if m > 0: mag_plot1[i, :] /= m
                        
                for i in range(mag_plot2.shape[0]):
                    if not np.all(np.isfinite(mag_plot2[i, :])): mag_plot2[i, :] = 0.0
                    else:
                        m = np.max(mag_plot2[i, :])
                        if m > 0: mag_plot2[i, :] /= m
                        
            else: # Nonlinear
                # Incoherent sum
                mag_plot2 = np.sum(np.abs(phi_omega_all_shifted), axis=0)
                title2_suffix = "(Sum all kx/Energy)"
                
                # Global normalization
                if not np.all(np.isfinite(mag_plot1)): mag_plot1[~np.isfinite(mag_plot1)] = 0.0
                m = np.max(mag_plot1)
                if m > 0: mag_plot1 /= m
                
                if not np.all(np.isfinite(mag_plot2)): mag_plot2[~np.isfinite(mag_plot2)] = 0.0
                m = np.max(mag_plot2)
                if m > 0: mag_plot2 /= m

            # Plotting
            ky_grid, omega_grid = np.meshgrid(data.ky, omega_shifted, indexing='ij')
            
            self.fig.clf()
            ax1 = self.fig.add_subplot(211)
            ax2 = self.fig.add_subplot(212)
            
            cs1 = ax1.contourf(ky_grid, omega_grid, mag_plot1, levels=50, cmap='jet')
            self.fig.colorbar(cs1, ax=ax1)
            ax1.set_title(f'Phi FFT (kx=0): {label}')
            ax1.set_ylabel(r'$\omega (c_s/a)$')
            
            cs2 = ax2.contourf(ky_grid, omega_grid, mag_plot2, levels=50, cmap='jet')
            self.fig.colorbar(cs2, ax=ax2)
            ax2.set_title(f'Phi FFT {title2_suffix}: {label}')
            ax2.set_ylabel(r'$\omega (c_s/a)$')
            ax2.set_xlabel(r'$k_y \rho_s$')
            
            # Overlay Linear Frequency
            if self.fft_overlay_freq_var.get() and hasattr(data, 'freq'):
                freq_ky = data.ky
                freq_omega = None
                if data.freq.ndim == 3:
                    shape = data.freq.shape
                    if shape[1] == len(freq_ky): freq_omega = data.freq[0, :, -1]
                    elif shape[2] == len(freq_ky): freq_omega = data.freq[0, -1, :]
                elif data.freq.ndim == 2:
                        freq_omega = data.freq[0, :]
                        
                if freq_omega is not None and len(freq_ky) == len(freq_omega):
                        freq_omega_plot = freq_omega * freq_mult
                        ax1.plot(freq_ky, freq_omega_plot, 'k--', linewidth=1.5, label='Linear Freq')
                        ax2.plot(freq_ky, freq_omega_plot, 'k--', linewidth=1.5, label='Linear Freq')
                        ax1.legend(loc='upper right')
                        ax2.legend(loc='upper right')
            
            self.ax = ax2

        elif view_mode == "Omega vs kx":
            # --- Logic for Omega vs kx ---
            # 1. ky=0 component
            # 2. Sum over all ky
            
            # phi_omega_all_shifted: [n_radial, n_ky, n_omega]
            # We need [n_radial, n_omega]
            
            # phi_omega_all_shifted is [n_radial, n_ky, n_omega]
            # The radial dimension is in "standard FFT order" (0, 1, ... -N/2 ...)
            # We need to shift it to match kx_shifted (which is centered at 0)
            
            # Shift radial dimension (axis 0) to center kx=0
            phi_shifted_radial = np.fft.fftshift(phi_omega_all_shifted, axes=0)
            
            # ky=0 component
            phi_ky0 = phi_shifted_radial[:, ky_idx_0, :]
            mag_plot1 = np.abs(phi_ky0)
            
            # Sum over all ky
            if analysis_mode == "Linear":
                # Coherent sum
                complex_sum = np.sum(phi_shifted_radial, axis=1) # Sum over ky
                mag_plot2 = np.abs(complex_sum)
                title2_suffix = "(Sum all ky/y=0)"
                
                # Per-kx normalization
                for i in range(mag_plot1.shape[0]):
                    if not np.all(np.isfinite(mag_plot1[i, :])): mag_plot1[i, :] = 0.0
                    else:
                        m = np.max(mag_plot1[i, :])
                        if m > 0: mag_plot1[i, :] /= m
                        
                for i in range(mag_plot2.shape[0]):
                    if not np.all(np.isfinite(mag_plot2[i, :])): mag_plot2[i, :] = 0.0
                    else:
                        m = np.max(mag_plot2[i, :])
                        if m > 0: mag_plot2[i, :] /= m
                        
            else: # Nonlinear
                # Incoherent sum
                mag_plot2 = np.sum(np.abs(phi_shifted_radial), axis=1)
                title2_suffix = "(Sum all ky/Energy)"
                
                # Global normalization
                if not np.all(np.isfinite(mag_plot1)): mag_plot1[~np.isfinite(mag_plot1)] = 0.0
                m = np.max(mag_plot1)
                if m > 0: mag_plot1 /= m
                
                if not np.all(np.isfinite(mag_plot2)): mag_plot2[~np.isfinite(mag_plot2)] = 0.0
                m = np.max(mag_plot2)
                if m > 0: mag_plot2 /= m

            # Plotting
            # Ensure dimensions match: kx_shifted (n_radial) vs mag_plot (n_radial, n_omega)
            if len(kx_shifted) != mag_plot1.shape[0]:
                # kx is defined on n_radial-1 grid points sometimes?
                # Or data.kx is not including the 0 point correctly?
                # Check size
                if len(kx_shifted) == mag_plot1.shape[0] - 1:
                        # Append or prepend?
                        # data.kx usually has n_radial-1 points.
                        # Let's try to regenerate kx based on shape of mag_plot1
                        n_r = mag_plot1.shape[0]
                        lx = data.length
                        dkx = 2 * np.pi / lx
                        kx = np.fft.fftfreq(n_r, d=1.0/n_r) * n_r * dkx
                        kx_shifted = np.fft.fftshift(kx)
                elif len(kx_shifted) == mag_plot1.shape[0] + 1:
                        kx_shifted = kx_shifted[:-1] # Truncate?

            # Filter for positive kx only
            pos_kx_indices = np.where(kx_shifted >= 0)[0]
            kx_shifted = kx_shifted[pos_kx_indices]
            mag_plot1 = mag_plot1[pos_kx_indices, :]
            mag_plot2 = mag_plot2[pos_kx_indices, :]

            kx_grid, omega_grid_kx = np.meshgrid(kx_shifted, omega_shifted, indexing='ij')
            
            self.fig.clf()
            ax1 = self.fig.add_subplot(211)
            ax2 = self.fig.add_subplot(212)
            
            cs1 = ax1.contourf(kx_grid, omega_grid_kx, mag_plot1, levels=50, cmap='jet')
            self.fig.colorbar(cs1, ax=ax1)
            ax1.set_title(f'Phi FFT (ky=0): {label}')
            ax1.set_ylabel(r'$\omega (c_s/a)$')
            ax1.invert_xaxis()
            
            cs2 = ax2.contourf(kx_grid, omega_grid_kx, mag_plot2, levels=50, cmap='jet')
            self.fig.colorbar(cs2, ax=ax2)
            ax2.set_title(f'Phi FFT {title2_suffix}: {label}')
            ax2.set_ylabel(r'$\omega (c_s/a)$')
            ax2.set_xlabel(r'$k_x \rho_s$')
            ax2.invert_xaxis()
            
            self.ax = ax2

        self.fig.tight_layout()
    def _plot_flux(self, data, label, plot_type, t_indices, t_start, t_end, species_override_index):
    # Ensure flux data is loaded
        if not hasattr(data, 'ky_flux'):
            try:
                print(f"Loading flux for {label}...")
                data.getflux()
            except Exception as e:
                print(f"Could not load flux for {label}: {e}")
        
        moment_idx = 1 if "Energy" in plot_type else 0 # 0=particle, 1=energy
        
        # Determine target species indices
        target_indices = []
        spec_label = ""
        
        if species_override_index is not None:
            target_indices = [species_override_index]
            # Get label
            try:
                specs = self._get_case_species(data)
                z, m = specs[species_override_index]
                spec_label = self._get_species_name(z, m)
            except:
                spec_label = f"s{species_override_index+1}"
        else:
            selected_spec_str = self.species_var.get()
            
            if selected_spec_str == "Main Ion (D+T)":
                target_indices = self._get_main_ion_indices(data)
                if not target_indices:
                    print(f"Warning: No Main Ion (D/T) found in {label}")
                    return
                spec_label = "Main Ion"
            else:
                # Parse "Z=1.0, M=2.0"
                match = re.search(r"Z=([-\d\.]+), M=([-\d\.]+)", selected_spec_str)
                if match:
                    target_z = float(match.group(1))
                    target_m = float(match.group(2))
                    
                    # Find index in current data
                    current_specs = self._get_case_species(data) # Returns [(z,m), ...]
                    # Find index matching target_z, target_m
                    for idx, (z, m) in enumerate(current_specs):
                            if abs(z - target_z) < 0.01 and abs(m - target_m) < 0.01:
                                target_indices = [idx]
                                spec_label = self._get_species_name(z, m)
                                break
                    if not target_indices:
                            print(f"Warning: Species {selected_spec_str} not found in {label}")
                            return
                else:
                    # Fallback
                    target_indices = [0]
                    spec_label = "s1"

        label = f"{label} - {spec_label}"

        x = []
        y = []

        if "vs ky" in plot_type:
            if hasattr(data, 'ky_flux'):
                # Native: [n_species, n_flux, n_field, n_n, nt]
                try:
                    if len(t_indices) > 0 and data.ky_flux.ndim == 5:
                        # Sum flux over target species
                        flux_sum = 0
                        for s_idx in target_indices:
                             # data.ky_flux[s_idx, moment_idx] -> [field, ky, t]
                             # Sum over field (axis 0) -> [ky, t]
                             flux_sum += np.sum(data.ky_flux[s_idx, moment_idx], axis=0)
                        
                        # Average over selected time indices
                        y = np.mean(flux_sum[:, t_indices], axis=1)
                        label = f"{label} (Time Avg: {t_start:.1f}-{t_end:.1f})"
                    else:
                        # Fallback to last time step if no range valid
                        y = 0
                        for s_idx in target_indices:
                             y += np.sum(data.ky_flux[s_idx, moment_idx, :, :, -1], axis=0)

                except Exception as e:
                    print(f"Error calculating time average for {label}: {e}, falling back to last step")
                    y = 0
                    for s_idx in target_indices:
                        y += np.sum(data.ky_flux[s_idx, moment_idx, :, :, -1], axis=0)
                    
                x = data.ky

        elif "vs Time" in plot_type:
            if hasattr(data, 'ky_flux'):
                # Native: [n_species, n_flux, n_field, n_n, nt]
                flux_sum_t = 0
                for s_idx in target_indices:
                    # flux_t = data.ky_flux[s_idx, moment_idx, :, :, :]
                    # Sum over field (axis 2) and ky (axis 3) -> [nt]
                    # Note: ky_flux is [species, moment, field, ky, t]
                    # data.ky_flux[s_idx, moment_idx] is [field, ky, t]
                    flux_sum_t += np.sum(data.ky_flux[s_idx, moment_idx], axis=(0, 1))

                y = flux_sum_t
                x = data.t
                
                # Calculate statistics for selected time window
                if len(t_indices) > 0:
                    y_subset = y[t_indices]
                    mean_val = np.mean(y_subset)
                    std_val = np.std(y_subset)
                    
                    label = f"{label} (Mean: {mean_val:.2e}, Std: {std_val:.2e})"
                    
                    # Plot dashed line for mean
                    # Use the time range from t_indices
                    t_mean_start = x[t_indices[0]]
                    t_mean_end = x[t_indices[-1]]
                    
                    # Plot main line first to get color
                    line, = self.ax.plot(x, y, label=label) # No marker
                    
                    self.ax.plot([t_mean_start, t_mean_end], [mean_val, mean_val],
                                    linestyle='--', color=line.get_color(), linewidth=1.5)
                    
                    # Skip standard plotting
                    return

        self._plot_1d(x, y, label, plot_type)

    def _plot_phi_1d(self, data, label, plot_type, t_indices, t_start, t_end):
        if not hasattr(data, 'kxky_phi'):
            try:
                print(f"Loading big field data for {label}...")
                data.getbigfield()
            except Exception as e:
                print(f"Could not load big field data for {label}: {e}")
                return

        if not hasattr(data, 'kxky_phi'):
            return

        # Extract midplane phi
        ndim = data.kxky_phi.ndim
        i_theta = data.theta_plot * 4 // 8
        
        phi = None
        if ndim == 5: # [2, nr, theta, ny, nt]
             phi = data.kxky_phi[0, :, i_theta, :, :] + 1j * data.kxky_phi[1, :, i_theta, :, :]
        elif ndim == 4: # [nr, theta, ny, nt]
             phi = data.kxky_phi[:, i_theta, :, :]
        else:
             print(f"Unknown kxky_phi shape: {data.kxky_phi.shape}")
             return
             
        # phi shape: [nr, ny, nt]
        # Calculate amplitude squared |phi|^2
        phi_sq = np.abs(phi)**2
        
        if plot_type == "Phi vs ky":
            # Sum over kx (axis 0) -> [ny, nt]
            phi_ky_t = np.sum(phi_sq, axis=0)
            
            # Time average
            if len(t_indices) > 0:
                y = np.mean(phi_ky_t[:, t_indices], axis=1)
                label = f"{label} (Avg: {t_start:.1f}-{t_end:.1f})"
            else:
                y = phi_ky_t[:, -1]
                
            # Plot RMS amplitude: sqrt(Sum |phi|^2)
            y = np.sqrt(y)
            
            x = data.ky
            self._plot_1d(x, y, label, plot_type)
            
        elif plot_type == "Phi vs Time":
            # Sum over kx (axis 0) and ky (axis 1) -> [nt]
            phi_t = np.sum(phi_sq, axis=(0, 1))
            
            # Sqrt for RMS amplitude
            y = np.sqrt(phi_t)
            x = data.t
            
            # Calculate statistics for selected time window
            if len(t_indices) > 0:
                y_subset = y[t_indices]
                mean_val = np.mean(y_subset)
                std_val = np.std(y_subset)
                
                label = f"{label} (Mean: {mean_val:.2e}, Std: {std_val:.2e})"
                
                # Plot main line
                line, = self.ax.plot(x, y, label=label)
                
                # Plot dashed line for mean
                t_mean_start = x[t_indices[0]]
                t_mean_end = x[t_indices[-1]]
                self.ax.plot([t_mean_start, t_mean_end], [mean_val, mean_val],
                                linestyle='--', color=line.get_color(), linewidth=1.5)
            else:
                self.ax.plot(x, y, label=label)

    def _plot_single_case(self, data, label, plot_type, species_override_index=None):
        # Get time indices if relevant
        t_indices = []
        t_start, t_end = 0, 0
        if hasattr(data, 't'):
             t_indices, t_start, t_end = self._get_time_indices(data.t)

        try:
            if plot_type in ["Frequency", "Growth Rate"]:
                self._plot_frequency_growth(data, label, plot_type, t_indices, t_start, t_end)
            elif "Flux" in plot_type:
                self._plot_flux(data, label, plot_type, t_indices, t_start, t_end, species_override_index)
            elif plot_type in ["Phi vs ky", "Phi vs Time"]:
                self._plot_phi_1d(data, label, plot_type, t_indices, t_start, t_end)
            elif plot_type == "Fluctuation 2D":
                self._plot_fluctuation_2d(data, label, plot_type, t_indices, t_start, t_end)
            elif plot_type == "Phi FFT":
                self._plot_phi_fft(data, label, t_indices)
            else:
                # Fallback for unknown types or standard plotting if logic was simple
                pass

        except Exception as e:
            print(f"Error processing data for {label}: {e}")
            traceback.print_exc()
    


if __name__ == "__main__":
    root = tk.Tk()
    app = CGYRO_Comparison(root)
    root.mainloop()
