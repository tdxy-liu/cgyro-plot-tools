"""
Data-export mixin for CGYRO comparison GUI.

Provides:
- Data -> transfer bin to readable
- Export all `bin.cgyro.*` files for selected cases to Origin-friendly text tables
- Special handling for `.cgyro.triad`: split by physical channel into Sheet1..Sheet8
"""

import os
import re
from tkinter import filedialog, messagebox

import numpy as np

try:
    from cgyro_comparison_bootstrap import default_share_dir
except ImportError:
    from .cgyro_comparison_bootstrap import default_share_dir


class CgyroDataExportMixin:
    """Standalone data-export tools mixed into the main GUI controller."""

    _TRIAD_CHANNEL_LABELS = [
        "T_a",
        "N_a_or_Ta_star",
        "d_Ta_dS_dt",
        "dW_em_dt",
        "delta_S_a",
        "D_r_a",
        "D_theta_a",
        "D_c_a",
    ]

    def _clear_current_plot_data(self):
        """Reset cached x-y-z datasets for the next plot action."""
        self._current_plot_xyz_datasets = []

    def _record_current_plot_xyz_dataset(self, label, x, y, z):
        """
        Cache one plotted dataset as x-y-z columns.

        For 2D maps, `x` and `y` may be 1D axes or full mesh grids. `z` is
        stored in the displayed orientation: z[row, column] maps to y, x.
        """
        # Some Matplotlib artists do not retain the original physical axes
        # after rendering (especially pcolormesh/contour/image paths).  Plotting
        # routines can call this hook while they still have clean numerical
        # axes, making Data -> Save current plot data reproduce the actual map.
        try:
            x_arr = np.asarray(x, dtype=float)
            y_arr = np.asarray(y, dtype=float)
            z_arr = np.ma.filled(np.asarray(z, dtype=float), np.nan)
        except Exception:
            return

        if z_arr.size == 0:
            return

        if not hasattr(self, "_current_plot_xyz_datasets"):
            self._current_plot_xyz_datasets = []

        self._current_plot_xyz_datasets.append({
            "label": str(label) if label else "dataset",
            "x": x_arr.copy(),
            "y": y_arr.copy(),
            "z": np.asarray(z_arr, dtype=float).copy(),
        })

    def save_current_plot_data(self):
        """Export currently displayed plot data as vertical x-y-z text columns."""
        datasets = self._collect_current_plot_xyz_datasets()
        if not datasets:
            messagebox.showwarning("Save current plot data", "No plotted data found to export.")
            return

        out_path = filedialog.asksaveasfilename(
            title="Save current plot data",
            defaultextension=".txt",
            filetypes=[
                ("Text files", "*.txt"),
                ("All files", "*.*"),
            ],
            initialdir=self._default_data_export_dir(),
        )
        if not out_path:
            return

        try:
            with open(out_path, "w", encoding="utf-8") as f:
                self._write_current_plot_xyz_file(f, datasets)
            messagebox.showinfo("Save current plot data", f"Plot data saved:\n{out_path}")
        except Exception as e:
            messagebox.showerror("Save current plot data", f"Failed to save plot data:\n{e}")

    def _collect_current_plot_xyz_datasets(self):
        """Collect cached 2D datasets plus visible 1D line/image artists."""
        datasets = []

        # Prefer explicitly cached datasets because they preserve the physical
        # axes chosen by the plotting code.  Artist fallback below is only for
        # older/simple plots that did not call `_record_current_plot_xyz_dataset`.
        for item in getattr(self, "_current_plot_xyz_datasets", []):
            flat = self._flatten_xyz_dataset(
                item.get("x"),
                item.get("y"),
                item.get("z"),
            )
            if flat is not None:
                label, x, y, z = item.get("label", "dataset"), flat[0], flat[1], flat[2]
                datasets.append((str(label), x, y, z))

        fig = getattr(self, "fig", None)
        if fig is None:
            return datasets

        for ax_idx, ax in enumerate(getattr(fig, "axes", [])):
            ax_title = ""
            try:
                ax_title = ax.get_title().strip()
            except Exception:
                ax_title = ""

            for line_idx, line in enumerate(ax.get_lines()):
                try:
                    x = np.asarray(line.get_xdata(), dtype=float).reshape(-1)
                    y = np.asarray(line.get_ydata(), dtype=float).reshape(-1)
                except Exception:
                    continue
                n = min(x.size, y.size)
                if n <= 0:
                    continue
                label = line.get_label()
                if not label or str(label).startswith("_"):
                    label = f"axis{ax_idx + 1}_line{line_idx + 1}"
                # Origin imports a uniform XYZ worksheet more easily than a
                # mixed XY/XYZ file.  For 1D curves, Z=0 is a harmless plane.
                z = np.zeros(n, dtype=float)
                datasets.append((str(label), x[:n], y[:n], z))

            for img_idx, img in enumerate(getattr(ax, "images", [])):
                try:
                    z_img = np.ma.filled(np.asarray(img.get_array(), dtype=float), np.nan)
                except Exception:
                    continue
                if z_img.ndim != 2 or z_img.size <= 0:
                    continue
                try:
                    xmin, xmax, ymin, ymax = [float(v) for v in img.get_extent()]
                except Exception:
                    xmin, xmax = 0.0, float(z_img.shape[1] - 1)
                    ymin, ymax = 0.0, float(z_img.shape[0] - 1)
                x = self._image_axis_centers(xmin, xmax, z_img.shape[1])
                y = self._image_axis_centers(ymin, ymax, z_img.shape[0])
                flat = self._flatten_xyz_dataset(x, y, z_img)
                if flat is None:
                    continue
                label = ax_title or f"axis{ax_idx + 1}_image{img_idx + 1}"
                datasets.append((str(label), flat[0], flat[1], flat[2]))

        return datasets

    @staticmethod
    def _image_axis_centers(vmin, vmax, n):
        """Return pixel-center coordinates for an imshow extent."""
        n = int(n)
        if n <= 1:
            return np.array([(float(vmin) + float(vmax)) * 0.5], dtype=float)
        step = (float(vmax) - float(vmin)) / float(n)
        return float(vmin) + (np.arange(n, dtype=float) + 0.5) * step

    @staticmethod
    def _flatten_xyz_dataset(x, y, z):
        """Normalize 1D or 2D x-y-z data to three equal-length 1D arrays."""
        try:
            x_arr = np.asarray(x, dtype=float)
            y_arr = np.asarray(y, dtype=float)
            z_arr = np.ma.filled(np.asarray(z, dtype=float), np.nan)
        except Exception:
            return None

        if z_arr.ndim == 1:
            x_flat = x_arr.reshape(-1)
            y_flat = y_arr.reshape(-1)
            z_flat = z_arr.reshape(-1)
            n = min(x_flat.size, y_flat.size, z_flat.size)
            if n <= 0:
                return None
            return x_flat[:n], y_flat[:n], z_flat[:n]

        if z_arr.ndim != 2:
            return None

        n_y, n_x = z_arr.shape
        if x_arr.shape == z_arr.shape and y_arr.shape == z_arr.shape:
            xx = x_arr
            yy = y_arr
        elif x_arr.ndim == 1 and y_arr.ndim == 1 and x_arr.size == n_x and y_arr.size == n_y:
            # Standard map convention: z[row, col] == z[y_index, x_index].
            xx, yy = np.meshgrid(x_arr, y_arr)
        elif x_arr.ndim == 1 and y_arr.ndim == 1 and x_arr.size == n_y and y_arr.size == n_x:
            # Accept transposed axis inputs to make export robust to older
            # plotting routines that passed axes in storage order.
            xx, yy = np.meshgrid(y_arr, x_arr)
        else:
            return None

        return xx.reshape(-1), yy.reshape(-1), z_arr.reshape(-1)

    @staticmethod
    def _format_xyz_value(value):
        """Format one numeric text-table value."""
        try:
            v = float(value)
        except Exception:
            return "nan"
        if not np.isfinite(v):
            return "nan"
        return f"{v:.16e}"

    def _write_current_plot_xyz_file(self, f, datasets):
        """Write an Origin-friendly tab-separated XYZ worksheet."""
        multi_dataset = len(datasets) > 1
        if multi_dataset:
            # A Dataset column lets Origin import multiple plotted cases into
            # one worksheet while still allowing filtering/grouping by case.
            f.write("X\tY\tZ\tDataset\n")
        else:
            f.write("X\tY\tZ\n")

        for label, x, y, z in datasets:
            label_txt = self._origin_safe_dataset_label(label)
            n = min(len(x), len(y), len(z))
            for i in range(n):
                row = [
                    self._format_xyz_value(x[i]),
                    self._format_xyz_value(y[i]),
                    self._format_xyz_value(z[i]),
                ]
                if multi_dataset:
                    f.write("\t".join(row) + f"\t{label_txt}\n")
                else:
                    f.write("\t".join(row) + "\n")

    @staticmethod
    def _origin_safe_dataset_label(label):
        """Keep dataset labels on one Origin-importable text token."""
        text = str(label) if label is not None else "dataset"
        text = re.sub(r"[\t\r\n]+", " ", text).strip()
        return text or "dataset"

    @staticmethod
    def _default_data_export_dir():
        """Return preferred data-export directory, with local fallbacks."""
        return default_share_dir()

    def transfer_bin_to_readable(self):
        """
        Export selected-case `bin.cgyro.*` files into readable text tables.

        Behavior:
        - Iterates all bin files in each case.
        - Writes one output folder per case.
        - For `.cgyro.triad`, writes Sheet1..Sheet8 files by physical channel.
        - Keeps all remaining dimensions unchanged (index columns + value columns).
        """
        if not hasattr(self, "cases") or len(self.cases) == 0:
            messagebox.showwarning("Warning", "No loaded cases to export.")
            return

        selected_cases = self._get_selected_case_names()
        if len(selected_cases) == 0:
            messagebox.showwarning("Warning", "No selected cases to export.")
            return

        out_root = filedialog.askdirectory(
            title="Select output directory for readable bin export",
            initialdir=self._default_data_export_dir(),
        )
        if not out_root:
            return

        total_cases = 0
        total_files = 0
        total_sheets = 0
        failed_cases = []

        for case_name in selected_cases:
            data = self.cases.get(case_name, None)
            if data is None:
                failed_cases.append(f"{case_name} (case object missing)")
                continue

            case_dir = self._resolve_case_dir_for_export(data)
            if (not case_dir) or (not os.path.isdir(case_dir)):
                failed_cases.append(f"{case_name} (invalid case directory)")
                continue

            bin_dir = os.path.join(case_dir, "bin")
            if not os.path.isdir(bin_dir):
                failed_cases.append(f"{case_name} (no bin directory)")
                continue

            case_out = os.path.join(out_root, f"{self._sanitize_name(case_name)}_readable")
            os.makedirs(case_out, exist_ok=True)

            n_files, n_sheets = self._export_case_bin_files(data, case_name, bin_dir, case_out)
            total_cases += 1
            total_files += int(n_files)
            total_sheets += int(n_sheets)

        summary = (
            f"Export completed.\n\n"
            f"Cases exported: {total_cases}\n"
            f"Bin files exported: {total_files}\n"
            f"Triad sheets exported: {total_sheets}"
        )
        if failed_cases:
            summary += "\n\nSkipped:\n- " + "\n- ".join(failed_cases)
        messagebox.showinfo("transfer bin to readable", summary)

    def _resolve_case_dir_for_export(self, data):
        """Resolve case directory using plotting helper when available."""
        try:
            if hasattr(self, "_resolve_case_dir"):
                return self._resolve_case_dir(data)
        except Exception:
            pass
        try:
            p = getattr(data, "dir", None) or getattr(data, "path", None)
            if not p:
                return None
            p = str(p)
            if os.path.isfile(p):
                return os.path.dirname(os.path.abspath(p))
            return os.path.abspath(p.rstrip("/\\"))
        except Exception:
            return None

    @staticmethod
    def _sanitize_name(text):
        """Sanitize one text token for safe file/folder names."""
        return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("._") or "case"

    @staticmethod
    def _iter_bin_suffixes(bin_dir):
        """Yield `.cgyro.*` suffixes for files named `bin.cgyro.*`."""
        for fname in sorted(os.listdir(bin_dir)):
            if not fname.startswith("bin.cgyro."):
                continue
            yield fname[len("bin"):]  # e.g. ".cgyro.triad"

    def _export_case_bin_files(self, data, case_name, bin_dir, case_out):
        """Export all bin files for one case. Returns `(n_files, n_triad_sheets)`."""
        n_files = 0
        n_sheets = 0

        for suffix in self._iter_bin_suffixes(bin_dir):
            if suffix == ".cgyro.triad":
                wrote = self._export_triad_by_channel(data, case_name, case_out)
                n_files += 1 if wrote > 0 else 0
                n_sheets += int(wrote)
                continue

            ok = self._export_generic_bin_file(data, suffix, case_out)
            if ok:
                n_files += 1

        return n_files, n_sheets

    def _export_generic_bin_file(self, data, suffix, case_out):
        """
        Export one non-triad bin file as tab-separated table.

        For generic files with unknown shape, this keeps original ordering as 1D:
        columns = `flat_index`, `value` (or `real`,`imag` for complex).
        """
        arr = self._extract_bin_array(data, suffix)
        if arr is None:
            return False

        out_name = self._sanitize_name(suffix.lstrip(".")) + ".txt"
        out_path = os.path.join(case_out, out_name)
        self._write_flat_table(arr, out_path, source_suffix=suffix)
        return True

    def _extract_bin_array(self, data, suffix):
        """
        Read one `bin/out` payload using `data.extract`.

        Strategy:
        - default float read (`cmplx=False`),
        - for `kxky_*` files, fallback to complex read when needed.
        """
        if not hasattr(data, "extract"):
            return None

        try:
            _tmsg, fmt, raw = data.extract(suffix, cmplx=False)
            if fmt != "null":
                arr = np.asarray(raw)
                if arr.size > 0:
                    return arr
        except Exception:
            pass

        # Complex fallback for kxky-like bins.
        if suffix.startswith(".cgyro.kxky_"):
            try:
                _tmsg, fmt, raw = data.extract(suffix, cmplx=True)
                if fmt != "null":
                    arr = np.asarray(raw)
                    if arr.size > 0:
                        return arr
            except Exception:
                pass
        return None

    def _export_triad_by_channel(self, data, case_name, case_out):
        """
        Export `.cgyro.triad` into per-channel "Sheet" files.

        Output:
        - `triad_channels/Sheet1_*.txt` ... `Sheet8_*.txt`
        - each sheet preserves dimensions except channel:
          `[ri, species, radial, n, time, value]`
        """
        triad = self._load_triad_array_for_export(data, case_name)
        if triad is None:
            return 0

        if triad.ndim != 6 or triad.shape[0] != 2:
            print(f"Skip triad export for {case_name}: unsupported shape {triad.shape}")
            return 0

        _ri, _n_species, _n_radial, n_chan, _n_n, _n_t = triad.shape
        n_chan_use = min(int(n_chan), len(self._TRIAD_CHANNEL_LABELS))
        if n_chan_use <= 0:
            return 0

        out_dir = os.path.join(case_out, "triad_channels")
        os.makedirs(out_dir, exist_ok=True)

        wrote = 0
        for ch in range(n_chan_use):
            label = self._TRIAD_CHANNEL_LABELS[ch]
            out_name = f"Sheet{ch+1}_{self._sanitize_name(label)}.txt"
            out_path = os.path.join(out_dir, out_name)
            arr_ch = triad[:, :, :, ch, :, :]  # [ri, species, radial, n, time]
            self._write_nd_table(
                arr_ch,
                out_path,
                dim_labels=["ri", "species", "radial", "n", "time"],
                source_suffix=".cgyro.triad",
                channel_index=ch + 1,
                channel_label=label,
            )
            wrote += 1

        self._write_triad_readme(out_dir, n_chan_use)
        return wrote

    def _load_triad_array_for_export(self, data, case_name):
        """Load triad array with same shape convention used in plotting backend."""
        try:
            if hasattr(self, "_load_triad_only_if_needed"):
                if self._load_triad_only_if_needed(data, case_name):
                    triad = np.asarray(getattr(data, "triad", None))
                    if isinstance(triad, np.ndarray) and triad.size > 0:
                        return triad
        except Exception:
            pass

        arr = self._extract_bin_array(data, ".cgyro.triad")
        if arr is None:
            return None

        # Fallback reshape inference: [2, species, radial, 8, n, time] (Fortran order).
        n_n = int(getattr(data, "n_n", 0))
        n_radial = int(getattr(data, "n_radial", 0))
        n_species = int(getattr(data, "n_species", 0))
        n_chan = 8
        if n_n <= 0 or n_radial <= 0 or n_species <= 0:
            return None

        block = 2 * n_species * n_radial * n_chan * n_n
        if block <= 0:
            return None
        n_time = int(arr.size // block)
        if n_time <= 0:
            return None
        nd = block * n_time
        if arr.size < nd:
            return None
        arr = np.asarray(arr[:nd], dtype=float)
        try:
            return np.reshape(arr, (2, n_species, n_radial, n_chan, n_n, n_time), order="F")
        except Exception:
            return None

    @staticmethod
    def _write_flat_table(arr, out_path, source_suffix=""):
        """Write one flat-index table for arbitrary 1D payload."""
        a = np.asarray(arr).reshape(-1)
        is_complex = np.iscomplexobj(a)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(f"# source: {source_suffix}\n")
            f.write(f"# original_ndim: 1 (raw flatten)\n")
            f.write(f"# original_size: {a.size}\n")
            if is_complex:
                f.write("flat_index\treal\timag\n")
                for i, v in enumerate(a):
                    f.write(f"{i}\t{float(np.real(v)):.16e}\t{float(np.imag(v)):.16e}\n")
            else:
                f.write("flat_index\tvalue\n")
                for i, v in enumerate(a):
                    f.write(f"{i}\t{float(v):.16e}\n")

    @staticmethod
    def _write_nd_table(arr, out_path, dim_labels, source_suffix="", channel_index=None, channel_label=None):
        """
        Write an N-D array to tab-separated rows with explicit index columns.

        This preserves all dimensions without any aggregation.
        """
        a = np.asarray(arr)
        is_complex = np.iscomplexobj(a)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(f"# source: {source_suffix}\n")
            f.write(f"# shape: {tuple(int(x) for x in a.shape)}\n")
            if channel_index is not None:
                f.write(f"# triad_channel_index: {int(channel_index)}\n")
            if channel_label is not None:
                f.write(f"# triad_channel_label: {channel_label}\n")

            cols = list(dim_labels)
            if is_complex:
                cols += ["real", "imag"]
            else:
                cols += ["value"]
            f.write("\t".join(cols) + "\n")

            for idx in np.ndindex(a.shape):
                idx_txt = "\t".join(str(int(i)) for i in idx)
                v = a[idx]
                if is_complex:
                    f.write(f"{idx_txt}\t{float(np.real(v)):.16e}\t{float(np.imag(v)):.16e}\n")
                else:
                    f.write(f"{idx_txt}\t{float(v):.16e}\n")

    def _write_triad_readme(self, out_dir, n_chan_use):
        """Write channel mapping note for triad sheet files."""
        readme = os.path.join(out_dir, "README_channels.txt")
        with open(readme, "w", encoding="utf-8") as f:
            f.write("triad sheet mapping (Origin import helper)\n")
            f.write("Each SheetN file preserves dims: [ri, species, radial, n, time]\n\n")
            for i in range(int(n_chan_use)):
                f.write(f"Sheet{i+1}: {self._TRIAD_CHANNEL_LABELS[i]}\n")
