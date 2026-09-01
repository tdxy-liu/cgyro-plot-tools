"""
Data-export mixin for CGYRO comparison GUI.

Provides:
- Data -> transfer bin to readable
- Export all `bin.cgyro.*` files for selected cases to Origin-friendly text tables
- Special handling for `.cgyro.triad`: split by physical channel into Sheet1..Sheet8
"""

import csv
import datetime
import json
import os
import re
import sys
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

    def batch_export_current_plot_data(self):
        """Export every dataset in the current figure into a new folder."""
        datasets = self._collect_current_plot_xyz_datasets()
        if not datasets:
            messagebox.showwarning(
                "Batch export",
                "No plotted data found to export.",
                parent=getattr(self, "root", None),
            )
            return

        selected_dir = filedialog.askdirectory(
            title="Select folder for current plot data export",
            initialdir=self._default_data_export_dir(),
            mustexist=False,
        )
        if not selected_dir:
            return

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        export_dir = os.path.join(
            selected_dir,
            "cgyro_plot_export_{}".format(timestamp),
        )

        try:
            os.makedirs(export_dir, exist_ok=False)
            dataset_records = []
            combined_path = os.path.join(export_dir, "all_datasets.csv")
            with open(combined_path, "w", encoding="utf-8-sig", newline="") as combined_file:
                combined_writer = csv.writer(combined_file)
                combined_writer.writerow(["Dataset", "X", "Y", "Z"])

                for index, (label, x, y, z) in enumerate(datasets, start=1):
                    safe_label = self._export_filename_label(label)
                    filename = "dataset_{:03d}_{}.csv".format(index, safe_label)
                    path = os.path.join(export_dir, filename)
                    n = min(len(x), len(y), len(z))
                    with open(path, "w", encoding="utf-8-sig", newline="") as data_file:
                        writer = csv.writer(data_file)
                        writer.writerow(["X", "Y", "Z"])
                        for row_index in range(n):
                            row = [x[row_index], y[row_index], z[row_index]]
                            writer.writerow(row)
                            combined_writer.writerow([label, *row])

                    dataset_records.append({
                        "label": str(label),
                        "file": filename,
                        "points": int(n),
                    })

            metadata = self._current_plot_export_metadata(dataset_records)
            metadata["combined_file"] = "all_datasets.csv"
            metadata_path = os.path.join(export_dir, "plot_metadata.json")
            with open(metadata_path, "w", encoding="utf-8") as metadata_file:
                json.dump(metadata, metadata_file, ensure_ascii=False, indent=2)

            messagebox.showinfo(
                "Batch export complete",
                "Exported {} dataset(s) to:\n{}".format(len(datasets), export_dir),
                parent=getattr(self, "root", None),
            )
        except Exception as exc:
            messagebox.showerror(
                "Batch export",
                "Failed to export current plot data:\n{}".format(exc),
                parent=getattr(self, "root", None),
            )

    @staticmethod
    def _export_filename_label(label):
        """Convert a plotted label into a safe, readable filename fragment."""
        text = str(label) if label is not None else "dataset"
        text = re.sub(r"[^A-Za-z0-9._+-]+", "_", text).strip("._")
        return text[:100] or "dataset"

    def _current_plot_export_metadata(self, dataset_records):
        """Collect lightweight context alongside a batch plot-data export."""
        plot_type = ""
        display_plot_type = ""
        try:
            built = self._build_effective_plot_type()
            if isinstance(built, tuple) and len(built) >= 3:
                plot_type = str(built[1])
                display_plot_type = str(built[2])
        except Exception:
            pass
        if not plot_type:
            try:
                plot_type = str(self.plot_type_var.get()).strip()
                display_plot_type = plot_type
            except Exception:
                plot_type = ""

        axes = []
        for index, axis in enumerate(getattr(getattr(self, "fig", None), "axes", []), start=1):
            try:
                axes.append({
                    "index": index,
                    "title": axis.get_title(),
                    "xlabel": axis.get_xlabel(),
                    "ylabel": axis.get_ylabel(),
                })
            except Exception:
                pass

        case_names = []
        try:
            case_names = [str(name) for name in self._get_selected_case_names()]
        except Exception:
            pass

        metadata = {
            "tool": "CGYRO Comparison Tool",
            "exported_at": datetime.datetime.now().astimezone().isoformat(),
            "plot_type": plot_type,
            "display_plot_type": display_plot_type,
            "cases": case_names,
            "axes": axes,
            "datasets": dataset_records,
        }
        try:
            metadata["workspace_state"] = self._get_workspace_state()
        except Exception:
            pass
        return metadata

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

            suffixes = self._collect_bin_suffixes(case_dir)
            if not suffixes:
                failed_cases.append(f"{case_name} (no bin.cgyro.* files)")
                continue

            case_out = os.path.join(out_root, f"{self._sanitize_name(case_name)}_readable")
            os.makedirs(case_out, exist_ok=True)

            n_files, n_sheets = self._export_case_bin_files(data, case_name, suffixes, case_out)
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
        if not os.path.isdir(bin_dir):
            return
        for fname in sorted(os.listdir(bin_dir)):
            if not fname.startswith("bin.cgyro."):
                continue
            yield fname[len("bin"):]  # e.g. ".cgyro.triad"

    def _collect_bin_suffixes(self, case_dir):
        """
        Return available binary diagnostic suffixes for one case.

        CGYRO normally writes flat files like `case/bin.cgyro.freq`.  A few
        local workflows have used `case/bin/bin.cgyro.freq`, so keep that
        fallback too, while de-duplicating suffixes.
        """
        suffixes = []
        seen = set()
        for search_dir in (case_dir, os.path.join(case_dir, "bin")):
            for suffix in self._iter_bin_suffixes(search_dir):
                if suffix in seen:
                    continue
                suffixes.append(suffix)
                seen.add(suffix)
        return suffixes

    def _export_case_bin_files(self, data, case_name, suffixes, case_out):
        """Export all bin files for one case. Returns `(n_files, n_triad_sheets)`."""
        n_files = 0
        n_sheets = 0

        for suffix in suffixes:
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

        Known CGYRO diagnostics are reshaped back to named dimensions so Origin
        can filter and plot by coordinates.  Unknown files still export as a
        flat, first-row-header table.
        """
        arr = self._extract_bin_array(data, suffix)
        if arr is None:
            return False

        out_name = self._sanitize_name(suffix.lstrip(".")) + ".txt"
        out_path = os.path.join(case_out, out_name)

        if self._write_structured_origin_table(data, suffix, arr, out_path):
            return True

        self._write_flat_table(arr, out_path)
        return True

    def _write_structured_origin_table(self, data, suffix, arr, out_path):
        """Write known CGYRO diagnostics as coordinate-rich Origin tables."""
        if suffix == ".cgyro.freq":
            return self._write_freq_origin_table(data, arr, out_path)

        if suffix in (".cgyro.fullt", ".cgyro.fullt_asym"):
            return self._write_fullt_origin_table(data, suffix, arr, out_path)

        if suffix in (".cgyro.kxky_phi", ".cgyro.kxky_apar", ".cgyro.kxky_bpar"):
            return self._write_kxky_origin_table(data, arr, out_path, species_dependent=False)

        if suffix in (".cgyro.kxky_n", ".cgyro.kxky_e", ".cgyro.kxky_v"):
            return self._write_kxky_origin_table(data, arr, out_path, species_dependent=True)

        if suffix in (".cgyro.phib", ".cgyro.aparb", ".cgyro.bparb"):
            return self._write_ballooning_origin_table(data, arr, out_path)

        if suffix in (".cgyro.lky_flux_n", ".cgyro.lky_flux_e", ".cgyro.lky_flux_v"):
            return self._write_lky_flux_origin_table(data, arr, out_path)

        if suffix in (".cgyro.ky_flux", ".cgyro.ky_cflux"):
            return self._write_ky_flux_origin_table(data, arr, out_path)

        return False

    def _write_freq_origin_table(self, data, arr, out_path):
        """Write `.cgyro.freq` as ky/time rows with omega and gamma columns."""
        raw = np.asarray(arr).reshape(-1)
        n_n = self._positive_int(getattr(data, "n_n", 0))
        n_time = self._positive_int(getattr(data, "n_time", 0))
        if n_n <= 0 or n_time <= 0:
            return False

        nd = 2 * n_n * n_time
        if raw.size < nd:
            return False
        try:
            freq = np.reshape(np.asarray(raw[:nd], dtype=float), (2, n_n, n_time), order="F")
        except Exception:
            return False

        ky_axis = self._ky_axis(data, n_n)
        time_axis = self._time_axis(data, n_time)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("ky_index\tky\ttime_index\ttime\tomega\tgamma\n")
            for iky in range(n_n):
                ky = self._axis_value(ky_axis, iky)
                for it in range(n_time):
                    row = [
                        str(iky),
                        self._format_origin_value(ky),
                        str(it),
                        self._format_origin_value(self._axis_value(time_axis, it)),
                        self._format_origin_value(freq[0, iky, it]),
                        self._format_origin_value(freq[1, iky, it]),
                    ]
                    f.write("\t".join(row) + "\n")
        return True

    def _write_fullt_origin_table(self, data, suffix, arr, out_path):
        """Write FULLT/FULLT_ASYM as source/target coordinate rows."""
        raw = np.asarray(arr)
        if raw.size <= 0:
            return False
        if np.iscomplexobj(raw):
            raw = np.column_stack((np.real(raw).reshape(-1), np.imag(raw).reshape(-1))).reshape(-1)
        raw = np.asarray(raw, dtype=float).reshape(-1)

        n_n = self._positive_int(getattr(data, "n_n", 0))
        n_radial = self._positive_int(getattr(data, "n_radial", 0))
        n_time = self._positive_int(getattr(data, "n_time", 0))
        if n_n <= 0 or n_radial <= 0:
            return False

        n_signed = 2 * n_n - 1
        spatial = n_radial * n_signed * n_n
        n_channel, n_time, n_source_kx, legacy_complex = self._infer_fullt_shape(
            raw.size,
            spatial,
            n_radial,
            n_time,
        )
        if n_channel <= 0 or n_time <= 0 or n_source_kx <= 0:
            return False

        scale = 2 if legacy_complex else 1
        nd = scale * n_radial * n_signed * n_channel * n_source_kx * n_n * n_time
        try:
            if legacy_complex:
                shaped_raw = np.reshape(
                    raw[:nd],
                    (2, n_radial, n_signed, n_channel, n_source_kx, n_n, n_time),
                    order="F",
                )
                shaped = np.transpose(shaped_raw, (5, 4, 1, 2, 3, 6, 0))
            else:
                shaped_raw = np.reshape(
                    raw[:nd],
                    (n_radial, n_signed, n_channel, n_source_kx, n_n, n_time),
                    order="F",
                )
                shaped = np.zeros((n_n, n_source_kx, n_radial, n_signed, n_channel, n_time, 2), dtype=float)
                shaped[:, :, :, :, :, :, 0] = np.transpose(shaped_raw, (4, 3, 0, 1, 2, 5))
        except Exception:
            return False

        dim_specs = [
            ("source_ky_index", "source_ky", self._fullt_source_ky_axis(data, n_n)),
            ("source_kx_index", "source_kx", self._fullt_source_kx_axis(data, n_source_kx)),
            ("target_kx_index", "target_kx", self._fullt_target_kx_axis(data, n_radial)),
            ("target_ky_index", "target_ky", self._fullt_target_ky_axis(data, n_n, n_signed)),
            ("channel", None, None),
            ("time_index", "time", self._time_axis(data, n_time)),
            ("real_imag", None, None),
        ]
        self._write_origin_nd_table(shaped, out_path, dim_specs)
        return True

    def _write_kxky_origin_table(self, data, arr, out_path, species_dependent=False):
        """Write kxky fields as radial/theta/species/ky/time rows."""
        payload = self._complex_payload_from_raw(arr)
        if payload is None:
            return False

        n_radial = self._positive_int(getattr(data, "n_radial", 0))
        theta_plot = self._positive_int(getattr(data, "theta_plot", 0))
        n_species = self._positive_int(getattr(data, "n_species", 1))
        n_n = self._positive_int(getattr(data, "n_n", 0))
        n_time = self._positive_int(getattr(data, "n_time", 0))
        if n_radial <= 0 or theta_plot <= 0 or n_n <= 0:
            return False

        spatial = n_radial * theta_plot * n_n
        if species_dependent:
            spatial *= max(1, n_species)
        n_channel, n_time = self._infer_channel_time(payload.size, spatial, n_time, allowed_channels=(1,))
        if n_channel <= 0 or n_time <= 0:
            return False

        if species_dependent:
            shape = (n_radial, theta_plot, n_species, n_n, n_time)
            dim_specs = [
                ("radial_index", "kx", self._full_kx_axis(data, n_radial)),
                ("theta_index", "theta", self._theta_plot_axis(data, theta_plot)),
                ("species", None, None),
                ("ky_index", "ky", self._ky_axis(data, n_n)),
                ("time_index", "time", self._time_axis(data, n_time)),
            ]
        else:
            shape = (n_radial, theta_plot, n_n, n_time)
            dim_specs = [
                ("radial_index", "kx", self._full_kx_axis(data, n_radial)),
                ("theta_index", "theta", self._theta_plot_axis(data, theta_plot)),
                ("ky_index", "ky", self._ky_axis(data, n_n)),
                ("time_index", "time", self._time_axis(data, n_time)),
            ]

        nd = int(np.prod(shape))
        try:
            shaped = np.reshape(payload[:nd], shape, order="F")
        except Exception:
            return False
        self._write_origin_nd_table(shaped, out_path, dim_specs)
        return True

    def _write_ballooning_origin_table(self, data, arr, out_path):
        """Write phib/aparb/bparb as ballooning-point/time rows."""
        payload = self._complex_payload_from_raw(arr)
        if payload is None:
            return False

        n_theta = self._positive_int(getattr(data, "n_theta", 0))
        n_radial = self._positive_int(getattr(data, "n_radial", 0))
        n_time = self._positive_int(getattr(data, "n_time", 0))
        n_point = n_theta * n_radial
        if n_point <= 0:
            return False

        n_channel, n_time = self._infer_channel_time(payload.size, n_point, n_time, allowed_channels=(1,))
        if n_channel <= 0 or n_time <= 0:
            return False
        nd = n_point * n_time
        try:
            shaped = np.reshape(payload[:nd], (n_point, n_time), order="F")
        except Exception:
            return False

        dim_specs = [
            ("ballooning_point", None, None),
            ("time_index", "time", self._time_axis(data, n_time)),
        ]
        self._write_origin_nd_table(shaped, out_path, dim_specs)
        return True

    def _write_lky_flux_origin_table(self, data, arr, out_path):
        """Write lky flux diagnostics as global/species/field/ky/time rows."""
        payload = self._complex_payload_from_raw(arr)
        if payload is None:
            return False

        n_global = self._positive_int(getattr(data, "n_global", 0)) + 1
        n_species = self._positive_int(getattr(data, "n_species", 0))
        n_field = self._positive_int(getattr(data, "n_field", 0))
        n_n = self._positive_int(getattr(data, "n_n", 0))
        n_time = self._positive_int(getattr(data, "n_time", 0))
        spatial = n_global * n_species * n_field * n_n
        if spatial <= 0:
            return False

        n_channel, n_time = self._infer_channel_time(payload.size, spatial, n_time, allowed_channels=(1,))
        if n_channel <= 0 or n_time <= 0:
            return False

        shape = (n_global, n_species, n_field, n_n, n_time)
        nd = int(np.prod(shape))
        try:
            shaped = np.reshape(payload[:nd], shape, order="F")
        except Exception:
            return False
        dim_specs = [
            ("global_index", None, None),
            ("species", None, None),
            ("field", None, None),
            ("ky_index", "ky", self._ky_axis(data, n_n)),
            ("time_index", "time", self._time_axis(data, n_time)),
        ]
        self._write_origin_nd_table(shaped, out_path, dim_specs)
        return True

    def _write_ky_flux_origin_table(self, data, arr, out_path):
        """Write ky flux diagnostics as species/moment/field/ky/time rows."""
        raw = np.asarray(arr).reshape(-1)
        n_species = self._positive_int(getattr(data, "n_species", 0))
        n_field = self._positive_int(getattr(data, "n_field", 0))
        n_n = self._positive_int(getattr(data, "n_n", 0))
        n_time = self._positive_int(getattr(data, "n_time", 0))
        base = n_species * n_field * n_n * n_time
        if base <= 0 or raw.size < base:
            return False

        n_moment = raw.size // base
        nd = base * n_moment
        if n_moment <= 0:
            return False
        try:
            shaped = np.reshape(np.asarray(raw[:nd], dtype=float), (n_species, n_moment, n_field, n_n, n_time), order="F")
        except Exception:
            return False

        dim_specs = [
            ("species", None, None),
            ("moment", None, None),
            ("field", None, None),
            ("ky_index", "ky", self._ky_axis(data, n_n)),
            ("time_index", "time", self._time_axis(data, n_time)),
        ]
        self._write_origin_nd_table(shaped, out_path, dim_specs)
        return True

    @staticmethod
    def _positive_int(value):
        """Return a positive int, or 0 when metadata is absent/invalid."""
        try:
            ivalue = int(value)
        except Exception:
            return 0
        return ivalue if ivalue > 0 else 0

    @staticmethod
    def _complex_payload_from_raw(arr):
        """
        Return a 1D complex payload from CGYRO binary data.

        Direct binary reads use complex dtype for known complex diagnostics.
        If an older reader has already supplied an interleaved real array,
        recover `real + i*imag` here.
        """
        a = np.asarray(arr).reshape(-1)
        if a.size <= 0:
            return None
        if np.iscomplexobj(a):
            return np.asarray(a, dtype=complex)
        if a.size % 2 != 0:
            return None
        real = np.asarray(a[0::2], dtype=float)
        imag = np.asarray(a[1::2], dtype=float)
        return real + 1j * imag

    @staticmethod
    def _infer_channel_time(total_size, spatial_size, metadata_time, allowed_channels=(1,)):
        """Infer `(n_channel, n_time)` from payload length and metadata."""
        total_size = int(total_size)
        spatial_size = int(spatial_size)
        metadata_time = int(metadata_time) if metadata_time is not None else 0
        if total_size <= 0 or spatial_size <= 0:
            return 0, 0

        for n_channel in allowed_channels:
            block = spatial_size * int(n_channel)
            if block > 0 and total_size % block == 0:
                n_time = total_size // block
                if metadata_time <= 0 or n_time == metadata_time:
                    return int(n_channel), int(n_time)

        # Some incomplete/running cases have a stale n_time in the metadata.
        # Prefer the cleanest divisibility over failing the export.
        for n_channel in allowed_channels:
            block = spatial_size * int(n_channel)
            if block > 0 and total_size % block == 0:
                return int(n_channel), int(total_size // block)

        return 0, 0

    @staticmethod
    def _infer_fullt_shape(total_size, spatial_size, n_radial, metadata_time):
        """Infer `(n_channel, n_time, n_source_kx, legacy_complex)` for FULLT payloads."""
        total_size = int(total_size)
        spatial_size = int(spatial_size)
        n_radial = int(n_radial)
        metadata_time = int(metadata_time) if metadata_time is not None else 0
        if total_size <= 0 or spatial_size <= 0 or n_radial <= 0:
            return 0, 0, 0, False

        source_kx_candidates = [n_radial, 1] if n_radial > 1 else [1]
        channel_candidates = [2, 1]

        for scale in (1, 2):
            for n_source_kx in source_kx_candidates:
                for n_channel in channel_candidates:
                    block = scale * spatial_size * n_source_kx * n_channel
                    if block <= 0:
                        continue
                    if metadata_time > 0 and total_size == block * metadata_time:
                        return int(n_channel), int(metadata_time), int(n_source_kx), bool(scale == 2)

        for scale in (1, 2):
            for n_source_kx in source_kx_candidates:
                for n_channel in channel_candidates:
                    block = scale * spatial_size * n_source_kx * n_channel
                    if block > 0 and total_size % block == 0:
                        return int(n_channel), int(total_size // block), int(n_source_kx), bool(scale == 2)

        return 0, 0, 0, False

    @staticmethod
    def _format_origin_value(value):
        """Format one Origin table numeric cell."""
        try:
            v = float(value)
        except Exception:
            return "nan"
        if not np.isfinite(v):
            return "nan"
        return f"{v:.16e}"

    @staticmethod
    def _axis_value(axis, index):
        """Return finite axis value at index, or NaN when unavailable."""
        if axis is None:
            return np.nan
        arr = np.asarray(axis, dtype=float).reshape(-1)
        if 0 <= int(index) < arr.size:
            return arr[int(index)]
        return np.nan

    def _time_axis(self, data, n_time):
        """Return time axis with a length compatible with exported data."""
        t = np.asarray(getattr(data, "t", []), dtype=float).reshape(-1)
        if t.size >= n_time:
            return t[:n_time]
        return np.arange(n_time, dtype=float)

    def _ky_axis(self, data, n_n):
        """Return stored source-ky axis."""
        ky = np.asarray(getattr(data, "kynorm", getattr(data, "ky", [])), dtype=float).reshape(-1)
        if ky.size >= n_n and np.all(np.isfinite(ky[:n_n])):
            return ky[:n_n]
        return np.arange(n_n, dtype=float)

    def _theta_plot_axis(self, data, theta_plot):
        """Return theta-plot axis."""
        theta = np.asarray(getattr(data, "thetap", []), dtype=float).reshape(-1)
        if theta.size >= theta_plot and np.all(np.isfinite(theta[:theta_plot])):
            return theta[:theta_plot]

        theta_full = np.asarray(getattr(data, "theta", []), dtype=float).reshape(-1)
        if theta_full.size >= theta_plot and theta_plot > 0:
            if theta_full.size == theta_plot:
                return theta_full
            step = max(1, theta_full.size // theta_plot)
            idx = np.arange(theta_plot, dtype=int) * step
            idx = np.clip(idx, 0, theta_full.size - 1)
            return theta_full[idx]

        return np.arange(theta_plot, dtype=float)

    def _full_kx_axis(self, data, n_radial):
        """Return kx axis including the leftmost CGYRO radial storage bin."""
        p = np.asarray(getattr(data, "p", []), dtype=float).reshape(-1)
        length = float(getattr(data, "length", 0.0) or 0.0)
        if p.size >= n_radial and np.isfinite(length) and abs(length) > 1.0e-12:
            return 2.0 * np.pi * p[:n_radial] / length

        kx = np.asarray(getattr(data, "kxnorm", getattr(data, "kx", [])), dtype=float).reshape(-1)
        if kx.size == n_radial:
            return kx
        if kx.size == n_radial - 1 and n_radial > 1:
            if np.isfinite(length) and abs(length) > 1.0e-12:
                dkx = 2.0 * np.pi / length
            else:
                dkx = 1.0
            p_index = np.arange(n_radial, dtype=float) - (n_radial // 2)
            return p_index * dkx
        return np.arange(n_radial, dtype=float)

    def _fullt_target_kx_axis(self, data, n_radial):
        """Return target kx axis for FULLT/FULLT_ASYM diagnostics."""
        return self._full_kx_axis(data, n_radial)

    def _fullt_target_ky_axis(self, data, n_n, n_signed):
        """Return signed target-ky axis for FULLT/FULLT_ASYM diagnostics."""
        ky_native = np.asarray(getattr(data, "full_t_target_ky", []), dtype=float).reshape(-1)
        if ky_native.size >= n_signed and np.all(np.isfinite(ky_native[:n_signed])):
            return ky_native[:n_signed]

        ky_source = self._ky_axis(data, n_n)
        if n_n > 1 and ky_source.size > 1:
            dky = float(ky_source[1])
            return np.arange(-(n_n - 1), n_n, dtype=float) * dky
        if n_signed == 1:
            return np.asarray([0.0], dtype=float)
        return np.arange(-(n_signed // 2), n_signed - (n_signed // 2), dtype=float)

    def _fullt_source_ky_axis(self, data, n_source):
        """Return source-ky axis for FULLT/FULLT_ASYM diagnostics."""
        ky_native = np.asarray(getattr(data, "full_t_source_ky", []), dtype=float).reshape(-1)
        if ky_native.size >= n_source and np.all(np.isfinite(ky_native[:n_source])):
            return ky_native[:n_source]
        return self._ky_axis(data, n_source)

    def _fullt_source_kx_axis(self, data, n_source):
        """Return source-kx axis for FULLT/FULLT_ASYM diagnostics."""
        if n_source == 1:
            return np.asarray([0.0], dtype=float)

        kx_native = np.asarray(getattr(data, "full_t_source_kx", []), dtype=float).reshape(-1)
        if kx_native.size >= n_source and np.all(np.isfinite(kx_native[:n_source])):
            return kx_native[:n_source]

        kx_index = np.asarray(getattr(data, "full_t_source_kx_index", []), dtype=float).reshape(-1)
        if kx_index.size >= n_source and np.all(np.isfinite(kx_index[:n_source])):
            return kx_index[:n_source]

        p = np.asarray(getattr(data, "p", []), dtype=float).reshape(-1)
        length = float(getattr(data, "length", 0.0) or 0.0)
        if p.size >= n_source and np.all(np.isfinite(p[:n_source])):
            if np.isfinite(length) and abs(length) > 1.0e-12:
                return 2.0 * np.pi * p[:n_source] / length
            return p[:n_source]

        return np.arange(n_source, dtype=float) - (n_source // 2)

    def _write_origin_nd_table(self, arr, out_path, dim_specs):
        """
        Write an N-D array as an Origin-friendly long table.

        `dim_specs` entries are `(index_column, coordinate_column, coordinate_axis)`.
        Coordinate columns are optional but make Origin filtering and plotting
        much less painful.
        """
        a = np.asarray(arr)
        is_complex = np.iscomplexobj(a)
        cols = []
        for index_col, coord_col, _axis in dim_specs:
            cols.append(index_col)
            if coord_col:
                cols.append(coord_col)
        if is_complex:
            cols += ["real", "imag", "abs", "phase"]
        else:
            cols += ["value"]

        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\t".join(cols) + "\n")
            for idx in np.ndindex(a.shape):
                row = []
                for dim, (index_col, coord_col, axis) in enumerate(dim_specs):
                    del index_col
                    row.append(str(int(idx[dim])))
                    if coord_col:
                        row.append(self._format_origin_value(self._axis_value(axis, idx[dim])))
                v = a[idx]
                if is_complex:
                    row.extend([
                        self._format_origin_value(np.real(v)),
                        self._format_origin_value(np.imag(v)),
                        self._format_origin_value(np.abs(v)),
                        self._format_origin_value(np.angle(v)),
                    ])
                else:
                    row.append(self._format_origin_value(v))
                f.write("\t".join(row) + "\n")

    def _extract_bin_array(self, data, suffix):
        """
        Read one `bin/out` payload using `data.extract`.

        pygacode.cgyro.data exposes `extract(f)` and returns the raw real
        stream.  Known complex diagnostics are converted later from CGYRO's
        packed real/imag layout, not by passing a `cmplx` keyword.
        """
        arr = self._extract_payload_compat(data, suffix)
        if arr is not None:
            return arr
        return None

    def _extract_payload_compat(self, data, suffix):
        """Read one raw payload using pygacode's `extract(f)` API."""
        if not hasattr(data, "extract"):
            return self._read_payload_direct(data, suffix)

        try:
            _tmsg, fmt, raw = data.extract(suffix)
            if fmt != "null":
                arr = np.asarray(raw)
                if arr.size > 0:
                    return arr
        except Exception:
            pass
        return self._read_payload_direct(data, suffix)

    def _read_payload_direct(self, data, suffix):
        """Directly read `bin/out` payload when `data.extract` cannot."""
        case_dir = self._resolve_case_dir_for_export(data)
        if not case_dir:
            return None

        bin_paths = [
            os.path.join(case_dir, "bin" + suffix),
            os.path.join(case_dir, "bin", "bin" + suffix),
        ]
        out_path = os.path.join(case_dir, "out" + suffix)
        for bin_path in bin_paths:
            if not os.path.isfile(bin_path):
                continue
            try:
                arr = np.fromfile(bin_path, dtype=self._real_dtype_for_data(data))
            except Exception:
                continue
            if arr.size > 0:
                return arr

        if os.path.isfile(out_path):
            try:
                arr = np.fromfile(out_path, dtype=float, sep=" ")
            except Exception:
                return None
            return arr if arr.size > 0 else None
        return None

    @staticmethod
    def _real_dtype_for_data(data):
        """Return the real dtype used by the loaded CGYRO data object.

        Recent pygacode exposes ``BYTE`` at module scope rather than on each
        ``cgyrodata`` instance.  Prefer the current case precision flag when
        available, then support both pygacode layouts, and finally use CGYRO's
        single-precision default instead of silently choosing float64.
        """
        case_dir = getattr(data, "dir", None) or getattr(data, "path", None) or ""
        equilibrium_path = os.path.join(str(case_dir), "out.cgyro.equilibrium")
        if os.path.isfile(equilibrium_path):
            try:
                equilibrium = np.loadtxt(equilibrium_path, dtype=float, ndmin=1)
                if equilibrium.size > 0:
                    flag = float(equilibrium[-1])
                    flag_int = int(round(flag))
                    if flag_int in (0, 1) and abs(flag - flag_int) <= 1.0e-8:
                        return np.dtype("float64" if flag_int == 1 else "float32")
            except (OSError, ValueError):
                pass

        dtype = getattr(data, "BYTE", None)
        try:
            dt = np.dtype(dtype)
            if dt.kind == "f" and dt.itemsize in (4, 8):
                return dt
        except (TypeError, ValueError):
            pass

        data_module = sys.modules.get(type(data).__module__)
        try:
            dt = np.dtype(getattr(data_module, "BYTE", None))
            if dt.kind == "f" and dt.itemsize in (4, 8):
                return dt
        except (TypeError, ValueError):
            pass

        return np.dtype("float32")

    @staticmethod
    def _complex_dtype_for_data(data):
        """Return complex dtype matching the loaded CGYRO precision."""
        dtype = getattr(data, "CBYTE", None)
        if dtype is not None:
            try:
                dt = np.dtype(dtype)
                if np.issubdtype(dt, np.complexfloating):
                    return dt
            except Exception:
                pass
        real_dtype = CgyroDataExportMixin._real_dtype_for_data(data)
        return np.dtype("complex64") if real_dtype == np.dtype("float32") else np.dtype("complex128")

    @staticmethod
    def _suffix_is_complex_payload(suffix):
        """True for CGYRO diagnostics written by bcomplex writers."""
        suffix = str(suffix)
        exact = {
            ".cgyro.triad",
            ".cgyro.phib",
            ".cgyro.aparb",
            ".cgyro.bparb",
        }
        if suffix in exact:
            return True
        prefixes = (
            ".cgyro.kxky_",
            ".cgyro.lky_flux_",
        )
        return any(suffix.startswith(prefix) for prefix in prefixes)

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
            if arr_ch.ndim == 5 and arr_ch.shape[0] == 2:
                arr_complex = arr_ch[0] + 1j * arr_ch[1]
                dim_specs = [
                    ("species", None, None),
                    ("radial_index", "kx", self._full_kx_axis(data, arr_complex.shape[1])),
                    ("ky_index", "ky", self._ky_axis(data, arr_complex.shape[2])),
                    ("time_index", "time", self._time_axis(data, arr_complex.shape[3])),
                ]
                self._write_origin_nd_table(arr_complex, out_path, dim_specs)
            else:
                self._write_nd_table(
                    arr_ch,
                    out_path,
                    dim_labels=["ri", "species", "radial", "n", "time"],
                )
            wrote += 1

        self._write_triad_readme(out_dir, n_chan_use)
        return wrote

    def _load_triad_array_for_export(self, data, case_name):
        """Load triad array with same shape convention used in plotting backend."""
        triad_raw = getattr(data, "triad", None)
        if triad_raw is not None:
            triad = np.asarray(triad_raw)
            if triad.ndim > 0 and triad.size > 0:
                return triad

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

        if np.iscomplexobj(arr):
            arr = np.column_stack(
                (np.real(arr).reshape(-1), np.imag(arr).reshape(-1))
            ).reshape(-1)

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
        del source_suffix
        a = np.asarray(arr).reshape(-1)
        is_complex = np.iscomplexobj(a)
        with open(out_path, "w", encoding="utf-8") as f:
            if is_complex:
                f.write("flat_index\treal\timag\tabs\tphase\n")
                for i, v in enumerate(a):
                    row = [
                        str(i),
                        CgyroDataExportMixin._format_origin_value(np.real(v)),
                        CgyroDataExportMixin._format_origin_value(np.imag(v)),
                        CgyroDataExportMixin._format_origin_value(np.abs(v)),
                        CgyroDataExportMixin._format_origin_value(np.angle(v)),
                    ]
                    f.write("\t".join(row) + "\n")
            else:
                f.write("flat_index\tvalue\n")
                for i, v in enumerate(a):
                    f.write(f"{i}\t{CgyroDataExportMixin._format_origin_value(v)}\n")

    @staticmethod
    def _write_nd_table(arr, out_path, dim_labels, source_suffix="", channel_index=None, channel_label=None):
        """
        Write an N-D array to tab-separated rows with explicit index columns.

        This preserves all dimensions without any aggregation.
        """
        a = np.asarray(arr)
        is_complex = np.iscomplexobj(a)
        with open(out_path, "w", encoding="utf-8") as f:
            del source_suffix, channel_index, channel_label

            cols = list(dim_labels)
            if is_complex:
                cols += ["real", "imag", "abs", "phase"]
            else:
                cols += ["value"]
            f.write("\t".join(cols) + "\n")

            for idx in np.ndindex(a.shape):
                idx_txt = "\t".join(str(int(i)) for i in idx)
                v = a[idx]
                if is_complex:
                    row = [
                        idx_txt,
                        CgyroDataExportMixin._format_origin_value(np.real(v)),
                        CgyroDataExportMixin._format_origin_value(np.imag(v)),
                        CgyroDataExportMixin._format_origin_value(np.abs(v)),
                        CgyroDataExportMixin._format_origin_value(np.angle(v)),
                    ]
                    f.write("\t".join(row) + "\n")
                else:
                    f.write(f"{idx_txt}\t{CgyroDataExportMixin._format_origin_value(v)}\n")

    def _write_triad_readme(self, out_dir, n_chan_use):
        """Write channel mapping note for triad sheet files."""
        readme = os.path.join(out_dir, "README_channels.txt")
        with open(readme, "w", encoding="utf-8") as f:
            f.write("triad sheet mapping (Origin import helper)\n")
            f.write("Each SheetN file uses first-row column headers and columns:\n")
            f.write("species, radial_index, kx, ky_index, ky, time_index, time, real, imag, abs, phase\n\n")
            for i in range(int(n_chan_use)):
                f.write(f"Sheet{i+1}: {self._TRIAD_CHANNEL_LABELS[i]}\n")
