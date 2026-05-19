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
            initialdir=os.getcwd(),
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
