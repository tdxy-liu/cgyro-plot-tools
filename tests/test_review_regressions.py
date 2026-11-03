"""Regression cases found during the 2026-10-29 code review (offline/headless)."""

import json
import queue
from pathlib import Path
import tempfile
import subprocess
import sys
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np
from matplotlib.figure import Figure

import cgyro_update as updater
from cgyro_comparison_ui import CgyroUiMixin
from cgyro_comparison_plotting import Plotting
from cgyro_data_export import CgyroDataExportMixin


def variable(value):
    return SimpleNamespace(get=lambda: value)


class Listbox:
    def __init__(self, names=()):
        self.names = list(names)
        self.selected = set()

    def size(self):
        return len(self.names)

    def get(self, index):
        return self.names[index]

    def delete(self, first, last=None):
        self.names.clear()
        self.selected.clear()

    def insert(self, index, name):
        self.names.append(name)

    def curselection(self):
        return tuple(sorted(self.selected))

    def selection_clear(self, first, last=None):
        self.selected.clear()

    def selection_set(self, index):
        self.selected.add(index)


class WorkspaceHarness(CgyroUiMixin):
    def __init__(self):
        self.original = SimpleNamespace(dir="original")
        self.cases = {"original": self.original}
        self.case_listbox = Listbox(self.cases)
        self._locked_case_names = []
        self._fluc_advanced_case_values = {}
        self._fluc_advanced_seeded = False
        self.root = Mock()
        for name in ("_update_species_list", "_restore_workspace_state",
                     "update_options", "_ensure_fluc_advanced_case_values",
                     "_refresh_case_summary", "_stop_animation", "clear_plot"):
            setattr(self, name, Mock())
        self._case_selection_lock_enabled = lambda: False


class WorkspaceReviewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "workspace.json"

    def test_failed_serialization_preserves_previous_workspace(self):
        self.path.write_text('{"cases": []}', encoding="utf-8")
        original = self.path.read_bytes()
        with self.assertRaises(TypeError):
            CgyroUiMixin._write_workspace_file(self.path, {"bad": object()})
        self.assertEqual(self.path.read_bytes(), original)
        self.assertEqual(list(self.path.parent.iterdir()), [self.path])

    def test_failed_replace_preserves_previous_workspace_and_cleans_temp(self):
        self.path.write_text('{"cases": []}', encoding="utf-8")
        with patch("cgyro_io.os.replace", side_effect=OSError("disk unavailable")):
            with self.assertRaises(OSError):
                CgyroUiMixin._write_workspace_file(self.path, {"cases": ["new"]})
        self.assertEqual(json.loads(self.path.read_text()), {"cases": []})
        self.assertEqual(list(self.path.parent.iterdir()), [self.path])

    def test_unavailable_workspace_preserves_loaded_cases(self):
        self.path.write_text(json.dumps({"cases": [{"path": "missing"}]}))
        app = WorkspaceHarness()
        with patch.object(app, "_instantiate_cgyrodata_light", side_effect=OSError("missing")), \
                patch("cgyro_comparison_ui.cgyrodata_plot", side_effect=OSError("missing")):
            result = app._load_workspace_file(self.path, show_message=False)
        self.assertIsNone(result)
        self.assertEqual(app.cases, {"original": app.original})
        self.assertEqual(app.case_listbox.names, ["original"])

    def test_restore_preserves_case_aliases_selection_and_relative_paths(self):
        self.path.write_text(json.dumps({
            "cases": [{"name": "custom-A", "path": "one/shared"},
                      {"name": "custom-B", "path": "two/shared"}],
            "selected_cases": ["custom-B"],
        }))
        app = WorkspaceHarness()
        with patch.object(app, "_instantiate_cgyrodata_light",
                          side_effect=lambda path: SimpleNamespace(dir=path)) as load:
            result = app._load_workspace_file(self.path, show_message=False)
        self.assertEqual(result["loaded_count"], 2)
        self.assertEqual(app.case_listbox.names, ["custom-A", "custom-B"])
        self.assertEqual(app._case_listbox_selected_names(), ["custom-B"])
        expected = (self.path.parent / "one/shared").resolve()
        self.assertEqual(Path(load.call_args_list[0].args[0]).resolve(), expected)

    def test_invalid_selection_schema_does_not_clear_workspace(self):
        self.path.write_text(json.dumps({"cases": [], "selected_cases": 42}))
        app = WorkspaceHarness()
        self.assertIsNone(app._load_workspace_file(self.path, show_message=False))
        self.assertEqual(list(app.cases), ["original"])

    def test_partial_workspace_load_keeps_valid_cases_and_reports_failure(self):
        self.path.write_text(json.dumps({"cases": [
            {"name": "valid", "path": "valid"}, {"name": "bad", "path": "missing"},
        ]}))
        app = WorkspaceHarness()
        def load(path):
            if "missing" in path:
                raise OSError("missing")
            return SimpleNamespace(dir=path)
        with patch.object(app, "_instantiate_cgyrodata_light", side_effect=load), \
                patch("cgyro_comparison_ui.cgyrodata_plot", side_effect=OSError("missing")):
            result = app._load_workspace_file(self.path, show_message=False)
        self.assertEqual(list(app.cases), ["valid"])
        self.assertEqual(result["loaded_count"], 1)
        self.assertEqual(len(result["failed"]), 1)

    def test_failed_auto_restore_retains_recovery_file(self):
        app = WorkspaceHarness()
        app._auto_workspace_path = str(self.path)
        app.status_var = Mock()
        self.path.write_text("invalid JSON")
        app._restore_auto_workspace()
        self.assertTrue(self.path.exists())


class UpdateReviewTests(unittest.TestCase):
    def release(self, version):
        return json.dumps({"tag_name": "v" + version, "name": version})

    def test_stale_release_does_not_hide_newer_main_version(self):
        with patch.object(updater, "_fetch_text",
                          side_effect=[self.release("0.2.4"), "0.2.23"]):
            result = updater.check_for_updates("0.2.22", proxy_config={"mode": "direct"})
        self.assertTrue(result.update_available)
        self.assertEqual(result.latest_version, "0.2.23")
        self.assertEqual(result.source, "GitHub VERSION file")

    def test_newer_release_wins_over_older_manifest(self):
        with patch.object(updater, "_fetch_text",
                          side_effect=[self.release("0.2.24"), "0.2.23"]):
            result = updater.check_for_updates("0.2.22", proxy_config={"mode": "direct"})
        self.assertEqual(result.latest_version, "0.2.24")

    def test_valid_release_survives_manifest_network_failure(self):
        with patch.object(updater, "_fetch_text", side_effect=[
                self.release("0.2.24"), updater.UpdateCheckError("raw offline"),
                updater.UpdateCheckError("API offline")]):
            result = updater.check_for_updates("0.2.22", proxy_config={"mode": "direct"})
        self.assertEqual(result.latest_version, "0.2.24")

    def test_ssh_option_injection_rejected_before_starting_process(self):
        with patch.object(updater.subprocess, "Popen") as process:
            with self.assertRaises(updater.UpdateConnectionError):
                updater._validate_update_proxy_config({
                    "mode": "ssh-socks", "ssh_host": "-oProxyCommand=unexpected",
                })
        process.assert_not_called()

    def test_git_option_injection_rejected_before_running_git(self):
        with patch.object(updater, "_run_git") as git:
            with self.assertRaises(updater.GitUpdateError):
                updater.update_from_git(remote="--upload-pack=unexpected")
        git.assert_not_called()

    def test_wrong_branch_and_dirty_worktree_never_reach_pull(self):
        repo = str(Path(__file__).resolve().parents[1])
        for replies in ((repo, "refs/heads/feature"), (repo, "refs/heads/main", " M file.py")):
            processes = [SimpleNamespace(returncode=0, stdout=text) for text in replies]
            with self.subTest(replies=replies), patch.object(updater, "_run_git", side_effect=processes) as git:
                with self.assertRaises(updater.GitUpdateError):
                    updater.update_from_git(proxy_config={"mode": "direct"})
                self.assertFalse(any("pull" in call.args[1] for call in git.call_args_list))

    def test_non_fast_forward_is_reported_as_failure(self):
        repo = str(Path(__file__).resolve().parents[1])
        replies = [SimpleNamespace(returncode=0, stdout=text)
                   for text in (repo, "refs/heads/main", "", "old-sha")]
        replies.append(SimpleNamespace(returncode=1, stdout="Not possible to fast-forward"))
        with patch.object(updater, "_run_git", side_effect=replies) as git:
            with self.assertRaisesRegex(updater.GitUpdateError, "fast-forward"):
                updater.update_from_git(proxy_config={"mode": "direct"})
        self.assertIn("--ff-only", git.call_args_list[-1].args[1])

    def test_invalid_socks_urls_rejected(self):
        for value in ("socks5h://localhost:0", "socks5h://localhost:1080/path",
                      "socks5h://localhost:1080?x=1", "socks5h://localhost:1080#frag"):
            with self.subTest(value=value), self.assertRaises(updater.UpdateConnectionError):
                updater._parse_socks_proxy(value)

    def test_shell_metacharacters_in_proxy_hosts_are_rejected(self):
        for host in ("$(unexpected)", "host;unexpected", "host`unexpected`", "host'unexpected"):
            for parse, scheme in ((updater._parse_socks_proxy, "socks5h"),
                                  (updater._parse_http_proxy, "http")):
                with self.subTest(host=host, scheme=scheme), self.assertRaises(updater.UpdateConnectionError):
                    parse(f"{scheme}://{host}:1080")

    def test_both_update_services_failing_produces_an_error(self):
        with patch.object(updater, "_fetch_text", side_effect=updater.UpdateCheckError("offline")):
            with self.assertRaises(updater.UpdateCheckError):
                updater.check_for_updates("0.2.22", proxy_config={"mode": "direct"})

    def test_ipv6_socks_url_roundtrip(self):
        url = updater._socks_url("::1", 1080)
        self.assertEqual(updater._parse_socks_proxy(url), ("::1", 1080))

    def test_cli_help_needs_no_gui_or_scientific_dependencies(self):
        script = Path(__file__).resolve().parents[1] / "cgyro_comparison.py"
        run = subprocess.run([sys.executable, "-S", str(script), "--help"],
                             capture_output=True, text=True, timeout=10)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn("--update", run.stdout)


class BackgroundUpdateReviewTests(unittest.TestCase):
    def application(self):
        app = WorkspaceHarness()
        app._update_events = queue.Queue()
        app._update_poll_id = None
        app._update_check_in_progress = False
        app._git_update_in_progress = False
        app._update_check_menu = None
        app.status_var = Mock()
        return app

    def test_network_worker_does_not_call_tk_from_background_thread(self):
        app = self.application()
        app.root.after.side_effect = AssertionError("worker touched Tk")
        result = SimpleNamespace(update_available=False)
        with patch("cgyro_comparison_ui.check_for_updates", return_value=result):
            worker = threading.Thread(target=app._update_check_worker, args=(True,))
            worker.start()
            worker.join(timeout=2)
        self.assertFalse(worker.is_alive())
        callback, args = app._update_events.get_nowait()
        self.assertEqual(args, (True, result, None))
        self.assertEqual(callback, app._finish_update_check)
        app.root.after.assert_not_called()

    def test_git_update_returns_to_event_loop_before_download_finishes(self):
        app = self.application()
        started, finish = threading.Event(), threading.Event()
        app._save_auto_workspace = Mock(return_value="recovery.json")
        app._restart_after_update = Mock()
        def download(**kwargs):
            started.set()
            if not finish.wait(timeout=5):
                raise RuntimeError("test download was not released")
            return SimpleNamespace(output="updated")
        with patch("cgyro_comparison_ui.update_from_git", side_effect=download):
            app._pull_update_from_git("0.2.23")
            try:
                self.assertTrue(started.wait(timeout=2))
                self.assertTrue(app._git_update_in_progress)
                app._restart_after_update.assert_not_called()
            finally:
                finish.set()
            callback, args = app._update_events.get(timeout=2)
            callback(*args)
        app._restart_after_update.assert_called_once_with("recovery.json", "0.2.23")
        self.assertFalse(app._git_update_in_progress)

    def test_unexpected_git_error_releases_busy_state(self):
        app = self.application()
        app._git_update_in_progress = True
        app._remove_auto_workspace = Mock()
        app._show_manual_git_update = Mock()
        with patch("cgyro_comparison_ui.update_from_git", side_effect=RuntimeError("unexpected")):
            app._git_update_worker("repo", "recovery.json", "0.2.23")
        callback, args = app._update_events.get_nowait()
        callback(*args)
        self.assertFalse(app._git_update_in_progress)
        app._show_manual_git_update.assert_called_once()


class ExportHarness(CgyroDataExportMixin):
    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot()
        self._clear_current_plot_data()


class ExportReviewTests(unittest.TestCase):
    def test_csv_case_labels_are_text_not_spreadsheet_formulas(self):
        for label in ("=1+1", "+SUM(A1)", "-SUM(A1)", "  @SUM(A1)"):
            self.assertTrue(CgyroDataExportMixin._spreadsheet_safe_label(label).startswith("'"))
        self.assertEqual(CgyroDataExportMixin._spreadsheet_safe_label("case-A"), "case-A")

    def test_repeated_or_colliding_case_exports_never_reuse_destination(self):
        with tempfile.TemporaryDirectory() as folder:
            first = Path(CgyroDataExportMixin._new_readable_export_dir(folder, "case A"))
            marker = first / "original.txt"
            marker.write_text("previous result")
            second = Path(CgyroDataExportMixin._new_readable_export_dir(folder, "case_A"))
            third = Path(CgyroDataExportMixin._new_readable_export_dir(folder, "case A"))
            self.assertEqual(len({first, second, third}), 3)
            self.assertEqual(marker.read_text(), "previous result")

    def test_missing_precision_metadata_defaults_to_single_precision(self):
        data = SimpleNamespace()
        self.assertEqual(CgyroDataExportMixin._real_dtype_for_data(data), np.dtype("float32"))
        data.hiprec_flag = 1
        self.assertEqual(CgyroDataExportMixin._real_dtype_for_data(data), np.dtype("float64"))

    def test_binary_fallback_preserves_actual_values_without_byte_attribute(self):
        with tempfile.TemporaryDirectory() as folder:
            values = np.array([1.25, -2.5, 9.75, 0.5], dtype=np.float32)
            values.tofile(Path(folder) / "bin.cgyro.test")
            app = ExportHarness()
            result = app._read_payload_direct(SimpleNamespace(dir=folder), ".cgyro.test")
            np.testing.assert_equal(result, values)

    def test_current_case_precision_beats_reader_defaults(self):
        with tempfile.TemporaryDirectory() as folder:
            data = SimpleNamespace(dir=folder, BYTE="float64", hiprec_flag=1)
            path = Path(folder) / "out.cgyro.equilibrium"
            path.write_text("9 8\n7 0\n")
            self.assertEqual(CgyroDataExportMixin._real_dtype_for_data(data), np.dtype("float32"))

    def test_cached_masked_values_remain_missing(self):
        app = ExportHarness()
        app._record_current_plot_xyz_dataset("map", [0, 1], [0],
                                              np.ma.array([[1., 999.]], mask=[[0, 1]]))
        self.assertTrue(np.isnan(app._collect_current_plot_xyz_datasets()[0][3][1]))

    def test_masked_line_values_remain_missing(self):
        app = ExportHarness()
        app.ax.plot([0, 1], np.ma.array([1., 999.], mask=[0, 1]), label="curve")
        self.assertTrue(np.isnan(app._collect_current_plot_xyz_datasets()[0][2][1]))

    def test_cached_image_is_not_exported_twice(self):
        app = ExportHarness()
        app.ax.imshow([[1., 2.]], extent=[-1, 3, -1, 1], origin="lower")
        app._record_current_plot_xyz_dataset("map", [0., 2.], [0.], [[1., 2.]])
        self.assertEqual(len(app._collect_current_plot_xyz_datasets()), 1)

    def test_upper_origin_image_keeps_coordinate_value_pairing(self):
        app = ExportHarness()
        app.ax.imshow([[10.], [20.]], extent=[0, 2, 0, 4], origin="upper")
        _, x, y, z = app._collect_current_plot_xyz_datasets()[0]
        np.testing.assert_equal(x, [1, 1])
        np.testing.assert_equal(y, [3, 1])
        np.testing.assert_equal(z, [10, 20])

    def test_hidden_artists_and_colorbar_not_exported(self):
        app = ExportHarness()
        app.ax.plot([0, 1], [4, 5], label="hidden", visible=False)
        bar = app.fig.add_axes([.8, .1, .1, .8], label="<colorbar>")
        bar.plot([0, 1], [1, 2])
        self.assertEqual(app._collect_current_plot_xyz_datasets(), [])

    def test_export_metadata_uses_plotted_not_later_selected_options(self):
        app = ExportHarness()
        state = {"t_start_var": "500"}
        app._get_workspace_state = lambda: state
        app._capture_current_plot_context("Phi vs theta", "Phi vs theta", ["case-A"])
        state["t_start_var"] = "900"
        metadata = app._current_plot_export_metadata([])
        self.assertEqual(metadata["workspace_state"]["t_start_var"], "500")
        self.assertEqual(metadata["cases"], ["case-A"])


class NumericReviewTests(unittest.TestCase):
    def test_small_nonzero_profiles_can_be_normalized(self):
        app = Plotting()
        for mode, expected in (("Max value", [.5, 1]), ("Min value", [1, 2]),
                               ("Min-Max [0, 1]", [0, 1])):
            with self.subTest(mode=mode):
                app.fluc_normalization_var = variable(mode)
                normalized, label = app._normalize_fluctuation_profile([1e-20, 2e-20])
                np.testing.assert_allclose(normalized, expected)
                self.assertTrue(label)

    def test_complex_field_is_reused_without_reconverting_whole_array(self):
        app = Plotting()
        data = SimpleNamespace(kxky_phi=np.ones((2, 4, 3, 2, 5), dtype=np.float32))
        first = app._load_kxky_complex(data, "case", "kxky_phi", ".cgyro.kxky_phi")
        second = app._load_kxky_complex(data, "case", "kxky_phi", ".cgyro.kxky_phi")
        self.assertIs(first, second)

    def test_custom_small_value_and_zero_normalization(self):
        app = Plotting()
        app.fluc_normalization_var = variable("Custom value")
        app.fluc_normalization_custom_value_var = variable("1e-30")
        values, label = app._normalize_fluctuation_profile([2e-30, 3e-30])
        np.testing.assert_allclose(values, [2, 3])
        self.assertTrue(label)
        app.fluc_normalization_var = variable("Max absolute value")
        values, label = app._normalize_fluctuation_profile([0, 0])
        np.testing.assert_equal(values, [0, 0])
        self.assertEqual(label, "")

    def test_complex_view_preserves_pygacode_array_and_reloads_on_new_source(self):
        app = Plotting()
        raw = np.asfortranarray(np.arange(120, dtype=np.float32).reshape(2, 3, 2, 2, 5))
        data = SimpleNamespace(kxky_phi=raw)
        first = app._load_kxky_complex(data, "case", "kxky_phi", ".cgyro.kxky_phi")
        self.assertIs(data.kxky_phi, raw)
        self.assertTrue(np.shares_memory(first, raw))
        np.testing.assert_equal(first, raw[0] + 1j * raw[1])
        data.kxky_phi = raw + 2
        second = app._load_kxky_complex(data, "case", "kxky_phi", ".cgyro.kxky_phi")
        self.assertIsNot(first, second)
        np.testing.assert_equal(second, data.kxky_phi[0] + 1j * data.kxky_phi[1])

    def test_fft_aligns_payload_to_time_and_rejects_invalid_sampling(self):
        app = Plotting()
        field = np.ones((3, 2, 4))
        sampled, times, dt = app._prepare_fft_samples(field, np.arange(6), np.arange(6))
        self.assertEqual(sampled.shape[-1], 4)
        np.testing.assert_equal(times, np.arange(4))
        self.assertEqual(dt, 1)
        for time_axis in ([0, 1, 3, 4], [0, 1, 1, 2], [3, 2, 1, 0], [0, np.nan, 2, 3]):
            with self.subTest(time_axis=time_axis), self.assertRaises(ValueError):
                app._prepare_fft_samples(field, time_axis, [])

    def test_rcorr_reads_only_requested_field(self):
        app = Plotting()
        app.ax = Figure().add_subplot()
        app.others_rcorr_theta_var = variable("-1")
        app.others_rcorr_field_var = variable("Apar")
        field = np.ones((5, 3, 2, 4), dtype=complex)
        app._load_named_kxky_complex = Mock(return_value=field)
        data = SimpleNamespace(t=np.arange(4), n_radial=5, kx=np.arange(4),
                               getbigfield=Mock(), kxky_select=Mock())
        with patch("cgyro_comparison_plotting_other.sp_signal", None):
            app._plot_other_rcorr_phi(data, "case", np.arange(4), 0, 3)
        app._load_named_kxky_complex.assert_called_once_with(data, "case", "apar")
        data.getbigfield.assert_not_called()
        data.kxky_select.assert_not_called()
        self.assertTrue(app.ax.lines)


class RefreshPerformanceReviewTests(unittest.TestCase):
    def test_formula_canvas_only_redraws_when_content_changes(self):
        app = WorkspaceHarness()
        figure = Figure()
        axis = figure.add_subplot()
        canvas = SimpleNamespace(draw_idle=Mock(), get_tk_widget=Mock())
        for _ in range(3):
            app._draw_formula_panel(figure, axis, canvas, ["Title", "Same explanation"])
        self.assertEqual(canvas.draw_idle.call_count, 1)
        app._draw_formula_panel(figure, axis, canvas, ["Title", "New explanation"])
        self.assertEqual(canvas.draw_idle.call_count, 2)

    def test_input_metadata_cached_and_invalidated_without_shared_mutable_state(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "input.cgyro"
            path.write_text("BETAE = 1D-2 ! comment\nN_RADIAL 8\n")
            with patch("builtins.open", wraps=open) as read:
                first = Plotting._read_input_cgyro_scalars(folder)
                second = Plotting._read_input_cgyro_scalars(folder)
                self.assertEqual(read.call_count, 1)
            self.assertEqual(first, {"BETAE": .01, "N_RADIAL": 8.})
            first["BETAE"] = 9
            self.assertEqual(second["BETAE"], .01)
            path.write_text("BETAE = 0.025\n")
            self.assertEqual(Plotting._read_input_cgyro_scalars(folder), {"BETAE": .025})


if __name__ == "__main__":
    unittest.main()
