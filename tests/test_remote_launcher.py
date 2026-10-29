"""Exercise launcher parsing locally without SSH, rsync, Slurm, or a cluster."""

import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "remote_cgyro"
BASH = (r"C:\Program Files\Git\bin\bash.exe" if os.name == "nt" else shutil.which("bash"))


def function_source(name):
    text = SCRIPT.read_text(encoding="utf-8")
    start = text.index(name + "() {\n")
    end = text.index("\n}\n", start) + 3
    return text[start:end] + "\n"


@unittest.skipUnless(BASH and Path(BASH).is_file(), "Bash required for isolated launcher tests")
class RemoteLauncherReviewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.folder = Path(self.tmp.name)
        self.info = self.folder / "run.env"
        self.functions = function_source("decode_run_info_value") + function_source("load_run_info")
        self.setup = 'REGION=TY; REMOTE_RUN_DIR_OVERRIDE=; RUN_INFO_FILE=$1; LOCAL_DIR=/old\n'

    def run_bash(self, code, *args):
        return subprocess.run([BASH, "--noprofile", "--norc", "-c", code,
                               "launcher-test", self.info.as_posix(), *args],
                              capture_output=True, encoding="utf-8", timeout=10)

    def test_launcher_shell_syntax(self):
        run = subprocess.run([BASH, "-n", SCRIPT.as_posix()], capture_output=True,
                             encoding="utf-8", timeout=10)
        self.assertEqual(run.returncode, 0, run.stderr)

    def test_run_info_cannot_execute_commands(self):
        marker = self.folder / "must-not-exist"
        self.info.write_text("SAVED_REMOTE_RUN_DIR=/safe/path\n"
                             + "printf exploited > " + shlex.quote(marker.as_posix()) + "\n")
        run = self.run_bash(self.functions + self.setup + "load_run_info")
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse(marker.exists())

    def test_command_substitution_in_value_is_rejected(self):
        marker = self.folder / "must-not-exist"
        self.info.write_text("SAVED_REMOTE_RUN_DIR=$(touch " + shlex.quote(marker.as_posix()) + ")\n")
        run = self.run_bash(self.functions + self.setup + "load_run_info")
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse(marker.exists())

    def test_printf_quoted_legacy_paths_remain_supported(self):
        path = "/remote/case with spaces;literal'quote"
        run = self.run_bash(
            self.functions + self.setup
            + "printf 'SAVED_REGION=%q\\nSAVED_REMOTE_RUN_DIR=%q\\n' TY \"$2\" > \"$RUN_INFO_FILE\"\n"
            + 'load_run_info && printf "%s" "$REMOTE_RUN_DIR"', path,
        )
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertEqual(run.stdout, path)

    def test_region_mismatch_is_rejected(self):
        self.info.write_text("SAVED_REGION=HF\nSAVED_REMOTE_RUN_DIR=/safe/path\n")
        run = self.run_bash(self.functions + self.setup + "load_run_info")
        self.assertNotEqual(run.returncode, 0)


if __name__ == "__main__":
    unittest.main()
