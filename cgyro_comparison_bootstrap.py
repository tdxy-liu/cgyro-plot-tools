"""
CGYRO comparison bootstrap helpers.

Centralizes pygacode import bootstrap, fallback mocks, and shared constants.
"""

import os
import sys
import getpass
import json
import sysconfig
import numpy as np

# Ensure pygacode can be imported.
_script_dir = os.path.dirname(__file__)
# Machine-local settings survive ordinary launches and remain outside Git.
# An explicit GACODE_ROOT still has priority, including over a stale local file.
_runtime_config_path = os.path.join(_script_dir, "cgyro_runtime.local.json")
_selected_gacode_root = os.environ.get("GACODE_ROOT")
_root_selection_source = "GACODE_ROOT" if _selected_gacode_root else None
if not _selected_gacode_root and os.path.isfile(_runtime_config_path):
    try:
        with open(_runtime_config_path, encoding="utf-8-sig") as handle:
            _runtime_config = json.load(handle)
    except (OSError, ValueError) as exc:
        raise RuntimeError("Cannot read CGYRO local runtime configuration: " + _runtime_config_path) from exc
    if not isinstance(_runtime_config, dict):
        raise RuntimeError("CGYRO local runtime configuration must be a JSON object")
    _configured_root = _runtime_config.get("gacode_root")
    if not isinstance(_configured_root, str) or not _configured_root.strip():
        raise RuntimeError("CGYRO local runtime configuration requires a nonempty gacode_root")
    _selected_gacode_root = os.path.expanduser(_configured_root.strip())
    if not os.path.isabs(_selected_gacode_root):
        _selected_gacode_root = os.path.join(_script_dir, _selected_gacode_root)
    _root_selection_source = "cgyro_runtime.local.json"

# Only the current interpreter/platform's private dependencies may be imported.
# In particular, never load a cp312 Windows decoder into a cp39/Linux process.
_runtime_abi = sys.implementation.cache_tag + "-" + sysconfig.get_platform()
_private_dependencies = os.path.join(_script_dir, ".runtime", "python", _runtime_abi)
if os.path.isfile(os.path.join(_private_dependencies, "zstandard", "__init__.py")):
    _private_dependencies = os.path.abspath(_private_dependencies)
    sys.path[:] = [p for p in sys.path if os.path.normcase(os.path.abspath(p)) !=
                  os.path.normcase(_private_dependencies)]
    sys.path.insert(0, _private_dependencies)

_pygacode_candidates = [
    os.path.abspath(os.path.join(_script_dir, rel_path))
    for rel_path in (
        '../f2py',
        '../gacode/f2py',
        '../gacode-master/f2py',
    )
]
if _selected_gacode_root:
    # A selected root must win even when it was already in PYTHONPATH.
    # Do not quietly select an older sibling repository if it is invalid.
    _pygacode_candidates = [os.path.abspath(os.path.join(_selected_gacode_root, "f2py"))]
    if not os.path.isfile(os.path.join(_pygacode_candidates[0], "pygacode", "cgyro", "data.py")):
        raise RuntimeError(_root_selection_source + " does not contain f2py/pygacode/cgyro/data.py")
for candidate in reversed(_pygacode_candidates):
    if os.path.isdir(candidate):
        sys.path[:] = [p for p in sys.path if os.path.normcase(os.path.abspath(p)) !=
                      os.path.normcase(os.path.abspath(candidate))]
        sys.path.insert(0, candidate)

if _selected_gacode_root and "pygacode" in sys.modules:
    loaded = getattr(sys.modules["pygacode"], "__file__", "") or ""
    expected = os.path.join(_pygacode_candidates[0], "pygacode", "__init__.py")
    if os.path.normcase(os.path.abspath(loaded)) != os.path.normcase(expected):
        raise RuntimeError("An older pygacode is already loaded; restart the GUI with the selected GACODE root")

try:
    from pygacode.cgyro.data import cgyrodata
    from pygacode.cgyro.data_plot import cgyrodata_plot
except ImportError as exc:
    if os.environ.get("CGYRO_ALLOW_MOCK_DATA") != "1":
        raise RuntimeError(
            "Cannot import pygacode. Set GACODE_ROOT to the intended repository "
            "or install pygacode. Mock data are disabled for real analysis."
        ) from exc
    print("Error: Could not import cgyrodata. Please ensure pygacode is available.")

    class cgyrodata:
        """Lightweight mock used only when pygacode imports fail."""

        def __init__(self, path):
            """Populate minimal attributes used by the GUI during local debugging."""
            self.path = path
            self.ky = np.linspace(0, 1, 5)
            self.freq = np.zeros((2, 5, 1))
            self.t = np.linspace(0, 10, 10)
            self.ky_flux = np.zeros((1, 2, 2, 5, 10))
            self.z = np.array([1.0])
            self.mass = np.array([2.0])
            self.n_radial = 1
            self.theta_plot = 1

        def getflux(self):
            """Mock no-op for flux loading."""
            return None

        def getbigfield(self):
            """Mock no-op for big-field loading."""
            return None

    class cgyrodata_plot:
        """Fallback mock placeholder for plot-wrapper loader."""

        pass


DEFAULT_APP_TITLE = "CGYRO Comparison Tool"
DEFAULT_WINDOW_GEOMETRY = "1200x800"


def _detect_current_user():
    """Return a best-effort username for shared-directory defaults."""
    user_name = os.environ.get("USER") or os.environ.get("USERNAME")
    if user_name:
        return str(user_name)
    try:
        user_name = getpass.getuser()
    except Exception:
        user_name = ""
    return str(user_name) if user_name else "unknown"


_current_user = _detect_current_user()
DEFAULT_CASE_PICKER_ROOT = f"/data/share/{_current_user}"
DEFAULT_LINEAR_GAMMA_FILE = "omega_gamma_vs_ky.txt"
DEFAULT_EXPORT_DIRNAME = "CGYRO_vs_CGYRO_exports"


def default_share_dir(*, fallback_to_cwd=True):
    """Return `/data/share/$USER` when available, then shared/local fallbacks."""
    # Most production CGYRO cases and exported comparison tables live on the
    # cluster shared filesystem.  Using this helper everywhere keeps file
    # dialogs from silently falling back to the local launch directory.
    candidates = []
    user_name = os.environ.get("USER") or os.environ.get("USERNAME") or _current_user
    if user_name:
        candidates.append(os.path.join("/data/share", str(user_name)))
    candidates.append("/data/share")
    if fallback_to_cwd:
        candidates.append(os.getcwd())

    for path in candidates:
        try:
            if path and os.path.isdir(path):
                return path
        except Exception:
            pass
    return os.getcwd() if fallback_to_cwd else ""
