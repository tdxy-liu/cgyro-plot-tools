"""
CGYRO comparison bootstrap helpers.

Centralizes pygacode import bootstrap, fallback mocks, and shared constants.
"""

import os
import sys
import getpass
import numpy as np

# Ensure pygacode can be imported.
_script_dir = os.path.dirname(__file__)
for rel_path in ('../gacode/f2py', '../gacode-master/f2py'):
    candidate = os.path.abspath(os.path.join(_script_dir, rel_path))
    if os.path.isdir(candidate) and candidate not in sys.path:
        sys.path.append(candidate)

try:
    from pygacode.cgyro.data import cgyrodata
    from pygacode.cgyro.data_plot import cgyrodata_plot
except ImportError:
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
