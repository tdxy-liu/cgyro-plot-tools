"""
CGYRO comparison plotting GUI.

This Tkinter tool loads one or more CGYRO case directories and provides
interactive comparison plots for frequencies, growth rates, flux channels,
fluctuations, FFT views, and zonal ExB diagnostics.
"""

import tkinter as tk
import sys

try:
    from cgyro_comparison_ui import CgyroUiMixin
    from cgyro_comparison_plotting import Plotting
    from cgyro_data_export import CgyroDataExportMixin
except ImportError:
    # Support package-style import, e.g. `from CGYRO.cgyro_comparison import ...`.
    from .cgyro_comparison_ui import CgyroUiMixin
    from .cgyro_comparison_plotting import Plotting
    from .cgyro_data_export import CgyroDataExportMixin


class CGYRO_Comparison(CgyroUiMixin, CgyroDataExportMixin, Plotting):
    """
    Main GUI/controller class for interactive CGYRO case comparison.

    The tool is split into mixins so UI state, binary/data export, and plotting
    backends can evolve independently while sharing one Tk controller instance.
    """


if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            from cgyro_update import _run_cli
        except ImportError:
            from .cgyro_update import _run_cli
        raise SystemExit(_run_cli(sys.argv[1:]))

    root = tk.Tk()
    app = CGYRO_Comparison(root)
    root.mainloop()
