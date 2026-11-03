"""Standalone CGYRO Transfer Explorer: GUI by default, batch export with --output."""

import argparse
from pathlib import Path
import sys


def parser():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--version", action="version", version=(Path(__file__).resolve().parent / "VERSION").read_text().strip())
    cli.add_argument("--gui", action="store_true", help="Open the independent GUI.")
    cli.add_argument("--case", action="append", default=[], help="Case folder; repeat for batch export.")
    cli.add_argument("--kx", type=float, default=0.)
    cli.add_argument("--ky", type=float, default=.2, help="Fixed physical ky magnitude; matched once, then held for every step.")
    cli.add_argument("--steps", type=int, default=4, help="Number of panels, including the initial mode (1-20).")
    cli.add_argument("--direction", choices=("positive", "negative"), default="positive")
    cli.add_argument("--diagnostic", choices=("asym", "fullt"), default="asym")
    cli.add_argument("--start", type=float, help="Blank: midpoint between first saved time and chosen end.")
    cli.add_argument("--end", type=float, help="Blank: last saved time.")
    cli.add_argument("--cycles", choices=("stop", "allow", "skip"), default="stop")
    cli.add_argument("--output", type=Path, help="Export parent folder. With --case, selects headless batch mode.")
    cli.add_argument("--formats", nargs="+", choices=("pdf", "svg", "png", "tiff"), default=["pdf", "svg", "png"])
    cli.add_argument("--width-mm", type=float, default=183.)
    cli.add_argument("--panel-height-mm", type=float, default=32.)
    cli.add_argument("--font-pt", type=float, default=9.)
    cli.add_argument("--dpi", type=int, default=600)
    cli.add_argument("--color-max", type=float, help="Fixed shared symmetric bound; blank uses the unmultiplied global max |T|.")
    cli.add_argument("--normalize", action="store_true", help="Divide ALL panels by one global max |T|.")
    cli.add_argument("--glyph", choices=("smooth", "diamonds", "cells"), default="smooth",
                     help="Smooth: classic FULLT red-white-blue lobes (default).")
    cli.add_argument("--panel-scaling", choices=("auto", "none", "custom"), default="auto",
                     help="Auto: amplify weak panels by labelled 1/2/5 decades. None: preserve relative amplitudes.")
    cli.add_argument("--panel-multipliers", nargs="+", type=float,
                     help="Custom display factors in panel order, e.g. 1 1 10; selects custom scaling.")
    cli.add_argument("--xlim", nargs=2, type=float)
    cli.add_argument("--ylim", nargs=2, type=float)
    cli.add_argument("--title", help="Custom context line; pass an empty string to hide it.")
    cli.add_argument("--colorbar-label", default="", help="Optional label; blank by default.")
    return cli


def main(argv=None):
    cli = parser()
    args = cli.parse_args(argv)
    try:
        from cgyro_fullt_trace_core import TraceSettings, trace_case
        from cgyro_fullt_trace_plot import FigureSettings, export_bundle, suggested_bundle_name
    except ImportError:
        from .cgyro_fullt_trace_core import TraceSettings, trace_case
        from .cgyro_fullt_trace_plot import FigureSettings, export_bundle, suggested_bundle_name
    trace = TraceSettings(args.kx, args.ky, args.steps, args.direction, args.diagnostic, args.start, args.end, args.cycles)
    figure = FigureSettings(width_mm=args.width_mm, panel_height_mm=args.panel_height_mm, font_pt=args.font_pt,
                            color_max=args.color_max, normalize=args.normalize, glyph=args.glyph,
                            xlim=tuple(args.xlim) if args.xlim else None, ylim=tuple(args.ylim) if args.ylim else None,
                            title=args.title, colorbar_label=args.colorbar_label, dpi=args.dpi,
                            panel_scaling="custom" if args.panel_multipliers is not None else args.panel_scaling,
                            panel_multipliers=tuple(args.panel_multipliers) if args.panel_multipliers is not None else None)
    try:
        trace.validate()
        figure.validate()
        if args.gui or args.output is None:
            try:
                from cgyro_fullt_trace_ui import launch
            except ImportError:
                from .cgyro_fullt_trace_ui import launch
            launch(args.case, trace, figure, args.output, args.formats)
            return 0
        if not args.case:
            cli.error("--output requires at least one --case (or add --gui).")
        failures = 0
        for case in args.case:
            try:
                result = trace_case(case, trace, progress=lambda text: print(text, flush=True))
                output = export_bundle(result, figure, args.output / suggested_bundle_name(result), args.formats)
                print(f"Exported {len(result.panels)} panels: {output}\n{result.stop_reason}", flush=True)
            except (ValueError, OSError, RuntimeError) as exc:
                failures += 1
                print(f"Failed: {case}: {exc}", file=sys.stderr)
        return 1 if failures else 0
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"FULLT Trace: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
