"""Shared offline font policy for Tk controls and scientific figures."""

import ctypes
from ctypes.util import find_library
from functools import lru_cache
import logging
import os
from pathlib import Path
import sys

import matplotlib as mpl
from matplotlib import font_manager

GUI_FONT_FAMILY = "Noto Sans CJK SC"
GUI_FONT_SIZE = 11
PLOT_FONT_FAMILY = "Arial"
CHINESE_PLOT_FONT_FAMILY = "SimSun"  # 宋体; reserved for Chinese figure text.
GUI_FONT_FALLBACKS = (GUI_FONT_FAMILY, "Noto Sans SC", "Source Han Sans SC",
                      "Microsoft YaHei UI", "Microsoft YaHei", "微软雅黑",
                      "WenQuanYi Micro Hei", "WenQuanYi Zen Hei", "SimHei", "SimSun")
PLOT_FONT_FALLBACKS = (PLOT_FONT_FAMILY, "Liberation Sans", "DejaVu Sans")
_FONT_DIR = Path(__file__).resolve().parent / "third_party" / "fonts" / "noto-sans-cjk-sc"
_GUI_FONT_FILES = tuple(_FONT_DIR / ("NotoSansCJKsc-" + weight + ".otf")
                        for weight in ("Regular", "Bold"))
_LOG = logging.getLogger(__name__)


@lru_cache(maxsize=1)
def _register_gui_fonts():
    """Register bundled fonts in this process, without installing/downloading."""
    paths = [path for path in _GUI_FONT_FILES if path.is_file()]
    if not paths:
        return
    try:
        if sys.platform == "win32":
            add_font = ctypes.WinDLL("gdi32", use_last_error=True).AddFontResourceExW
            add_font.argtypes = [ctypes.c_wchar_p, ctypes.c_uint32, ctypes.c_void_p]
            add_font.restype = ctypes.c_int
            for path in paths:
                if not add_font(str(path), 0x10, None):  # FR_PRIVATE; removed on exit.
                    raise OSError("Cannot load bundled GUI font: " + path.name)
        elif sys.platform.startswith("linux"):
            library = find_library("fontconfig")
            if not library:
                raise OSError("Fontconfig is unavailable")
            config = ctypes.CDLL(library)
            config.FcConfigGetCurrent.argtypes = []
            config.FcConfigGetCurrent.restype = ctypes.c_void_p
            config.FcConfigAppFontAddFile.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
            config.FcConfigAppFontAddFile.restype = ctypes.c_int
            config.FcConfigBuildFonts.argtypes = [ctypes.c_void_p]
            config.FcConfigBuildFonts.restype = ctypes.c_int
            current = config.FcConfigGetCurrent()
            if not current:
                raise OSError("Cannot get the Fontconfig configuration")
            for path in paths:
                if not config.FcConfigAppFontAddFile(current, os.fsencode(path)):
                    raise OSError("Cannot load bundled GUI font: " + path.name)
            if not config.FcConfigBuildFonts(current):
                raise OSError("Cannot refresh application fonts")
    except (OSError, AttributeError) as exc:
        _LOG.warning("Bundled Noto GUI font unavailable; using an installed fallback: %s", exc)


@lru_cache(maxsize=1)
def _register_matplotlib_gui_fonts():
    # The GUI's empty-canvas prompt is interface text, not a scientific label.
    if any(font.name == GUI_FONT_FAMILY for font in font_manager.fontManager.ttflist):
        return
    for path in _GUI_FONT_FILES:
        if path.is_file():
            font_manager.fontManager.addfont(str(path))


def gui_font_family(root):
    from tkinter import font as tkfont
    families = set(tkfont.families(root))
    if GUI_FONT_FAMILY not in families:
        _register_gui_fonts()
        families = set(tkfont.families(root))
    _register_matplotlib_gui_fonts()
    family = next((name for name in GUI_FONT_FALLBACKS if name in families), None)
    return (family or tkfont.nametofont("TkDefaultFont", root=root).actual("family"),
            family is not None)


def interface_plot_font_family():
    """Font for bilingual UI instructions drawn inside a Matplotlib canvas."""
    _register_matplotlib_gui_fonts()
    installed = {font.name for font in font_manager.fontManager.ttflist}
    return next((name for name in GUI_FONT_FALLBACKS if name in installed), "DejaVu Sans")


@lru_cache(maxsize=1)
def plot_font_family():
    installed = {font.name for font in font_manager.fontManager.ttflist}
    return next((name for name in PLOT_FONT_FALLBACKS if name in installed), "DejaVu Sans")


def plot_font_families():
    families = [plot_font_family()]
    if any(font.name == CHINESE_PLOT_FONT_FAMILY for font in font_manager.fontManager.ttflist):
        families.append(CHINESE_PLOT_FONT_FAMILY)
    return families + ["STIXGeneral"]


def plot_font_style():
    """Arial Latin/Greek glyphs first; STIX supplies missing math symbols."""
    family = plot_font_family()
    style = {
        "font.family": plot_font_families(),
        "font.sans-serif": plot_font_families(),
        "mathtext.fontset": "custom", "mathtext.default": "it",
        "mathtext.rm": family, "mathtext.it": family + ":italic",
        "mathtext.bf": family + ":bold", "mathtext.bfit": family + ":bold:italic",
        "mathtext.sf": family, "mathtext.tt": family,
        "mathtext.cal": "STIXGeneral:italic", "mathtext.fallback": "stix",
        "axes.unicode_minus": True, "text.usetex": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    }
    # Older supported Matplotlib versions do not expose mathtext.bfit.
    return {key: value for key, value in style.items() if key in mpl.rcParams}


def configure_plot_fonts():
    # Tk preview redraws happen after the figure-building rc_context exits.
    # Keep the same math font mapping active for later canvas draws/exports.
    mpl.rcParams.update(plot_font_style())
