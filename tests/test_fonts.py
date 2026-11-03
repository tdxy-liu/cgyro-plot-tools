"""Check real glyph selection and the shared figure typography policy."""

import unittest
from io import StringIO
from unittest.mock import patch

import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.mathtext import MathTextParser
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from cgyro_fonts import (CHINESE_PLOT_FONT_FAMILY, configure_plot_fonts,
                         gui_font_family, plot_font_family, plot_font_style)


class GuiFontTests(unittest.TestCase):
    def test_legacy_and_lowercase_tk_cjk_names_are_recognized(self):
        # Hefei's core-X11 Tk lists these legacy fonts, not Fontconfig families.
        for families, expected in (
                (("helvetica", "fangsong ti", "song ti"), "song ti"),
                (("helvetica", "fangsong ti"), "fangsong ti"),
                (("helvetica", "wenquanyi zen hei"), "wenquanyi zen hei"),
                (("song ti", "noto sans cjk sc"), "noto sans cjk sc")):
            with self.subTest(families=families), \
                    patch("tkinter.font.families", return_value=families), \
                    patch("cgyro_fonts._register_gui_fonts"), \
                    patch("cgyro_fonts._register_matplotlib_gui_fonts"):
                self.assertEqual(gui_font_family(object()), (expected, True))

    def test_bundled_noto_becomes_preferred_after_registration(self):
        with patch("tkinter.font.families", side_effect=[("song ti",), ("song ti", "Noto Sans CJK SC")]), \
                patch("cgyro_fonts._register_gui_fonts"), \
                patch("cgyro_fonts._register_matplotlib_gui_fonts"):
            self.assertEqual(gui_font_family(object()), ("Noto Sans CJK SC", True))

    def test_missing_cjk_family_keeps_available_default_font(self):
        with patch("tkinter.font.families", return_value=("helvetica",)), \
                patch("tkinter.font.Font") as default, \
                patch("cgyro_fonts._register_gui_fonts"), \
                patch("cgyro_fonts._register_matplotlib_gui_fonts"):
            default.return_value.actual.return_value = "helvetica"
            self.assertEqual(gui_font_family(object()), ("helvetica", False))


class FigureFontTests(unittest.TestCase):
    def test_math_letters_use_plot_font_and_missing_symbols_use_stix(self):
        with mpl.rc_context(plot_font_style()):
            parsed = MathTextParser("path").parse(r"$k_x\rho_s=0.2\quad\mathbb{R}$")
        # Newer Matplotlib adds glyph metadata after the font/size/position.
        families = {glyph[0].family_name for glyph in parsed.glyphs}
        self.assertIn(plot_font_family(), families)
        self.assertTrue(any(name.startswith("STIX") for name in families), families)
        self.assertNotIn("DejaVu Sans", families - {plot_font_family()})

    def test_titles_ticks_and_legend_share_font_in_export(self):
        with mpl.rc_context():
            configure_plot_fonts()
            fig = Figure(figsize=(4, 3))
            ax = fig.subplots()
            ax.plot([0, 1], [0, 1], label="Case 1")
            ax.set(title="Transfer", xlabel=r"$k_x\rho_s$", ylabel=r"$T^\Phi$")
            legend = ax.legend()
            FigureCanvasAgg(fig).draw()
            texts = [ax.title, ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(),
                     *ax.get_yticklabels(), *legend.get_texts()]
            for text in texts:
                self.assertEqual(text.get_fontfamily()[0], plot_font_family())
                self.assertEqual(text.get_math_fontfamily(), "custom")
            svg = StringIO()
            fig.savefig(svg, format="svg")
            self.assertIn(plot_font_family(), svg.getvalue())
        self.assertEqual(CHINESE_PLOT_FONT_FAMILY, "SimSun")

    def test_arial_wins_when_installed(self):
        if any(font.name == "Arial" for font in font_manager.fontManager.ttflist):
            self.assertEqual(plot_font_family(), "Arial")


if __name__ == "__main__":
    unittest.main()
