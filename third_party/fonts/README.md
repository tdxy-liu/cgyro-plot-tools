# Application fonts

The shared policy lives in `cgyro_fonts.py`.

| Use | Font |
| --- | --- |
| GUI Chinese, English, digits | Noto Sans CJK SC, default 10 pt |
| Figure titles, axes, legends, numbers | Arial |
| Mathematics | Arial regular/italic/bold; missing glyphs from STIX |
| Optional Chinese figure text | SimSun (宋体) |

The two unmodified Noto Sans CJK SC OpenType fonts are from the
[official Noto CJK repository](https://github.com/notofonts/noto-cjk/tree/f8d157532fbfaeda587e826d4cd5b21a49186f7c/Sans/OTF/SimplifiedChinese).
They are licensed under the [SIL Open Font License 1.1](noto-sans-cjk-sc/LICENSE).
Their original copyright records are retained in the font metadata. The source
commit, paths and SHA-256 checksums are recorded in `noto-sans-cjk-sc/manifest.json`.

No fonts are downloaded during startup. If Noto is not already installed,
Windows uses private process registration and Linux uses Fontconfig's application
font API. Registration is cached; option and language changes do not repeat it.
No administrator access or persistent OS font installation is required.
Other platforms can install Noto normally; a readable installed GUI font remains
available if application-font registration is unsupported.

Older core-X11 Tk builds cannot use Fontconfig's OpenType application fonts.
For these builds the GUI also recognizes legacy CJK families such as `song ti`
and `fangsong ti`, including lowercase family names. Bitmap fonts may round the
requested point size to the nearest available size. Language selection remains
independent of font detection; a missing font never silently forces English.

For an affected Linux installation, `cgyro_tk_runtime.py` supports an optional
**private Xft-enabled Tk 8.6** under `.runtime/tk-xft`. Both GUI entry points load
it before importing Tk. Its manifest pins the Python prefix, `_tkinter` extension
checksum and library checksum; another Python environment keeps its own Tk.
The updater and headless exports do not need this runtime. No compilation or
network access happens during normal application startup.

To prepare it explicitly, obtain and unpack the
[official Tk 8.6 source](https://github.com/tcltk/tk/tree/core-8-6-12)
and run with the Python interpreter used to start the GUI:

```sh
python cgyro_tk_runtime.py --install-from /path/to/tk-core-8-6-12
```

This requires gcc, make, X11/Xft development headers, and the matching Tcl
headers plus `lib/tclConfig.sh` in the Python prefix. It leaves shared Conda and
system libraries untouched and refuses to overwrite an existing private runtime.
The Hefei Python 3.11 / Tcl 8.6.12 installation was verified with Tk 8.6.12:
Latin, Chinese and digit glyphs all resolve to Noto Sans CJK SC at the configured GUI size.

Arial and SimSun are system fonts and are **not redistributed** here. STIX is
provided by Matplotlib. On machines without Arial the Latin fallback is
Liberation Sans, then DejaVu Sans. For identical Arial output on a Linux server,
install a licensed Arial font using that system's normal font installation.
Installing SimSun is only needed for optional Chinese plot text, not the GUI.
Scientific labels are English by default and do not change with GUI language.

Math uses [Matplotlib's custom font mapping and STIX fallback](https://matplotlib.org/stable/users/explain/text/mathtext.html#custom-fonts).
PDF/PS use Type 42 fonts; SVG retains text. Changing the GUI font does not change
the figure font size or any plotted values.
