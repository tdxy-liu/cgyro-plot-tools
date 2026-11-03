# Interface and language guide

Start the comparison app with `python cgyro_comparison.py`.

Start **CGYRO Transfer Explorer** with `python cgyro_fullt_trace.py`, or use **Plot → CGYRO Transfer Explorer...** in the comparison app. Existing command-line options and saved JSON sessions remain compatible.

## Workflow

1. **Select cases.** The case list stays visible. Add and Remove are direct actions; the comparison app groups multiple-case loading, reloading, case information and input comparison under More. No selection uses all cases. Lock selection preserves the comparison selection when editing options. Long case names have horizontal scrolling.
2. **Configure the diagnostic or trace.** Choose the diagnostic first, then its quantities and spectral selection. Transfer Explorer has a Trace settings tab with initial kx, fixed ky, panel count and extremum direction.
3. **Set time and display options.** The comparison app's Time and display section follows diagnostic settings; scroll down in a shorter window. Automatic time modes disable manual start/end fields. Switching back to Manual retains the earlier values. Transfer Explorer has a Figure style tab for dimensions, rendering, multipliers and the shared color scale. Titles, resolution and repeat handling can be expanded when needed.
4. **Run.** The blue Plot or Trace cases button stays visible. Apply figure style reuses cached trace results. Changes to the initial mode or averaging interval require tracing again.
5. **Inspect and export.** Use the comparison result area's data export, animation/page controls and Matplotlib image toolbar. Transfer Explorer separates Stacked figure from Trace coordinates and values. Its export controls sit below the results; expand Export folder and formats before exporting figures and raw trace data.

Definitions and calculation notes are collapsed by default and typeset only when opened. Editing a diagnostic does not reload all cases.

## Language

Use **Settings → Language → English / 中文** in either app. Changes apply immediately, preserving parameters, case selection and existing figures. The preference is saved in `~/.cgyro_comparison_tool/ui.json` and used the next time either app starts. Another already-open window can be switched independently.

Scientific symbols, exported figures/data, case names, paths and underlying technical error messages are not translated. GUI Chinese, English and digits use **Noto Sans CJK SC, 11 pt** by default; secondary text is 10 pt, section headings 11 pt bold and the app heading 13 pt bold. Regular and bold Noto fonts are bundled and loaded into the current Windows/Linux process when needed, without downloads or system font changes. Tk/system DPI controls scaling.

Some servers use core-X11 Tk builds that cannot load the bundled OpenType fonts. The app falls back to CJK families visible to Tk, such as `song ti` or `fangsong ti`; available bitmap sizes determine the actual size. Font detection never overrides the selected language. If characters remain missing, install Noto Sans CJK SC and use a Fontconfig/Xft-enabled Tk build.

Figure titles, axes, legends and numbers use **Arial**. Math uses the matching regular/italic/bold faces, with **STIX** for unavailable symbols. **SimSun (宋体)** is reserved for Chinese figure text; scientific labels remain English by default. Preview and export share this policy. Arial and SimSun are read from the local system. Without Arial, Linux falls back to Liberation Sans, then DejaVu Sans; exact Arial reproduction requires a licensed local Arial installation. See [fonts and licenses](third_party/fonts/README.md).

## Workspaces and sessions

Use **File → Save Workspace... / Open Workspace...** in the comparison app to save and restore case lists and plotting settings. Transfer Explorer saves and opens sessions from its **File** menu. Successful saves are reported in the status bar.

`CGYRO_UI_CONFIG` selects an independent language preference file; `CGYRO_UI_LANGUAGE=en` overrides the startup language for the current process.

## Verification

Frequency and Growth Rate load only `bin.cgyro.freq` or `out.cgyro.freq` when
first requested, then average the selected time samples. Missing or malformed
diagnostics are reported with the case name; valid cases still plot. If the case
is still running, use Reload after its output files finish writing.

On Linux with an old core-X11 Tk, an optional private Xft Tk enables the bundled
Noto font without changing the system Python. See [font setup](third_party/fonts/README.md).

Layout checks cover 920×640, 1280×900 and 1600×1000 windows at multiple Tk scaling factors. They check complete-field reflow, control boundaries, retained inputs, independent scrolling and language round trips. Font regression checks cover legacy X11 names; missing font detection does not override the selected language.
