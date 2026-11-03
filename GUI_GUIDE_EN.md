# Interface and language guide · v0.2.32

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

Figure titles, axes, legends and numbers use **Arial**. Math uses the matching regular/italic/bold faces, with **STIX** for unavailable symbols. **SimSun (宋体)** is reserved for Chinese figure text; scientific labels remain English by default. Preview and export share this policy. Arial and SimSun are read from the local system. Without Arial, Linux falls back to Liberation Sans, then DejaVu Sans; exact Arial reproduction requires a licensed local Arial installation. See [fonts and licenses](third_party/fonts/README.md).

## Presets and sessions

Choose preset applies a saved parameter set; its fields remain editable. Save preset stores the current settings. Saving the same name replaces that preset.

- Comparison presets contain general plotting parameters, excluding case directories, selection locks, per-case kx/ky, species selection and the linear-data file path.
- Transfer Explorer presets contain trace and figure settings, excluding case directories, current results and export directories.
- Workspaces and sessions continue to store case lists and session settings. Successful saves are reported in the status bar.

Preset files are `~/.cgyro_comparison_tool/comparison_presets.json` and `fullt_trace_presets.json`. Invalid files produce an error and are not silently overwritten. Isolated configurations can use `CGYRO_PRESET_DIR` and `CGYRO_UI_CONFIG`; `CGYRO_UI_LANGUAGE=en` overrides the startup language for the current process.

## Verification

Layout checks cover 920×640, 1280×900 and 1600×1000 windows at multiple Tk scaling factors. They check complete-field reflow, control boundaries, retained inputs, independent scrolling, language round trips and preset persistence. Linux registers bundled Noto fonts using Fontconfig; installed-font fallback remains available if registration fails.
