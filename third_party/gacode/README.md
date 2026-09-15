# Bundled FULLT reader provenance

`../../cgyro_fullt_reader.py` is the FTZ v1 reader/writer from the user's
customized GACODE tree, `f2py/pygacode/cgyro/fullt_compressed.py`, snapshotted
on 2026-09-15. That tree is based on GACODE commit
`6f0cb880e09e4d48aa72a6871ed9fb96b9f7ed4a`; the custom compression module is
an addition, **not** a file present in that upstream commit.

Snapshot SHA-256 (original file bytes):
`d55908d58075a39d6424f1f91124f15d7b44c35189814397800783dbbb95670e`.

It is shipped here under Apache-2.0, with the source tree's LICENSE and NOTICE
retained in this directory. This does not imply endorsement by General Atomics.

Local adaptations: module rename, provenance comments, and actionable optional
dependency errors (missing zstandard is not reported as corrupt data). No file
format, numeric algorithm, pairing, precision, sign or normalization change.
The Python encoder is retained for bitwise regression tests; the GUI opens
existing data read-only unless the user explicitly requests a separate export.

All Comparison FULLT paths use this pinned module, **not** an older codec found
on PYTHONPATH. Ordinary fields/fluxes still use the user's pygacode. Updating
this tool neither replaces that installation nor changes the CGYRO executable.
Future codec synchronization must retain the cross-reader compatibility tests.
The Zstandard Python extension is installed separately per platform; no binary
extension, simulation data or machine-specific paths are distributed here.
