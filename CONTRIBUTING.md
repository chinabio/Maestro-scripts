# Contributing

Thanks for your interest in improving Maestro Exporters!

## Ground rules

These are **Maestro user scripts** — they run inside Schrödinger's bundled
Python interpreter and depend on `schrodinger.structure`,
`schrodinger.maestro`, and (for the 2D converter) Schrödinger's RDKit adapter.
You generally cannot run them outside Maestro, so test changes by installing the
script via **Scripts → Manage Scripts… → Install…** and exercising the dialog.

## Style & conventions

* Keep the three scripts consistent with each other — they intentionally share
  the same shape: `DEFAULTS`, `load_settings` / `save_settings`,
  `get_selected_structures`, `is_protein`, `group_by_protein`, writer helpers,
  launch/reveal helpers, a Qt `ExportDialog` + `SuccessDialog`, and a `panel()`
  entry point.
* Support **both PyQt6 and PyQt5** (the `QT = 6/5` probe at the top).
* Fail soft: a single bad structure should never abort an export. Print a
  `[ScriptName] ...` diagnostic and continue.
* When launching an external GUI, scrub Schrödinger paths from `PATH` and clear
  the Qt/Python environment variables so the target's own runtime is used.

## Screenshots

If your change alters the UI, please refresh the relevant image in
[`screenshots/`](screenshots/) (keep the existing filenames so README links keep
working) and confirm the README renders correctly.

## Pull requests

1. Fork and create a feature branch.
2. Describe what you changed and how you tested it inside Maestro.
3. Update `README.md` / `docs/USAGE.md` and `CHANGELOG.md` if behaviour changes.

## Reporting bugs

Open an issue with your Schrödinger version, OS, the target tool + version, and
the console output printed in the Maestro terminal (lines prefixed with
`[Mae2Flare]`, `[Mae2Moe]`, `[Mae2SD2D]`, or `[pyflare]`).
