# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- `screenshots/` folder with a placeholder guide and a **Screenshots** section
  in the README wired to the expected filenames.

## [1.0.0] - 2026-06-02

### Added
- `Mae2Flare.py` — export selected Project Table entries to Cresset Flare as
  Maestro + SDF (recommended) or a `.flrp` project via `pyflare`, with automatic
  fallback when `pyflare` hits the SQLite project-cache bug.
- `Mae2Moe.py` — export selected entries to CCG MOE as PDB + SDF (recommended)
  or a `.moe` session built by an auto-generated SVL startup script.
- `Mae2SD2D.py` — convert selected ligands/poses to a single clean 2D V2000 SD
  file, with optional aggressive/balanced charge neutralization, hydrogen
  stripping, SMILES tagging, and fault-tolerant per-structure handling.
- Shared niceties across all scripts: protein-vs-ligand grouping, persisted
  `*.settings.json`, PyQt6/PyQt5 support, and Reveal/Open-folder shortcuts.
