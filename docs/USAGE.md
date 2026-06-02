# Usage Guide

This guide walks through a typical workflow shared by all three scripts.

## 1. Prepare your selection

In the **Project Table**, select the entries you want to export. The expected
shape is a **receptor followed by its docking poses**, optionally repeated for
several targets:

```text
[x] Receptor_A           (protein)
[x] Receptor_A_pose_1    (ligand)
[x] Receptor_A_pose_2    (ligand)
[x] Receptor_B           (protein)
[x] Receptor_B_pose_1    (ligand)
```

The scripts walk the selection **in Project Table order**. Each protein-like
entry starts a new group; following non-protein entries become that group's
ligands. If you select nothing, the **included** Workspace entries are used.

### How "protein" is decided

An entry is treated as a protein if **either**:

- it contains at least one standard amino-acid residue
  (`ALA, ARG, ASN, ... VAL`, plus protonation variants like `HID/HIE/HIP`,
  `CYX`, `ASH`, `GLH`, `LYN`), **or**
- its atom count exceeds the **protein atom threshold** (default `500`,
  adjustable in the dialog).

## 2. Run the script

**Scripts → User Scripts →** *MAE to Flare... / MAE to MOE... / Selected to 2D SD...*

## 3. Choose options in the dialog

Common fields: output folder, output name (no extension), protein atom
threshold, and "open in target when finished".

| Script | Format choices | Extra options |
| --- | --- | --- |
| Flare | Maestro+SDF (recommended) - Flare project (`.flrp`) | Schrödinger & Flare folders |
| MOE | PDB+SDF (recommended) - MOE session (`.moe`) | Schrödinger & MOE folders |
| 2D SD | Force V2000 (recommended) | skip proteins - strip H - neutralize (aggressive/balanced) - add SMILES tag |

Your choices are saved to a `*.settings.json` next to the script and restored
next time.

## 4. Review the result

A summary dialog reports what was written and offers **Reveal in Explorer /
Finder** and **Open folder** buttons. For the 2D converter it also reports how
many structures were written, skipped, failed, and neutralized — watch for a
warning if any molecules were dropped.

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| "No entries are selected and nothing is included" | Select rows in the Project Table or include entries in the Workspace. |
| Target tool not found | Correct the tool's install folder in the dialog. |
| `.flrp` export fails | Known `pyflare` SQLite bug — the script auto-falls back to Maestro+SDF. |
| Poses imported as a protein | Raise the protein atom threshold, or check the receptor sits *before* its poses in the table. |
| 2D conversion warns about dropped molecules | Check the Maestro terminal for `[Mae2SD2D]` lines naming the failing structures. |

## Console diagnostics

Each script prints progress to the Maestro terminal, prefixed by its name:
`[Mae2Flare]`, `[Mae2Moe]`, `[Mae2SD2D]`, and (for the Flare driver) `[pyflare]`.
Include these when filing an issue.
