#Name: User Scripts:MAE to MOE...
#Command: pythonrun Mae2Moe.panel
"""
Mae2Moe.py  -  Maestro user script

Export the SELECTED Project Table entries (protein + poses, repeated)
to Chemical Computing Group's MOE (Molecular Operating Environment).

Strategy (mirrors the Flare exporter):
  * Receptors are written as PDB (.pdb)  -> MOE imports as a system
    into the Sequence Editor / MOE window, preserving chain/residue info.
  * Ligands/poses are written as SDF (.sdf) -> MOE's SDF reader treats
    each record as a small molecule (loaded into a MOE database .mdb,
    or into the MOE window depending on the chosen action).
  * MOE is launched with an auto-generated SVL startup script that:
        - opens the receptor PDB into the MOE window,
        - imports the ligand SDF into a new .mdb database,
        - saves a combined .moe session file next to the inputs.

Output modes:
  * "PDB + SDF (recommended)" - always works, MOE opens both natively.
  * "MOE session (.moe)"      - same files, but also runs an SVL script
    so MOE saves a .moe session file on disk you can re-open later.

Install:  Maestro  ->  Scripts  ->  Manage Scripts...  ->  Install...
Run:      Maestro  ->  Scripts  ->  User Scripts  ->  MAE to MOE...
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from schrodinger import structure
from schrodinger.maestro import maestro

try:
    from PyQt6 import QtWidgets, QtCore
    QT = 6
except ImportError:
    from PyQt5 import QtWidgets, QtCore
    QT = 5


# ---------------------------------------------------------------------------
# Defaults & persisted settings
# ---------------------------------------------------------------------------

DEFAULTS = {
    "schrodinger_dir": r"C:\Program Files\Schrodinger2026-2",
    "moe_dir":         r"C:\Program Files\Chemical Computing Group\MOE",
    "output_dir":      r"C:\temp",
    "output_name":     "maestro_export",      # extension chosen by format
    "output_format":   "split",               # "split" (default) or "moe"
    "protein_atom_threshold": 500,
    "open_in_moe":     True,
}

SETTINGS_FILE = Path(__file__).with_suffix(".settings.json")


def load_settings() -> dict:
    s = dict(DEFAULTS)
    if SETTINGS_FILE.exists():
        try:
            s.update(json.loads(SETTINGS_FILE.read_text()))
        except Exception:
            pass
    return s


def save_settings(s: dict) -> None:
    try:
        SETTINGS_FILE.write_text(json.dumps(s, indent=2))
    except Exception as e:
        print(f"[Mae2Moe] could not save settings: {e}")


# ---------------------------------------------------------------------------
# Selection -> ordered list of structures
# ---------------------------------------------------------------------------

def get_selected_structures():
    pt = maestro.project_table_get()

    rows = [r for r in pt.selected_rows]
    if not rows:
        rows = [r for r in pt.included_rows]
    if not rows:
        raise RuntimeError(
            "No entries are selected and nothing is included in the Workspace."
        )

    out = []
    for r in rows:
        st = r.getStructure()
        title = st.title or f"entry_{r.index}"
        out.append((r.index, title, st))
    out.sort(key=lambda x: x[0])
    return out


def is_protein(st, atom_threshold: int) -> bool:
    AA = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
    }
    try:
        for res in st.residue:
            if res.pdbres.strip().upper() in AA:
                return True
    except Exception:
        pass
    return st.atom_total > atom_threshold


def group_by_protein(entries, atom_threshold: int):
    groups = []
    current = None
    for idx, title, st in entries:
        if is_protein(st, atom_threshold):
            current = {"protein": st, "ligands": [], "name": title}
            groups.append(current)
        else:
            if current is None:
                current = {"protein": None, "ligands": [],
                           "name": "ligands_only"}
                groups.append(current)
            current["ligands"].append(st)
    return groups


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_split_native(groups, out_dir: Path, base_name: str):
    """Write receptors as .pdb and ligands as .sdf.

    Returns (proteins_path or None, ligands_path or None).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    proteins_path = out_dir / f"{base_name}_proteins.pdb"
    ligands_path  = out_dir / f"{base_name}_ligands.sdf"

    have_p = False
    have_l = False

    pw = structure.StructureWriter(str(proteins_path))
    lw = structure.StructureWriter(str(ligands_path))
    try:
        for gi, g in enumerate(groups, 1):
            if g["protein"] is not None:
                pw.append(g["protein"])
                have_p = True
            for lig in g["ligands"]:
                lw.append(lig)
                have_l = True
    finally:
        pw.close()
        lw.close()

    if not have_p:
        try: proteins_path.unlink()
        except OSError: pass
        proteins_path = None
    if not have_l:
        try: ligands_path.unlink()
        except OSError: pass
        ligands_path = None

    return proteins_path, ligands_path


# ---------------------------------------------------------------------------
# SVL startup script (only used when output_format == "moe")
# ---------------------------------------------------------------------------

SVL_TEMPLATE = "\n".join([
    '// Mae2Moe.svl  --  auto-generated by Mae2Moe.py',
    '//',
    '// Opens the exported PDB receptor in the MOE window, imports the SDF',
    '// ligands into a fresh .mdb database, and saves a combined .moe session',
    '// file on disk so the user can re-open the whole assembly later.',
    '',
    'local function main []',
    "    local protein_pdb = '__PROTEIN_PDB__';",
    "    local ligand_sdf  = '__LIGAND_SDF__';",
    "    local out_moe     = '__OUT_MOE__';",
    "    local out_mdb     = '__OUT_MDB__';",
    '',
    '    // Load receptor into the MOE window (if provided).',
    "    if length protein_pdb > 0 and ftype protein_pdb == 'file' then",
    '        ReadAuto protein_pdb;',
    '    endif',
    '',
    '    // Import ligand SDF into a new MOE database (.mdb).',
    "    if length ligand_sdf > 0 and ftype ligand_sdf == 'file' then",
    "        if ftype out_mdb == 'file' then fdelete out_mdb; endif",
    '        db_ImportSD [out_mdb, ligand_sdf, []];',
    '        // Also load ligands into the MOE window alongside the protein',
    '        // so the user immediately sees the complex.',
    '        ReadAuto ligand_sdf;',
    '    endif',
    '',
    '    // Save the combined session.',
    '    SaveAs out_moe;',
    '',
    "    write ['[Mae2Moe] saved session: {}\\n', out_moe];",
    'endfunction',
    '',
    'main [];'
])


def write_svl_startup(svl_path: Path, protein_pdb, ligand_sdf,
                      out_moe: Path, out_mdb: Path) -> None:
    def svl_str(p):
        if p is None:
            return ""
        # SVL strings use single quotes; forward slashes are safe on Windows.
        return str(p).replace("\\", "/")

    body = (SVL_TEMPLATE
            .replace("__PROTEIN_PDB__", svl_str(protein_pdb))
            .replace("__LIGAND_SDF__",  svl_str(ligand_sdf))
            .replace("__OUT_MOE__",     svl_str(out_moe))
            .replace("__OUT_MDB__",     svl_str(out_mdb)))
    svl_path.write_text(body, encoding="utf-8")


# ---------------------------------------------------------------------------
# Locate MOE executable
# ---------------------------------------------------------------------------

def find_executable(moe_dir: Path, names) -> Path:
    if isinstance(names, str):
        names = [names]
    for sub in ("", "bin", "bin-win64", "bin-lnx64"):
        for n in names:
            p = (moe_dir / sub / n) if sub else (moe_dir / n)
            if p.exists():
                return p
    for n in names:
        for p in moe_dir.rglob(n):
            return p
    raise FileNotFoundError(f"Could not find any of {names} inside {moe_dir}.")


def find_moe_gui(moe_dir: Path) -> Path:
    if os.name == "nt":
        return find_executable(moe_dir, ["moe.exe", "MOE.exe"])
    return find_executable(moe_dir, ["moe", "MOE"])


# ---------------------------------------------------------------------------
# Launching MOE
# ---------------------------------------------------------------------------

def open_in_moe(moe_dir: Path, paths, svl_script: Path = None) -> None:
    """Launch MOE on one or more file paths.

    MOE's GUI accepts file arguments on its command line and will dispatch
    them to the correct reader (.pdb -> system, .sdf -> molecule, .moe ->
    session, .mdb -> database browser).

    If svl_script is provided, MOE is invoked with '-run <script>' so the
    script runs at startup -- this is how we build a .moe session.
    """
    if isinstance(paths, (str, Path)):
        paths = [paths]
    paths = [Path(p) for p in paths if p is not None]

    moe_exe = find_moe_gui(moe_dir)
    moe_exe_dir = moe_exe.parent

    cmd = [str(moe_exe)]
    if svl_script is not None:
        cmd += ["-run", str(svl_script)]
    else:
        cmd += [str(p) for p in paths]

    print("[Mae2Moe] launching: " + " ".join(f'"{c}"' for c in cmd))

    env = os.environ.copy()

    # Strip Schrödinger dirs from PATH so they don't shadow MOE's own Qt/
    # Python DLLs (MOE ships its own runtime, same hazard as Flare).
    schro_root = (os.environ.get("SCHRODINGER") or "").strip()
    keep = []
    for part in env.get("PATH", "").split(os.pathsep):
        if not part:
            continue
        p_lower = part.lower()
        if schro_root and p_lower.startswith(schro_root.lower()):
            continue
        if "schrodinger" in p_lower:
            continue
        keep.append(part)
    keep.insert(0, str(moe_exe_dir))
    env["PATH"] = os.pathsep.join(keep)

    for var in ("QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH",
                "QT_QPA_PLATFORM", "PYTHONPATH",
                "QML2_IMPORT_PATH", "QML_IMPORT_PATH",
                "PYTHONHOME"):
        env.pop(var, None)

    if os.name == "nt":
        DETACHED_PROCESS = 0x00000008
        CREATE_NEW_PROCESS_GROUP = 0x00000200
        subprocess.Popen(
            cmd,
            cwd=str(moe_exe_dir),
            env=env,
            creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
            close_fds=True,
        )
    else:
        subprocess.Popen(
            cmd,
            cwd=str(moe_exe_dir),
            env=env,
            start_new_session=True,
            close_fds=True,
        )


# ---------------------------------------------------------------------------
# Reveal helpers
# ---------------------------------------------------------------------------

def reveal_in_file_manager(path: Path) -> None:
    path = Path(path)
    if os.name == "nt":
        subprocess.Popen(["explorer", f"/select,{str(path)}"], close_fds=True)
    elif sys.platform == "darwin":
        subprocess.Popen(["open", "-R", str(path)], close_fds=True)
    else:
        subprocess.Popen(["xdg-open", str(path.parent)], close_fds=True)


def open_folder(path: Path) -> None:
    folder = Path(path)
    if folder.is_file():
        folder = folder.parent
    if os.name == "nt":
        os.startfile(str(folder))                          # noqa: S606
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(folder)], close_fds=True)
    else:
        subprocess.Popen(["xdg-open", str(folder)], close_fds=True)


# ---------------------------------------------------------------------------
# Dialogs
# ---------------------------------------------------------------------------

class ExportDialog(QtWidgets.QDialog):
    def __init__(self, settings: dict, n_selected: int, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export selected entries to MOE")
        self.settings = dict(settings)

        form = QtWidgets.QFormLayout()
        form.addRow(QtWidgets.QLabel(
            f"<b>{n_selected}</b> entries will be exported "
            f"(in Project Table order)."))

        self.out_dir_edit = QtWidgets.QLineEdit(self.settings["output_dir"])
        b1 = QtWidgets.QPushButton("Browse...")
        b1.clicked.connect(self._pick_dir)
        r1 = QtWidgets.QHBoxLayout(); r1.addWidget(self.out_dir_edit); r1.addWidget(b1)
        form.addRow("Output folder:", self._wrap(r1))

        self.out_name_edit = QtWidgets.QLineEdit(self.settings["output_name"])
        form.addRow("Output name (no extension):", self.out_name_edit)

        self.fmt_combo = QtWidgets.QComboBox()
        self.fmt_combo.addItem(
            "PDB proteins + SDF ligands  -- recommended", "split")
        self.fmt_combo.addItem(
            "MOE session (.moe via SVL startup script)", "moe")
        cur = settings.get("output_format", "split")
        idx = self.fmt_combo.findData(cur)
        if idx >= 0:
            self.fmt_combo.setCurrentIndex(idx)
        form.addRow("Output format:", self.fmt_combo)

        self.schro_edit = QtWidgets.QLineEdit(self.settings["schrodinger_dir"])
        b2 = QtWidgets.QPushButton("Browse...")
        b2.clicked.connect(lambda: self._pick_dir_into(self.schro_edit))
        r2 = QtWidgets.QHBoxLayout(); r2.addWidget(self.schro_edit); r2.addWidget(b2)
        form.addRow("Schrödinger folder:", self._wrap(r2))

        self.moe_edit = QtWidgets.QLineEdit(self.settings["moe_dir"])
        b3 = QtWidgets.QPushButton("Browse...")
        b3.clicked.connect(lambda: self._pick_dir_into(self.moe_edit))
        r3 = QtWidgets.QHBoxLayout(); r3.addWidget(self.moe_edit); r3.addWidget(b3)
        form.addRow("MOE folder:", self._wrap(r3))

        self.thr_spin = QtWidgets.QSpinBox()
        self.thr_spin.setRange(50, 100000)
        self.thr_spin.setValue(int(self.settings["protein_atom_threshold"]))
        form.addRow("Protein atom threshold:", self.thr_spin)

        self.open_chk = QtWidgets.QCheckBox("Open result in MOE when finished")
        self.open_chk.setChecked(bool(self.settings.get("open_in_moe", True)))
        form.addRow("", self.open_chk)

        if QT == 6:
            bb = QtWidgets.QDialogButtonBox(
                QtWidgets.QDialogButtonBox.StandardButton.Ok |
                QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        else:
            bb = QtWidgets.QDialogButtonBox(
                QtWidgets.QDialogButtonBox.Ok |
                QtWidgets.QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)

        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(form); lay.addWidget(bb)

    def _wrap(self, sublayout):
        w = QtWidgets.QWidget(); w.setLayout(sublayout); return w

    def _pick_dir(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Choose output folder",
            self.out_dir_edit.text() or r"C:\temp")
        if d:
            self.out_dir_edit.setText(d)

    def _pick_dir_into(self, edit):
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Choose folder", edit.text() or "")
        if d:
            edit.setText(d)

    def values(self) -> dict:
        name = self.out_name_edit.text().strip() or "maestro_export"
        for ext in (".moe", ".mdb", ".pdb", ".maegz", ".mae", ".sdf"):
            if name.lower().endswith(ext):
                name = name[: -len(ext)]
                break
        return {
            "output_dir":             self.out_dir_edit.text().strip() or r"C:\temp",
            "output_name":            name,
            "output_format":          self.fmt_combo.currentData(),
            "schrodinger_dir":        self.schro_edit.text().strip(),
            "moe_dir":                self.moe_edit.text().strip(),
            "protein_atom_threshold": int(self.thr_spin.value()),
            "open_in_moe":            bool(self.open_chk.isChecked()),
        }


class SuccessDialog(QtWidgets.QDialog):
    def __init__(self, paths, opened_in_moe: bool,
                 note: str = "", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export complete")
        if isinstance(paths, (str, Path)):
            paths = [paths]
        self.paths = [Path(p) for p in paths]
        self.primary_path = self.paths[0]

        list_html = "<br>".join(f"<code>{p}</code>" for p in self.paths)
        body = f"<b>Saved:</b><br>{list_html}"
        if note:
            body += f"<br><br><i>{note}</i>"
        if opened_in_moe:
            body += "<br><br>Opening in MOE..."

        msg = QtWidgets.QLabel(body)
        if QT == 6:
            msg.setTextInteractionFlags(
                QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        else:
            msg.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)

        btn_reveal = QtWidgets.QPushButton("Reveal in Explorer")
        btn_reveal.clicked.connect(self._reveal)
        btn_folder = QtWidgets.QPushButton("Open folder")
        btn_folder.clicked.connect(self._open_folder)
        btn_ok = QtWidgets.QPushButton("OK")
        btn_ok.setDefault(True); btn_ok.clicked.connect(self.accept)

        if sys.platform == "darwin":
            btn_reveal.setText("Reveal in Finder")
        elif os.name != "nt":
            btn_reveal.setText("Show in file manager")

        row = QtWidgets.QHBoxLayout()
        row.addWidget(btn_reveal); row.addWidget(btn_folder)
        row.addStretch(1); row.addWidget(btn_ok)

        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(msg); lay.addSpacing(8); lay.addLayout(row)

    def _reveal(self):
        try:
            reveal_in_file_manager(self.primary_path)
        except Exception as e:
            QtWidgets.QMessageBox.warning(
                self, "Mae2Moe", f"Could not reveal file:\n{e}")

    def _open_folder(self):
        try:
            open_folder(self.primary_path)
        except Exception as e:
            QtWidgets.QMessageBox.warning(
                self, "Mae2Moe", f"Could not open folder:\n{e}")


# ---------------------------------------------------------------------------
# Entry point Maestro calls
# ---------------------------------------------------------------------------

def panel():
    settings = load_settings()

    try:
        entries = get_selected_structures()
    except Exception as e:
        QtWidgets.QMessageBox.critical(None, "Mae2Moe", str(e))
        return

    dlg = ExportDialog(settings, n_selected=len(entries))
    accepted_code = (QtWidgets.QDialog.DialogCode.Accepted if QT == 6
                     else QtWidgets.QDialog.Accepted)
    if dlg.exec() != accepted_code:
        return

    new_settings = dlg.values()
    settings.update(new_settings)
    save_settings(settings)

    out_dir = Path(settings["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    moe_dir = Path(settings["moe_dir"])
    if not moe_dir.exists():
        QtWidgets.QMessageBox.critical(
            None, "Mae2Moe",
            f"MOE folder does not exist:\n{moe_dir}")
        return

    groups = group_by_protein(
        entries, atom_threshold=settings["protein_atom_threshold"])
    if not groups:
        QtWidgets.QMessageBox.warning(None, "Mae2Moe", "Nothing to export.")
        return

    summary = "\n".join(
        f"  group {i}: protein={'yes' if g['protein'] else 'no'}, "
        f"ligands={len(g['ligands'])}  ({g['name']})"
        for i, g in enumerate(groups, 1)
    )
    print(f"[Mae2Moe] {len(groups)} group(s):\n{summary}")

    # Always write the split native files; .moe mode just adds an SVL
    # startup script and a target .moe session path on top of them.
    prot, lig = write_split_native(
        groups, out_dir, settings["output_name"])
    final_paths = [p for p in (prot, lig) if p is not None]
    if not final_paths:
        QtWidgets.QMessageBox.warning(
            None, "Mae2Moe", "Nothing to export.")
        return

    svl_script = None
    note_to_user = ""
    fmt = settings.get("output_format", "split")
    if fmt == "moe":
        out_moe = out_dir / f"{settings['output_name']}.moe"
        out_mdb = out_dir / f"{settings['output_name']}_ligands.mdb"
        svl_script = out_dir / f"{settings['output_name']}_open.svl"
        write_svl_startup(svl_script, prot, lig, out_moe, out_mdb)
        # The .moe / .mdb won't exist until MOE has actually run the SVL.
        final_paths = [p for p in (prot, lig, svl_script) if p is not None]
        note_to_user = (
            "MOE will run the generated SVL script at startup to build "
            f"<code>{out_moe.name}</code> and "
            f"<code>{out_mdb.name}</code> in the same folder."
        )

    opened = False
    if settings.get("open_in_moe", True):
        try:
            open_in_moe(moe_dir, final_paths, svl_script=svl_script)
            opened = True
        except Exception as e:
            QtWidgets.QMessageBox.warning(
                None, "Mae2Moe",
                f"Saved {final_paths},\nbut could not launch MOE:\n{e}")

    SuccessDialog(final_paths, opened_in_moe=opened, note=note_to_user).exec()


if __name__ == "__main__":
    panel()
