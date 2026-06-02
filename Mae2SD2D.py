#Name: User Scripts:Selected to 2D SD...
#Command: pythonrun Mae2SD2D.panel
"""
Mae2SD2D.py  -  Maestro user script  (Schrodinger 2026-1 ready)

Convert the SELECTED ligand / docking-pose entries in the Project Table
to a single 2D SD file (V2000), with clean heavy-atom-only depictions
and an optional charge-neutralization pass (aggressive or balanced).

This revision is robust against per-structure RDKit failures:
  * neutralize_mol: non-strict sanitize; reverts to original mol on error
  * remove_hs_mol : retries with sanitize=False if strict RemoveHs fails
  * compute_2d    : retries after non-strict sanitize
  * panel()       : pops up a loud warning if any molecules were dropped
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from schrodinger import structure
from schrodinger.maestro import maestro

# ---------------------------------------------------------------------------
# Adapter probe (Schrodinger 2026-1)
# ---------------------------------------------------------------------------

_TO_RDKIT      = None
_RDKIT_OPTIONS = None
_GEN2D_ENABLE  = None
_ADAPTER_PATH  = None

for _modname in ("schrodinger.adapter", "schrodinger.rdkit_extensions"):
    try:
        _m = __import__(_modname, fromlist=["*"])
    except Exception:
        continue
    _to_rdkit = getattr(_m, "to_rdkit", None)
    if callable(_to_rdkit):
        _TO_RDKIT      = _to_rdkit
        _RDKIT_OPTIONS = getattr(_m, "RDKitOptions", None)
        _GEN2D_ENABLE  = getattr(_m, "Generate2DCoordinates_Enable", None)
        _ADAPTER_PATH  = _modname
        break

# RDKit
_RDKIT = None
try:
    from rdkit import Chem as _Chem
    from rdkit.Chem import AllChem
    _RDKIT = _Chem
except Exception:
    pass

# RDKit MolStandardize Uncharger (preferred neutralizer)
_rdMS      = None
_UNCHARGER = None
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as _rdMS
    _UNCHARGER = _rdMS.Uncharger()
except Exception:
    _UNCHARGER = None

_HAS_PIPELINE = bool(_RDKIT and _TO_RDKIT)

print(f"[Mae2SD2D] adapter           : {_ADAPTER_PATH or '(not available)'}")
print(f"[Mae2SD2D] RDKit             : "
      f"{'available' if _RDKIT else 'NOT available (script cannot run)'}")
print(f"[Mae2SD2D] neutralizer       : "
      f"{'rdMolStandardize.Uncharger' if _UNCHARGER else ('SMARTS fallback' if _HAS_PIPELINE else 'NO')}")

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
    "output_dir":      r"C:\temp",
    "output_name":     "ligands_2d",
    "skip_proteins":   True,
    "protein_atom_threshold": 500,
    "strip_hydrogens": True,
    "neutralize":      False,
    "neutralize_mode": "aggressive",     # "aggressive" or "balanced"
    "force_v2000":     True,
    "add_smiles_tag":  True,
}

SETTINGS_FILE = Path(__file__).with_suffix(".settings.json")


def load_settings() -> dict:
    s = dict(DEFAULTS)
    if SETTINGS_FILE.exists():
        try:
            s.update(json.loads(SETTINGS_FILE.read_text()))
        except Exception:
            pass
    if s.get("neutralize_mode") not in ("aggressive", "balanced"):
        s["neutralize_mode"] = "aggressive"
    return s


def save_settings(s: dict) -> None:
    try:
        SETTINGS_FILE.write_text(json.dumps(s, indent=2))
    except Exception as e:
        print(f"[Mae2SD2D] could not save settings: {e}")


# ---------------------------------------------------------------------------
# Selection
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


# ---------------------------------------------------------------------------
# RDKit-level helpers
# ---------------------------------------------------------------------------

def st_to_mol(st):
    if _TO_RDKIT is None:
        return None
    try:
        return _TO_RDKIT(st)
    except Exception as e:
        print(f"[Mae2SD2D] to_rdkit failed for "
              f"'{getattr(st, 'title', '?')}': {e}")
        return None


# SMARTS-based neutralization patterns
_NEUTRALIZE_PATTERNS = None
def _compile_neutralize_patterns():
    global _NEUTRALIZE_PATTERNS
    if _NEUTRALIZE_PATTERNS is not None or _RDKIT is None:
        return _NEUTRALIZE_PATTERNS
    raw = [
        # Deprotonated acids -> reprotonate
        ("[O-][C,S,P,N]",                                0),
        ("[n-]",                                         0),
        ("[N-;X2]=[C,N]",                                0),
        ("[N-;X1]#[C]",                                  0),
        # Protonated bases -> deprotonate
        ("[$([N+;H,H2,H3,H4]);!$([N+]=*);!$([N+]#*)]",   0),
        # Aromatic N+: match either H1/H2 explicit OR any trivalent n+
        # that isn't part of an n-oxide.
        ("[$([n+;H1,H2]);!$([n+][O-])]",                 0),
        ("[$([n+;X3]);!$([n+][O-])]",                    0),
        ("[$([P+;H,H2,H3,H4])]",                         0),
        ("[$([S+;H,H2,H3])]",                            0),
    ]
    pats = []
    for smarts, target in raw:
        try:
            patt = _RDKIT.MolFromSmarts(smarts)
            if patt is not None:
                pats.append((patt, target))
        except Exception:
            continue
    _NEUTRALIZE_PATTERNS = pats
    return pats


def _neutralize_with_smarts(mol):
    """Adjust formal charges + explicit Hs on protonatable atoms.
    Returns (mol, changed_bool)."""
    if _RDKIT is None:
        return mol, False
    pats = _compile_neutralize_patterns()
    if not pats:
        return mol, False
    rw = _RDKIT.RWMol(mol)
    changed = False
    for patt, target in pats:
        for match in rw.GetSubstructMatches(patt):
            a = rw.GetAtomWithIdx(match[0])
            current = a.GetFormalCharge()
            if current == target:
                continue
            delta = target - current
            a.SetFormalCharge(target)
            try:
                nh = a.GetNumExplicitHs()
                a.SetNumExplicitHs(max(0, nh - delta))
                a.SetNoImplicit(False)
                changed = True
            except Exception:
                a.SetFormalCharge(current)
    if changed:
        try:
            _RDKIT.SanitizeMol(rw)
        except Exception as e:
            print(f"[Mae2SD2D] SMARTS-neutralize sanitize warning: {e}")
    return rw.GetMol(), changed


def _make_forced_uncharger():
    """Build an aggressive Uncharger across RDKit API revisions."""
    if _rdMS is None:
        return None
    for kwargs in (
        {"canonicalOrdering": True, "force": True},
        {"force": True},
    ):
        try:
            return _rdMS.Uncharger(**kwargs)
        except TypeError:
            continue
        except Exception:
            continue
    try:
        return _rdMS.Uncharger(True)
    except Exception:
        return None


def neutralize_mol(mol, mode: str = "aggressive"):
    """Return (mol, was_changed).
       mode == "aggressive" -> force-Uncharger + SMARTS sweep
       mode == "balanced"   -> default Uncharger only (preserves net charge)
    Robust: if anything blows up, returns the original mol untouched
    instead of raising, so downstream steps still get a valid Mol.
    """
    if _RDKIT is None or mol is None:
        return mol, False

    original = _RDKIT.Mol(mol)  # deep copy for fallback
    before   = [a.GetFormalCharge() for a in mol.GetAtoms()]

    try:
        if mode == "balanced":
            if _UNCHARGER is not None:
                mol = _UNCHARGER.uncharge(mol)
        else:
            forced = _make_forced_uncharger()
            if forced is not None:
                try:
                    mol = forced.uncharge(mol)
                except Exception as e:
                    print(f"[Mae2SD2D] forced Uncharger failed ({e}); "
                          "trying default Uncharger")
                    if _UNCHARGER is not None:
                        try:
                            mol = _UNCHARGER.uncharge(mol)
                        except Exception:
                            pass
            elif _UNCHARGER is not None:
                try:
                    mol = _UNCHARGER.uncharge(mol)
                except Exception:
                    pass
            mol, _ = _neutralize_with_smarts(mol)

        # Non-strict valence / sanitize: don't let one bad atom kill the Mol.
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception as e:
            print(f"[Mae2SD2D] UpdatePropertyCache(strict=False) warning: {e}")
        try:
            _RDKIT.SanitizeMol(
                mol,
                sanitizeOps=(_RDKIT.SanitizeFlags.SANITIZE_ALL
                             ^ _RDKIT.SanitizeFlags.SANITIZE_PROPERTIES),
            )
        except Exception as e:
            print(f"[Mae2SD2D] post-neutralize SanitizeMol warning: {e}; "
                  "reverting this structure to its un-neutralized form")
            return original, False

    except Exception as e:
        print(f"[Mae2SD2D] neutralize_mol unexpected error ({e}); "
              "reverting to original mol")
        return original, False

    after = [a.GetFormalCharge() for a in mol.GetAtoms()]
    return mol, after != before


def remove_hs_mol(mol):
    """Strip hydrogens robustly. RDKit's default RemoveHs() will throw if
    valence/charge state is inconsistent (often happens right after the
    SMARTS neutralization pass). We retry with sanitize=False and then
    re-sanitize separately so a bad atom doesn't drop the molecule.
    """
    if _RDKIT is None or mol is None:
        return mol
    try:
        return _RDKIT.RemoveHs(mol)
    except Exception as e:
        print(f"[Mae2SD2D] RemoveHs strict failed ({e}); retrying "
              "with sanitize=False")
    try:
        out = _RDKIT.RemoveHs(mol, sanitize=False)
        try:
            out.UpdatePropertyCache(strict=False)
            _RDKIT.SanitizeMol(
                out,
                sanitizeOps=(_RDKIT.SanitizeFlags.SANITIZE_ALL
                             ^ _RDKIT.SanitizeFlags.SANITIZE_PROPERTIES),
            )
        except Exception as e2:
            print(f"[Mae2SD2D] post-RemoveHs sanitize warning: {e2}")
        return out
    except Exception as e:
        print(f"[Mae2SD2D] RemoveHs(sanitize=False) also failed ({e}); "
              "returning Mol with Hs intact")
        return mol


def compute_2d(mol):
    """Lay out the Mol in 2D. Belt-and-braces against sanitize failures."""
    if _RDKIT is None or mol is None:
        return mol
    try:
        mol.RemoveAllConformers()
    except Exception:
        pass
    try:
        AllChem.Compute2DCoords(mol)
        return mol
    except Exception as e:
        print(f"[Mae2SD2D] Compute2DCoords strict failed ({e}); "
              "retrying after non-strict sanitize")
    try:
        mol.UpdatePropertyCache(strict=False)
        _RDKIT.SanitizeMol(
            mol,
            sanitizeOps=(_RDKIT.SanitizeFlags.SANITIZE_ALL
                         ^ _RDKIT.SanitizeFlags.SANITIZE_PROPERTIES),
        )
        AllChem.Compute2DCoords(mol)
    except Exception as e2:
        print(f"[Mae2SD2D] Compute2DCoords retry also failed ({e2}); "
              "molecule will be written with whatever coords it has")
    return mol


def copy_props_st_to_mol(st, mol):
    if mol is None:
        return mol
    try:
        for k, v in st.property.items():
            try:
                mol.SetProp(str(k), str(v))
            except Exception:
                pass
    except Exception:
        pass
    title = (st.title or "").strip()
    if title:
        try:
            mol.SetProp("_Name", title)
        except Exception:
            pass
    return mol


def add_smiles_prop(mol, key="s_user_SMILES"):
    if _RDKIT is None or mol is None:
        return
    try:
        smi = _RDKIT.MolToSmiles(mol)
        if smi:
            mol.SetProp(key, smi)
    except Exception as e:
        print(f"[Mae2SD2D] could not compute SMILES: {e}")


# ---------------------------------------------------------------------------
# SD writing (direct from RDKit Mol; coords + charges preserved)
# ---------------------------------------------------------------------------

def write_mols_to_sdf(mols, out_path: Path, force_v2000: bool = True) -> int:
    writer = _RDKIT.SDWriter(str(out_path))
    if force_v2000:
        try:
            writer.SetForceV3000(False)
        except Exception:
            pass
    n = 0
    for m in mols:
        if m is None:
            continue
        try:
            writer.write(m)
            n += 1
        except Exception as e:
            name = "?"
            try:
                name = m.GetPropsAsDict().get("_Name", "?")
            except Exception:
                pass
            print(f"[Mae2SD2D] SDWriter.write failed for '{name}': {e}")
    writer.close()
    return n


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
        self.setWindowTitle("Export selected entries as 2D SD")
        self.settings = dict(settings)

        if _HAS_PIPELINE:
            gen_desc = ("to_rdkit -> (optional neutralize -> RemoveHs) -> "
                        "Compute2DCoords -> SDWriter")
        else:
            gen_desc = ("<font color=red>NOT AVAILABLE - needs RDKit and "
                        "schrodinger.adapter.to_rdkit</font>")

        if _UNCHARGER is not None:
            neut_desc = "rdMolStandardize.Uncharger"
        elif _HAS_PIPELINE:
            neut_desc = "SMARTS fallback (RDKit MolStandardize not built)"
        else:
            neut_desc = "<font color=red>not available</font>"

        form = QtWidgets.QFormLayout()
        form.addRow(QtWidgets.QLabel(
            f"<b>{n_selected}</b> entries will be processed "
            f"(in Project Table order)."))
        form.addRow(QtWidgets.QLabel(f"<i>Pipeline:</i> {gen_desc}"))
        form.addRow(QtWidgets.QLabel(f"<i>Neutralizer:</i> {neut_desc}"))

        self.out_dir_edit = QtWidgets.QLineEdit(self.settings["output_dir"])
        b1 = QtWidgets.QPushButton("Browse...")
        b1.clicked.connect(self._pick_dir)
        r1 = QtWidgets.QHBoxLayout(); r1.addWidget(self.out_dir_edit); r1.addWidget(b1)
        form.addRow("Output folder:", self._wrap(r1))

        self.out_name_edit = QtWidgets.QLineEdit(self.settings["output_name"])
        form.addRow("Output name (no extension):", self.out_name_edit)

        # --- Neutralization group ----------------------------------------
        self.neutralize_chk = QtWidgets.QCheckBox(
            "Neutralize charges (acids/bases) before 2D layout")
        self.neutralize_chk.setChecked(bool(self.settings.get("neutralize", False)))
        self.neutralize_chk.setEnabled(_HAS_PIPELINE)
        form.addRow("", self.neutralize_chk)

        self.mode_aggressive = QtWidgets.QRadioButton(
            "Aggressive - remove every protonatable charge "
            "(pyridinium -> pyridine, ammonium -> amine)")
        self.mode_balanced = QtWidgets.QRadioButton(
            "Balanced - only neutralize if net charge is preserved "
            "(carboxylate <-> ammonium pairs)")
        mode_group = QtWidgets.QButtonGroup(self)
        mode_group.addButton(self.mode_aggressive)
        mode_group.addButton(self.mode_balanced)

        cur_mode = self.settings.get("neutralize_mode", "aggressive")
        if cur_mode == "balanced":
            self.mode_balanced.setChecked(True)
        else:
            self.mode_aggressive.setChecked(True)

        mode_box = QtWidgets.QVBoxLayout()
        mode_box.setContentsMargins(24, 0, 0, 0)
        mode_box.addWidget(self.mode_aggressive)
        mode_box.addWidget(self.mode_balanced)
        mode_wrap = QtWidgets.QWidget()
        mode_wrap.setLayout(mode_box)
        form.addRow("", mode_wrap)

        def _sync_mode_enabled():
            on = self.neutralize_chk.isChecked() and _HAS_PIPELINE
            self.mode_aggressive.setEnabled(on)
            self.mode_balanced.setEnabled(on)
        self.neutralize_chk.toggled.connect(lambda _=None: _sync_mode_enabled())
        _sync_mode_enabled()
        # -----------------------------------------------------------------

        self.strip_chk = QtWidgets.QCheckBox(
            "Strip hydrogens before 2D layout (clean depiction)")
        self.strip_chk.setChecked(bool(self.settings.get("strip_hydrogens", True)))
        form.addRow("", self.strip_chk)

        self.v2000_chk = QtWidgets.QCheckBox(
            "Force V2000 SDF (recommended)")
        self.v2000_chk.setChecked(bool(self.settings.get("force_v2000", True)))
        form.addRow("", self.v2000_chk)

        self.smiles_chk = QtWidgets.QCheckBox(
            "Add canonical SMILES as 's_user_SMILES' tag")
        self.smiles_chk.setChecked(bool(self.settings.get("add_smiles_tag", True)))
        form.addRow("", self.smiles_chk)

        self.skip_chk = QtWidgets.QCheckBox(
            "Skip proteins / receptors in the selection")
        self.skip_chk.setChecked(bool(self.settings.get("skip_proteins", True)))
        form.addRow("", self.skip_chk)

        self.thr_spin = QtWidgets.QSpinBox()
        self.thr_spin.setRange(50, 100000)
        self.thr_spin.setValue(int(self.settings["protein_atom_threshold"]))
        form.addRow("Protein atom threshold:", self.thr_spin)

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

    def values(self) -> dict:
        name = self.out_name_edit.text().strip() or "ligands_2d"
        for ext in (".sd", ".sdf"):
            if name.lower().endswith(ext):
                name = name[: -len(ext)]
                break
        mode = "balanced" if self.mode_balanced.isChecked() else "aggressive"
        return {
            "output_dir":             self.out_dir_edit.text().strip() or r"C:\temp",
            "output_name":            name,
            "skip_proteins":          bool(self.skip_chk.isChecked()),
            "protein_atom_threshold": int(self.thr_spin.value()),
            "strip_hydrogens":        bool(self.strip_chk.isChecked()),
            "neutralize":             bool(self.neutralize_chk.isChecked()),
            "neutralize_mode":        mode,
            "force_v2000":            bool(self.v2000_chk.isChecked()),
            "add_smiles_tag":         bool(self.smiles_chk.isChecked()),
        }


class SuccessDialog(QtWidgets.QDialog):
    def __init__(self, out_path: Path, n_written: int, n_skipped: int,
                 n_failed: int, n_neutralized: int = 0,
                 mode_used: str = "",
                 format_note: str = "", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export complete")
        self.out_path = Path(out_path)

        body = (f"<b>Wrote:</b> {n_written} ligand(s) "
                f"to<br><code>{self.out_path}</code>")
        if format_note:
            body += f"<br><br>{format_note}"
        if n_neutralized:
            extra = f" ({mode_used} mode)" if mode_used else ""
            body += (f"<br><br>Neutralized charges on "
                     f"<b>{n_neutralized}</b> structure(s){extra} "
                     f"(tagged <code>s_user_Neutralized = True</code>).")
        if n_skipped:
            body += f"<br><br>Skipped {n_skipped} protein/receptor entry/entries."
        if n_failed:
            body += (f"<br><br><i>{n_failed} structure(s) could not be "
                     f"converted; details in the Maestro Python Shell.</i>")

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
            reveal_in_file_manager(self.out_path)
        except Exception as e:
            QtWidgets.QMessageBox.warning(
                self, "Mae2SD2D", f"Could not reveal file:\n{e}")

    def _open_folder(self):
        try:
            open_folder(self.out_path)
        except Exception as e:
            QtWidgets.QMessageBox.warning(
                self, "Mae2SD2D", f"Could not open folder:\n{e}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def panel():
    settings = load_settings()

    try:
        entries = get_selected_structures()
    except Exception as e:
        QtWidgets.QMessageBox.critical(None, "Mae2SD2D", str(e))
        return

    if not _HAS_PIPELINE:
        QtWidgets.QMessageBox.critical(
            None, "Mae2SD2D",
            "This script needs both RDKit and schrodinger.adapter.to_rdkit, "
            "neither of which was found in this Maestro. Cannot continue.")
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
    final_path = out_dir / f"{settings['output_name']}.sdf"

    threshold      = int(settings["protein_atom_threshold"])
    skip_proteins  = bool(settings["skip_proteins"])
    strip_h        = bool(settings["strip_hydrogens"])
    neutralize     = bool(settings["neutralize"])
    neut_mode      = settings.get("neutralize_mode", "aggressive")
    force_v2000    = bool(settings["force_v2000"])
    add_smiles     = bool(settings["add_smiles_tag"])

    targets = []
    n_skipped = 0
    for idx, title, st in entries:
        if skip_proteins and is_protein(st, threshold):
            n_skipped += 1
            print(f"[Mae2SD2D] skip protein: row {idx}, '{title}', "
                  f"{st.atom_total} atoms")
            continue
        targets.append((idx, title, st))

    if not targets:
        QtWidgets.QMessageBox.warning(
            None, "Mae2SD2D",
            "Nothing to write - all selected entries were skipped "
            "(treated as proteins). Untick 'Skip proteins / receptors' "
            "if you want them written too.")
        return

    print(f"[Mae2SD2D] preparing {len(targets)} ligand(s) for {final_path}")
    print(f"[Mae2SD2D] strip H = {strip_h}, neutralize = {neutralize} "
          f"(mode={neut_mode}), "
          f"force V2000 = {force_v2000}, add SMILES tag = {add_smiles}")

    prepared = []
    n_failed       = 0
    n_neutralized  = 0
    for idx, title, st in targets:
        try:
            mol = st_to_mol(st)
            if mol is None:
                raise RuntimeError("to_rdkit returned None")

            was_neutralized = False
            if neutralize:
                mol, was_neutralized = neutralize_mol(mol, mode=neut_mode)
                if was_neutralized:
                    n_neutralized += 1

            if strip_h:
                mol = remove_hs_mol(mol)

            mol = compute_2d(mol)

            copy_props_st_to_mol(st, mol)
            if not mol.HasProp("_Name") and title:
                try:
                    mol.SetProp("_Name", title)
                except Exception:
                    pass
            if was_neutralized:
                try:
                    mol.SetProp("s_user_Neutralized", "True")
                except Exception:
                    pass
            if add_smiles:
                add_smiles_prop(mol)

            prepared.append(mol)
        except Exception as e:
            n_failed += 1
            print(f"[Mae2SD2D] FAILED row {idx} '{title}': {e}")

    n_written = write_mols_to_sdf(prepared, final_path,
                                  force_v2000=force_v2000)
    format_note = ("Format: <b>V2000 SDF</b>" if force_v2000
                   else "Format: SDF (V3000 may appear for large molecules)")

    # Loud warning if anything was dropped
    if n_failed and n_written < len(targets):
        QtWidgets.QMessageBox.warning(
            None, "Mae2SD2D - some structures were dropped",
            f"{n_written} of {len(targets)} ligand(s) were written to:\n"
            f"{final_path}\n\n"
            f"{n_failed} ligand(s) failed somewhere in the pipeline.\n"
            f"Open Maestro's Python Shell to see the per-structure "
            f"'[Mae2SD2D] FAILED row N ...' lines and the underlying RDKit "
            f"error messages."
        )

    SuccessDialog(final_path, n_written=n_written, n_skipped=n_skipped,
                  n_failed=n_failed, n_neutralized=n_neutralized,
                  mode_used=(neut_mode if neutralize else ""),
                  format_note=format_note).exec()


if __name__ == "__main__":
    panel()
