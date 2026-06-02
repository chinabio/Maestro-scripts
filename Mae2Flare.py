#Name: User Scripts:MAE to Flare...
#Command: pythonrun Mae2Flare.panel
"""
Mae2Flare.py  -  Maestro user script

Author:  Pulan Yu
Email:   chinabio@gmail.com

Export the SELECTED Project Table entries (protein + poses, repeated)
to Cresset Flare, ensuring ligands are imported as LIGANDS (not proteins).

Strategy:
  * Receptors are written as Maestro (.maegz) -> preserves bond orders,
    formal charges, and metadata Schrödinger has set on the protein.
  * Ligands/poses are written as SDF (.sdf)   -> Flare's SDF reader
    unconditionally treats each record as a Ligand, so we don't rely
    on its protein-vs-ligand heuristic.
  * Flare is launched with BOTH files on its command line; it imports
    them correctly into the Proteins and Ligands tables respectively.

Output modes:
  * "Maestro+SDF (recommended)" - the path above, always works.
  * "Flare project (.flrp)"     - calls pyflare to build a .flrp,
    following Cresset's fepcreate.py pattern (load-if-exists, save to
    .tmp + atomic rename). Falls back to Maestro+SDF if pyflare hits
    the known SQLite project-cache bug.

Install:  Maestro  ->  Scripts  ->  Manage Scripts...  ->  Install...
Run:      Maestro  ->  Scripts  ->  User Scripts  ->  MAE to Flare...
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
    "flare_dir":       r"C:\Program Files\Cresset-BMD\Flare",
    "output_dir":      r"C:\temp",
    "output_name":     "maestro_export",      # extension chosen by format
    "output_format":   "split",               # "split" (default) or "flrp"
    "protein_atom_threshold": 500,
    "open_in_flare":   True,
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
        print(f"[Mae2Flare] could not save settings: {e}")


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
    """Write receptors as .maegz and ligands as .sdf.

    Returns (proteins_path or None, ligands_path or None).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    proteins_path = out_dir / f"{base_name}_proteins.maegz"
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


def write_split_for_pyflare(groups, tmpdir: Path):
    """Write a proteins .maegz and a ligands .sdf for the pyflare driver.

    Each structure's title is tagged with " _gNN" so the pyflare side can
    reconstruct per-group Role membership after a single read of each file.
    Returns (proteins_path or None, ligands_path or None, group_names dict).
    """
    proteins_path = tmpdir / "all_proteins.maegz"
    ligands_path  = tmpdir / "all_ligands.sdf"
    group_names = {}

    have_p = False
    have_l = False

    pw = structure.StructureWriter(str(proteins_path))
    lw = structure.StructureWriter(str(ligands_path))
    try:
        for gi, g in enumerate(groups, 1):
            group_names[gi] = g["name"] or f"Group {gi}"
            if g["protein"] is not None:
                st = g["protein"]
                st.title = f"{(st.title or '').strip()} _g{gi:02d}"
                pw.append(st)
                have_p = True
            for lig in g["ligands"]:
                lig.title = f"{(lig.title or '').strip()} _g{gi:02d}"
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

    return proteins_path, ligands_path, group_names


# ---------------------------------------------------------------------------
# pyflare driver (only used when output_format == "flrp")
# ---------------------------------------------------------------------------

PYFLARE_SCRIPT = r'''
"""pyflare driver for Mae2Flare.py.

Follows Cresset's fepcreate.py pattern: load-if-exists, save to .tmp,
atomic-replace. Proteins come in via project.proteins.extend, ligands
via project.ligands.extend -- this is the documented API and matches
fepcreate.py exactly, so Flare classifies them correctly.
"""
import json
import os
import re
import sys
import traceback

try:
    from cresset import flare
except Exception as e:
    sys.stderr.write("Could not import cresset.flare: %s\n" % e)
    sys.exit(2)


GROUP_TAG_RE = re.compile(r"\s*_g(\d{2,})\s*$")


def log(msg):
    sys.stderr.write("[pyflare] " + msg + "\n")
    sys.stderr.flush()


def parse_group(title):
    if not title:
        return None
    m = GROUP_TAG_RE.search(title)
    return int(m.group(1)) if m else None


def strip_group_tag(title):
    return GROUP_TAG_RE.sub("", title).strip() if title else title


def safe_role_name(base, used):
    name = (base or "Group").strip() or "Group"
    n = name
    k = 2
    while n in used:
        n = "%s (%d)" % (name, k)
        k += 1
    used.add(n)
    return n


def assign_role(role, item):
    for meth in ("add", "append", "include"):
        fn = getattr(role, meth, None)
        if callable(fn):
            try:
                fn(item); return True
            except Exception:
                pass
    try:
        item.role = role; return True
    except Exception:
        pass
    roles_attr = getattr(item, "roles", None)
    if roles_attr is not None:
        for meth in ("add", "append", "include"):
            fn = getattr(roles_attr, meth, None)
            if callable(fn):
                try:
                    fn(role); return True
                except Exception:
                    pass
    return False


def main():
    if len(sys.argv) < 3:
        log("usage: pyflare driver.py <out.flrp> <group_names.json> "
            "[--proteins PATH] [--ligands PATH]")
        sys.exit(64)

    out_path   = sys.argv[1]
    names_path = sys.argv[2]

    proteins_path = ligands_path = None
    rest = sys.argv[3:]
    i = 0
    while i < len(rest):
        if rest[i] == "--proteins" and i + 1 < len(rest):
            proteins_path = rest[i + 1]; i += 2
        elif rest[i] == "--ligands" and i + 1 < len(rest):
            ligands_path = rest[i + 1]; i += 2
        else:
            log("unknown arg: %r" % rest[i]); sys.exit(64)

    with open(names_path, "r", encoding="utf-8") as f:
        group_names = {int(k): v for k, v in json.load(f).items()}

    log("output: %s" % out_path)
    log("proteins: %s" % proteins_path)
    log("ligands : %s" % ligands_path)

    if os.path.isfile(out_path):
        log("loading existing project")
        project = flare.Project.load(out_path)
    else:
        log("creating new project")
        project = flare.Project()

    used = set()
    role_by_group = {}
    for gi in sorted(group_names):
        rname = safe_role_name(group_names[gi], used)
        try:
            role_by_group[gi] = project.roles.append(rname)
            log("role %d -> %r" % (gi, rname))
        except Exception as e:
            role_by_group[gi] = None
            log("could not create role %r: %s" % (rname, e))

    proteins_added = []
    if proteins_path:
        proteins_added = list(
            project.proteins.extend(flare.read_file(proteins_path, None))
        )
        log("imported %d protein(s)" % len(proteins_added))

    ligands_added = []
    if ligands_path:
        ligands_added = list(
            project.ligands.extend(flare.read_file(ligands_path, None))
        )
        log("imported %d ligand(s)" % len(ligands_added))

    failures = 0
    for items in (proteins_added, ligands_added):
        for it in items:
            t = getattr(it, "title", "") or ""
            gi = parse_group(t)
            try:
                it.title = strip_group_tag(t)
            except Exception:
                pass
            if gi is None:
                continue
            role = role_by_group.get(gi)
            if role is None:
                continue
            if not assign_role(role, it):
                failures += 1
    if failures:
        log("note: %d item(s) not added to a Role "
            "(unsupported on this Flare version)." % failures)

    tmp_path = out_path + ".tmp"
    try:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        log("saving to %s" % tmp_path)
        project.save(tmp_path)
    except Exception as e:
        log("project.save(tmp) failed: %s" % e)
        sys.exit(10)

    try:
        os.replace(tmp_path, out_path)
        log("renamed %s -> %s" % (tmp_path, out_path))
    except Exception as e:
        log("rename failed: %s" % e)
        sys.exit(11)

    log("done -> %s" % out_path)


try:
    main()
except SystemExit:
    raise
except Exception:
    sys.stderr.write("Unhandled exception in pyflare driver:\n")
    traceback.print_exc()
    sys.exit(1)
'''


# ---------------------------------------------------------------------------
# Locate Flare executables and Qt plugin tree
# ---------------------------------------------------------------------------

def find_executable(flare_dir: Path, names) -> Path:
    if isinstance(names, str):
        names = [names]
    for sub in ("", "bin", "python", "scripts"):
        for n in names:
            p = (flare_dir / sub / n) if sub else (flare_dir / n)
            if p.exists():
                return p
    for n in names:
        for p in flare_dir.rglob(n):
            return p
    raise FileNotFoundError(f"Could not find any of {names} inside {flare_dir}.")


def find_pyflare(flare_dir: Path) -> Path:
    return find_executable(flare_dir, "pyflare.exe")


def find_flare_gui(flare_dir: Path) -> Path:
    return find_executable(flare_dir, ["Flare.exe", "flare.exe"])


def find_flare_qt_dirs(flare_dir: Path, flare_exe: Path):
    """Locate Flare's Qt 'plugins' dir and 'plugins\\platforms' dir.

    Handles both layouts:
      * plugins\\platforms\\qwindows.dll  (current Flare installs)
      * platforms\\qwindows.dll           (older / repackaged installs)
    """
    candidates = [
        flare_exe.parent / "plugins" / "platforms",
        flare_dir / "plugins" / "platforms",
        flare_dir / "bin" / "plugins" / "platforms",
        flare_exe.parent / "platforms",
        flare_dir / "platforms",
        flare_dir / "bin" / "platforms",
    ]
    for cand in candidates:
        if (cand / "qwindows.dll").exists():
            return (cand.parent, cand)
    for p in flare_dir.rglob("qwindows.dll"):
        if p.parent.name.lower() == "platforms":
            return (p.parent.parent, p.parent)
    return (None, None)


def run_pyflare(flare_dir: Path, out_flrp: Path, proteins_path,
                ligands_path, group_names: dict):
    """Returns (success: bool, message: str)."""
    pyflare = find_pyflare(flare_dir)
    driver_py = Path(tempfile.gettempdir()) / "mae_to_flare_driver.py"
    driver_py.write_text(PYFLARE_SCRIPT, encoding="utf-8")
    names_json = Path(tempfile.gettempdir()) / "mae_to_flare_group_names.json"
    names_json.write_text(
        json.dumps({str(k): v for k, v in group_names.items()}),
        encoding="utf-8",
    )
    args = [str(pyflare), str(driver_py), str(out_flrp), str(names_json)]
    if proteins_path is not None:
        args += ["--proteins", str(proteins_path)]
    if ligands_path is not None:
        args += ["--ligands", str(ligands_path)]

    print("[Mae2Flare] running:", " ".join(f'"{c}"' for c in args))
    result = subprocess.run(
        args, capture_output=True, text=True,
        creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
    )
    if result.stdout:
        print("[pyflare stdout]\n" + result.stdout)
    if result.stderr:
        print("[pyflare stderr]\n" + result.stderr)

    if result.returncode == 0 and out_flrp.exists():
        return True, "saved"
    return False, (
        f"pyflare exit {result.returncode}\n\n"
        f"STDERR:\n{result.stderr or '(empty)'}\n\n"
        f"STDOUT:\n{result.stdout or '(empty)'}"
    )


# ---------------------------------------------------------------------------
# Launching Flare (with a Qt-clean environment, multi-file capable)
# ---------------------------------------------------------------------------

def open_in_flare(flare_dir: Path, paths) -> None:
    """Launch Flare.exe on one or more file paths with an environment that
    won't let Maestro's Qt DLLs shadow Flare's own Qt plugin.

    `paths` may be a single Path or a list of Paths. Empty entries are
    skipped. Flare's GUI accepts multiple files on its command line.
    """
    if isinstance(paths, (str, Path)):
        paths = [paths]
    paths = [Path(p) for p in paths if p is not None]
    if not paths:
        raise ValueError("open_in_flare: no paths to open")

    flare_exe = find_flare_gui(flare_dir)
    flare_exe_dir = flare_exe.parent
    print(f"[Mae2Flare] launching: \"{flare_exe}\" " +
          " ".join(f'"{p}"' for p in paths))

    env = os.environ.copy()

    # 1) Remove Schrödinger dirs from PATH so they don't shadow Flare's Qt.
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
    # 2) Put Flare's own dir FIRST on PATH so its Qt DLLs win.
    keep.insert(0, str(flare_exe_dir))
    env["PATH"] = os.pathsep.join(keep)

    # 3) Drop env vars that override Qt's plugin search or Python.
    for var in ("QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH",
                "QT_QPA_PLATFORM", "PYTHONPATH",
                "QML2_IMPORT_PATH", "QML_IMPORT_PATH",
                "PYTHONHOME"):
        env.pop(var, None)

    # 4) Point Qt at Flare's own plugin tree.
    plugins_dir, platforms_dir = find_flare_qt_dirs(flare_dir, flare_exe)
    if platforms_dir is not None:
        env["QT_QPA_PLATFORM_PLUGIN_PATH"] = str(platforms_dir)
        env["QT_PLUGIN_PATH"] = str(plugins_dir)
        print(f"[Mae2Flare] Qt plugins dir : {plugins_dir}")
        print(f"[Mae2Flare] Qt platforms dir: {platforms_dir}")
    else:
        print(f"[Mae2Flare] WARNING: could not locate qwindows.dll under "
              f"{flare_dir}. Flare may fail to start.")

    cmd = [str(flare_exe)] + [str(p) for p in paths]

    if os.name == "nt":
        DETACHED_PROCESS = 0x00000008
        CREATE_NEW_PROCESS_GROUP = 0x00000200
        subprocess.Popen(
            cmd,
            cwd=str(flare_exe_dir),
            env=env,
            creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
            close_fds=True,
        )
    else:
        subprocess.Popen(
            cmd,
            cwd=str(flare_exe_dir),
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
        self.setWindowTitle("Export selected entries to Flare")
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
            "Maestro proteins + SDF ligands  -- recommended", "split")
        self.fmt_combo.addItem(
            "Flare project (.flrp)  -- may fail on buggy Flare builds", "flrp")
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

        self.flare_edit = QtWidgets.QLineEdit(self.settings["flare_dir"])
        b3 = QtWidgets.QPushButton("Browse...")
        b3.clicked.connect(lambda: self._pick_dir_into(self.flare_edit))
        r3 = QtWidgets.QHBoxLayout(); r3.addWidget(self.flare_edit); r3.addWidget(b3)
        form.addRow("Flare folder:", self._wrap(r3))

        self.thr_spin = QtWidgets.QSpinBox()
        self.thr_spin.setRange(50, 100000)
        self.thr_spin.setValue(int(self.settings["protein_atom_threshold"]))
        form.addRow("Protein atom threshold:", self.thr_spin)

        self.open_chk = QtWidgets.QCheckBox("Open result in Flare when finished")
        self.open_chk.setChecked(bool(self.settings.get("open_in_flare", True)))
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
        for ext in (".flrp", ".flr", ".maegz", ".mae", ".sdf"):
            if name.lower().endswith(ext):
                name = name[: -len(ext)]
                break
        return {
            "output_dir":             self.out_dir_edit.text().strip() or r"C:\temp",
            "output_name":            name,
            "output_format":          self.fmt_combo.currentData(),
            "schrodinger_dir":        self.schro_edit.text().strip(),
            "flare_dir":              self.flare_edit.text().strip(),
            "protein_atom_threshold": int(self.thr_spin.value()),
            "open_in_flare":          bool(self.open_chk.isChecked()),
        }


class SuccessDialog(QtWidgets.QDialog):
    def __init__(self, paths, opened_in_flare: bool,
                 note: str = "", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export complete")
        if isinstance(paths, (str, Path)):
            paths = [paths]
        self.paths = [Path(p) for p in paths]
        # Used by Reveal/Open-folder buttons:
        self.primary_path = self.paths[0]

        list_html = "<br>".join(f"<code>{p}</code>" for p in self.paths)
        body = f"<b>Saved:</b><br>{list_html}"
        if note:
            body += f"<br><br><i>{note}</i>"
        if opened_in_flare:
            body += "<br><br>Opening in Flare..."

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
                self, "Mae2Flare", f"Could not reveal file:\n{e}")

    def _open_folder(self):
        try:
            open_folder(self.primary_path)
        except Exception as e:
            QtWidgets.QMessageBox.warning(
                self, "Mae2Flare", f"Could not open folder:\n{e}")


# ---------------------------------------------------------------------------
# Entry point Maestro calls
# ---------------------------------------------------------------------------

def panel():
    settings = load_settings()

    try:
        entries = get_selected_structures()
    except Exception as e:
        QtWidgets.QMessageBox.critical(None, "Mae2Flare", str(e))
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

    flare_dir = Path(settings["flare_dir"])
    if not flare_dir.exists():
        QtWidgets.QMessageBox.critical(
            None, "Mae2Flare",
            f"Flare folder does not exist:\n{flare_dir}")
        return

    groups = group_by_protein(
        entries, atom_threshold=settings["protein_atom_threshold"])
    if not groups:
        QtWidgets.QMessageBox.warning(None, "Mae2Flare", "Nothing to export.")
        return

    summary = "\n".join(
        f"  group {i}: protein={'yes' if g['protein'] else 'no'}, "
        f"ligands={len(g['ligands'])}  ({g['name']})"
        for i, g in enumerate(groups, 1)
    )
    print(f"[Mae2Flare] {len(groups)} group(s):\n{summary}")

    fmt = settings.get("output_format", "split")
    note_to_user = ""
    final_paths = []

    if fmt == "flrp":
        out_flrp = out_dir / f"{settings['output_name']}.flrp"
        with tempfile.TemporaryDirectory(prefix="mae2flare_") as td:
            tmpdir = Path(td)
            proteins_tmp, ligands_tmp, group_names = write_split_for_pyflare(
                groups, tmpdir)
            ok, msg = run_pyflare(flare_dir, out_flrp,
                                  proteins_tmp, ligands_tmp, group_names)
        if ok:
            final_paths = [out_flrp]
        else:
            # Fall back to the split-native path
            print("[Mae2Flare] pyflare failed -- falling back to split files")
            print(msg)
            prot, lig = write_split_native(
                groups, out_dir, settings["output_name"])
            final_paths = [p for p in (prot, lig) if p is not None]
            note_to_user = (
                "Could not write a .flrp on this Flare build "
                "(internal SQLite project-cache bug). "
                "Saved as Maestro proteins + SDF ligands instead — "
                "Flare opens both natively and classifies them correctly."
            )
    else:
        # Default: Maestro proteins + SDF ligands
        prot, lig = write_split_native(
            groups, out_dir, settings["output_name"])
        final_paths = [p for p in (prot, lig) if p is not None]
        if not final_paths:
            QtWidgets.QMessageBox.warning(
                None, "Mae2Flare", "Nothing to export.")
            return

    opened = False
    if settings.get("open_in_flare", True):
        try:
            open_in_flare(flare_dir, final_paths)
            opened = True
        except Exception as e:
            QtWidgets.QMessageBox.warning(
                None, "Mae2Flare",
                f"Saved {final_paths},\nbut could not launch Flare:\n{e}")

    SuccessDialog(final_paths, opened_in_flare=opened, note=note_to_user).exec()


if __name__ == "__main__":
    panel()