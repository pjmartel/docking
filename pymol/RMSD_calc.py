from pymol import cmd
import numpy as np
import os
from collections import Counter, defaultdict

try:
    from scipy.optimize import linear_sum_assignment
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


def _get_element(atom):
    """Extract element symbol reliably from either PDB or PDBQT atoms."""
    elem = getattr(atom, 'elem', '') or getattr(atom, 'symbol', '') or ''
    elem = elem.strip().upper()
    ad_type_map = {'A': 'C', 'OA': 'O', 'NA': 'N', 'SA': 'S', 'HD': 'H',
                   'HS': 'H', 'HA': 'H', 'NS': 'N', 'OS': 'O'}
    if elem in ad_type_map:
        elem = ad_type_map[elem]
    real_2char = {'CL', 'BR', 'FE', 'ZN', 'MG', 'MN', 'CO', 'NI', 'CU',
                  'SE', 'SI', 'CA', 'NA', 'AL'}
    if len(elem) > 1 and elem not in real_2char:
        elem = elem[0]
    return elem


def _unique_obj_name(base):
    """Return base name that does not clash with any object already in the session."""
    existing = set(cmd.get_object_list())
    if base not in existing:
        return base
    k = 2
    while f"{base}_{k}" in existing:
        k += 1
    return f"{base}_{k}"


def _build_atom_map(dock_coords, dock_elems, crys_coords, crys_elems):
    """
    Build an optimal bijective mapping dock_idx -> crys_idx constrained by
    element type. Uses Hungarian assignment (scipy) or greedy fallback.
    Returns a dict, or None if per-element atom counts do not match.
    """
    crys_by_elem = defaultdict(list)
    for j, e in enumerate(crys_elems):
        crys_by_elem[e].append(j)

    dock_by_elem = defaultdict(list)
    for i, e in enumerate(dock_elems):
        dock_by_elem[e].append(i)

    atom_map = {}
    for elem, di in dock_by_elem.items():
        ci = crys_by_elem.get(elem, [])
        if len(di) != len(ci):
            return None         # element count mismatch - caller handles error
        if len(di) == 1:
            atom_map[di[0]] = ci[0]
            continue

        cost = np.array([[np.sum((dock_coords[i] - crys_coords[j]) ** 2)
                          for j in ci] for i in di])

        if _HAS_SCIPY:
            row_ind, col_ind = linear_sum_assignment(cost)
            for a, b in zip(row_ind, col_ind):
                atom_map[di[a]] = ci[b]
        else:
            pairs = sorted((cost[a, b], a, b)
                           for a in range(len(di)) for b in range(len(ci)))
            used_a, used_b = set(), set()
            for _, a, b in pairs:
                if a not in used_a and b not in used_b:
                    atom_map[di[a]] = ci[b]
                    used_a.add(a)
                    used_b.add(b)

    return atom_map


def vina_rmsd(pdbqt_file, ligand_resn, pdb_file=None, chain="A", show=1, title_states=1):
    """
    Calculate heavy-atom RMSD of Vina docking poses vs a crystal ligand.

    USAGE:
        vina_rmsd pdbqt_file, ligand_resn [, pdb_file [, chain [, show [, title_states]]]]

    ARGUMENTS:
        pdbqt_file  = Vina output PDBQT file with multiple poses
        ligand_resn = residue name of the ligand (e.g. STI, P17)
        pdb_file    = PDB file to load as crystal reference (optional).
                      If omitted, ligand_resn must already be present in the
                      current session (e.g. as part of an already-loaded receptor).
        chain       = chain containing the crystal ligand (default: A)
        show        = 1 to update visualization, 0 to suppress (default: 1)
        title_states= 1 to set per-state titles with RMSD values (default: 1)

    NOTES:
        - RMSD is a pure placement RMSD — no fitting is performed. Vina output
          must share the same coordinate frame as the crystal structure.
        - Atom correspondence is built once from pose 1 using Hungarian assignment
          per element group (falls back to greedy if scipy is not installed).
        - All poses are validated for atom count and element consistency before
          RMSD is computed; any inconsistency aborts with an error.
        - Visualization is additive: nothing is hidden or removed from the session.

    EXAMPLES:
        # Reference from a PDB file to load:
        vina_rmsd imatinib_docked.pdbqt, STI, 1IEP.pdb, B

        # Reference from a ligand already in the session (e.g. inside the receptor):
        vina_rmsd imatinib_docked.pdbqt, STI, chain=B
        vina_rmsd imatinib_docked.pdbqt, STI
        vina_rmsd imatinib_docked.pdbqt, STI, chain=B, title_states=1
    """
    show = bool(int(show))
    title_states = bool(int(title_states))

    print("=== RMSD Calculation: Vina Poses vs Crystal Structure ===\n")

    # --- Crystal ligand: load from file or locate in current session ---
    pdb_obj = None
    if pdb_file:
        pdb_obj = _unique_obj_name(
            os.path.splitext(os.path.basename(pdb_file))[0]
        )
        cmd.load(pdb_file, pdb_obj)
        crystal_sel = (f"{pdb_obj} and chain {chain} "
                       f"and resn {ligand_resn} and not hydrogen")
    else:
        # Locate the ligand among already-loaded objects
        crystal_sel = f"resn {ligand_resn} and chain {chain} and not hydrogen"
        if cmd.count_atoms(crystal_sel) == 0:
            # Relax chain constraint and warn
            crystal_sel_noc = f"resn {ligand_resn} and not hydrogen"
            if cmd.count_atoms(crystal_sel_noc) > 0:
                crystal_sel = crystal_sel_noc
                print(f"Warning: no atoms found for chain '{chain}'; "
                      f"using all chains for resn {ligand_resn}.")
            else:
                print(f"ERROR: resn '{ligand_resn}' not found in the current "
                      f"session. Load a PDB file or provide pdb_file=.")
                return

    cmd.select("crystal_lig", crystal_sel)
    n_crystal = cmd.count_atoms("crystal_lig")
    if n_crystal == 0:
        print(f"ERROR: No atoms found for '{crystal_sel}'. Check resn/chain.")
        return
    print(f"Crystal ligand ({ligand_resn}, chain {chain}): {n_crystal} heavy atoms")

    crystal_model = cmd.get_model("crystal_lig")
    crystal_coords = np.array([a.coord for a in crystal_model.atom])
    crystal_elems = [_get_element(a) for a in crystal_model.atom]

    # --- Load Vina poses ---
    vina_obj = _unique_obj_name("vina_poses")
    cmd.load(pdbqt_file, vina_obj)
    n_poses = cmd.count_states(vina_obj)
    if n_poses == 0:
        print(f"ERROR: No poses loaded from '{pdbqt_file}'.")
        return
    print(f"Loaded {n_poses} docking poses from '{pdbqt_file}'\n")

    # --- Build atom map from pose 1 ---
    pose1_model = cmd.get_model(f"{vina_obj} and not hydrogen", state=1)
    pose1_coords = np.array([a.coord for a in pose1_model.atom])
    pose1_elems = [_get_element(a) for a in pose1_model.atom]

    n_dock = len(pose1_elems)
    n_crys = len(crystal_elems)

    dock_counts = Counter(pose1_elems)
    crys_counts = Counter(crystal_elems)
    print(f"Docked elements:  {dict(sorted(dock_counts.items()))}")
    print(f"Crystal elements: {dict(sorted(crys_counts.items()))}")

    if n_dock != n_crys:
        print(f"ERROR: Heavy atom count mismatch — docked={n_dock}, crystal={n_crys}.")
        return
    if dock_counts != crys_counts:
        print("ERROR: Element distributions differ between docked and crystal ligand. Aborting.")
        return

    atom_map = _build_atom_map(pose1_coords, pose1_elems, crystal_coords, crystal_elems)
    if atom_map is None or len(atom_map) != n_dock:
        mapped = len(atom_map) if atom_map else 0
        print(f"ERROR: Could not map all atoms ({mapped}/{n_dock}). "
              f"Check element types.")
        return

    method = "Hungarian (optimal)" if _HAS_SCIPY else "greedy (install scipy for optimal)"
    print(f"Atom mapping: {method}\n")

    crystal_reordered = np.array([crystal_coords[atom_map[i]] for i in range(n_dock)])

    # --- Validate all poses for consistency before computing RMSD ---
    print(f"Validating {n_poses} poses for consistency...")
    for s in range(1, n_poses + 1):
        m = cmd.get_model(f"{vina_obj} and not hydrogen", state=s)
        s_elems = [_get_element(a) for a in m.atom]
        if len(s_elems) != n_dock:
            print(f"ERROR: Pose {s} has {len(s_elems)} heavy atoms, "
                  f"expected {n_dock}. Aborting.")
            return
        if Counter(s_elems) != dock_counts:
            print(f"ERROR: Pose {s} has a different element distribution. Aborting.")
            return
    print("All poses consistent.\n")

    # --- Compute RMSD ---
    print("=== Heavy-Atom RMSD Results ===")
    print(f"  {'Pose':>4}   {'RMSD (Å)':>10}")
    print("  " + "-" * 20)
    rmsd_by_state = {}
    for s in range(1, n_poses + 1):
        model = cmd.get_model(f"{vina_obj} and not hydrogen", state=s)
        coords = np.array([a.coord for a in model.atom])
        if coords.shape != crystal_reordered.shape:
            print(f"  {s:4d}   {'shape mismatch':>10}")
            continue
        rmsd = np.sqrt(np.mean(np.sum((coords - crystal_reordered) ** 2, axis=1)))
        rmsd_by_state[s] = rmsd
        print(f"  {s:4d}   {rmsd:10.3f}")

    if title_states:
        for s in range(1, n_poses + 1):
            if s in rmsd_by_state:
                cmd.set_title(vina_obj, s, f"Pose {s} | RMSD {rmsd_by_state[s]:.3f} Å")
            else:
                cmd.set_title(vina_obj, s, f"Pose {s} | RMSD n/a")
        print("\nState titles updated with RMSD values.")

    # --- Visualization (additive — nothing is hidden or removed) ---
    if show:
        if pdb_obj:
            cmd.show("cartoon", f"{pdb_obj} and polymer")
        cmd.show("sticks", "crystal_lig")
        cmd.color("forest", "crystal_lig")
        cmd.set("all_states", 0, vina_obj)
        cmd.show("sticks", vina_obj)
        cmd.set("stick_transparency", 0.0, vina_obj)
        cmd.color("gray70", vina_obj)
        cmd.zoom("crystal_lig", 5)
        cmd.center("crystal_lig")
        print("\nGreen sticks = Crystal ligand (reference)")
        print("Gray sticks = Docked pose (current state)")


cmd.extend("vina_rmsd", vina_rmsd)