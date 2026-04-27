from pymol import cmd
import numpy as np
from collections import Counter


def _get_element(atom):
    """Extract element symbol reliably from either PDB or PDBQT atoms."""
    elem = getattr(atom, 'elem', '') or atom.symbol or ''
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


def vina_rmsd(pdbqt_file, ligand_resn, pdb_file, chain="A"):
    """
    Calculate heavy-atom RMSD of Vina docking poses vs a crystal ligand.

    Usage in PyMOL:
        vina_rmsd imatinib_docked_e20.pdbqt, STI, 1IEP.pdb, B
        vina_rmsd P17_docked_e20_1m52.pdbqt, P17, 1M52.pdb, B
    """
    pdb_obj = pdb_file.split('.')[0]

    print("=== RMSD Calculation: Vina Poses vs Crystal Structure ===\n")

    # Load crystal structure
    cmd.load(pdb_file, pdb_obj)
    cmd.remove(f"{pdb_obj} and not (polymer or resn {ligand_resn})")
    crystal_sel = f"{pdb_obj} and chain {chain} and resn {ligand_resn} and not hydrogen"
    cmd.select("crystal_lig", crystal_sel)
    n_crystal = cmd.count_atoms('crystal_lig')
    if n_crystal == 0:
        print(f"ERROR: No atoms found for '{crystal_sel}'. Check resn/chain.")
        return
    print(f"Crystal ligand ({ligand_resn} chain {chain}): {n_crystal} heavy atoms")

    # Get crystal coordinates and elements
    crystal_model = cmd.get_model(crystal_sel)
    crystal_coords = np.array([a.coord for a in crystal_model.atom])
    crystal_elems = [_get_element(a) for a in crystal_model.atom]

    # Load Vina poses
    cmd.load(pdbqt_file, "vina_poses")
    n_poses = cmd.count_states("vina_poses")
    print(f"Loaded {n_poses} docking poses\n")

    # Build atom mapping using pose 1 (best Vina result, closest to crystal).
    # Vina PDBQT uses generic atom names so PyMOL matchers fail;
    # we map by element + nearest-neighbor distance instead.
    pose1_model = cmd.get_model("vina_poses and not hydrogen", state=1)
    pose1_coords = np.array([a.coord for a in pose1_model.atom])
    pose1_elems = [_get_element(a) for a in pose1_model.atom]

    n_dock = len(pose1_elems)
    n_crys = len(crystal_elems)

    dock_counts = Counter(pose1_elems)
    crys_counts = Counter(crystal_elems)
    print(f"Docked elements:  {dict(dock_counts)}")
    print(f"Crystal elements: {dict(crys_counts)}")

    if n_dock != n_crys:
        print(f"ERROR: Heavy atom count mismatch: docked={n_dock}, crystal={n_crys}")
        return
    if dock_counts != crys_counts:
        print("WARNING: element distributions differ!")

    # Globally-sorted greedy mapping constrained by element type
    pairs = []
    for i in range(n_dock):
        for j in range(n_crys):
            if pose1_elems[i] == crystal_elems[j]:
                d2 = np.sum((pose1_coords[i] - crystal_coords[j]) ** 2)
                pairs.append((d2, i, j))
    pairs.sort()

    used_dock = set()
    used_crys = set()
    atom_map = {}
    for d2, i, j in pairs:
        if i not in used_dock and j not in used_crys:
            atom_map[i] = j
            used_dock.add(i)
            used_crys.add(j)

    if len(atom_map) != n_dock:
        print(f"ERROR: Could not map all atoms: {len(atom_map)}/{n_dock}.")
        return

    # Reorder crystal coords to match docked atom order
    crystal_reordered = np.array([crystal_coords[atom_map[i]] for i in range(n_dock)])

    # Compute RMSD for each pose
    print("\n=== Heavy-Atom RMSD Results ===")
    print(f"{'Pose':>4}   {'RMSD (A)':>10}")
    print("-" * 28)

    for i in range(1, n_poses + 1):
        model = cmd.get_model("vina_poses and not hydrogen", state=i)
        coords = np.array([a.coord for a in model.atom])
        diffs = coords - crystal_reordered
        rmsd = np.sqrt(np.mean(np.sum(diffs ** 2, axis=1)))
        print(f"{i:4d}   {rmsd:10.3f}")

    # Visualization
    cmd.hide("everything")
    cmd.show("cartoon", f"{pdb_obj} and polymer")
    cmd.show("sticks", "crystal_lig")
    cmd.color("forest", "crystal_lig")

    cmd.set("all_states", 1)
    cmd.show("sticks", "vina_poses")
    cmd.set("stick_transparency", 0.6, "vina_poses")
    cmd.color("gray70", "vina_poses")

    cmd.zoom("crystal_lig", 5)
    cmd.center("crystal_lig")

    print("\nGreen sticks = Crystal ligand (reference)")
    print("Gray transparent sticks = All docked poses")


cmd.extend("vina_rmsd", vina_rmsd)