# PyMOL Scripts

This folder contains PyMOL helper scripts for setting up and analyzing docking runs.

## Scripts

| Script | Description |
|--------|-------------|
| [make_box.py](make_box.py) | Draws a docking box in PyMOL as a pseudoatom wireframe. Supports explicit dimensions/center, bounding-box derivation from any named selection (`around=`), optional padding (`pad=`), and `quiet=1` to print only Vina-ready center/size lines. |
| [RMSD_calc.py](RMSD_calc.py) | Computes per-pose heavy-atom placement RMSD (no fitting) of Vina PDBQT output vs a crystal ligand. The crystal reference can be a PDB file to load or a residue already present in the session. Uses Hungarian assignment for optimal atom mapping (greedy fallback if scipy absent). Validates all poses for consistency before computing. Visualization is additive (`show=0` to suppress). |
