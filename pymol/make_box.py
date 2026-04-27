from pymol import cmd


def _bbox_from_selection(sel):
    """Return (cx, cy, cz, ext_x, ext_y, ext_z) for a PyMOL selection string."""
    model = cmd.get_model(sel)
    atoms = model.atom
    xs = [a.coord[0] for a in atoms]
    ys = [a.coord[1] for a in atoms]
    zs = [a.coord[2] for a in atoms]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin, zmax = min(zs), max(zs)
    cx = (xmin + xmax) / 2.0
    cy = (ymin + ymax) / 2.0
    cz = (zmin + zmax) / 2.0
    return cx, cy, cz, xmax - xmin, ymax - ymin, zmax - zmin


def _as_bool(value):
    """Parse PyMOL-style truthy/falsy values from bool/int/float/string."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if value is None:
        return False
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off", ""}:
        return False
    # Fallback: any non-empty string is treated as true
    return True


def box(dx=None, dy=None, dz=None, cx=None, cy=None, cz=None,
        pad=None, around=None, name="box", quiet=0):
    """
    Create a box using PyMOL pseudoatoms connected by bonds.

    USAGE:
        box pad=<p>, around=<sel>           # box around a named selection with padding
        box dx, dy, dz [, cx, cy, cz]       # explicit dimensions and optional center
        box dx, dy, dz [, pad=<p>]          # dimensions, center from current selection
        box ..., quiet=1                    # print only Vina center/size lines

    ARGUMENTS:
        dx, dy, dz  = box dimensions (Å); derived from selection + pad when around is used
        cx, cy, cz  = box center coordinates (optional; derived from selection if not given)
        pad         = padding (Å) added to each face of the bounding box
        around      = selection string to build the box around (e.g. "STI/", "chain A")
        name        = name for the box object (default: "box")
        quiet       = when true, output only the six Vina config lines

    NOTES:
        - When `around` is given, center and dimensions are always derived from the
          bounding box of that selection plus `pad` on every face.
        - When no `around` and no explicit center, the current PyMOL selection ("sele")
          is used if it exists; otherwise the origin (0, 0, 0) is used.
        - `pad` has no effect when explicit cx, cy, cz are given without `around`.

    EXAMPLES:
        box pad=6.0, around=STI/            # box around residue STI with 6 Å padding
        box pad=4, around=chain A           # box around chain A with 4 Å padding
        box pad=6, around=STI/, quiet=1     # Vina-ready output only
        box 20, 15, 10, 5, 5, 5            # 20x15x10 box centered at (5,5,5)
        box 10, 10, 10                      # 10x10x10 box centered on current selection
        box 10, 10, 10, name=mybox          # custom object name
    """

    quiet = _as_bool(quiet)

    # Convert optional numeric arguments
    if pad is not None:
        pad = float(pad)
    if dx is not None:
        dx = float(dx)
    if dy is not None:
        dy = float(dy)
    if dz is not None:
        dz = float(dz)
    if cx is not None:
        cx = float(cx)
    if cy is not None:
        cy = float(cy)
    if cz is not None:
        cz = float(cz)

    p = pad if pad is not None else 0.0

    # --- Determine center and dimensions ---

    if around is not None:
        # Named selection: always derive center and dimensions from its bounding box
        sel = str(around)
        if cmd.count_atoms(sel) == 0:
            print(f"Error: selection '{sel}' is empty or does not exist.")
            return
        scx, scy, scz, ext_x, ext_y, ext_z = _bbox_from_selection(sel)
        cx = scx if cx is None else cx
        cy = scy if cy is None else cy
        cz = scz if cz is None else cz
        dx = ext_x + 2 * p
        dy = ext_y + 2 * p
        dz = ext_z + 2 * p
        if not quiet:
            print(f"Box derived from selection '{sel}': center ({cx:.3f}, {cy:.3f}, {cz:.3f})")

    elif cx is None or cy is None or cz is None:
        # No explicit center: fall back to current selection "sele"
        if cmd.count_atoms('sele') > 0:
            scx, scy, scz, ext_x, ext_y, ext_z = _bbox_from_selection('sele')
            cx = scx if cx is None else cx
            cy = scy if cy is None else cy
            cz = scz if cz is None else cz
            if dx is None:
                dx = ext_x + 2 * p
            if dy is None:
                dy = ext_y + 2 * p
            if dz is None:
                dz = ext_z + 2 * p
            if not quiet:
                print(f"Box derived from current selection: center ({cx:.3f}, {cy:.3f}, {cz:.3f})")
        else:
            # No selection: fall back to origin
            cx = cx if cx is not None else 0.0
            cy = cy if cy is not None else 0.0
            cz = cz if cz is not None else 0.0
            if not quiet:
                print(f"No selection found. Using origin (0.000, 0.000, 0.000)")

    else:
        # Fully explicit center
        if pad is not None and not quiet:
            print("Warning: 'pad' is ignored when explicit center coordinates are provided.")
        if not quiet:
            print(f"Box centered at: ({cx:.3f}, {cy:.3f}, {cz:.3f})")

    # Final check: dimensions must be known by now
    if dx is None or dy is None or dz is None:
        print("Error: box dimensions (dx, dy, dz) could not be determined. "
              "Please provide them explicitly.")
        return

    # --- Build the box ---

    hx, hy, hz = dx / 2.0, dy / 2.0, dz / 2.0

    vertices = [
        (cx - hx, cy - hy, cz - hz),  # 0: ---
        (cx + hx, cy - hy, cz - hz),  # 1: +--
        (cx + hx, cy + hy, cz - hz),  # 2: ++-
        (cx - hx, cy + hy, cz - hz),  # 3: -+-
        (cx - hx, cy - hy, cz + hz),  # 4: --+
        (cx + hx, cy - hy, cz + hz),  # 5: +-+
        (cx + hx, cy + hy, cz + hz),  # 6: +++
        (cx - hx, cy + hy, cz + hz),  # 7: -++
    ]

    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),  # bottom face
        (4, 5), (5, 6), (6, 7), (7, 4),  # top face
        (0, 4), (1, 5), (2, 6), (3, 7),  # vertical edges
    ]

    cmd.delete(name)
    for i, (x, y, z) in enumerate(vertices):
        cmd.pseudoatom(name, pos=[x, y, z], name=f"V{i}")

    for v1, v2 in edges:
        cmd.bond(f"{name} and name V{v1}", f"{name} and name V{v2}")

    cmd.show("sticks", name)
    cmd.color("white", name)
    cmd.set("stick_radius", 0.1, name)

    if not quiet:
        print("\n" + "=" * 50)
        print("BOX INFORMATION")
        print("=" * 50)
        print(f"Center:      ({cx:.3f}, {cy:.3f}, {cz:.3f})")
        print(f"Dimensions:  {dx:.3f} x {dy:.3f} x {dz:.3f}")
        print(f"X extent:    {cx - hx:.3f} to {cx + hx:.3f}")
        print(f"Y extent:    {cy - hy:.3f} to {cy + hy:.3f}")
        print(f"Z extent:    {cz - hz:.3f} to {cz + hz:.3f}")
        print("=" * 50 + "\n")
        print("Autodock Vina config snippet:")

    print(f"center_x = {cx:.3f}")
    print(f"center_y = {cy:.3f}")
    print(f"center_z = {cz:.3f}")
    print(f"size_x = {dx:g}")
    print(f"size_y = {dy:g}")
    print(f"size_z = {dz:g}\n")


# Extend PyMOL command interface
cmd.extend("box", box)

print("Box command loaded successfully!")
print("Usage: box pad=<p>, around=<sel> [, quiet=1]  |  box dx, dy, dz [, cx, cy, cz] [, pad=<p>] [, name=<n>] [, quiet=1]")
