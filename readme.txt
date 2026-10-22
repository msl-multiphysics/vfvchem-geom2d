__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
2D Geometry Utility
Copyright (c) 2026 VFVChem Development Team

Overview
--------
This is a geometry utility for VFVChem, a vertex-centered finite volume method (FVM) solver for chemical engineering applications.
It provides a high-level wrapper for Gmsh geometry that writes `.obj` files that can be readily used with VFVChem.
VFVChem places its own nodes and points, so the `.obj` file holds the polygonal geometry and its named groups, not a mesh.

Note: The project is under active development. The public API and numerical models may change between releases.

Features
--------
- High-level 2D geometry primitives: rectangles, ellipses, and polygons
- Boolean operations on surfaces: union, intersection, and difference
- Automatic splitting of overlapping surfaces into domains, boundaries, and interfaces
- Approximation of curved edges by straight segments with a maximum segment length
- Named domain, boundary, and interface groups for VFVChem
- Interactive geometry visualization with PyVista, including ids and group names
- Export to `.obj` files for use with VFVChem

Requirements
------------
- Python 3.9 or newer
- Dependencies: gmsh, pyvista, matplotlib, and numpy (installed with the package)

Quick Start
-----------
Clone or download the project, then install it in editable mode from the project root:

    pip install -e .

Run an example from the `examples` folder, for instance:

    python examples/obj_step.py

This builds a channel with a backward-facing step, shows the geometry with domain and boundary ids and names, and writes `obj_step.obj`.

Workflow
--------
1. Define the geometry with `add_rectangle`, `add_ellipse`, `add_polygon`, `add_union`, `add_intersection`, and `add_difference`.
2. Call `finalize_geometry(size_curve)`. Overlapping surfaces are split into domains, and every curve is classified as a
   boundary (one domain) or an interface (two domains). Curved edges are split into straight segments no longer than
   `size_curve` (default: 1% of the bounding box diagonal). The Gmsh session is closed afterward.
3. Call `show_geometry()` to see the automatically generated labels (`dom_0`, `bnd_0`, `itf_8`, ...).
   The number in each label is the id passed to `group_dom`, `group_bnd`, and `group_itf`.
4. Name and join the ids with `group_dom`, `group_bnd`, and `group_itf`, then call `update_group()` to check the groups.
5. Call `export_geometry(file_name, h_min)` to write the `.obj` file.

Group Rules
-----------
- Group names must begin with `dom_`, `bnd_`, or `itf_` (case-sensitive), followed by at least one character.
  Names cannot contain whitespace or `#`.
- Every domain must be in one domain group. A domain group may hold several domains; curves between two domains of the
  same group lie inside it, are not written, and cannot be grouped.
- Every other curve must be in a boundary group, an interface group, or both. An interface curve may be in one interface
  group and one boundary group (e.g., a wall that is an interface for heat transfer and a no-slip boundary for momentum transfer).
- A boundary group must lie on the edge of a single domain group, since VFVChem takes its normals from that domain.
- An interface group must separate a single pair of domain groups.

OBJ Format
----------
The `.obj` file uses the subset of the Wavefront format that VFVChem reads:
- `v x y z` - vertex; only x and y are used
- `g name` - group; the `l` records that follow belong to it
- `l i j k ...` - polyline through 1-based vertex indices
- `#` - comment

Boundary and interface groups are open or closed polylines. Each domain group line is one closed loop with the domain on
its left: outer loops are counterclockwise and holes are clockwise. Curves shared by groups use the same vertex indices,
so VFVChem places the same nodes on them.

VFVChem insets its boundary nodes from each vertex by 0.1 * h_min, where h_min is the boundary node spacing in
`Variables::new_obj`. If `h_min` is given to `export_geometry`, every segment is checked to be longer than 0.2 * h_min.
Note that `size_curve` only sets how closely curves are followed; the node spacing is still set by h_min in VFVChem.

The groups are selected by name in VFVChem, for instance:

    let mut vars = Variables::new_obj("obj_fourcyl.obj", 0.001, 0.001, 0.001, 100, 42)?;
    let dom_c = vars.add_dom("dom_c")?;
    let dom_h = vars.add_dom("dom_h")?;
    let bnd_i = vars.add_bnd(dom_c, "bnd_i")?;
    let itf_h = vars.add_itf(dom_c, dom_h, "itf_h")?;

Project Structure
-----------------
- `src/vfvchem_geom2d/` — `Geom2D` API and visualization helpers
- `examples/` — example scripts and corresponding `.obj` outputs

Contributing
------------
Due to limited maintenance capacity, we are not currently accepting external contributions.
However, users are welcome to fork and modify the repository in accordance with the MIT License.

License
-------
VFVChem geom2d is available under the MIT License. See license.txt for details.
