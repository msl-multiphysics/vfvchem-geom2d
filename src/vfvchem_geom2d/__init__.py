"""2D geometry helpers built on gmsh for VFVChem."""

import math
import re
from collections import Counter

import gmsh

__version__ = "0.0.1"

__all__ = ["Geom2D"]

# vfvm2d insets boundary nodes from each vertex by OFFSET_RATIO * h_min
# so every segment must be longer than 2 * OFFSET_RATIO * h_min
OFFSET_RATIO = 0.1


class Geom2D:
    def __init__(self, model_name: str, z: float = 0.0):
        # store variables
        self.model_name = model_name
        self.z = z

        # states
        self.final_geom = False
        self.dom_groups = {}  # name -> [dom id, ...]
        self.bnd_groups = {}  # name -> [curve id, ...]
        self.itf_groups = {}  # name -> [curve id, ...]
        self._groups_dirty = True

        # topology after finalize_geometry
        # curve ids list the boundary curves first, then the interface curves
        self.scale = 1.0  # bounding box diagonal
        self.num_dom = 0
        self.num_bnd = 0
        self.num_curve = 0
        self.vert_pos = {}  # vertex key -> (x, y); key is a gmsh point tag or (curve id, i)
        self.curve_vids = []  # curve id -> [vertex key, ...] along the curve
        self.curve_doms = []  # curve id -> [dom id, ...]
        self.curve_fwd = []  # curve id -> direction of the curve in the loop of its lowest dom id
        self.dom_loops = []  # dom id -> [[(curve id, forward), ...], ...] with the domain on the left

        # initialize gmsh
        gmsh.initialize()
        gmsh.model.add(self.model_name)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        # finalize gmsh if the geometry was not finalized (e.g., after an error)
        if gmsh.isInitialized():
            gmsh.finalize()

    def add_rectangle(self, x_min: float, y_min: float, x_max: float, y_max: float):
        # error if geometry is finalized
        if self.final_geom:
            raise ValueError("Cannot edit geometry after finalizing geometry.")

        # error if x_min >= x_max or y_min >= y_max
        if x_min >= x_max:
            raise ValueError("x_min must be less than x_max.")
        if y_min >= y_max:
            raise ValueError("y_min must be less than y_max.")

        # add rectangle to model
        gmsh.model.occ.addRectangle(x_min, y_min, self.z, x_max - x_min, y_max - y_min)
        gmsh.model.occ.synchronize()

        # return (dim, tag) of added surface
        surf_dt = gmsh.model.getEntities(2)[-1]
        return surf_dt

    def add_ellipse(self, x_cen: float, y_cen: float, x_rad: float, y_rad: float):
        # error if geometry is finalized
        if self.final_geom:
            raise ValueError("Cannot edit geometry after finalizing geometry.")

        # gmsh requires the major radius along x
        if x_rad <= 0 or y_rad <= 0:
            raise ValueError("x_rad and y_rad must be positive.")
        if x_rad < y_rad:
            raise ValueError("x_rad must be greater than or equal to y_rad.")

        # add ellipse to model
        gmsh.model.occ.addDisk(x_cen, y_cen, self.z, x_rad, y_rad)
        gmsh.model.occ.synchronize()

        # return (dim, tag) of added surface
        surf_dt = gmsh.model.getEntities(2)[-1]
        return surf_dt

    def add_polygon(self, points: list):
        # error if geometry is finalized
        if self.final_geom:
            raise ValueError("Cannot edit geometry after finalizing geometry.")

        if len(points) < 3:
            raise ValueError("Polygon requires at least 3 points.")

        # drop duplicate closing vertex if first == last
        pts_xy = list(points)
        if pts_xy[0][0] == pts_xy[-1][0] and pts_xy[0][1] == pts_xy[-1][1]:
            pts_xy = pts_xy[:-1]
            if len(pts_xy) < 3:
                raise ValueError("Polygon requires at least 3 distinct points.")

        # build closed wire from consecutive points
        pts = [gmsh.model.occ.addPoint(x, y, self.z) for x, y in pts_xy]
        lines = [
            gmsh.model.occ.addLine(pts[i], pts[(i + 1) % len(pts)])
            for i in range(len(pts))
        ]
        loop = gmsh.model.occ.addCurveLoop(lines)
        gmsh.model.occ.addPlaneSurface([loop])
        gmsh.model.occ.synchronize()

        # return (dim, tag) of added surface
        surf_dt = gmsh.model.getEntities(2)[-1]
        return surf_dt

    def add_union(self, surf1: tuple, surf2: tuple, del1: bool = True, del2: bool = True):
        # error if geometry is finalized
        if self.final_geom:
            raise ValueError("Cannot edit geometry after finalizing geometry.")

        # union two surfaces
        # assumes that the inputs overlap
        out, _ = gmsh.model.occ.fuse([surf1], [surf2], -1, del1, del2)
        gmsh.model.occ.synchronize()

        # error if inputs do not overlap
        if len(out) != 1:
            raise ValueError("Inputs do not overlap.")

        # return (dim, tag) of unioned surface
        return out[0]

    def add_intersection(self, surf1: tuple, surf2: tuple, del1: bool = True, del2: bool = True):
        # error if geometry is finalized
        if self.final_geom:
            raise ValueError("Cannot edit geometry after finalizing geometry.")

        # intersect two surfaces
        out, _ = gmsh.model.occ.intersect([surf1], [surf2], -1, del1, del2)
        gmsh.model.occ.synchronize()

        # error if inputs do not overlap
        if len(out) != 1:
            raise ValueError("Inputs do not overlap")

        # return (dim, tag) of intersected surface
        return out[0]

    def add_difference(self, surf1: tuple, surf2: tuple, del1: bool = True, del2: bool = True):
        # error if geometry is finalized
        if self.final_geom:
            raise ValueError("Cannot edit geometry after finalizing geometry.")

        # difference of two surfaces (surf1 minus surf2)
        out, _ = gmsh.model.occ.cut([surf1], [surf2], -1, del1, del2)
        gmsh.model.occ.synchronize()

        # error if result is not a single surface
        if len(out) != 1:
            raise ValueError("Difference did not produce a single surface")

        # return (dim, tag) of remaining surface
        return out[0]

    def finalize_geometry(self, size_curve: float = None):
        """Split the geometry into domains and curves, then close the gmsh session.

        Curved edges are approximated by straight segments no longer than size_curve
        (default: 1% of the bounding box diagonal). Straight edges are kept as one segment.
        """
        # error if geometry is finalized
        if self.final_geom:
            raise ValueError("Geometry is already finalized.")
        if size_curve is not None and size_curve <= 0:
            raise ValueError("size_curve must be positive.")

        # get list of surfaces
        surf_list = gmsh.model.getEntities(2)
        if not surf_list:
            raise ValueError("No surfaces in geometry.")

        # fragment the surfaces so overlapping surfaces become conformal domains
        if len(surf_list) > 1:
            gmsh.model.occ.fragment([surf_list[0]], surf_list[1:])
            gmsh.model.occ.synchronize()

        # model scale
        x_min, y_min, _, x_max, y_max, _ = gmsh.model.getBoundingBox(-1, -1)
        self.scale = math.hypot(x_max - x_min, y_max - y_min)
        if size_curve is None:
            size_curve = 0.01 * self.scale

        # curves bounding each domain
        surf_curves = []
        curve_doms = {}  # curve tag -> [dom id, ...]
        for dom, (_, surf) in enumerate(gmsh.model.getEntities(2)):
            curves = sorted({abs(c) for _, c in gmsh.model.getBoundary([(2, surf)], oriented=False)})
            surf_curves.append(curves)
            for c in curves:
                curve_doms.setdefault(c, []).append(dom)

        # boundary curves touch one domain and interface curves touch two
        bnd_list = sorted(c for c, doms in curve_doms.items() if len(doms) == 1)
        itf_list = sorted(c for c, doms in curve_doms.items() if len(doms) == 2)
        if len(bnd_list) + len(itf_list) != len(curve_doms):
            raise ValueError("A curve is shared by more than two domains.")
        curve_list = bnd_list + itf_list
        curve_id = {c: i for i, c in enumerate(curve_list)}
        self.num_dom = len(surf_curves)
        self.num_bnd = len(bnd_list)
        self.num_curve = len(curve_list)
        self.curve_doms = [curve_doms[c] for c in curve_list]

        # approximate each curve by straight segments
        # points shared by curves keep their gmsh tag so the curves share vertices
        for cid, c in enumerate(curve_list):
            (t0,), (t1,) = gmsh.model.getParametrizationBounds(1, c)
            x0, y0, _ = gmsh.model.getValue(1, c, [t0])
            x1, y1, _ = gmsh.model.getValue(1, c, [t1])
            pts = gmsh.model.getAdjacencies(1, c)[1]
            for p in pts:
                x, y, _ = gmsh.model.getValue(0, p, [])
                self.vert_pos[p] = (float(x), float(y))
            if len(pts) > 0:
                start = min(pts, key=lambda p: math.dist(self.vert_pos[p], (x0, y0)))
                end = min(pts, key=lambda p: math.dist(self.vert_pos[p], (x1, y1)))
            else:
                start = end = (cid, 0)
                self.vert_pos[start] = (float(x0), float(y0))

            # straight curves are one segment; curved ones need at least two (three if closed)
            length = gmsh.model.occ.getMass(1, c)
            if start != end and length - math.dist(self.vert_pos[start], self.vert_pos[end]) < 1e-9 * self.scale:
                num_seg = 1
            else:
                num_seg = max(math.ceil(length / size_curve), 2 if start != end else 3)
            ts = [t0 + (t1 - t0) * i / num_seg for i in range(1, num_seg)]
            xyz = gmsh.model.getValue(1, c, ts) if ts else []
            vids = [start]
            for i in range(1, num_seg):
                self.vert_pos[(cid, i)] = (float(xyz[3 * i - 3]), float(xyz[3 * i - 2]))
                vids.append((cid, i))
            vids.append(end)
            self.curve_vids.append(vids)

        # join the curves of each domain into loops
        # the outer loop is ccw and holes are cw so the domain is on the left
        for curves in surf_curves:
            cids = [curve_id[c] for c in curves]
            loops = self._chain([(cid, True) for cid in cids] + [(cid, False) for cid in cids])
            areas = []
            for lp in loops:
                vids = self._join(lp)
                if vids[0] != vids[-1]:
                    raise ValueError("Domain boundary does not form closed loops.")
                areas.append(self._area(vids))
            outer = max(range(len(loops)), key=lambda i: abs(areas[i]))
            for i, lp in enumerate(loops):
                if (areas[i] > 0) != (i == outer):
                    loops[i] = [(cid, not fwd) for cid, fwd in reversed(lp)]
            self.dom_loops.append(loops)

        # direction of each curve in the loop of its lowest dom id
        self.curve_fwd = [None] * self.num_curve
        for loops in reversed(self.dom_loops):
            for lp in loops:
                for cid, fwd in lp:
                    self.curve_fwd[cid] = fwd

        # the geometry is cached, so gmsh is no longer needed
        self.final_geom = True
        gmsh.finalize()

    def _vids(self, cid: int, fwd: bool):
        # vertex keys of a curve in the given direction
        return self.curve_vids[cid] if fwd else self.curve_vids[cid][::-1]

    def _chain(self, dcurves: list):
        # join (curve id, forward) pairs end to end into chains
        # a curve given in both directions is used once
        # chains start where no other curve ends, so open polylines are not split
        ends = {self._vids(*dc)[-1] for dc in dcurves}
        todo = sorted(dcurves, key=lambda dc: self._vids(*dc)[0] in ends)
        used = set()
        chains = []
        for dc in todo:
            if dc[0] in used:
                continue
            chain = [dc]
            used.add(dc[0])
            while self._vids(*chain[-1])[-1] != self._vids(*chain[0])[0]:
                end = self._vids(*chain[-1])[-1]
                nxt = next((dc for dc in todo if dc[0] not in used and self._vids(*dc)[0] == end), None)
                if nxt is None:
                    break
                chain.append(nxt)
                used.add(nxt[0])
            chains.append(chain)
        return chains

    def _join(self, chain: list):
        # vertex keys of a chain without repeating the shared end points
        vids = list(self._vids(*chain[0]))
        for dc in chain[1:]:
            vids.extend(self._vids(*dc)[1:])
        return vids

    def _area(self, vids: list):
        # signed area of a closed polyline; positive if ccw
        pos = [self.vert_pos[v] for v in vids]
        return 0.5 * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pos, pos[1:]))

    def _curve_label(self, cid: int):
        return f"bnd_{cid}" if cid < self.num_bnd else f"itf_{cid}"

    def _add_group(self, groups: dict, kind: str, name: str, ids: list, valid_ids: range):
        if not self.final_geom:
            raise ValueError(f"Cannot call group_{kind} before finalizing geometry.")

        # vfvm2d reads one token per group name and treats '#' as a comment
        if not re.fullmatch(rf"{kind}_[^\s#]+", name):
            raise ValueError(
                f"group_{kind}: name '{name}' must begin with '{kind}_' "
                f"followed by at least one character, with no whitespace or '#'."
            )
        if len(ids) < 1:
            raise ValueError(f"group_{kind} requires at least one id.")
        if len(set(ids)) != len(ids):
            raise ValueError(f"group_{kind} ids must be unique.")

        # each id belongs to at most one group of each kind
        for i in ids:
            if i not in valid_ids:
                raise ValueError(f"group_{kind}: id {i} is not a valid {kind} id.")
            for other, other_ids in groups.items():
                if i in other_ids:
                    raise ValueError(f"group_{kind}: id {i} is already in group '{other}'.")

        # repeated names continue the same group
        groups.setdefault(name, []).extend(ids)
        self._groups_dirty = True

    def group_dom(self, name: str, ids: list):
        """Join domains into a named group. ids are the numbers of the dom_# labels from show_geometry()."""
        self._add_group(self.dom_groups, "dom", name, ids, range(self.num_dom))

    def group_bnd(self, name: str, ids: list):
        """Join curves into a named boundary group. ids are the numbers of the bnd_# or itf_# labels from show_geometry().

        An interface curve may also be in a boundary group (e.g., a no-slip wall between a fluid and a solid).
        """
        self._add_group(self.bnd_groups, "bnd", name, ids, range(self.num_curve))

    def group_itf(self, name: str, ids: list):
        """Join curves into a named interface group. ids are the numbers of the itf_# labels from show_geometry()."""
        self._add_group(self.itf_groups, "itf", name, ids, range(self.num_bnd, self.num_curve))

    def update_group(self):
        """Check that the named groups describe a valid vfvm2d geometry.

        Every domain must be in a domain group, and every curve on the edge of a domain
        group must be in a boundary or interface group. Curves between two domains of
        the same domain group lie inside that group and must not be grouped.
        """
        if not self.final_geom:
            raise ValueError("Cannot update groups before finalizing geometry.")

        # every domain is in a domain group
        dom_name = {d: name for name, ids in self.dom_groups.items() for d in ids}
        for d in range(self.num_dom):
            if d not in dom_name:
                raise ValueError(f"update_group: dom_{d} is not in a domain group.")

        # domain groups on the sides of each curve
        sides = [{dom_name[d] for d in doms} for doms in self.curve_doms]
        line_name = {}  # curve id -> [group name, ...]
        for name, ids in [*self.bnd_groups.items(), *self.itf_groups.items()]:
            for cid in ids:
                line_name.setdefault(cid, []).append(name)

        # curves inside a domain group are not grouped; all other curves are
        for cid, doms in enumerate(self.curve_doms):
            inside = len(doms) == 2 and len(sides[cid]) == 1
            if inside and cid in line_name:
                raise ValueError(
                    f"update_group: {self._curve_label(cid)} lies inside '{dom_name[doms[0]]}' "
                    f"and cannot be in '{line_name[cid][0]}'."
                )
            if not inside and cid not in line_name:
                raise ValueError(
                    f"update_group: {self._curve_label(cid)} is not in a boundary or interface group."
                )

        # vfvm2d takes boundary normals from one domain, so a boundary group lies on one domain group
        for name, ids in self.bnd_groups.items():
            if not set.intersection(*(sides[cid] for cid in ids)):
                raise ValueError(f"update_group: '{name}' does not lie on a single domain group.")

        # an interface group separates one pair of domain groups
        for name, ids in self.itf_groups.items():
            if len({frozenset(sides[cid]) for cid in ids}) != 1:
                raise ValueError(f"update_group: '{name}' separates more than one pair of domain groups.")

        self._groups_dirty = False

    def _ensure_groups_synced(self):
        if self._groups_dirty:
            self.update_group()

    def show_geometry(
        self,
        *,
        show_domains: bool = True,
        show_boundaries: bool = True,
        show_interfaces: bool = True,
        show_labels: bool = True,
    ):
        """Plot the domains, boundaries, and interfaces with their group names or ids."""
        if not self.final_geom:
            raise ValueError("Cannot show geometry before finalizing geometry.")
        try:
            from .viz import render_geometry
        except ImportError as exc:
            raise ImportError(
                "show_geometry() requires pyvista, matplotlib, and numpy. "
                "Install with: pip install pyvista matplotlib numpy"
            ) from exc

        # group names, or ids for entities that are not grouped yet
        dom_name = {d: name for name, ids in self.dom_groups.items() for d in ids}
        line_name = {}  # curve id -> [group name, ...]
        for name, ids in [*self.itf_groups.items(), *self.bnd_groups.items()]:
            for cid in ids:
                line_name.setdefault(cid, []).append(name)

        pos = self.vert_pos
        domains = [
            (dom_name.get(d, f"dom_{d}"), [[pos[v] for v in self._join(lp)] for lp in loops])
            for d, loops in enumerate(self.dom_loops)
        ]
        boundaries = []
        interfaces = []
        for cid, doms in enumerate(self.curve_doms):
            # skip curves inside a domain group
            if len(doms) == 2 and all(d in dom_name for d in doms) and dom_name[doms[0]] == dom_name[doms[1]]:
                continue
            label = " / ".join(line_name.get(cid, [self._curve_label(cid)]))
            region = (label, [[pos[v] for v in self.curve_vids[cid]]])
            if cid < self.num_bnd:
                boundaries.append(region)
            else:
                interfaces.append(region)

        render_geometry(
            self.model_name,
            self.z,
            domains,
            boundaries,
            interfaces,
            show_domains=show_domains,
            show_boundaries=show_boundaries,
            show_interfaces=show_interfaces,
            show_labels=show_labels,
        )

    def _fmt(self, x: float):
        # trim floating-point noise from the cad kernel
        x = float(f"{x:.12g}") if abs(x) > 1e-12 * self.scale else 0.0
        return repr(x)

    def export_geometry(self, file_path: str, h_min: float = None):
        """Write the groups to an obj file for vfvm2d.

        h_min is the optional boundary node spacing used in vfvm2d (Variables::new_obj).
        If given, every segment is checked to be longer than 2 * 0.1 * h_min.
        """
        # error if geometry is not finalized
        if not self.final_geom:
            raise ValueError("Cannot export geometry before finalizing geometry.")
        if h_min is not None and h_min <= 0:
            raise ValueError("h_min must be positive.")
        self._ensure_groups_synced()

        # polylines of the 1d groups follow the loop direction of their lowest dom id
        records = []  # (group name, [[vertex key, ...], ...])
        for name, ids in [*self.bnd_groups.items(), *self.itf_groups.items()]:
            chains = self._chain([(cid, self.curve_fwd[cid]) for cid in ids])
            records.append((name, [self._join(ch) for ch in chains]))
        num_1d = len(records)

        # loops of the 2d groups; curves shared by two domains of a group are inside it
        for name, ids in self.dom_groups.items():
            dcurves = [dc for d in ids for lp in self.dom_loops[d] for dc in lp]
            count = Counter(cid for cid, _ in dcurves)
            chains = self._chain([dc for dc in dcurves if count[dc[0]] == 1])
            loops = [self._join(ch) for ch in chains]
            if any(lp[0] != lp[-1] for lp in loops):
                raise ValueError(f"export_geometry: '{name}' does not form closed loops.")
            if sum(self._area(lp) for lp in loops) <= 0:
                raise ValueError(f"export_geometry: '{name}' loops are not ccw around the domain.")
            records.append((name, loops))

        # segment lengths
        min_len = 2 * OFFSET_RATIO * h_min if h_min is not None else 0.0
        for name, lines in records:
            for line in lines:
                for a, b in zip(line, line[1:]):
                    seg_len = math.dist(self.vert_pos[a], self.vert_pos[b])
                    if seg_len <= min_len:
                        raise ValueError(
                            f"export_geometry: segment in '{name}' has length {seg_len:.6g}; "
                            f"vfvm2d requires a length greater than {min_len:.6g}."
                        )

        # number the vertices (1-based) in order of first use
        obj_vid = {}
        for _, lines in records:
            for line in lines:
                for v in line:
                    obj_vid.setdefault(v, len(obj_vid) + 1)

        # write obj file
        with open(file_path, "w") as f:
            f.write("# vertices\n")
            for v in obj_vid:
                x, y = self.vert_pos[v]
                if not (math.isfinite(x) and math.isfinite(y)):
                    raise ValueError(f"export_geometry: vertex ({x}, {y}) is not finite.")
                f.write(f"v {self._fmt(x)} {self._fmt(y)} {self._fmt(self.z)}\n")
            for i, (name, lines) in enumerate(records):
                if i == 0:
                    f.write("\n# 1d groups\n")
                if i == num_1d:
                    f.write("\n# 2d groups\n")
                f.write(f"g {name}\n")
                for line in lines:
                    f.write("l " + " ".join(str(obj_vid[v]) for v in line) + "\n")
