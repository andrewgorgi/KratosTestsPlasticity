import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
from typing import List, Tuple

import KratosMultiphysics
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.path import Path

PARAM_SPACE_COLOR = 'gray'
SURROGATE_COLOR = 'darkred'
TRUE_BOUNDARY_COLOR = 'blue'
ESBM_GAUSS_COLOR = 'red'
ELEMENT_GAUSS_COLOR = 'deepskyblue'
ACTIVE_ELEMENT_COLOR = '#FFFACD'
GAP_FILL_COLOR = "#ffa2a2"  # light red
POLYNOMIAL_NURBS_COLOR = 'deepskyblue'

INTERPOLATION_LEVEL = 3



def trim_loop(xc, yc):
    """Remove duplicate closing point if present."""
    if len(xc) > 1 and np.isclose(xc[0], xc[-1]) and np.isclose(yc[0], yc[-1]):
        return xc[:-1], yc[:-1]
    return xc, yc


def extract_boundary_coords(mp):
    """Return closed boundary coordinates from a skin model part."""
    if len(mp.Conditions) == 0:
        return [], []

    x_coords, y_coords = [], []
    for i, cond in enumerate(mp.Conditions):
        geom = cond.GetGeometry()
        if i == 0:
            x_coords.append(geom[0].X)
            y_coords.append(geom[0].Y)
        x_coords.append(geom[1].X)
        y_coords.append(geom[1].Y)

    x_coords.append(x_coords[0])
    y_coords.append(y_coords[0])
    return x_coords, y_coords


def add_gap_quads(ax, surrogate_mp, skin_mp, color):
    """Fill gap as quads between surrogate segments and projected skin nodes."""
    gap_added = False
    for cond in surrogate_mp.Conditions:
        geom = cond.GetGeometry()
        p0 = (geom[0].X, geom[0].Y)
        p1 = (geom[1].X, geom[1].Y)

        proj_pts = []
        for node in geom:
            try:
                proj_id = node.GetValue(KratosMultiphysics.IgaApplication.PROJECTION_NODE_ID)
            except Exception:
                proj_id = 0
            if proj_id and skin_mp.HasNode(proj_id):
                proj_node = skin_mp.GetNode(proj_id)
                proj_pts.append((proj_node.X, proj_node.Y))

        if len(proj_pts) == 2:
            quad = [p0, p1, proj_pts[1], proj_pts[0]]
            patch = patches.Polygon(quad, closed=True, facecolor=color, edgecolor='none', lw=0, alpha=1.0, zorder=1)
            ax.add_patch(patch)
            gap_added = True
    return gap_added


def build_ring_clip(outer_loop, inner_loop=None):
    """Build a PathPatch to use as clip (outer minus inner)."""
    ox, oy = trim_loop(*outer_loop)
    if len(ox) < 3:
        return None

    verts = []
    codes = []

    for idx, (x, y) in enumerate(zip(ox, oy)):
        codes.append(Path.MOVETO if idx == 0 else Path.LINETO)
        verts.append((x, y))
    codes.append(Path.CLOSEPOLY)
    verts.append((ox[0], oy[0]))

    if inner_loop:
        ix, iy = trim_loop(*inner_loop)
        if len(ix) >= 3:
            for idx, (x, y) in enumerate(zip(reversed(ix), reversed(iy))):
                codes.append(Path.MOVETO if idx == 0 else Path.LINETO)
                verts.append((x, y))
            codes.append(Path.CLOSEPOLY)
            verts.append((ix[-1], iy[-1]))

    clip_path = Path(verts, codes)
    return patches.PathPatch(clip_path, facecolor='none', edgecolor='none')


def build_ring_patch(inner_loop, outer_loop, color, alpha=1.0, zorder=1):
    """Return a filled patch between outer and inner loops without painting back over the inner loop."""
    ox, oy = trim_loop(*outer_loop)
    ix, iy = trim_loop(*inner_loop)
    if len(ox) < 3 or len(ix) < 3:
        return None

    def ensure_ccw(xs, ys):
        if polygon_area(xs, ys) < 0:
            xs = list(reversed(xs))
            ys = list(reversed(ys))
        return xs, ys

    def ensure_cw(xs, ys):
        if polygon_area(xs, ys) > 0:
            xs = list(reversed(xs))
            ys = list(reversed(ys))
        return xs, ys

    ox, oy = ensure_ccw(list(ox), list(oy))
    ix, iy = ensure_cw(list(ix), list(iy))

    outer_pts = list(zip(ox, oy))
    inner_pts = list(zip(ix, iy))
    if outer_pts[0] != outer_pts[-1]:
        outer_pts.append(outer_pts[0])
    if inner_pts[0] != inner_pts[-1]:
        inner_pts.append(inner_pts[0])

    ring_path = Path.make_compound_path(Path(outer_pts), Path(inner_pts))
    return patches.PathPatch(
        ring_path,
        facecolor=color,
        edgecolor='none',
        lw=0,
        alpha=alpha,
        zorder=zorder,
        fill=True,
        joinstyle='round',
        capstyle='round',
        antialiased=True,
    )


def point_key(pt, tol=1e-6):
    """Hashable point key with tolerance."""
    return (round(pt[0] / tol) * tol, round(pt[1] / tol) * tol)


def build_closed_polyline(segments, tol=1e-6):
    """Chain segments by matching start/end keys into a closed polyline."""
    if not segments:
        return []
    from collections import defaultdict

    adj = defaultdict(list)
    for idx, seg in enumerate(segments):
        if not seg:
            continue
        skey = point_key(seg[0], tol)
        ekey = point_key(seg[-1], tol)
        adj[skey].append((idx, True))   # forward uses seg as-is
        adj[ekey].append((idx, False))  # backward means reverse seg

    start_key = next(iter(adj))
    polyline = []
    used = set()
    current_key = start_key

    while True:
        candidates = [(i, forward) for (i, forward) in adj[current_key] if i not in used]
        if not candidates:
            break
        seg_idx, forward = candidates[0]
        used.add(seg_idx)
        seg = segments[seg_idx]
        if forward:
            points = seg
            next_key = point_key(seg[-1], tol)
        else:
            points = list(reversed(seg))
            next_key = point_key(seg[0], tol)

        if not polyline:
            polyline.extend(points)
        else:
            if point_key(polyline[-1], tol) == point_key(points[0], tol):
                polyline.extend(points[1:])
            else:
                polyline.extend(points)

        current_key = next_key
        if current_key == start_key:
            break

    if polyline and point_key(polyline[0], tol) != point_key(polyline[-1], tol):
        polyline.append(polyline[0])
    return polyline


def polygon_area(xs, ys):
    """Signed polygon area (shoelace); positive if CCW."""
    if len(xs) < 3:
        return 0.0
    area = 0.0
    for i in range(len(xs) - 1):
        area += xs[i] * ys[i + 1] - xs[i + 1] * ys[i]
    return 0.5 * area


def closest_index(point, polyline):
    """Index of polyline point closest to given point."""
    if not polyline:
        return None
    arr = np.array(polyline)
    dists = np.linalg.norm(arr - np.array(point), axis=1)
    return int(np.argmin(dists))


def shortest_path(polyline, idx_start, idx_end):
    """Return shortest cyclic path from idx_start to idx_end (inclusive)."""
    if idx_start is None or idx_end is None or not polyline:
        return []
    n = len(polyline)
    if n == 0:
        return []
    # ensure open polyline (avoid duplicate last point if present)
    if np.allclose(polyline[0], polyline[-1]):
        polyline = polyline[:-1]
        n -= 1
    forward = []
    i = idx_start
    while True:
        forward.append(polyline[i])
        if i == idx_end:
            break
        i = (i + 1) % n

    backward = []
    i = idx_start
    while True:
        backward.append(polyline[i])
        if i == idx_end:
            break
        i = (i - 1 + n) % n

    return forward if len(forward) <= len(backward) else backward


def ray_segment_intersection(p, r, a, b, t_min=1e-8):
    """Return intersection point of ray p+t*r (t>=t_min) with segment ab, or None."""
    r = np.array(r, dtype=float)
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    s = b - a
    cross_rs = r[0] * s[1] - r[1] * s[0]
    if abs(cross_rs) < 1e-12:
        return None
    q_p = a - p
    t = (q_p[0] * s[1] - q_p[1] * s[0]) / cross_rs
    u = (q_p[0] * r[1] - q_p[1] * r[0]) / cross_rs
    if t >= t_min and 0.0 <= u <= 1.0:
        return p + t * r
    return None


def find_skin_intersection(origin, direction, skin_loop):
    """Find first intersection of ray with a skin loop polyline."""
    if not skin_loop:
        return None
    xs, ys = skin_loop
    pts = list(zip(xs, ys))
    if len(pts) < 2:
        return None
    best_pt = None
    best_t = float("inf")
    p = np.array(origin, dtype=float)
    r = np.array(direction, dtype=float)
    for i in range(len(pts) - 1):
        inter = ray_segment_intersection(p, r, pts[i], pts[i + 1])
        if inter is not None:
            t_val = np.linalg.norm(inter - p)
            if t_val < best_t:
                best_t = t_val
                best_pt = inter
    return best_pt


def oriented_normal(p0, p1, ref_point):
    """Unit normal to segment p0->p1 pointing toward ref_point."""
    v = np.array(p1) - np.array(p0)
    n = np.array([-v[1], v[0]])
    ref_vec = np.array(ref_point) - (np.array(p0) + np.array(p1)) * 0.5
    if np.dot(n, ref_vec) < 0:
        n = -n
    norm = np.linalg.norm(n)
    if norm == 0:
        return np.array([0.0, 0.0])
    return n / norm


def point_on_polyline(polyline, fraction):
    """Return point at given fraction (0-1) of polyline arc length."""
    if not polyline or len(polyline) < 2:
        return polyline[0] if polyline else (0.0, 0.0)
    seg_lengths = []
    total = 0.0
    for a, b in zip(polyline[:-1], polyline[1:]):
        l = np.hypot(b[0] - a[0], b[1] - a[1])
        seg_lengths.append(l)
        total += l
    if total == 0:
        return polyline[0]
    target = fraction * total
    acc = 0.0
    for (a, b), l in zip(zip(polyline[:-1], polyline[1:]), seg_lengths):
        if acc + l >= target:
            t = (target - acc) / l if l != 0 else 0
            return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))
        acc += l
    return polyline[-1]


def point_and_tangent_on_polyline(polyline, fraction):
    """Return point and unit tangent at a given fraction of arc length."""
    if not polyline or len(polyline) < 2:
        return (polyline[0] if polyline else (0.0, 0.0)), np.array([1.0, 0.0])
    seg_lengths = []
    total = 0.0
    for a, b in zip(polyline[:-1], polyline[1:]):
        l = np.hypot(b[0] - a[0], b[1] - a[1])
        seg_lengths.append(l)
        total += l
    if total == 0:
        return polyline[0], np.array([1.0, 0.0])
    target = fraction * total
    acc = 0.0
    for (a, b), l in zip(zip(polyline[:-1], polyline[1:]), seg_lengths):
        if acc + l >= target:
            t = (target - acc) / l if l != 0 else 0
            point = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))
            tangent = np.array([b[0] - a[0], b[1] - a[1]])
            norm = np.linalg.norm(tangent)
            if norm == 0:
                tangent = np.array([1.0, 0.0])
            else:
                tangent = tangent / norm
            return point, tangent
        acc += l
    point = polyline[-1]
    last_tangent = np.array([polyline[-1][0] - polyline[-2][0], polyline[-1][1] - polyline[-2][1]])
    norm = np.linalg.norm(last_tangent)
    if norm == 0:
        last_tangent = np.array([1.0, 0.0])
    else:
        last_tangent = last_tangent / norm
    return point, last_tangent


def find_condition_intersection(origin, direction, segments):
    """Return closest intersection of ray with any segment in list."""
    if not segments:
        return None
    p = np.array(origin, dtype=float)
    r = np.array(direction, dtype=float)
    best_pt = None
    best_t = float("inf")
    for a, b in segments:
        inter = ray_segment_intersection(p, r, np.array(a), np.array(b))
        if inter is not None:
            t_val = np.linalg.norm(inter - p)
            if t_val < best_t:
                best_t = t_val
                best_pt = inter
    return best_pt


def ensure_knot_vector(knots: List[float], degree: int, num_ctrl_pts: int) -> List[float]:
    """Return a knot vector of length num_ctrl_pts + degree + 1, falling back to open uniform."""
    expected_len = num_ctrl_pts + degree + 1
    if len(knots) == expected_len:
        return knots
    if len(knots) >= 2:
        start, end = knots[0], knots[-1]
        interior_count = max(expected_len - 2 * (degree + 1), 0)
        interior = []
        if interior_count > 0:
            interior = list(np.linspace(start, end, interior_count + 2)[1:-1])
        return [start] * (degree + 1) + interior + [end] * (degree + 1)
    return knots


def bspline_basis(i: int, degree: int, knots: List[float], u: float) -> float:
    """Cox-de Boor recursive definition."""
    if degree == 0:
        if (knots[i] <= u < knots[i + 1]) or (u == knots[-1] and knots[i] <= u <= knots[i + 1]):
            return 1.0
        return 0.0

    left_denom = knots[i + degree] - knots[i]
    right_denom = knots[i + degree + 1] - knots[i + 1]
    left = 0.0
    right = 0.0

    if left_denom > 0:
        left = (u - knots[i]) / left_denom * bspline_basis(i, degree - 1, knots, u)
    if right_denom > 0:
        right = (knots[i + degree + 1] - u) / right_denom * bspline_basis(i + 1, degree - 1, knots, u)
    return left + right


def evaluate_nurbs_curve(degree: int, knots: List[float], weights: List[float],
                         ctrl_pts: List[Tuple[float, float]], n_points: int = 100):
    """Compute sampled points of a NURBS curve."""
    knots = ensure_knot_vector(knots, degree, len(ctrl_pts))
    u_start = knots[degree]
    u_end = knots[-degree - 1]
    us = np.linspace(u_start, u_end, n_points)

    xs, ys = [], []
    for u in us:
        num_x = num_y = denom = 0.0
        for i, (x, y) in enumerate(ctrl_pts):
            N = bspline_basis(i, degree, knots, u)
            wN = weights[i] * N
            num_x += wN * x
            num_y += wN * y
            denom += wN
        if denom != 0:
            xs.append(num_x / denom)
            ys.append(num_y / denom)
    return xs, ys


def read_polynomial_nurbs(filepath):
    """Return list of curves with degree, knots, weights, and control points."""
    if not os.path.exists(filepath):
        return []

    curves = []
    with open(filepath, "r", encoding="utf-8") as infile:
        for line in infile:
            line = line.strip()
            if not line or not line.startswith("curve"):
                continue

            degree_match = re.search(r"degree=([0-9]+)", line)
            knots_match = re.search(r"knots=\[([^\]]+)\]", line)
            weights_match = re.search(r"weights=\[([^\]]+)\]", line)
            point_strings = re.findall(r"\(([^)]*)\)", line)

            if not (degree_match and knots_match and weights_match and point_strings):
                continue

            try:
                degree = int(degree_match.group(1))
                knots = [float(k) for k in knots_match.group(1).split() if k]
                weights = [float(w) for w in weights_match.group(1).split() if w]
            except ValueError:
                continue

            control_points = []
            for pt_str in point_strings:
                try:
                    parts = [float(v) for v in pt_str.split(",") if v]
                    if len(parts) >= 3:
                        control_points.append((parts[1], parts[2]))
                except ValueError:
                    control_points = []
                    break

            if control_points and len(weights) == len(control_points):
                curves.append({
                    "degree": degree,
                    "knots": knots,
                    "weights": weights,
                    "control_points": control_points
                })

    return curves


if __name__ == "__main__":

    # Clean old artifacts
    for fname in [
        "centers.txt",
        "centers_load.txt",
        "txt_files/Id_active_control_points.txt",
        "txt_files/Id_active_control_points_condition.txt",
        "txt_files/polynomial_curves_inner.txt",
        "txt_files/polynomial_curves_outer.txt",
    ]:
        if os.path.exists(fname):
            os.remove(fname)

    with open("ProjectParametersSolidGap.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Initialize()
    simulation.InitializeSolutionStep()

    index_param_space = 0
    while not (parameters["modelers"][index_param_space]["Parameters"].Has("number_of_knot_spans")):
        index_param_space += 1

    n_knot_span_u = parameters["modelers"][index_param_space]["Parameters"]["number_of_knot_spans"][0].GetInt()
    n_knot_span_v = parameters["modelers"][index_param_space]["Parameters"]["number_of_knot_spans"][1].GetInt()
    initial_u = parameters["modelers"][index_param_space]["Parameters"]["lower_point_uvw"][0].GetDouble()
    initial_v = parameters["modelers"][index_param_space]["Parameters"]["lower_point_uvw"][1].GetDouble()
    final_u = parameters["modelers"][index_param_space]["Parameters"]["upper_point_uvw"][0].GetDouble()
    final_v = parameters["modelers"][index_param_space]["Parameters"]["upper_point_uvw"][1].GetDouble()

    knots_u_1 = [initial_u + (final_u - initial_u) / n_knot_span_u * j for j in range(n_knot_span_u + 1)]
    knots_v_1 = [initial_v + (final_v - initial_v) / n_knot_span_v * j for j in range(n_knot_span_v + 1)]
    
    basis_functions_order = parameters["modelers"][index_param_space]["Parameters"]["polynomial_order"][0].GetInt()

    fig = plt.gcf()
    ax = fig.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.grid(False)
    label_used = set()
    param_grid_lines = []
    skin_loops = {}
    skin_segments = {}
    gap_label_added = False
    surrogate_loops = {}
    polynomial_segments = []
    polynomial_segments_outer = []
    surrogate_segments = {}

    for side in ["outer"]:
        surrogate_part_name = f"IgaModelPart.surrogate_{side}"
        skin_part_name = f"skin_model_part.{side}"

        # Surrogate boundary (segments)
        if model.HasModelPart(surrogate_part_name):
            mp = model.GetModelPart(surrogate_part_name)
            mp_skin = model.GetModelPart(skin_part_name)
            surrogate_segments_list = []
            if len(mp.Conditions) > 0:
                for cond in mp.Conditions:
                    geom = cond.GetGeometry()
                    x0, y0 = geom[0].X, geom[0].Y
                    x1, y1 = geom[1].X, geom[1].Y
                    surrogate_segments_list.append([(x0, y0), (x1, y1)])
                    

                    ax.plot([x0, x1], [y0, y1], color=SURROGATE_COLOR, linewidth=3,
                            label="Surrogate boundary" if "Surrogate boundary" not in label_used else "")

                    for i_vertex in range(0, 2):
                        projection_node_id = geom[i_vertex].GetValue(KratosMultiphysics.IgaApplication.PROJECTION_NODE_ID)
                        projection_node = mp_skin.GetNode(projection_node_id) 
                        x_skin = projection_node.X
                        y_skin = projection_node.Y
                        
                        ax.plot([geom[i_vertex].X, x_skin], [geom[i_vertex].Y, y_skin], color="k", linewidth=1.5,
                                label="Gap Interfaces" if "Gap Interfaces" not in label_used else "")
                        
                label_used.add("Surrogate boundary")
                label_used.add("Gap Interfaces")
            polyline = build_closed_polyline(surrogate_segments_list)
            if polyline:
                sx, sy = zip(*polyline)
                surrogate_loops[side] = (sx, sy)

        # Skin boundary
        if model.HasModelPart(skin_part_name):
            mp = model.GetModelPart(skin_part_name)
            segs = []
            for cond in mp.Conditions:
                geom = cond.GetGeometry()
                if len(geom) >= 2:
                    segs.append(((geom[0].X, geom[0].Y), (geom[1].X, geom[1].Y)))
            if segs:
                skin_segments[side] = segs
            x_skin, y_skin = extract_boundary_coords(mp)
            if x_skin:
                skin_loops[side] = (x_skin, y_skin)
                ax.plot(x_skin, y_skin, color=TRUE_BOUNDARY_COLOR, linewidth=2,
                        label="True boundary" if "True boundary" not in label_used else "")
                label_used.add("True boundary")

        # Gauss points support (if needed later)
        gp_part_name = f"IgaModelPart.SBM_Support_{side}"
        if model.HasModelPart(gp_part_name):
            mp = model.GetModelPart(gp_part_name)
            if len(mp.Conditions) > 0:
                for cond in mp.Conditions:
                    geom = cond.GetGeometry()
                    N = geom.ShapeFunctionsValues()
                    gp = np.zeros(2)
                    for i, node in enumerate(geom):
                        gp[0] += N[0, i] * node.X
                        gp[1] += N[0, i] * node.Y

        # Fill gap pieces for this side if both parts exist
        # if model.HasModelPart(surrogate_part_name) and model.HasModelPart(skin_part_name):
        #     mp_surr = model.GetModelPart(surrogate_part_name)
        #     mp_skin = model.GetModelPart(skin_part_name)
        #     gap_added_now = add_gap_quads(ax, mp_surr, mp_skin, GAP_FILL_COLOR)
        #     if gap_added_now and not gap_label_added:
        #         label_used.add("Gap elements")
        #         gap_label_added = True

    
    polynomial_curves_outer_path = os.path.join("txt_files", "polynomial_curves_outer.txt")
    for curve in read_polynomial_nurbs(polynomial_curves_outer_path):
        xs, ys = evaluate_nurbs_curve(
            degree=curve["degree"],
            knots=curve["knots"],
            weights=curve["weights"],
            ctrl_pts=curve["control_points"],
            n_points=200
        )
        if xs and ys:
            pts = list(zip(xs, ys))
            polynomial_segments_outer.append(pts)
            ax.plot(xs, ys, color=POLYNOMIAL_NURBS_COLOR, linewidth=2,
                    label="Polynomial NURBS" if "Polynomial NURBS" not in label_used else "")
            label_used.add("Polynomial NURBS")

    combined_outer_polyline = build_closed_polyline(polynomial_segments_outer)

    if not combined_outer_polyline and polynomial_segments_outer:
        flat_pts = [p for seg in polynomial_segments_outer for p in seg]
        cx = np.mean([p[0] for p in flat_pts])
        cy = np.mean([p[1] for p in flat_pts])
        combined_outer_polyline = sorted(flat_pts, key=lambda p: np.arctan2(p[1] - cy, p[0] - cx))
        combined_outer_polyline.append(combined_outer_polyline[0])


    # Fill the region between surrogate outer and outer polynomial polyline
    if combined_outer_polyline and "outer" in surrogate_loops:
        sox, soy = surrogate_loops["outer"]
        px, py = zip(*combined_outer_polyline)
        sur_area = abs(polygon_area(list(sox), list(soy)))
        poly_area = abs(polygon_area(list(px), list(py)))
        outer_loop = (sox, soy)
        inner_loop = (px, py)
        if poly_area > sur_area:
            outer_loop, inner_loop = (px, py), (sox, soy)

        outer_gap_patch = build_ring_patch(inner_loop, outer_loop, GAP_FILL_COLOR, alpha=1.0, zorder=0)
        if outer_gap_patch:
            ax.add_patch(outer_gap_patch)
            if not gap_label_added:
                label_used.add("Gap elements")
                gap_label_added = True
        

    # === PARAMETER SPACE GRID (clipped to surrogate area if available) ===
    for u in knots_u_1:
        lines = ax.plot([u, u], [knots_v_1[0], knots_v_1[-1]], linestyle=':', color=PARAM_SPACE_COLOR, linewidth=1.5,
                        label='Parameter space' if u == knots_u_1[0] else "")
        for ln in lines:
            ln.set_zorder(-0.5)
        param_grid_lines.extend(lines)
    for v in knots_v_1:
        lines = ax.plot([knots_u_1[0], knots_u_1[-1]], [v, v], linestyle=':', color=PARAM_SPACE_COLOR, linewidth=1.5)
        for ln in lines:
            ln.set_zorder(-0.5)
        param_grid_lines.extend(lines)
    label_used.add("Parameter space")
    if "outer" in skin_loops:
        clip_patch = build_ring_clip(
            skin_loops.get("outer"),
            skin_loops.get("inner")
        )
        if clip_patch:
            ax.add_patch(clip_patch)
            for ln in param_grid_lines:
                ln.set_clip_path(clip_patch)

    try:
        is_plot_all = int(input("Plot all the active knot spans and gauss points? 1: yes, 0: no "))
    except ValueError:
        is_plot_all = 0

    area = 0
    active_gauss_points = []

    if is_plot_all == 1:
        mp = model["IgaModelPart.StructuralAnalysisDomain"]
        for elem in mp.Elements:
            center = elem.GetGeometry().Center()
            x = center.X
            y = center.Y
            active_gauss_points.append((x, y))
            # ax.plot(x, y, "o", color=ELEMENT_GAUSS_COLOR, markersize=3,
            #         label="Elemental Gauss Points" if "Elemental Gauss Points" not in label_used else "")
            label_used.add("Elemental Gauss Points")

            weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
            area += weight

    print("Total area of the model:", area)

    # mp = model["IgaModelPart.GapElements"]
    # for elem in mp.Elements:
    #     center = elem.GetGeometry().Center()
    #     x = center.X
    #     y = center.Y
    #     ax.plot(x, y, "ro", markersize=2, zorder =5, label="Elemental Gauss Points" if "Elemental Gauss Points" not in label_used else "")
    #     weight = elem.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
    #     label_used.add("Gap Gauss Points")
    #     area += weight
        
        
    # mp = model["IgaModelPart.left"]
    # mp_skin = model["skin_model_part.inner"]
    perimeter = 0.0
    # for cond in mp.Conditions:
    #     center = cond.GetGeometry().Center()
    #     x = center.X
    #     y = center.Y
    #     #ax.plot(x, y, "o", color="blue", markersize=3, label="ESBM Gauss Points" if "ESBM Gauss Points" not in label_used else "")
        
    #     weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
    #     normal = cond.GetValue(KratosMultiphysics.NORMAL)
    #     # ax.arrow(x, y, normal[0] * 0.02, normal[1] * 0.02, head_width=0.001, head_length=0.02, fc='blue', ec='blue')
        
    #     projection = cond.GetValue(KratosMultiphysics.IgaApplication.PROJECTION_NODE_COORDINATES)
        
    #     if (len(projection) == 0): continue

    #     # ax.plot([x, projection[0]], [y, projection[1]], "g-o", markersize=5, linewidth = 3)
        
    #     # interpolation_node_ids = cond.GetValue(KratosMultiphysics.IgaApplication.TEMP_INTERPOLATION_NODES_ID)
    #     # for id in interpolation_node_ids:
    #     #     id =int(id)
    #     #     node = mp_skin.GetNode(id)
            
    #     #     ax.plot(node.X, node.Y, "ko", markersize=3)
    #     # perimeter += weight

    mp = model["IgaModelPart.GapInterfaces"]
    # perimeter = 0.0
    for cond in mp.Conditions:
        center = cond.GetGeometry().Center()
        x = center.X
        y = center.Y
        #ax.plot(x, y, "ko", markersize=2, label="Gap Interface Gauss Points" if "Gap Interface Gauss Points" not in label_used else "")
        weight = cond.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT)
        normal = cond.GetValue(KratosMultiphysics.NORMAL)
        # ax.arrow(x, y, normal[0] * 0.1, normal[1] * 0.1, head_width=0.05, head_length=0.1, fc='blue', ec='blue')
        # perimeter += weight

    true_perimeter = 2 * np.pi * 0.17
    true_area = 1 - np.pi * 0.17 * 0.17
    print("Total perimeter of the model:", perimeter)
    error_area = abs(area - true_area) / true_area if true_area != 0 else 0.0
    error_perimeter = abs(perimeter - true_perimeter) / true_perimeter if true_perimeter != 0 else 0.0
    print("Area error:", error_area)
    print("Perimeter error:", error_perimeter)

    # === ACTIVE CELLS SHADE ===
    if active_gauss_points:
        active_cells = set()
        for (x_gp, y_gp) in active_gauss_points:
            for i in range(len(knots_u_1) - 1):
                for j in range(len(knots_v_1) - 1):
                    u_min, u_max = knots_u_1[i], knots_u_1[i + 1]
                    v_min, v_max = knots_v_1[j], knots_v_1[j + 1]
                    if u_min <= x_gp <= u_max and v_min <= y_gp <= v_max:
                        active_cells.add((i, j))

        for (i, j) in active_cells:
            u_min, u_max = knots_u_1[i], knots_u_1[i + 1]
            v_min, v_max = knots_v_1[j], knots_v_1[j + 1]
            rect = patches.Rectangle((u_min, v_min), u_max - u_min, v_max - v_min,
                                     linewidth=0, edgecolor=None, facecolor=ACTIVE_ELEMENT_COLOR, zorder=-1,
                                     label="Active elements" if "Active elements" not in label_used else "")
            ax.add_patch(rect)
            label_used.add("Active elements")

    ax.set_aspect('equal')
    ax.set_xlabel("")
    ax.set_ylabel("")

    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)

    fig.tight_layout()
    # fig.savefig("explanation_sbm_zoom_gap_fill.pdf", format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig("bunny_GAP_linear_20x20.pdf", format="pdf", dpi=1000, bbox_inches="tight")

    custom_handles = [
        Patch(facecolor=ACTIVE_ELEMENT_COLOR, edgecolor='lightgray', linestyle=':', linewidth=1.0, label='Active elements'),
        Patch(facecolor=GAP_FILL_COLOR, edgecolor='black', linewidth=1.5, label='Gap elements'),
        Line2D([0], [0], color=PARAM_SPACE_COLOR, linestyle=':', label='Parameter space'),
        Line2D([0], [0], color=SURROGATE_COLOR, linewidth=2, label='Surrogate boundary'),
        Line2D([0], [0], color=TRUE_BOUNDARY_COLOR, linewidth=2, label='True boundary'),
        Line2D([0], [0], marker='o', color=ELEMENT_GAUSS_COLOR, linestyle='None', markersize=3, label='Elemental Gauss Points'),
        Line2D([0], [0], marker='o', color="red", linestyle='None', markersize=3, label='Gap Gauss Points'),
        Line2D([0], [0], marker='o', color=ESBM_GAUSS_COLOR, linestyle='None', markersize=3, label='Boundary Gauss Points'),
        Line2D([0], [0], marker='o', color="black", linestyle='None', markersize=3, label='Gap Interface Gauss Points'),
        Line2D([0], [0], color=POLYNOMIAL_NURBS_COLOR, linewidth=2, label='Polynomial NURBS'),
    ]

    legend_fig, legend_ax = plt.subplots(figsize=(4, 1.8))
    legend_ax.axis('off')
    legend_ax.legend(
        handles=custom_handles,
        ncol=3,
        frameon=False,
        columnspacing=1.5,
        handletextpad=1.2
    )

    legend_fig.tight_layout()
    plt.show()
