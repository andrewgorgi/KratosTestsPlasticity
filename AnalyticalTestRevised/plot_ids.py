"""
Plot curves from a geometry JSON in your "Lines" format.

Features:
- Each Line (piece) is plotted with a different color (matplotlib default cycle).
- The Line id is drawn next to the curve piece (at its mid-parameter point).
- Supports NURBS degree 1/2/3 with knotVector + weights + CPCoordinates.
- Plots in XY (ignores Z). You can switch to XZ/YZ if needed.

Usage:
- Set INPUT_JSON at the top.
- Run: python plot_lines_with_ids.py
"""

import json
import math
import numpy as np
import matplotlib.pyplot as plt


# ================== USER PARAMETERS ==================

# INPUT_JSON = "Gear_moved_rotated.json"
# INPUT_JSON = "Pinion_rotated.json"
INPUT_JSON = "geometry/rectangle.nurbs.json"
SAMPLES_PER_CURVE = 250
PLOT_CONTROL_POLYGON = False
LABEL_OFFSET = 0.0  # in plot units; set >0 to shift labels a bit

# =====================================================


def cox_de_boor(i: int, p: int, u: float, U: list) -> float:
    if p == 0:
        last = abs(u - U[-1]) < 1e-12
        return 1.0 if (U[i] <= u < U[i + 1] or (last and abs(u - U[i + 1]) < 1e-12)) else 0.0

    a = 0.0
    b = 0.0
    d1 = U[i + p] - U[i]
    d2 = U[i + p + 1] - U[i + 1]
    if d1 > 1e-14:
        a = (u - U[i]) / d1 * cox_de_boor(i, p - 1, u, U)
    if d2 > 1e-14:
        b = (U[i + p + 1] - u) / d2 * cox_de_boor(i + 1, p - 1, u, U)
    return a + b


def eval_nurbs_point(cps: list, wts: list, U: list, p: int, u: float) -> tuple:
    n = len(cps) - 1
    N = [cox_de_boor(i, p, u, U) for i in range(n + 1)]
    w = sum(N[i] * wts[i] for i in range(n + 1))
    if abs(w) < 1e-14:
        return (math.nan, math.nan)

    x = sum(N[i] * cps[i][0] * wts[i] for i in range(n + 1)) / w
    y = sum(N[i] * cps[i][1] * wts[i] for i in range(n + 1)) / w
    return (x, y)


def sample_curve(line: dict, samples: int) -> np.ndarray:
    cps = line["CPCoordinates"]
    wts = line["Weights"]
    U = line["knotVector"]
    p = int(line["pDegree"])

    # Parameter domain for clamped B-splines / NURBS: [U[p], U[-p-1]]
    u0 = float(U[p])
    u1 = float(U[-p - 1])

    us = np.linspace(u0, u1, samples)
    pts = [eval_nurbs_point(cps, wts, U, p, float(u)) for u in us]
    return np.array(pts, dtype=float)


def main():
    with open(INPUT_JSON, "r", encoding="utf-8") as f:
        obj = json.load(f)

    lines = obj.get("Lines", [])
    if not lines:
        raise RuntimeError("No 'Lines' found in the input JSON.")

    # Stable ordering by id
    lines = sorted(lines, key=lambda L: int(L.get("id", 0)))

    plt.figure(figsize=(8, 6))

    for line in lines:
        pts = sample_curve(line, SAMPLES_PER_CURVE)
        plt.plot(pts[:, 0], pts[:, 1])

        # Label at middle sample
        mid = pts[len(pts) // 2]
        lx = float(mid[0]) + LABEL_OFFSET
        ly = float(mid[1]) + LABEL_OFFSET
        plt.text(lx, ly, str(line["id"]), fontsize=9)

        if PLOT_CONTROL_POLYGON:
            cps = np.array(line["CPCoordinates"], dtype=float)
            plt.plot(cps[:, 0], cps[:, 1], linestyle="--", marker="o", alpha=0.6)

    plt.axis("equal")
    plt.title("Lines with per-piece color and id labels")
    plt.show()


if __name__ == "__main__":
    main()
