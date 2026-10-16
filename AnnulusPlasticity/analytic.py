"""Small-strain, plane-strain J2 annular torsion; monotone positive loading."""
import numpy as np

RI, RO = 1.0, 2.0
E, NU, YIELD, H = 67000.0, 0.0, 201.0, 10000.0
G = E / (2 * (1 + NU))
END_TIME = 0.002
FINAL_PLASTIC_RADIUS = 1.5
C_FINAL = YIELD * FINAL_PLASTIC_RADIUS**2 / np.sqrt(3)


def evaluate(points, time):
    """Return u, flattened grad(u), [sxx, syy, sxy], equivalent plastic strain.

    Outer boundary is fixed; gamma = w' - w/r, tau = C/r**2.
    The analytic extension outside the physical annulus is also used on SBM bulk.
    """
    if not 0 <= time <= END_TIME * (1 + 1e-10):
        raise ValueError('The reference assumes monotone loading on [0, END_TIME].')
    points = np.asarray(points, dtype=float)
    x, y = points.T
    r = np.hypot(x, y)
    if np.any(r <= 0):
        raise ValueError('The annular reference is singular at the origin.')
    c = C_FINAL * time / END_TIME
    rp = np.sqrt(np.sqrt(3) * c / YIELD)
    tau = c / r**2
    alpha = np.maximum((np.sqrt(3) * tau - YIELD) / H, 0)
    gamma = tau / G + np.sqrt(3) * alpha
    # f = w/r; integrate f'=gamma/r inward from f(RO)=0.
    f = -c / (2 * G) * (1/r**2 - 1/RO**2)
    plastic = r < rp
    if np.any(plastic):
        r_pl = r[plastic]
        f[plastic] += (-3*c/(2*H)*(1/r_pl**2 - 1/rp**2)
                       + np.sqrt(3)*YIELD/H*np.log(rp/r_pl))
    u = np.column_stack((-f*y, f*x))
    a = gamma / r**2
    grad = np.column_stack((-a*x*y, -f-a*y*y, f+a*x*x, a*x*y))
    stress = np.column_stack((-2*c*x*y/r**4, 2*c*x*y/r**4,
                              c*(x*x-y*y)/r**4))
    return u, grad, stress, alpha


def traction(points, normals, time):
    stress = evaluate(points, time)[2]
    normals = np.asarray(normals)
    return np.column_stack((stress[:,0]*normals[:,0] + stress[:,2]*normals[:,1],
                            stress[:,2]*normals[:,0] + stress[:,1]*normals[:,1]))


def boundary_traction(points, time):
    """Prescribed traction using the material-outward normal of the annulus.

    The linearized skin can store curve-oriented normals (+e_r on both
    supplied CCW circles). Those are not material-outward on the inner loop.
    Classify points by their nearest circular boundary, including the small
    geometric deviations of the discretized skin.
    """
    points = np.asarray(points,dtype=float)
    radii = np.linalg.norm(points,axis=1)
    if np.any(radii <= 0):
        raise ValueError('Boundary normal is undefined at the origin.')
    sign = np.where(radii < (RI+RO)/2,-1.,1.)
    normals = sign[:,None]*points/radii[:,None]
    return traction(points,normals,time)
