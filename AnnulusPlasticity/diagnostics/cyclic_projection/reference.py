"""Exact stress-controlled cyclic torsion with linear isotropic hardening.

Signed plastic shear follows stress direction; equivalent plastic strain follows
the maximum absolute torque. Unloading retains plastic strain. The outer circle
is fixed throughout. This reference is separate from the monotone production one.
"""
import numpy as np

E, NU, H, YIELD = 67000., 0., 10000., 201.
G = E/(2*(1+NU))
RO = 2.
C_FINAL = YIELD*1.5**2/np.sqrt(3)
TIMES = np.arange(7)*.001
LEVELS = np.array([0., 1., 0., -1.2, 0., 1.4, 0.])
END_TIME = float(TIMES[-1])


def plastic(r, c):
    alpha = np.maximum((np.sqrt(3)*c/r**2-YIELD)/H, 0.)
    f = np.zeros_like(r)
    rp = np.sqrt(np.sqrt(3)*c/YIELD)
    active = r < rp
    if active.any():
        rr = r[active]
        f[active] = -3*c/(2*H)*(1/rr**2-1/rp**2)+np.sqrt(3)*YIELD/H*np.log(rp/rr)
    return alpha, f


def evaluate(points, time):
    points = np.asarray(points, dtype=float)
    x, y = points.T
    r = np.hypot(x, y)
    assert np.all(r > 0.)
    current = float(np.interp(time, TIMES, LEVELS))
    history = list(LEVELS[TIMES < time-1e-12])+[current]
    maximum = 0.
    plastic_gamma, plastic_f = np.zeros_like(r), np.zeros_like(r)
    for level in history:
        if abs(level) > maximum:
            alpha_new, f_new = plastic(r, abs(level)*C_FINAL)
            alpha_old, f_old = plastic(r, maximum*C_FINAL)
            plastic_gamma += np.sign(level)*np.sqrt(3)*(alpha_new-alpha_old)
            plastic_f += np.sign(level)*(f_new-f_old)
            maximum = abs(level)
    alpha = plastic(r, maximum*C_FINAL)[0]
    c = current*C_FINAL
    f = -c/(2*G)*(1/r**2-1/RO**2)+plastic_f
    gamma = c/(G*r**2)+plastic_gamma
    u = np.column_stack((-f*y, f*x))
    a = gamma/r**2
    grad = np.column_stack((-a*x*y, -f-a*y*y, f+a*x*x, a*x*y))
    stress = np.column_stack((-2*c*x*y/r**4, 2*c*x*y/r**4, c*(x*x-y*y)/r**4))
    return u, grad, stress, alpha


def boundary_traction(points, time):
    points = np.asarray(points, dtype=float)
    stress = evaluate(points, time)[2]
    r = np.linalg.norm(points, axis=1)
    normal = np.where(r < 1.5, -1., 1.)[:, None]*points/r[:, None]
    return np.column_stack((stress[:, 0]*normal[:, 0]+stress[:, 2]*normal[:, 1],
                            stress[:, 2]*normal[:, 0]+stress[:, 1]*normal[:, 1]))


def verify():
    """Check closed forms against an independent stress-driven scalar return map."""
    r = np.array([.95, 1., 1.1, 1.35, 1.6, 1.9, 2.])
    alpha = np.zeros_like(r)
    gamma_p = np.zeros_like(r)
    maximum_error = 0.
    for time in np.linspace(0., END_TIME, 121):
        tau = C_FINAL*np.interp(time, TIMES, LEVELS)/r**2
        delta = np.maximum((np.sqrt(3)*abs(tau)-YIELD)/H-alpha, 0.)
        gamma_p += np.sqrt(3)*np.sign(tau)*delta
        alpha += delta
        _, grad, _, reference_alpha = evaluate(np.column_stack((r, np.zeros_like(r))), time)
        gamma = grad[:, 2]+grad[:, 1]
        maximum_error = max(maximum_error, float(np.max(abs(gamma-(tau/G+gamma_p)))),
                            float(np.max(abs(reference_alpha-alpha))))
    assert maximum_error < 1e-13, maximum_error
    # Verify displacement derivatives through all reversals, including residual u.
    points = np.array([[1.1, .3], [1.4, .2], [1.8, .1]])
    max_gradient_error = 0.
    for time in (0.0007, .001, .002, .0028, .003, .004, .0049, .006):
        exact = evaluate(points, time)[1]
        for axis in range(2):
            delta = np.zeros_like(points); delta[:, axis] = 1e-6
            fd = (evaluate(points+delta, time)[0]-evaluate(points-delta, time)[0])/2e-6
            max_gradient_error = max(max_gradient_error,
                float(np.max(abs(fd-exact[:, [axis, 2+axis]]))))
    assert max_gradient_error < 1e-8, max_gradient_error
    return dict(scalar_history_error=maximum_error, displacement_gradient_error=max_gradient_error)
