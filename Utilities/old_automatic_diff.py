import re
import sympy as sp

# Costanti materiali fissate
E_VALUE = 67000
NU_VALUE = 0.0
SIG_Y_VALUE = 201

# Symbols
x, y = sp.symbols("x y", real=True)
t = sp.symbols("t", real=True)
E = sp.Rational(67000)
nu = sp.Rational(NU_VALUE)
sig_y = sp.Rational(201)

# 1) Displacement field (your choice)
# u = sp.cos(x) * sp.sin(y) * t
# v = sp.cos(y) * sp.sin(x) * t

# u = (x*x+x*y)*t
# v = 0

u = -2*sp.cos(sp.pi*x/4)*sp.sin(sp.pi*y/2)*t
v = 4*sp.sin(sp.pi*x/4)*sp.cos(sp.pi*y/2)*t

# u = (x**2*(x-4)**2*y**2*(y-2)**2)*t
# v = 0

# 2) Strain tensor (small strains, plane strain: eps_zz = 0)
ux = sp.diff(u, x)
uy = sp.diff(u, y)
vx = sp.diff(v, x)
vy = sp.diff(v, y)

eps_xx = ux
eps_yy = vy
eps_xy = (uy + vx) / 2
eps_zz = sp.Integer(0)  # plane strain

# 3) Elastic 3D constitutive law (plane strain)
#    sigma = lambda * tr(eps) * I + 2G * eps
lam = E * nu / ((1 + nu) * (1 - 2 * nu))
G   = E / (2 * (1 + nu))

trace_eps = eps_xx + eps_yy + eps_zz

sigma_xx_el = lam * trace_eps + 2 * G * eps_xx
sigma_yy_el = lam * trace_eps + 2 * G * eps_yy
sigma_zz_el = lam * trace_eps + 2 * G * eps_zz  # here eps_zz = 0, but sigma_zz_el != 0
sigma_xy_el = 2 * G * eps_xy
# sigma_xz_el = 0, sigma_yz_el = 0 in this kinematics

# 4) Mean stress (pressure) and deviatoric part (3D)
trace_sigma_el = sigma_xx_el + sigma_yy_el + sigma_zz_el
p = trace_sigma_el / 3

s_xx_tr = sigma_xx_el - p
s_yy_tr = sigma_yy_el - p
s_zz_tr = sigma_zz_el - p
s_xy_tr = sigma_xy_el         # s_xz_tr = 0, s_yz_tr = 0

# 5) Equivalent von Mises stress (3D J2)
sigma_eq_el_sq = sp.Rational(3, 2) * (
    s_xx_tr**2 + s_yy_tr**2 + s_zz_tr**2
    + 2 * (s_xy_tr**2)   # s_xz_tr = 0, s_yz_tr = 0
)
sigma_eq_el = sp.sqrt(sigma_eq_el_sq)

# 6) Perfect plasticity J2 (radial return, single step, from zero stress)
#    If sigma_eq_el <= sig_y -> purely elastic: sigma = sigma_el
#    If sigma_eq_el > sig_y  -> plastic:
#       s = (sig_y / sigma_eq_el) * s_tr
#       sigma = p * I + s
ratio = sig_y / sigma_eq_el  # only used where sigma_eq_el > sig_y

# Elastic branch
sigma_xx_elastic_branch = sigma_xx_el
sigma_yy_elastic_branch = sigma_yy_el
sigma_xy_elastic_branch = sigma_xy_el
sigma_zz_elastic_branch = sigma_zz_el

# Plastic branch
sigma_xx_plastic_branch = p + ratio * s_xx_tr
sigma_yy_plastic_branch = p + ratio * s_yy_tr
sigma_xy_plastic_branch =       ratio * s_xy_tr
sigma_zz_plastic_branch = p + ratio * s_zz_tr
# (you could also compute sigma_zz_plastic_branch = p + ratio * s_zz_tr if needed)

# Final stress components with perfect plasticity (plane strain J2)
sigma_xx = sp.Piecewise(
    (sigma_xx_elastic_branch, sigma_eq_el <= sig_y),
    (sigma_xx_plastic_branch, True)
)

sigma_yy = sp.Piecewise(
    (sigma_yy_elastic_branch, sigma_eq_el <= sig_y),
    (sigma_yy_plastic_branch, True)
)

sigma_xy = sp.Piecewise(
    (sigma_xy_elastic_branch, sigma_eq_el <= sig_y),
    (sigma_xy_plastic_branch, True)
)

sigma_zz = sp.Piecewise(
    (sigma_zz_elastic_branch, sigma_eq_el <= sig_y),
    (sigma_zz_plastic_branch, True)
)

trace_sigma_pl = sigma_xx + sigma_yy + sigma_zz
p_pl = trace_sigma_pl / 3
s_xx = sigma_xx - p_pl
s_yy = sigma_yy - p_pl
s_zz = sigma_zz - p_pl
s_xy = sigma_xy

sigma_eq_pl_sq = sp.Rational(3, 2) * (
    s_xx**2 + s_yy**2 + s_zz**2 + 2 * s_xy**2
)
sigma_eq_pl = sp.sqrt(sigma_eq_pl_sq)

s_vm = sp.lambdify((x, y, t), sigma_eq_pl, "numpy")
s_xx = sp.lambdify((x, y, t), sigma_xx, "numpy")
s_yy = sp.lambdify((x, y, t), sigma_yy, "numpy")
s_zz = sp.lambdify((x, y, t), sigma_zz, "numpy")
s_xy = sp.lambdify((x, y, t), sigma_xy, "numpy")

# print(s_xx(1.5, 1.9, 0.03))
# print(s_yy(1.5, 1.9, 0.03))
# print(s_zz(1.5, 1.9, 0.03))
# print(s_xy(1.5, 1.9, 0.03))
# print(s_vm(1.5, 1.9, 0.03))
# ef

# 7) Body forces via automatic differentiation
#    b = -div(sigma) in 2D:
#    b_x = -(d sigma_xx/dx + d sigma_xy/dy)
#    b_y = -(d sigma_xy/dx + d sigma_yy/dy)

bx = - (sp.diff(sigma_xx, x) + sp.diff(sigma_xy, y))
by = - (sp.diff(sigma_xy, x) + sp.diff(sigma_yy, y))

# At this point:
# - sigma_xx, sigma_yy, sigma_xy are the J2 perfect-plastic stresses
# - bx, by are the manufactured body forces (symbolic expressions)

# (Optional) Funzioni numeriche con parametri già fissati
bx_fun = sp.lambdify((x, y), bx.subs({}), "numpy")
by_fun = sp.lambdify((x, y), by.subs({}), "numpy")



# ---------------------------
# Tangent: C_tan = d sigma / d eps  (4x4)
# ---------------------------
# ---------------------------
# Voigt stress vector and strain vector
# ---------------------------
sigma_vec = sp.Matrix([
    sigma_xx,
    sigma_yy,
    sigma_xy  # this is tau_xy
])

eps_vec = sp.Matrix([
    eps_xx,
    eps_yy,
    eps_zz
])

sig_vec = sp.Matrix([sigma_xx, sigma_yy, sigma_xy])


# Optional: export as Python code strings
def to_code(expr):
    return sp.printing.pycode(sp.simplify(expr))


def _to_expression_str(expr: sp.Expr) -> str:
    """Return a string expression compatible con il parser di Kratos."""
    expr_str = sp.printing.pycode(sp.simplify(expr))
    return expr_str.replace("math.", "")


def _to_python_code(expr: sp.Expr) -> str:
    """Return a Python expression string (mantiene i prefissi math.)."""
    return sp.printing.pycode(sp.simplify(expr))


def _to_piecewise_str(expr: sp.Expr) -> str:
    """Return a string representation keeping explicit Piecewise, with sp. prefixes."""
    expr_str = sp.sstr(sp.simplify(expr))
    expr_str = re.sub(r"(?<!\.)\bPiecewise\b", "sp.Piecewise", expr_str)
    expr_str = re.sub(r"(?<!\.)\bsin\(", "sp.sin(", expr_str)
    expr_str = re.sub(r"(?<!\.)\bcos\(", "sp.cos(", expr_str)
    expr_str = re.sub(r"(?<!\.)\bAbs\(", "sp.Abs(", expr_str)
    return expr_str

def cpp(expr):
    return sp.cxxcode(sp.simplify(expr))


if __name__ == "__main__":
    bx_str = _to_expression_str(bx)
    by_str = _to_expression_str(by)
    print(f'"value"           : ["{bx_str}", "{by_str}", "0.0"]')
    sigma_xx_str = _to_piecewise_str(sigma_xx)
    sigma_xy_str = _to_piecewise_str(sigma_xy)
    sigma_yy_str = _to_piecewise_str(sigma_yy)
    sigma_zz_str = _to_piecewise_str(sigma_zz)
    print(f"sigma_xx = {sigma_xx_str}")
    print(f"sigma_xy = {sigma_xy_str}")
    print(f"sigma_yy = {sigma_yy_str}")
    print(f"sigma_zz = {sigma_zz_str}")
    
    
    bx_str = _to_piecewise_str(bx)
    by_str = _to_piecewise_str(by)
    
    print(f"bx = {bx_str}")
    print(f"by = {by_str}")
    
    
    for i in range(3):
        print(f"rStressVector[{i}] = {cpp(sig_vec[i])};")
    
    
    
