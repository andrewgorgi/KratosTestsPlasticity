import sympy as sp

# Material parameters
E, nu, sig_y = sp.symbols("E nu sig_y", real=True, positive=True)

# Strains (plane strain, Voigt-3: [xx, yy, xy])
eps_xx, eps_yy, gamma_xy = sp.symbols("eps_xx eps_yy gamma_xy", real=True)
eps_xy = gamma_xy / 2   # tensorial shear strain
eps_zz = sp.Integer(0)  # plane strain

# Elastic 3D parameters
lam = E*nu/((1+nu)*(1-2*nu))
G   = E/(2*(1+nu))

# Elastic stresses
trace_eps = eps_xx + eps_yy + eps_zz

sigma_xx_el = lam*trace_eps + 2*G*eps_xx
sigma_yy_el = lam*trace_eps + 2*G*eps_yy
sigma_zz_el = lam*trace_eps + 2*G*eps_zz
sigma_xy_el = 2*G*eps_xy

# Deviatoric trial
trace_sigma = sigma_xx_el + sigma_yy_el + sigma_zz_el
p = trace_sigma/3

s_xx = sigma_xx_el - p
s_yy = sigma_yy_el - p
s_zz = sigma_zz_el - p
s_xy = sigma_xy_el

sigma_eq = sp.sqrt(sp.Rational(3,2)*(s_xx**2 + s_yy**2 + s_zz**2 + 2*s_xy**2))

ratio = sig_y/sigma_eq

# Final stresses (piecewise)
sigma_xx = sp.Piecewise((sigma_xx_el, sigma_eq <= sig_y),
                        (p + ratio*s_xx, True))

sigma_yy = sp.Piecewise((sigma_yy_el, sigma_eq <= sig_y),
                        (p + ratio*s_yy, True))

sigma_xy = sp.Piecewise((sigma_xy_el, sigma_eq <= sig_y),
                        (ratio*s_xy, True))

# Voigt vectors
sig_vec = sp.Matrix([sigma_xx, sigma_yy, sigma_xy])
eps_vec = sp.Matrix([eps_xx, eps_yy, gamma_xy])

# Tangente consistente (3x3)
C_tan = sig_vec.jacobian(eps_vec)

# --- Stampa C++ ---
def cpp(expr):
    return sp.cxxcode(sp.simplify(expr))

print("Matrix C(3,3);")
for i in range(3):
    for j in range(3):
        print(f"rD({i},{j}) = {cpp(C_tan[i,j])};")
        
        
        
for i in range(3):
    print(f"rStressVector[{i}] = {cpp(sig_vec[i])};")
