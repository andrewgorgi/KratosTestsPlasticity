"""Campi fabbricati coerenti: spostamento, gradiente, stress e forza di volume."""
from functools import lru_cache
import sympy as sp
import numpy as np

FIELDS = ('linear', 'polynomial', 'cosine_middle', 'cosine_boundary')
PHYSICS = ('plastic', 'elastic', 'plastic_hardening')
x, y, t = sp.symbols('x y t', real=True)
E, NU, YIELD = 67000, 0, 201
HARDENING = 10000  # ISOTROPIC_HARDENING_MODULUS for plastic_hardening

@lru_cache(None)
def expressions(field, physics):
    if field == 'linear':
        u = sp.Matrix([x*t, 0])
    elif field == 'polynomial':
        u = sp.Matrix([x*(x-4)*y*(y-2)*t, 0])
    elif field in ('cosine_middle', 'cosine_boundary'):
        phase = y if field == 'cosine_middle' else y+1
        u = sp.Matrix([-2*sp.cos(sp.pi*x/4)*sp.sin(sp.pi*phase/2)*t,
                       4*sp.sin(sp.pi*x/4)*sp.cos(sp.pi*phase/2)*t])
    else:
        raise ValueError(field)
    if physics not in PHYSICS:
        raise ValueError(physics)
    gradient = u.jacobian([x, y])
    # Piccole deformazioni, plane strain, nu=0: sigma_zz elastico = 0.
    elastic = sp.Matrix([E*gradient[0, 0], E*gradient[1, 1],
                         E*(gradient[0, 1]+gradient[1, 0])/2, 0])
    stress = elastic
    if physics in ('plastic', 'plastic_hardening'):
        # J2 con hardening isotropo lineare, carico proporzionale monotono t>=0.
        hardening = HARDENING if physics == 'plastic_hardening' else 0
        three_g = sp.Rational(3 * E, 2 * (1 + NU))
        linear_ratio = hardening / (three_g + hardening)
        yield_factor = three_g * YIELD / (three_g + hardening)
        trial = elastic.diff(t)
        mean = (trial[0]+trial[1]+trial[3])/3
        dev = trial-sp.Matrix([mean, mean, 0, mean])
        q_squared = sp.Rational(3, 2)*(dev[0]**2+dev[1]**2+dev[3]**2+2*dev[2]**2)
        q = sp.Pow(q_squared, sp.Rational(1,2), evaluate=False)
        plastic = (sp.Matrix([mean*t, mean*t, 0, mean*t])
                   + linear_ratio*t*dev + yield_factor*dev/q)
        stress = sp.Matrix([sp.Piecewise((a, t*q <= YIELD), (b, True))
                            for a, b in zip(elastic, plastic)])
    body = sp.Matrix([-sp.diff(elastic[0], x)-sp.diff(elastic[2], y),
                      -sp.diff(elastic[2], x)-sp.diff(elastic[1], y), 0])
    if physics in ('plastic', 'plastic_hardening'):
        # Derivata analitica di dev/sqrt(q_squared), senza derivare Abs/sign
        # introdotti dalle semplificazioni di SymPy. Il ramo plastico non vale per q=0.
        hydro = sp.Matrix([mean,mean,0,mean])
        ds_dx = (hydro.diff(x)*t + linear_ratio*t*dev.diff(x)
                 + yield_factor*(dev.diff(x)/q-dev*q_squared.diff(x)/(2*q**3)))
        ds_dy = (hydro.diff(y)*t + linear_ratio*t*dev.diff(y)
                 + yield_factor*(dev.diff(y)/q-dev*q_squared.diff(y)/(2*q**3)))
        bp = sp.Matrix([-ds_dx[0]-ds_dy[2], -ds_dx[2]-ds_dy[1],0])
        body = sp.Matrix([sp.Piecewise((a,t*q<=YIELD),(b,True)) for a,b in zip(body,bp)])
    return u, gradient, stress, body

def expression_string(expr):
    from sympy.printing.pycode import pycode
    expr = sp.piecewise_fold(expr)
    if isinstance(expr, sp.Piecewise):
        if len(expr.args) != 2 or expr.args[-1][1] != True:
            raise ValueError('Il parser Kratos richiede una condizione con due rami.')
        (first, condition), (second, _) = expr.args
        comparison = f'({pycode(condition.lhs)}) {condition.rel_op} ({pycode(condition.rhs)})'
        text = f'({pycode(first)}) if {comparison} else ({pycode(second)})'
    else:
        text = pycode(expr)
    return text.replace('math.', '')

@lru_cache(None)
def evaluator(field, physics):
    u, grad, stress, _ = expressions(field, physics)
    funcs = [sp.lambdify((x, y, t), list(values), 'numpy')
             for values in (u, grad, stress)]
    def evaluate(points, time):
        points = np.asarray(points)
        with np.errstate(divide='ignore', invalid='ignore'):
            return tuple(np.array([np.broadcast_to(v, (len(points),)) for v in fn(points[:,0], points[:,1], time)]).T
                         for fn in funcs)
    return evaluate
