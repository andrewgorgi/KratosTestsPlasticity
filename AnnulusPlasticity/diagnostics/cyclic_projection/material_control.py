"""Independent virgin pure-shear control of the installed J2 material law."""
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import annulus_study as study
from check_projection_history import MaterialPoint
from local import material
import reference


def main():
    study.KM.ParallelUtilities.SetNumThreads(1)
    rows = []
    reference_error = 0.
    for law in ('linear', 'saturation'):
        settings = SimpleNamespace(parameters=str(HERE/'parameters_input.json'),
                                   boundary='mixed', degree=1, dt=None)
        parameters = study.load_parameters(study.ROOT/'sbm', settings, [20]*2)
        path = material(law)
        parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(str(path))
        v = json.loads(path.read_text())['properties'][0]['Material']['Variables']
        G = v['YOUNG_MODULUS']/(2*(1+v['POISSON_RATIO']))
        sigma0, H = v['YIELD_STRESS'], v['ISOTROPIC_HARDENING_MODULUS']
        delta, b = v['EXPONENTIAL_SATURATION_YIELD_STRESS']-sigma0, v['HARDENING_EXPONENT']
        yield_stress = lambda alpha: sigma0+H*alpha+delta*(1-np.exp(-b*alpha))
        with study.shared.working_directory(study.ROOT/'sbm'):
            model = study.KM.Model()
            analysis = study.Analysis(model, parameters)
            analysis.Initialize()
            part = model['IgaModelPart']
            condition = next(c for c in part.Conditions if 'SbmLoadSolidCondition' in study.condition_name(c))
            for gamma in (.01, .02, .04):
                point = MaterialPoint(condition, part.ProcessInfo)
                strain = np.array([0., 0., gamma])
                stress, D = point.response(strain)
                alpha = brentq(lambda a: np.sqrt(3)*G*gamma-3*G*a-yield_stress(a),
                               0., gamma/np.sqrt(3), xtol=1e-15)
                tau_exact = G*(gamma-np.sqrt(3)*alpha)
                alpha_actual = (G*gamma-stress[2])/(np.sqrt(3)*G)
                fd_errors = []
                for h in (1e-6, 1e-7, 1e-8):
                    plus, minus = strain.copy(), strain.copy()
                    plus[2] += h; minus[2] -= h
                    fd = (point.response(plus)[0]-point.response(minus)[0])/(2*h)
                    fd_errors.append(dict(h=h, derivative=float(fd[2]),
                        tangent_relative_error=float(np.linalg.norm(D[:, 2]-fd)/max(1., np.linalg.norm(fd)))))
                rows.append(dict(law=law, gamma=gamma, tau_actual=float(stress[2]),
                    tau_exact=float(tau_exact), alpha_actual=float(alpha_actual), alpha_exact=float(alpha),
                    yield_residual=float(np.sqrt(3)*abs(stress[2])-yield_stress(alpha_actual)),
                    shear_tangent=float(D[2, 2]), finite_differences=fd_errors))
            if law == 'linear':
                # Full installed constitutive law, driven by the independently
                # constructed cyclic reference at six radii, including reversal.
                radii = np.array([.95, 1., 1.3, 1.6, 1.9, 2.1])
                points = [MaterialPoint(condition, part.ProcessInfo) for r in radii]
                xy = np.column_stack((radii, np.zeros_like(radii)))
                for time in np.linspace(0., reference.END_TIME, 121):
                    _, gradient, stress, _ = reference.evaluate(xy, time)
                    for index, point in enumerate(points):
                        strain = np.array([0., 0., gradient[index, 1]+gradient[index, 2]])
                        actual = point.response(strain)[0]
                        reference_error = max(reference_error, float(np.max(abs(actual-stress[index]))))
                        point.commit(strain)
    assert reference_error < 1e-8, reference_error
    (HERE/'material_control.json').write_text(json.dumps(rows, indent=2)+'\n')
    (HERE/'reference_validation.json').write_text(json.dumps({
        **reference.verify(), 'kratos_cyclic_stress_absolute_error': reference_error,
        'radii': radii.tolist(), 'steps': 121}, indent=2)+'\n')
    print('MATERIAL_CONTROL', json.dumps(rows), flush=True)


if __name__ == '__main__':
    main()
