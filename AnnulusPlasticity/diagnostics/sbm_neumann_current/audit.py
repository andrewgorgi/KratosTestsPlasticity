"""Independent traction and tensor-product spline checks of installed SBM Neumann."""
import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import BSpline, PPoly

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import cyclic_study as study
from check_projection_history import MaterialPoint

KM = study.KM


def strain_matrix(g):
    b = np.zeros((3, 2*len(g)))
    b[0, 0::2] = g[:, 0]
    b[1, 1::2] = g[:, 1]
    b[2, 0::2] = g[:, 1]
    b[2, 1::2] = g[:, 0]
    return b


def traction(n):
    return np.array([[n[0], 0., n[1]], [0., n[1], n[0]]])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, default=3)
    parser.add_argument('--spans', type=int, default=20)
    args = parser.parse_args()
    p, spans = args.degree, args.spans
    out = HERE/f'audit_p{p}_n{spans}'
    out.mkdir(exist_ok=True)
    KM.ParallelUtilities.SetNumThreads(1)
    parameters = KM.Parameters((HERE/'parameters_input.json').read_text())
    for modeler in parameters['modelers'].values():
        settings = modeler['Parameters']
        if settings.Has('input_filename'):
            settings['input_filename'].SetString(str((ROOT/'sbm'/settings['input_filename'].GetString()).resolve()))
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(str(HERE/'materials_input.json'))
    (out/'input.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    settings = SimpleNamespace(parameters='input.json', boundary='mixed', degree=p, dt=None, echo_level=0)
    parameters, _ = study.load_parameters(out, settings, [spans]*2)
    with study.base.shared.working_directory(out):
        model = KM.Model()
        analysis = study.base.Analysis(model, parameters)
        analysis.EvaluateSolution = lambda points, time: study.exact.evaluate(points, time)
        analysis.BoundaryTraction = lambda points, time: study.exact.boundary_traction(points, time)
        analysis.Initialize()
        part, skin = model['IgaModelPart'], model['skin_model_part']
        analysis.time = .00025
        part.CloneTimeStep(.00025)
        analysis.InitializeSolutionStep()
        conditions = [c for c in part.Conditions if 'SbmLoadSolidCondition' in study.base.condition_name(c)]
        def projection_nodes(part):
            result = list(part.Nodes) if part.Name=='projection_data' else []
            for child in part.SubModelParts:
                result += projection_nodes(child)
            return result
        nodes = sorted(projection_nodes(skin.GetSubModelPart('inner')), key=lambda n:n.Id)
        conditions.sort(key=lambda c:c.Id)
        assert len(nodes)==len(conditions), (len(nodes), len(conditions))
        breaks = np.linspace(-2.2, 2.2, spans+1)
        knots = np.r_[[breaks[0]]*p, breaks, [breaks[-1]]*p]
        count = spans+p
        greville = np.array([knots[i+1:i+p+1].mean() for i in range(count)])
        polynomials = [PPoly.from_spline(BSpline(knots, np.eye(count)[i], p)) for i in range(count)]
        D = np.diag([67000., 67000., 33500.])
        cache, rows = [], []
        rng = np.random.default_rng(471)
        # The modeler creates one projection node per condition, in matching ID order.
        # Agreement of the independently reconstructed K/R checks this correspondence.
        for condition, node in zip(conditions, nodes):
            geometry = condition.GetGeometry()
            coordinates = np.array([[n.X0, n.Y0] for n in geometry])
            center = np.array([geometry.Center().X, geometry.Center().Y])
            true = center/np.linalg.norm(center)
            projection = np.array([node.X0, node.Y0])
            distance = np.linalg.norm(projection-true)
            nt = -projection/np.linalg.norm(projection)
            ns = np.array(condition.GetValue(KM.NORMAL))[:2]
            weight = float(condition.GetValue(KM.INTEGRATION_WEIGHT))
            N = np.array(geometry.ShapeFunctionsValues())[0]
            g = np.array(geometry.ShapeFunctionDerivatives(1, 0))[:, :2]
            np.testing.assert_allclose(np.array(geometry.Jacobian(0))[:2, :2], np.eye(2), atol=1e-12)
            ij = np.argmin(abs(coordinates[:, :, None]-greville), axis=2)
            np.testing.assert_allclose(greville[ij], coordinates, atol=1e-12)
            span_x, span_y = ij.min(axis=0)+p

            def direct(point):
                values, gradients = [], []
                for i, j in ij:
                    x = point[0]-knots[span_x]
                    y = point[1]-knots[span_y]
                    a, b = polynomials[i].c[:, span_x], polynomials[j].c[:, span_y]
                    vx, vy = np.polyval(a, x), np.polyval(b, y)
                    values.append(vx*vy)
                    gradients.append([np.polyval(np.polyder(a), x)*vy, vx*np.polyval(np.polyder(b), y)])
                return np.array(values), np.array(gradients)

            Ns, gs = direct(center)
            _, gt = direct(projection)
            shift = projection-center
            expansions = {}
            for order in (p, 2*p):
                gradient = g.copy()
                for n in range(2, order+1):
                    derivative = np.array(geometry.ShapeFunctionDerivatives(n, 0))
                    for k in range(n):
                        coefficient = shift[0]**(n-1-k)*shift[1]**k/(factorial(n-1-k)*factorial(k))
                        gradient[:, 0] += derivative[:, k]*coefficient
                        gradient[:, 1] += derivative[:, k+1]*coefficient
                expansions[order] = gradient
            Bs, Bt = strain_matrix(g), strain_matrix(gt)
            dot = float(ns@nt)
            expectedK = (N[:, None, None]*(-traction(ns)@D@Bs+dot*traction(nt)@D@Bt)[None, :, :]*weight).reshape(2*len(N), -1)
            K, R = KM.Matrix(), KM.Vector()
            condition.CalculateLocalSystem(K, R, part.ProcessInfo)
            error = float(np.linalg.norm(np.array(K)-expectedK)/max(np.linalg.norm(expectedK), 1e-30))
            row = dict(id=condition.Id, projection_distance=float(distance), normal_dot=float(dot),
                       true_normal_error=float(np.linalg.norm(np.array(node.GetValue(KM.NORMAL))[:2]-nt)),
                       surrogate_basis_error=float(np.max(abs(N-Ns))),
                       surrogate_gradient_error=float(np.linalg.norm(g-gs)/np.linalg.norm(gs)),
                       truncated_gradient_error=float(np.linalg.norm(expansions[p]-gt)/np.linalg.norm(gt)),
                       complete_gradient_error=float(np.linalg.norm(expansions[2*p]-gt)/np.linalg.norm(gt)),
                       stiffness_direct_error=error)
            cache.append((condition, node, N, Bs, Bt, ns, nt, weight))
            rows.append(row)
        load_rows = []
        for time in (.0005, .001, .0015, .003, .004, .005, .006):
            analysis.time = time
            part.CloneTimeStep(time)
            analysis.InitializeSolutionStep()
            C = study.exact.C_FINAL*float(np.interp(time, study.exact.TIMES, study.exact.LEVELS))
            force_error, rhs_error, affine_error, random_error = 0., 0., 0., 0.
            for condition, node, N, Bs, Bt, ns, nt, weight in cache:
                point = np.array([node.X0, node.Y0])
                r = np.linalg.norm(point)
                # Independent polar traction: n=-e_r, sigma_rtheta=C/r^2.
                force = C*np.array([point[1], -point[0]])/r**3
                actual_force = np.array(node.GetValue(KM.FORCE))[:2]
                force_error = max(force_error, float(np.linalg.norm(force-actual_force)))
                rhs = KM.Vector()
                condition.CalculateRightHandSide(rhs, part.ProcessInfo)
                expected = (N[:, None]*force[None, :]*(ns@nt)*weight).ravel()
                rhs_error = max(rhs_error, float(np.linalg.norm(np.array(rhs)-expected)))
                if time == .001:
                    geometry = condition.GetGeometry()
                    coords = np.array([[n.X0, n.Y0] for n in geometry])
                    affine = 1e-6*np.column_stack((coords[:, 0]+.3*coords[:, 1], -.2*coords[:, 0]+.4*coords[:, 1])).ravel()
                    for kind, u in [('affine', affine), ('random', rng.normal(size=2*len(N))*1e-8)]:
                        for n, uv in zip(geometry, u.reshape(-1, 2)):
                            n.SetSolutionStepValue(KM.DISPLACEMENT_X, float(uv[0]))
                            n.SetSolutionStepValue(KM.DISPLACEMENT_Y, float(uv[1]))
                        condition.CalculateRightHandSide(rhs, part.ProcessInfo)
                        stress_s, stress_t = D@Bs@u, D@Bt@u
                        expected = (N[:, None]*(traction(ns)@stress_s+(ns@nt)*(force-traction(nt)@stress_t))[None, :]*weight).ravel()
                        err = float(np.linalg.norm(np.array(rhs)-expected))
                        if kind=='affine': affine_error=max(affine_error, err)
                        else: random_error=max(random_error, err)
                    for n in geometry:
                        n.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.)
                        n.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.)
            load_rows.append(dict(time=time, load_factor=float(C/study.exact.C_FINAL),
                max_force_error=force_error, max_zero_rhs_error=rhs_error,
                max_affine_rhs_error=affine_error, max_random_rhs_error=random_error))
        summary = dict(degree=p, spans=spans, conditions=len(rows), load_checks=load_rows,
                       maxima={key:max(row[key] for row in rows) for key in rows[0] if key!='id'},
                       minimum_normal_dot=min(row['normal_dot'] for row in rows))
        # Isolate history handling from Newton and from Taylor reconstruction:
        # independent material points receive direct polynomial gradients.
        chosen = [cache[i] for i in np.linspace(0, len(cache)-1, 12, dtype=int)]
        samples = [(item, MaterialPoint(item[0], part.ProcessInfo),
                    MaterialPoint(item[0], part.ProcessInfo)) for item in chosen]
        for node in skin.Nodes:
            node.SetValue(KM.FORCE, [0., 0., 0.])
        history = []
        scheme = analysis._GetSolver()._GetScheme()
        system, increment, system_rhs = KM.CompressedMatrix(), KM.Vector(), KM.Vector()
        for step, (scale, cross) in enumerate(((.008, 0.), (.014, .004), (.004, .001),
                (0., 0.), (-.008, .006), (-.016, -.004), (0., 0.), (.016, .002)), 1):
            part.CloneTimeStep(.006+step*.001)
            scheme.InitializeSolutionStep(part, system, increment, system_rhs)
            for node in part.Nodes:
                x, y = node.X0, node.Y0
                node.SetSolutionStepValue(KM.DISPLACEMENT_X, scale*(x+.35*x*y+.2*y*y)+cross*y)
                node.SetSolutionStepValue(KM.DISPLACEMENT_Y, scale*(-.4*y+.25*x*x)+cross*x*y)
            for item, surrogate, projected in samples:
                condition, _, N, Bs, Bt, ns, nt, weight = item
                geometry = condition.GetGeometry()
                u = np.array([[n.GetSolutionStepValue(KM.DISPLACEMENT_X),
                               n.GetSolutionStepValue(KM.DISPLACEMENT_Y)] for n in geometry]).ravel()
                ss, Ds = surrogate.response(Bs@u)
                st, Dt = projected.response(Bt@u)
                expectedR = (N[:, None]*(traction(ns)@ss-(ns@nt)*traction(nt)@st)[None, :]*weight).ravel()
                expectedK = (N[:, None, None]*(-traction(ns)@Ds@Bs+(ns@nt)*traction(nt)@Dt@Bt)[None, :, :]*weight).reshape(len(u), len(u))
                K, R = KM.Matrix(), KM.Vector()
                condition.CalculateLocalSystem(K, R, part.ProcessInfo)
                K, R = np.array(K), np.array(R)
                direction = rng.normal(size=len(u)); direction /= np.linalg.norm(direction)
                def response(v):
                    for node, uv in zip(geometry, v.reshape(-1, 2)):
                        node.SetSolutionStepValue(KM.DISPLACEMENT_X, float(uv[0]))
                        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, float(uv[1]))
                    rhs = KM.Vector()
                    condition.CalculateRightHandSide(rhs, part.ProcessInfo)
                    return np.array(rhs)
                for factor in (1., -2., .5):
                    response(u+factor*.03*direction)
                repeated = response(u)
                fd = -(response(u+1e-8*direction)-response(u-1e-8*direction))/2e-8
                response(u)
                relative = lambda difference, base:float(np.linalg.norm(difference)/max(1., np.linalg.norm(base)))
                history.append(dict(step=step, condition=condition.Id,
                    residual_error=relative(R-expectedR, expectedR),
                    tangent_error=relative(K-expectedK, expectedK),
                    trial_repeat_error=relative(repeated-R, R),
                    finite_difference_error=relative(K@direction-fd, fd),
                    true_surrogate_strain_distance=float(np.linalg.norm((Bt-Bs)@u))))
                surrogate.commit(Bs@u)
                projected.commit(Bt@u)
            scheme.FinalizeSolutionStep(part, system, increment, system_rhs)
        summary['history_checks'] = dict(samples=len(history),
            maxima={key:max(row[key] for row in history) for key in history[0] if key not in ('step', 'condition')})
        (out/'history.json').write_text(json.dumps(history, indent=2)+'\n')
        (out/'conditions.json').write_text(json.dumps(rows, indent=2)+'\n')
        (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
        print('NEUMANN_AUDIT', json.dumps(summary), flush=True)
        assert summary['maxima']['stiffness_direct_error'] < 1e-8
        assert max(r['max_force_error'] for r in load_rows) < 1e-8
        assert max(r['max_random_rhs_error'] for r in load_rows) < 1e-8
        assert summary['history_checks']['maxima']['residual_error'] < 1e-8
        assert summary['history_checks']['maxima']['tangent_error'] < 1e-8
        assert summary['history_checks']['maxima']['trial_repeat_error'] < 1e-10
        assert summary['history_checks']['maxima']['finite_difference_error'] < 1e-4


if __name__=='__main__':
    main()
