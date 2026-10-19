"""Cyclic/non-proportional audit of the two unmodified SBM macro branches.

Launch in separate processes with LD_PRELOAD set to the corresponding library.
Independent material points receive exactly B_s*u and B_t*u, isolating history
errors from Taylor truncation and from the nonlinear solver.
"""
import argparse
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import annulus_study as study
from check_projection_history import MaterialPoint, strain_matrix, traction

KM = study.KM


def material(law):
    data = json.loads(KM.Parameters((HERE/'materials_input.json').read_text()).PrettyPrintJsonString())
    if law == 'saturation':
        variables = data['properties'][0]['Material']['Variables']
        variables['EXPONENTIAL_SATURATION_YIELD_STRESS'] = 350.
        variables['HARDENING_EXPONENT'] = 80.
    path = HERE/f'material_{law}.json'
    temporary = path.with_suffix(f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    temporary.replace(path)
    return path


def operators(condition, degree, skin):
    geometry = condition.GetGeometry()
    # NEIGHBOUR_NODES is not Python-bound. Identify the stored projection using
    # the exact radial position for this circle. Agreement with the assembled
    # response below independently checks both the selected node and its normal.
    center = np.array([geometry.Center().X, geometry.Center().Y])
    radial = center/np.linalg.norm(center)
    # Use the stored projected node, including the modeler's numerical normal.
    # Its distance to the exact radial position is checked, and the condition
    # response check detects any wrong node selected among nearby candidates.
    node = min(skin.Nodes, key=lambda n: (n.X-radial[0])**2+(n.Y-radial[1])**2)
    projection = np.array([node.X, node.Y])
    assert np.linalg.norm(projection-radial) < 1e-4
    shift = projection-center
    ns = np.array(condition.GetValue(KM.NORMAL))[:2]
    nt = np.array(node.GetValue(KM.NORMAL))[:2]
    gradients = np.array(geometry.ShapeFunctionDerivatives(1, 0))
    np.testing.assert_allclose(np.array(geometry.Jacobian(0))[:2, :2], np.eye(2), atol=2e-13)
    extended = gradients.copy()
    if degree == 2:
        second = np.array(geometry.ShapeFunctionDerivatives(2, 0))
        extended[:, 0] += second[:, 0]*shift[0] + second[:, 1]*shift[1]
        extended[:, 1] += second[:, 1]*shift[0] + second[:, 2]*shift[1]
    return (strain_matrix(gradients), strain_matrix(extended),
            np.array(geometry.ShapeFunctionsValues())[0], ns, nt,
            float(condition.GetValue(KM.INTEGRATION_WEIGHT)))


def set_u(nodes, displacement):
    for node, uv in zip(nodes, displacement.reshape(-1, 2)):
        node.SetSolutionStepValue(KM.DISPLACEMENT_X, float(uv[0]))
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, float(uv[1]))


def run(mode, degree, spans, law):
    settings = SimpleNamespace(parameters=str(HERE/'parameters_input.json'),
                               boundary='mixed', degree=degree, dt=None)
    parameters = study.load_parameters(ROOT/'sbm', settings, [spans]*2)
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(str(material(law)))
    with study.shared.working_directory(ROOT/'sbm'):
        model = KM.Model()
        analysis = study.Analysis(model, parameters)
        analysis.Initialize()
        part = model['IgaModelPart']
        conditions = [c for c in part.Conditions
                      if 'SbmLoadSolidCondition' in study.condition_name(c)]
        conditions = [conditions[i] for i in np.linspace(0, len(conditions)-1, 8, dtype=int)]
        skin = model['skin_model_part']
        samples = [(c, operators(c, degree, skin), MaterialPoint(c, part.ProcessInfo),
                    MaterialPoint(c, part.ProcessInfo), MaterialPoint(c, part.ProcessInfo))
                   for c in conditions]
        scheme = analysis._GetSolver()._GetScheme()
        matrix, increment, rhs = KM.CompressedMatrix(), KM.Vector(), KM.Vector()
        rng = np.random.default_rng(635)
        rows = []
        # Initial plastic loading, unloading, rotation, reversal, reloading.
        path = ((.008, 0.), (.014, .004), (.004, .001), (0., 0.),
                (-.008, .006), (-.016, -.004), (0., 0.), (.016, .002))
        for step, (scale, cross) in enumerate(path, 1):
            part.CloneTimeStep(step*.001)
            scheme.InitializeSolutionStep(part, matrix, increment, rhs)
            for node in part.Nodes:
                x, y = node.X0, node.Y0
                node.SetSolutionStepValue(KM.DISPLACEMENT_X,
                    scale*(x+.35*x*y+.2*y*y)+cross*y)
                node.SetSolutionStepValue(KM.DISPLACEMENT_Y,
                    scale*(-.4*y+.25*x*x)+cross*x*y)
            for condition, data, surrogate, projected, reused in samples:
                B, Bt, N, ns, nt, weight = data
                nodes = list(condition.GetGeometry())
                u = np.array([[n.GetSolutionStepValue(KM.DISPLACEMENT_X),
                               n.GetSolutionStepValue(KM.DISPLACEMENT_Y)] for n in nodes]).ravel()
                ss, Ds = surrogate.response(B@u)
                st, Dt = projected.response(Bt@u)
                sw, Dw = reused.response(Bt@u)
                dot = ns@nt
                # Zero imposed traction suffices to audit the constitutive term.
                expected = (N[:, None]*(traction(ss, ns)-dot*traction(st, nt))*weight).ravel()
                reused_R = (N[:, None]*(traction(ss, ns)-dot*traction(sw, nt))*weight).ravel()
                actual_stress = st if mode == 'on' else sw
                actual_D = Dt if mode == 'on' else Ds
                branch_R = expected if mode == 'on' else reused_R

                def jacobian(D):
                    return np.column_stack([(N[:, None]*(-traction(Ds@B[:, j], ns)
                        + dot*traction(D@Bt[:, j], nt))*weight).ravel() for j in range(len(u))])

                expected_K, branch_K = jacobian(Dt), jacobian(actual_D)
                K, R = KM.Matrix(), KM.Vector()
                condition.CalculateLocalSystem(K, R, part.ProcessInfo)
                K, R = np.array(K), np.array(R)
                normalize = lambda difference, base: float(np.linalg.norm(difference)/max(1., np.linalg.norm(base)))
                # Normalize physical differences by the individual traction terms,
                # avoiding artificially large ratios from their near cancellation.
                flux_scale = max(1., np.linalg.norm((N[:, None]*traction(ss, ns)*weight).ravel()),
                                 np.linalg.norm((N[:, None]*traction(st, nt)*weight).ravel()))
                direction = rng.normal(size=len(u)); direction /= np.linalg.norm(direction)
                for factor in (1., -2., .5):
                    set_u(nodes, u+factor*.03*direction)
                    trial = KM.Vector()
                    condition.CalculateRightHandSide(trial, part.ProcessInfo)
                set_u(nodes, u)
                repeated = KM.Vector()
                condition.CalculateRightHandSide(repeated, part.ProcessInfo)
                perturbed = []
                for sign in (1, -1):
                    set_u(nodes, u+sign*1e-8*direction)
                    trial = KM.Vector()
                    condition.CalculateRightHandSide(trial, part.ProcessInfo)
                    perturbed.append(np.array(trial))
                set_u(nodes, u)
                fd = -(perturbed[0]-perturbed[1])/2e-8
                # R_s and R_t are also checked independently through their
                # tangents, not by assuming the off-branch Jacobian is correct.
                derivative_K = jacobian(Dt if mode == 'on' else Dw)
                rows.append(dict(mode=mode, degree=degree, spans=spans, law=law,
                    step=step, condition=condition.Id,
                    branch_residual_error=normalize(R-branch_R, branch_R),
                    branch_tangent_error=normalize(K-branch_K, branch_K),
                    physical_residual_error=normalize(R-expected, expected),
                    physical_flux_error=float(np.linalg.norm(R-expected)/flux_scale),
                    physical_tangent_error=normalize(K-expected_K, expected_K),
                    finite_difference_error=normalize(K@direction-fd, fd),
                    exact_derivative_fd_error=normalize(derivative_K@direction-fd, fd),
                    repeat_error=normalize(np.array(repeated)-R, R),
                    stress_history_relative_error=normalize(actual_stress-st, st),
                    stress_history_absolute_error=float(np.linalg.norm(actual_stress-st)),
                    true_surrogate_strain_distance=float(np.linalg.norm((Bt-B)@u))))
                surrogate.commit(B@u)
                projected.commit(Bt@u)
                reused.commit(B@u)
            scheme.FinalizeSolutionStep(part, matrix, increment, rhs)
        return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('on', 'off'), required=True)
    parser.add_argument('--spans', nargs='+', type=int, default=[20, 40])
    args = parser.parse_args()
    KM.ParallelUtilities.SetNumThreads(1)
    rows = []
    for law in ('linear', 'saturation'):
        for degree in (1, 2):
            for spans in args.spans:
                data = run(args.mode, degree, spans, law)
                rows.extend(data)
                print('LOCAL_CYCLIC', json.dumps({k: data[0][k] for k in ('mode', 'law', 'degree', 'spans')}), flush=True)
                (HERE/f'local_{args.mode}.json').write_text(json.dumps(rows, indent=2)+'\n')
                for row in data:
                    assert row['branch_residual_error'] < 1e-8, row
                    assert row['branch_tangent_error'] < 1e-8, row
                    assert row['repeat_error'] < 1e-10, row
                    # The nonlinear material itself is audited separately; a
                    # failed FD check is recorded rather than hiding later cases.
                    if law == 'linear':
                        assert row['exact_derivative_fd_error'] < 1e-4, row
                    if args.mode == 'on':
                        assert row['physical_tangent_error'] < 1e-8, row
                        assert row['physical_residual_error'] < 1e-8, row


if __name__ == '__main__':
    main()
