"""Audit SBM Neumann responses against two independent material-point histories.

Run from any folder with the installed Kratos build. No Newton solve is used:
prescribed trial coefficients isolate material history from global stability.
The current material/mesh parameters are copied; production files are untouched.
"""
import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
import annulus_study as study

KM = study.KM


def strain_matrix(gradients):
    result = np.zeros((3, 2 * len(gradients)))
    result[0, ::2] = gradients[:, 0]
    result[1, 1::2] = gradients[:, 1]
    result[2, ::2] = gradients[:, 1]
    result[2, 1::2] = gradients[:, 0]
    return result


class MaterialPoint:
    def __init__(self, condition, info):
        self.law = condition.Properties[KM.CONSTITUTIVE_LAW].Clone()
        self.values = KM.ConstitutiveLawParameters(
            condition.GetGeometry(), condition.Properties, info)
        self.strain = KM.Vector(3)
        self.stress = KM.Vector(3)
        self.tangent = KM.Matrix(3, 3)
        self.values.SetStrainVector(self.strain)
        self.values.SetStressVector(self.stress)
        self.values.SetConstitutiveMatrix(self.tangent)
        flags = self.values.GetOptions()
        for flag in (KM.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN,
                     KM.ConstitutiveLaw.COMPUTE_STRESS,
                     KM.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR):
            flags.Set(flag, True)
        shape = KM.Vector(np.array(condition.GetGeometry().ShapeFunctionsValues())[0])
        self.law.InitializeMaterial(condition.Properties, condition.GetGeometry(), shape)

    def response(self, strain):
        for i, value in enumerate(strain):
            self.strain[i] = float(value)
        self.law.CalculateMaterialResponseCauchy(self.values)
        return np.array(self.stress), np.array(self.tangent)

    def commit(self, strain):
        self.response(strain)
        self.law.FinalizeMaterialResponseCauchy(self.values)


def operators(condition, skin, degree):
    geometry = condition.GetGeometry()
    center = np.array([geometry.Center().X, geometry.Center().Y])
    # The linealized modeler selects the closest segment CENTER, then takes
    # its first node, rather than selecting the closest skin node directly.
    segments = list(skin.Conditions)
    segment = min(segments, key=lambda c: np.linalg.norm(
        np.array([c.GetGeometry().Center().X, c.GetGeometry().Center().Y]) - center))
    previous = skin.GetCondition(segment.Id - 1 if segment.Id != segments[0].Id
                                 else segments[-1].Id)
    projection = segment.GetGeometry()[0]
    shift = np.array([projection.X, projection.Y]) - center
    normals = []
    for item in (segment, previous):
        first, second = item.GetGeometry()
        tangent = np.array([second.X - first.X, second.Y - first.Y])
        normal = np.array([-tangent[1], tangent[0]])
        normals.append(normal / np.linalg.norm(normal))
    true_normal = sum(normals)
    true_normal /= np.linalg.norm(true_normal)  # inner boundary
    surrogate_normal = np.array(condition.GetValue(KM.NORMAL))[:2]
    gradients = np.array(geometry.ShapeFunctionDerivatives(1, 0))
    jacobian = np.array(geometry.Jacobian(0))[:2, :2]
    # Higher derivatives in the condition are used as physical derivatives.
    # This benchmark has an identity background mapping; assert that premise.
    np.testing.assert_allclose(jacobian, np.eye(2), atol=2e-13)
    extended = gradients.copy()
    if degree == 2:
        second = np.array(geometry.ShapeFunctionDerivatives(2, 0))
        extended[:, 0] += second[:, 0] * shift[0] + second[:, 1] * shift[1]
        extended[:, 1] += second[:, 1] * shift[0] + second[:, 2] * shift[1]
    return (strain_matrix(gradients), strain_matrix(extended),
            np.array(geometry.ShapeFunctionsValues())[0], surrogate_normal,
            true_normal, projection, float(condition.GetValue(KM.INTEGRATION_WEIGHT)))


def traction(stress, normal):
    return np.array([stress[0]*normal[0] + stress[2]*normal[1],
                     stress[2]*normal[0] + stress[1]*normal[1]])


def run(degree, spans, hardening, material_path):
    folder = study.ROOT / 'sbm'
    args = SimpleNamespace(parameters=None, boundary='mixed', degree=degree, dt=None)
    parameters = study.load_parameters(folder, args, [spans, spans])
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(str(material_path))
    with study.shared.working_directory(folder):
        model = KM.Model()
        analysis = study.Analysis(model, parameters)
        analysis.Initialize()
        part = model['IgaModelPart']
        skin = model['skin_model_part'].GetSubModelPart('inner')
        conditions = [c for c in part.Conditions
                      if 'SbmLoadSolidCondition' in study.condition_name(c)]
        # Spread samples around the hole, including different surrogate edges.
        conditions = [conditions[i] for i in np.linspace(0, len(conditions)-1, 8, dtype=int)]
        samples = [(c, operators(c, skin, degree),
                    MaterialPoint(c, part.ProcessInfo), MaterialPoint(c, part.ProcessInfo),
                    MaterialPoint(c, part.ProcessInfo)) for c in conditions]
        scheme = analysis._GetSolver()._GetScheme()
        matrix, increment, rhs = KM.CompressedMatrix(), KM.Vector(), KM.Vector()
        rng = np.random.default_rng(635)
        rows = []
        # Non-proportional loading, unloading, reversal, then reloading.
        for step, (scale, cross) in enumerate(((.008, 0.), (.014, .004),
                                             (.010, .003), (-.008, .006),
                                             (.004, -.009), (.016, .002)), 1):
            part.CloneTimeStep(step * .001)
            scheme.InitializeSolutionStep(part, matrix, increment, rhs)
            for node in part.Nodes:
                x, y = node.X0, node.Y0
                node.SetSolutionStepValue(KM.DISPLACEMENT_X,
                    scale * (x + .35*x*y + .2*y*y) + cross*y)
                node.SetSolutionStepValue(KM.DISPLACEMENT_Y,
                    scale * (-.4*y + .25*x*x) + cross*x*y)
            for condition, data, surrogate, projected, wrong_history in samples:
                B, Bt, N, ns, nt, projection, weight = data
                nodes = list(condition.GetGeometry())
                u = np.array([[n.GetSolutionStepValue(KM.DISPLACEMENT_X),
                               n.GetSolutionStepValue(KM.DISPLACEMENT_Y)] for n in nodes]).ravel()
                stress_s, Ds = surrogate.response(B @ u)
                stress_t, Dt = projected.response(Bt @ u)
                wrong_stress, _ = wrong_history.response(Bt @ u)
                dot = ns @ nt
                force = np.array(projection.GetValue(KM.FORCE))[:2]
                expected = (N[:, None] * (traction(stress_s, ns)
                            + dot*(force-traction(stress_t, nt)))*weight).ravel()
                expected_K = np.empty((len(u), len(u)))
                for j in range(len(u)):
                    column = -traction(Ds @ B[:, j], ns) + dot*traction(Dt @ Bt[:, j], nt)
                    expected_K[:, j] = (N[:, None]*column*weight).ravel()
                K, R = KM.Matrix(), KM.Vector()
                condition.CalculateLocalSystem(K, R, part.ProcessInfo)
                error_R = np.linalg.norm(np.array(R)-expected)/max(1., np.linalg.norm(expected))
                error_K = np.linalg.norm(np.array(K)-expected_K)/max(1., np.linalg.norm(expected_K))
                # Simulate rejected Newton/line-search trial states before
                # returning to the same coefficients. These must not commit.
                direction = rng.normal(size=len(u)); direction /= np.linalg.norm(direction)
                for factor in (1., -2., .5):
                    for n, uv in zip(nodes, (u + factor*.03*direction).reshape(-1, 2)):
                        n.SetSolutionStepValue(KM.DISPLACEMENT_X, float(uv[0]))
                        n.SetSolutionStepValue(KM.DISPLACEMENT_Y, float(uv[1]))
                    trial = KM.Vector()
                    condition.CalculateRightHandSide(trial, part.ProcessInfo)
                for n, uv in zip(nodes, u.reshape(-1, 2)):
                    n.SetSolutionStepValue(KM.DISPLACEMENT_X, float(uv[0]))
                    n.SetSolutionStepValue(KM.DISPLACEMENT_Y, float(uv[1]))
                repeated = KM.Vector()
                condition.CalculateRightHandSide(repeated, part.ProcessInfo)
                repeat_error = np.linalg.norm(np.array(repeated)-np.array(R))/max(1., np.linalg.norm(np.array(R)))
                # Jacobian check after committed multistep plastic history.
                perturbed = []
                for sign in (1, -1):
                    for n, uv in zip(nodes, (u+sign*1e-8*direction).reshape(-1, 2)):
                        n.SetSolutionStepValue(KM.DISPLACEMENT_X, float(uv[0]))
                        n.SetSolutionStepValue(KM.DISPLACEMENT_Y, float(uv[1]))
                    value = KM.Vector()
                    condition.CalculateRightHandSide(value, part.ProcessInfo)
                    perturbed.append(np.array(value))
                for n, uv in zip(nodes, u.reshape(-1, 2)):
                    n.SetSolutionStepValue(KM.DISPLACEMENT_X, float(uv[0]))
                    n.SetSolutionStepValue(KM.DISPLACEMENT_Y, float(uv[1]))
                fd = -(perturbed[0]-perturbed[1]) / 2e-8
                fd_error = np.linalg.norm(fd-np.array(K)@direction)/max(1., np.linalg.norm(fd))
                rows.append(dict(degree=degree, spans=spans, hardening=hardening,
                    step=step, condition=condition.Id, residual_error=float(error_R),
                    tangent_error=float(error_K), repeat_error=float(repeat_error),
                    finite_difference_error=float(fd_error),
                    true_surrogate_strain_distance=float(np.linalg.norm((Bt-B)@u)),
                    reused_history_stress_error=float(np.linalg.norm(wrong_stress-stress_t))))
                surrogate.commit(B @ u)
                projected.commit(Bt @ u)
                wrong_history.commit(B @ u)
            # The real Kratos scheme commits each actual condition, exactly
            # as in a converged solve. The next step tests that committed state.
            scheme.FinalizeSolutionStep(part, matrix, increment, rhs)
        return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=study.ROOT/'diagnostics/projection_history.json')
    args = parser.parse_args()
    KM.ParallelUtilities.SetNumThreads(1)
    rows = []
    material = json.loads(KM.Parameters((study.ROOT/'materials.json').read_text()).PrettyPrintJsonString())
    with TemporaryDirectory() as temporary:
        for hardening in (0., 10000.):
            material['properties'][0]['Material']['Variables']['ISOTROPIC_HARDENING_MODULUS'] = hardening
            path = Path(temporary)/'material.json'
            path.write_text(json.dumps(material))
            for degree in (1, 2):
                for spans in (10, 20):
                    result = run(degree, spans, hardening, path)
                    rows.extend(result)
                    print('HISTORY_AUDIT', json.dumps({
                        'hardening': hardening, 'degree': degree, 'spans': spans,
                        **{key: max(r[key] for r in result) for key in
                           ('residual_error', 'tangent_error', 'repeat_error',
                            'finite_difference_error', 'true_surrogate_strain_distance',
                            'reused_history_stress_error')}}), flush=True)
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(rows, indent=2)+'\n')
    for row in rows:
        assert row['residual_error'] < 1e-9, row
        assert row['tangent_error'] < 1e-9, row
        assert row['repeat_error'] < 1e-12, row
        assert row['finite_difference_error'] < 1e-5, row


if __name__ == '__main__':
    main()
