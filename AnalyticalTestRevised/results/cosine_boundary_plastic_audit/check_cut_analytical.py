"""Verifica nativa di stress, tangente e output analitici degli operatori Cut-SBM."""
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import study

KM, IGA = study.KM, study.IGA
KM.ParallelUtilities.SetNumThreads(4)


def radial_stress(strain, plastic):
    """Proiezione radiale 3D indipendente, con epsilon_zz=0 e nu=0."""
    trial = 67000. * np.array([strain[0], strain[1], strain[2] / 2., 0.])
    mean = (trial[0] + trial[1] + trial[3]) / 3.
    dev = trial - np.array([mean, mean, 0., mean])
    q = np.sqrt(1.5 * (dev[0]**2 + dev[1]**2 + dev[3]**2 + 2*dev[2]**2))
    if plastic and q > 201.:
        trial = np.array([mean, mean, 0., mean]) + 201. * dev / q
    return trial[:3]


def local_system(entity, info):
    matrix, rhs = KM.Matrix(), KM.Vector()
    entity.CalculateLocalSystem(matrix, rhs, info)
    return np.array(matrix), np.array(rhs)


def set_affine(nodes, strain):
    for node in nodes:
        node.SetSolutionStepValue(KM.DISPLACEMENT_X, strain[0]*node.X0 + strain[2]*node.Y0/2.)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, strain[1]*node.Y0 + strain[2]*node.X0/2.)


rows = []
with tempfile.TemporaryDirectory(prefix='cut_analytical_check_', dir='/tmp') as temporary:
    root = Path(temporary)
    (root / 'geometry').symlink_to(ROOT / 'geometry', target_is_directory=True)
    for physics in ('plastic', 'elastic'):
        (root / f'materials_{physics}.json').symlink_to(ROOT / f'materials_{physics}.json')
        folder = root / physics
        folder.mkdir()
        (folder / 'txt_files').mkdir()
        parameters = study.load_parameters(ROOT / 'gap_sbm', SimpleNamespace(
            parameters=None, field='cosine_boundary' if physics == 'plastic' else 'linear',
            physics=physics, spans=[8,4], degree=2))
        # La prova locale non richiede la quadratura fitta del bulk.
        bulk = parameters['modelers'][2]['Parameters']['element_condition_list'][0]
        if bulk.Has('number_of_integration_points_per_span'):
            bulk['number_of_integration_points_per_span'].SetInt(2)
        else:
            bulk.AddInt('number_of_integration_points_per_span',2)
        with study.working_directory(folder):
            model = KM.Model()
            analysis = study.Analysis(model, parameters)
            analysis.Initialize()
            part = model['IgaModelPart']
            info = part.ProcessInfo
            nodes = list(part.Nodes)
            node_map = {node.Id:node for node in nodes}
            cuts = list(part.GetSubModelPart('GapElements').Elements)
            entity_groups = {'CutSbmSolidElement':cuts}
            for name in ('CutSbmSolidCondition', 'CutSbmSolidInterfaceCondition'):
                entity_groups[name] = [c for c in part.Conditions if str(c).split(' #')[0].strip('"') == name]
            representatives = {name:max(entities, key=lambda e:abs(e.GetValue(KM.INTEGRATION_WEIGHT)))
                               for name,entities in entity_groups.items()}
            states = ([1.e-5,-4.e-5,2.e-5], [.001963495408493621,-.007853981633974483,0.],
                      [.003926990816987242,-.015707963267948967,0.], [.002,-.006,.0013])
            for step, strain in enumerate(states):
                strain = np.array(strain)
                if physics == 'elastic' and step == 1:
                    # Simulare il default nullo che può essere inserito da
                    # accessi non const degli operatori bulk già esistenti.
                    part.Properties[2].SetValue(KM.YIELD_STRESS,0.)
                analysis.time = analysis._GetSolver().AdvanceInTime(analysis.time)
                set_affine(nodes, strain)
                # Inizializzare gli operatori senza l'assemblaggio del residuo
                # globale: questa è una verifica locale dei soli operatori Cut.
                analysis._GetSolver()._GetScheme().InitializeSolutionStep(
                    analysis._GetSolver().GetComputingModelPart(),KM.CompressedMatrix(),KM.Vector(),KM.Vector())
                row = dict(physics=physics, strain=strain.tolist(), entities={})
                expected = radial_stress(strain, physics == 'plastic')
                element = representatives['CutSbmSolidElement']
                element.SetValue(KM.BODY_FORCE, [0.,0.,0.])
                K, rhs = local_system(element, info)
                dofs = element.GetDofList(info)
                virtual = np.zeros((len(dofs),3))
                for i,dof in enumerate(dofs):
                    node = node_map[dof.Id()]
                    if dof.GetVariable() == KM.DISPLACEMENT_X:
                        virtual[i] = [node.X0,0.,node.Y0/2.]
                    else:
                        virtual[i] = [0.,node.Y0,node.X0/2.]
                weight = element.GetValue(KM.INTEGRATION_WEIGHT)
                measured = -virtual.T @ rhs / weight
                stress_error = np.linalg.norm(measured-expected) / max(1.,np.linalg.norm(expected))
                assert stress_error < 1.e-10, (physics,strain,measured,expected)
                delta = 1.e-8
                D = np.column_stack([(radial_stress(strain+delta*np.eye(3)[j],physics=='plastic')
                                      -radial_stress(strain-delta*np.eye(3)[j],physics=='plastic'))/(2*delta)
                                     for j in range(3)])
                tangent_error = np.linalg.norm(virtual.T @ K @ virtual/weight-D) / np.linalg.norm(D)
                assert tangent_error < 1.e-7, (physics,strain,tangent_error)
                row.update(stress_measured=measured.tolist(),stress_expected=expected.tolist(),
                           stress_relative_error=float(stress_error),tangent_relative_error=float(tangent_error))
                for name,entity in representatives.items():
                    K, rhs = local_system(entity, info)
                    local_dofs = entity.GetDofList(info)
                    unique = {(d.Id(),d.GetVariable().Name()):d for d in local_dofs}
                    rng = np.random.default_rng(31)
                    direction = {key:float(rng.normal()) for key in unique}
                    local_direction = np.array([direction[(d.Id(),d.GetVariable().Name())] for d in local_dofs])
                    original = {key:d.GetSolutionStepValue() for key,d in unique.items()}
                    finite_rhs = []
                    for sign in (1.,-1.):
                        for key,dof in unique.items():dof.SetSolutionStepValue(original[key]+sign*delta*direction[key])
                        finite_rhs.append(local_system(entity, info)[1])
                    for key,dof in unique.items():dof.SetSolutionStepValue(original[key])
                    derivative = -(finite_rhs[0]-finite_rhs[1])/(2*delta)
                    expected_derivative = K @ local_direction
                    error = np.linalg.norm(derivative-expected_derivative)/max(1.,np.linalg.norm(expected_derivative))
                    assert error < 1.e-6, (physics,name,strain,error)
                    row['entities'][name] = dict(id=entity.Id,jacobian_relative_error=float(error))
                if physics == 'plastic':
                    analysis._GetSolver()._GetScheme().FinalizeSolutionStep(
                        analysis._GetSolver().GetComputingModelPart(),KM.CompressedMatrix(),KM.Vector(),KM.Vector())
                    for name in ('CutSbmSolidElement','CutSbmSolidCondition'):
                        entity = representatives[name]
                        published = np.array([entity.GetValue(v) for v in (IGA.CAUCHY_STRESS_XX,IGA.CAUCHY_STRESS_YY,IGA.CAUCHY_STRESS_XY)])
                        assert np.allclose(published,expected,rtol=1.e-10,atol=1.e-10),(name,published,expected)
                        row['entities'][name]['stress_published'] = published.tolist()
                else:
                    if step == 0:
                        assert not part.Properties[2].Has(KM.YIELD_STRESS), 'Inserimento spurio di YIELD_STRESS'
                    else:
                        assert part.Properties[2].GetValue(KM.YIELD_STRESS) == 0.
                rows.append(row)
                print('CUT_ANALYTICAL_CHECK',json.dumps(row),flush=True)

(Path(__file__).resolve().parent / 'cut_analytical_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
