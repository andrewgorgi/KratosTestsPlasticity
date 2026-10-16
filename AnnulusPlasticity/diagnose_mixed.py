"""Isolated SBM controls; preserve the annular production parameters/material.

Perfect plasticity uses an affine patch or a nonaffine cosine manufactured
field with its own body force and boundary traction, not H=0 annular torsion.
"""
import argparse
import json
from pathlib import Path
import traceback
import numpy as np
import sympy as sp
import annulus_study as study
from analytic import evaluate as annulus_evaluate, END_TIME, E, NU, H, YIELD
from analytic_fields import evaluator, expressions, x, y, t


def reference(case):
    if case in ('annulus_elastic','annulus_hardening'):
        return lambda points,time: annulus_evaluate(points,time)[:3], lambda points,time: np.zeros((len(points),3))
    field = 'linear' if case.startswith('patch_') else 'cosine_middle'
    scale = 10. if case in ('patch_perfect','patch_hardening') else 1.
    physics = ('elastic' if case == 'patch_elastic' else
               'plastic_hardening' if case == 'patch_hardening' else 'plastic')
    evaluate = evaluator(field,physics)
    force = sp.lambdify((x,y,t),list(expressions(field,physics)[3]),'numpy')
    def evaluate_scaled(points,time):
        return evaluate(points,scale*time)
    def body(points,time):
        points = np.asarray(points)
        with np.errstate(divide='ignore',invalid='ignore'):
            values = force(points[:,0],points[:,1],scale*time)
        return np.array([np.broadcast_to(v,(len(points),)) for v in values]).T
    return evaluate_scaled,body


def data(points,time,evaluate):
    points = np.asarray(points)
    u,_,stress = evaluate(points,time)
    radii = np.linalg.norm(points,axis=1)
    normals = np.where(radii < 1.5,-1.,1.)[:,None]*points/radii[:,None]
    force = np.column_stack((stress[:,0]*normals[:,0]+stress[:,2]*normals[:,1],
                             stress[:,2]*normals[:,0]+stress[:,1]*normals[:,1]))
    return u,force


def metrics(samples,evaluate):
    exact,gradient,_ = evaluate(samples['points'],samples['time'])
    u = samples['displacement']; w = samples['weights']
    bulk = np.isfinite(samples['gradient']).all(axis=1)
    l2 = np.sum(w[:,None]*(u-exact)**2)/np.sum(w[:,None]*exact**2)
    h1 = (np.sum(w[bulk,None]*(u[bulk]-exact[bulk])**2)
          +np.sum(w[bulk,None]*(samples['gradient'][bulk]-gradient[bulk])**2))
    h1 /= np.sum(w[bulk,None]*exact[bulk]**2)+np.sum(w[bulk,None]*gradient[bulk]**2)
    if l2 < 0 or h1 < 0:
        raise RuntimeError('Negative squared norm.')
    return dict(time=samples['time'],L2_relative=float(np.sqrt(l2)),H1_bulk_relative=float(np.sqrt(h1)),
                max_equivalent_plastic_strain=float(samples['equivalent_plastic_strain'].max()))


def analysis_type(evaluate,body):
    class ControlAnalysis(study.Analysis):
        def InitializeSolutionStep(self):
            time = float(self._GetSolver().GetComputingModelPart().ProcessInfo[study.KM.TIME])
            skin = self.model['skin_model_part']
            skin.ProcessInfo[study.KM.TIME] = time
            nodes = list(skin.Nodes)
            points = [[n.X0,n.Y0] for n in nodes]
            u,forces = data(points,time,evaluate)
            for n,value,force in zip(nodes,u,forces):
                n.SetValue(study.KM.DISPLACEMENT,[float(value[0]),float(value[1]),0.])
                n.SetValue(study.KM.FORCE,[float(force[0]),float(force[1]),0.])
            part = self.model['IgaModelPart']
            elements = list(part.Elements)
            points = [[e.GetGeometry().Center().X,e.GetGeometry().Center().Y] for e in elements]
            for element,force in zip(elements,body(points,time)):
                element.SetValue(study.KM.BODY_FORCE,[float(v) for v in force])
            study.shared.StructuralMechanicsAnalysis.InitializeSolutionStep(self)
    return ControlAnalysis


def run(case,degree,spans,output,args):
    evaluate,body = reference(case)
    ControlAnalysis = analysis_type(evaluate,body)
    folder = study.ROOT/'sbm'
    parameters = study.load_parameters(folder,args,[spans,spans])
    if args.projection_type is not None:
        for modeler in parameters['modelers'].values():
            settings = modeler['Parameters']
            if settings.Has('element_condition_list'):
                for entry in settings['element_condition_list'].values():
                    if entry.Has('sbm_parameters'):
                        entry['sbm_parameters']['projection_type'].SetString(args.projection_type)
    for i in range(2):
        study.shared.geometry_settings(parameters)['polynomial_order'][i].SetInt(degree)
    material = json.loads(study.KM.Parameters((study.ROOT/'materials.json').read_text()).PrettyPrintJsonString())
    if args.penalty_factor is not None:
        material['properties'][0]['Material']['Variables']['PENALTY_FACTOR'] = args.penalty_factor
    if case in ('patch_perfect','cosine_perfect'):
        material['properties'][0]['Material']['Variables']['ISOTROPIC_HARDENING_MODULUS'] = 0.
    expected = dict(YOUNG_MODULUS=E,POISSON_RATIO=NU,YIELD_STRESS=YIELD,
                    ISOTROPIC_HARDENING_MODULUS=0. if case in ('patch_perfect','cosine_perfect') else H)
    for key,value in expected.items():
        actual = material['properties'][0]['Material']['Variables'][key]
        if not np.isclose(actual,value,rtol=1e-12,atol=1e-12):
            raise ValueError(f'{key}={actual} disagrees with diagnostic reference {value}.')
    end = .0002 if case == 'annulus_elastic' else END_TIME
    dt = min(args.dt or .001,end/2)
    parameters['problem_data']['end_time'].SetDouble(end)
    parameters['solver_settings']['time_stepping']['time_step'].SetDouble(dt)
    output.mkdir(parents=True,exist_ok=True)
    (output/'materials.json').write_text(json.dumps(material,indent=2)+'\n')
    parameters['solver_settings']['material_import_settings']['materials_filename'].SetString(str(output/'materials.json'))
    (output/'parameters.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    result = dict(case=case,degree=degree,spans=spans,dt=dt,
                  material=material['properties'][0]['Material']['Variables'])
    with study.shared.working_directory(folder):
        model = study.KM.Model()
        analysis = ControlAnalysis(model,parameters)
        try:
            analysis.Run()
            result.update(status='converged',**metrics(analysis.samples,evaluate))
            np.savez(output/'solution.npz',**analysis.samples)
        except Exception as error:
            result.update(status='failed',error=str(error),time=float(getattr(analysis,'time',0.)))
            (output/'traceback.txt').write_text(traceback.format_exc())
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('CONTROL_RESULT '+json.dumps(result),flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',nargs='+',choices=('annulus_elastic','annulus_hardening','patch_elastic','patch_hardening','patch_perfect','cosine_perfect'),
                        default=['annulus_elastic','patch_elastic','patch_hardening','patch_perfect','cosine_perfect','annulus_hardening'])
    parser.add_argument('--degrees',type=int,nargs='+',default=[1,2])
    parser.add_argument('--meshes',type=int,nargs='+',default=[10,20])
    parser.add_argument('--output',type=Path,default=study.ROOT/'diagnostics/mixed_controls')
    parser.add_argument('--dt',type=float)
    parser.add_argument('--echo-level',type=int,default=0)
    parser.add_argument('--strategy',choices=('line_search','newton_raphson'))
    parser.add_argument('--penalty-factor',type=float,help='Override only the copied diagnostic material')
    parser.add_argument('--projection-type',choices=('linealized','analytical'),help='Override only the diagnostic projection')
    args = parser.parse_args()
    if min(*args.degrees,*args.meshes) < 1 or (args.dt is not None and args.dt <= 0):
        parser.error('Degrees, meshes and dt must be positive.')
    args.output = args.output.resolve()
    args.parameters = None; args.boundary = 'mixed'; args.degree = 1
    study.KM.ParallelUtilities.SetNumThreads(1)
    results = []
    for case in args.cases:
        for degree in args.degrees:
            for spans in args.meshes:
                results.append(run(case,degree,spans,args.output/case/f'p{degree}'/f'n{spans}',args))
                (args.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__ == '__main__':
    main()
