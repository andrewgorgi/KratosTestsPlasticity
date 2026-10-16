"""Check actual SBM mixed-boundary nodal data without solving the load steps."""
import argparse
from collections import Counter
import json
import numpy as np
import annulus_study as study
from analytic import C_FINAL, END_TIME


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree',type=int,default=1)
    parser.add_argument('--spans',type=int,nargs=2,default=[20,20])
    parser.add_argument('--time',type=float,default=END_TIME)
    args = parser.parse_args()
    if min(args.degree,*args.spans) < 1 or not 0 <= args.time <= END_TIME:
        parser.error('Use positive degree/spans and a time in [0, END_TIME].')
    args.parameters = None
    args.boundary = 'mixed'
    args.dt = None
    study.KM.ParallelUtilities.SetNumThreads(1)
    folder = study.ROOT/'sbm'
    parameters = study.load_parameters(folder,args,args.spans)
    with study.shared.working_directory(folder):
        model = study.KM.Model()
        analysis = study.Analysis(model,parameters)
        analysis.Initialize()
        analysis.time = args.time
        model['IgaModelPart'].CloneTimeStep(args.time)
        analysis.InitializeSolutionStep()
        counts = Counter(study.condition_name(c).strip('"') for c in model['IgaModelPart'].Conditions)
        report = dict(time=args.time,degree=args.degree,spans=args.spans,conditions=dict(counts),boundaries={})
        c = C_FINAL*args.time/END_TIME
        # All displacement DOFs are still zero: the Neumann RHS contains only
        # the prescribed external load. Check its work on a rigid rotation.
        moment = 0.
        resultant = np.zeros(2)
        part = model['IgaModelPart']
        for condition in part.Conditions:
            if study.condition_name(condition).strip('"') != 'SbmLoadSolidCondition':
                continue
            rhs = study.KM.Vector()
            condition.CalculateRightHandSide(rhs,part.ProcessInfo)
            for dof,value in zip(condition.GetDofList(part.ProcessInfo),rhs):
                node = part.GetNode(dof.Id())
                if dof.GetVariable().Name() == 'DISPLACEMENT_X':
                    resultant[0] += value
                    moment -= node.Y0*value
                elif dof.GetVariable().Name() == 'DISPLACEMENT_Y':
                    resultant[1] += value
                    moment += node.X0*value
        report['assembled_neumann_at_zero_displacement'] = dict(
            resultant=resultant.tolist(),moment=float(moment),
            physical_circle_moment_for_context=float(-2*np.pi*c),
            note='Raw RHS uses surrogate test functions and weights; compare the full shifted formulation before asserting physical torque balance.')
        for side in ('inner','outer'):
            rows = []
            for node in model['skin_model_part'].GetSubModelPart(side).Nodes:
                x = np.array([node.X,node.Y])
                r = np.linalg.norm(x)
                er = x/r
                et = np.array([-er[1],er[0]])
                force = np.array(node.GetValue(study.KM.FORCE))[:2]
                normal = np.array(node.GetValue(study.KM.NORMAL))[:2]
                expected = (-1 if side == 'inner' else 1)*c/r**2*et
                rows.append((r,np.linalg.norm(force-expected),force@et,force@er,normal@er))
            rows = np.array(rows)
            if not len(rows):
                raise RuntimeError(f'No nodes found on {side} boundary.')
            report['boundaries'][side] = dict(
                nodes=len(rows),radius_range=[float(rows[:,0].min()),float(rows[:,0].max())],
                max_force_error=float(rows[:,1].max()),
                tangential_force_range=[float(rows[:,2].min()),float(rows[:,2].max())],
                max_radial_force=float(np.abs(rows[:,3]).max()),
                stored_normal_dot_er_range=[float(rows[:,4].min()),float(rows[:,4].max())])
    output = folder/'results/mixed'/f'p{args.degree}'
    output.mkdir(parents=True,exist_ok=True)
    (output/'boundary_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    if any(data['max_force_error'] > 1e-9 for data in report['boundaries'].values()):
        raise RuntimeError('Prescribed nodal traction disagrees with the independent annular formula.')


if __name__ == '__main__':
    main()
