"""Frozen-history local branch audit, or native solve with the actual tangent utility."""
import argparse,csv,importlib.util,json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'sbm_cyclic_old_vs_new'))
import worker as native
sys.path.insert(0,json.loads((HERE/'bridge.json').read_text())['directory'])
import sbm_tangent_bridge as bridge
spec=importlib.util.spec_from_file_location('current_points',HERE.parent/'sbm_p3_taylor_iterations/worker.py')
current=importlib.util.module_from_spec(spec);spec.loader.exec_module(current)
from check_projection_history import MaterialPoint,strain_matrix
study=native.study;KM=study.KM

def write_csv(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def install_local(n,out):
    initialize=study.base.Analysis.Initialize
    def wrapped(self):
        initialize(self);solver=self._GetSolver();solve=solver.SolveSolutionStep;part=self.model['IgaModelPart'];info=part.ProcessInfo
        boundary=current.CurrentPoints(study,self.model);points=[]
        for c,Bs,Bt,s,t in boundary.points:
            points.extend([(c,'surrogate',Bs,s),(c,'true',Bt,t)])
        for e in self.model['IgaModelPart.StructuralAnalysisDomain'].Elements:
            if e.GetValue(KM.ACTIVATION_LEVEL)==2:continue
            g=e.GetGeometry();np.testing.assert_allclose(np.array(g.Jacobian(0))[:2,:2],np.eye(2),atol=1e-12)
            points.append((e,'bulk',strain_matrix(np.array(g.ShapeFunctionDerivatives(1,0))[:,:2]),MaterialPoint(e,info)))
        histories=[[np.zeros(6),0.] for _ in points]
        elastic=np.diag([67000.,67000.,33500.]);committed_validation=[]
        def coefficients(entity):return np.array([[v.GetSolutionStepValue(KM.DISPLACEMENT_X),v.GetSolutionStepValue(KM.DISPLACEMENT_Y)] for v in entity.GetGeometry()]).ravel()
        def calc(eps,hist):return current.phase.Audit.response(None,eps,hist)
        def checked():
            if float(self.time)<.00075-1e-12:
                result=solve();err=0.
                for j,(e,group,B,mp) in enumerate(points):
                    eps=B@coefficients(e);_,ref,new=calc(eps,histories[j]);stress,_=mp.response(eps)
                    err=max(err,float(np.linalg.norm(stress-ref)));mp.commit(eps);histories[j]=new
                assert err<1e-7,err
                committed_validation.append(dict(time=float(self.time),max_stress_error=err))
                return result
            target=HERE.parent/'sbm_p3_directional_warm_start'/f'audit_n{n}/native_iterates.npz'
            saved=np.load(target);values=saved['coefficients']
            builder=solver._GetBuilderAndSolver();dofs=sorted([d for d in builder.GetDofSet() if d.EquationId<values.shape[0]],key=lambda d:d.EquationId)
            assert [d.Id() for d in dofs]==list(saved['node_ids'])
            initial=[d.GetSolutionStepValue() for d in dofs];rows=[];validation=0.
            try:
                for index in (0,2):
                    for d,v in zip(dofs,values[:,index]):d.SetSolutionStepValue(float(v))
                    direction=values[:,index+1]-values[:,index]
                    for j,(e,group,B,mp) in enumerate(points):
                        eps=B@coefficients(e);local=np.array([direction[d.EquationId] for d in e.GetDofList(info)]);de=B@local
                        f0,ref,_=calc(eps,histories[j]);sigma,D=mp.response(eps)
                        validation=max(validation,float(np.linalg.norm(ref-sigma)))
                        tol=1e-8*max(201.,201.+10000.*histories[j][1]);near=abs(f0)<=tol
                        active=np.linalg.norm(D-elastic)/np.linalg.norm(elastic)>1e-8
                        predictions={};fds={};branches={}
                        scale=max(float(np.linalg.norm(D@de)),1e-12)
                        for eta in (1e-6,1e-8):
                            plus,Dp=mp.response(eps+eta*de);minus,Dm=mp.response(eps-eta*de)
                            fp=calc(eps+eta*de,histories[j])[0];fm=calc(eps-eta*de,histories[j])[0]
                            fds[eta]=((plus-sigma)/eta,(sigma-minus)/eta)
                            branches[eta]=(np.linalg.norm(Dp-elastic)/np.linalg.norm(elastic)>1e-8,np.linalg.norm(Dm-elastic)/np.linalg.norm(elastic)>1e-8,fp,fm)
                        mp.response(eps)
                        for order in (1,2):
                            numeric=np.array(bridge.numerical(mp.values,mp.law,order,-1.))
                            predictions[order]=numeric@de
                        small=None
                        if j%20==0:
                            mp.response(eps);small=np.array(bridge.numerical(mp.values,mp.law,2,1e-10))@de
                        repeat,_=mp.response(eps);repeat_error=float(np.linalg.norm(repeat-sigma));assert repeat_error<1e-10
                        for eta in (1e-6,1e-8):
                            forward,backward=fds[eta];ap,am,fp,fm=branches[eta]
                            rows.append(dict(iteration=index,entity=e.Id,group=group,eta=eta,
                                f_trial=float(f0),f_plus=float(fp),f_minus=float(fm),yield_tolerance=float(tol),near_yield=bool(near),
                                active_base=bool(active),active_plus=bool(ap),active_minus=bool(am),
                                switches_forward=bool(active!=ap),switches_backward=bool(active!=am),
                                robust_elastic_same_branch=bool(f0 < -tol and fp < -tol),robust_plastic_same_branch=bool(f0>tol and fp>tol),
                                direction_response_norm=scale,
                                analytic_forward_error=float(np.linalg.norm(forward-D@de)/scale),analytic_backward_error=float(np.linalg.norm(backward-D@de)/scale),
                                utility_forward_order1_error=float(np.linalg.norm(forward-predictions[1])/scale),utility_forward_order2_error=float(np.linalg.norm(forward-predictions[2])/scale),
                                utility_central_small_h_error=float(np.linalg.norm(forward-small)/scale) if small is not None else '',
                                repeat_stress_error=repeat_error))
            finally:
                for d,v in zip(dofs,initial):d.SetSolutionStepValue(float(v))
            assert validation<1e-7,validation
            write_csv(out/f'points_n{n}.csv',rows)
            (out/f'validation_n{n}.json').write_text(json.dumps(dict(points=len(points),committed_steps=committed_validation,max_reference_stress_error=validation,history_committed_at_failed_step=False),indent=2)+'\n')
            raise RuntimeError('Local material audit complete; no solve or commit at target time.')
        solver.SolveSolutionStep=checked
    study.base.Analysis.Initialize=wrapped

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--spans',type=int,required=True);parser.add_argument('--variant',required=True,choices=['local','utility_central','utility_forward']);parser.add_argument('--full',action='store_true');parser.add_argument('--h',type=float,default=-1.);parser.add_argument('--tag');args=parser.parse_args()
    out=HERE/(args.tag or (args.variant+'_full' if args.full else args.variant))
    if args.variant=='local':install_local(args.spans,out)
    else:
        modify=study.base.Analysis.ModifyInitialProperties
        def wrapped(self):
            modify(self);changed=[]
            for prop in self.model['IgaModelPart'].Properties:
                if prop.Has(KM.CONSTITUTIVE_LAW):
                    prop.SetValue(KM.CONSTITUTIVE_LAW,bridge.NumericalJ2(2 if args.variant=='utility_central' else 1,args.h));changed.append(prop.Id)
            assert changed,changed
            (out/f'law_replacement_n{args.spans}.json').write_text(json.dumps(dict(properties=changed,scope='volume and condition clones',order=2 if args.variant=='utility_central' else 1),indent=2))
        study.base.Analysis.ModifyInitialProperties=wrapped
    load=study.load_parameters
    def limited(*a,**kw):
        parameters,material=load(*a,**kw);
        if not args.full:parameters['problem_data']['end_time'].SetDouble(.00075)
        return parameters,material
    study.load_parameters=limited
    sys.argv=[__file__,'--mode','new','--degree','3','--spans',str(args.spans),'--output',str(out)];native.main()
if __name__=='__main__':main()
