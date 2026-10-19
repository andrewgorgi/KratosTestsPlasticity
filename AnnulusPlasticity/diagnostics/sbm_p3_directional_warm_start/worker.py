"""Native failed iterates, grouped one-sided derivatives, and a true L2 projection."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from KratosMultiphysics.scipy_conversion_tools import to_csr
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'sbm_cyclic_old_vs_new'))
import worker as native
study=native.study
KM=study.KM

def write_csv(path,rows):
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)

def install(n,fresh=False):
    initialize=study.base.Analysis.Initialize
    output=(HERE/'fresh' if fresh else HERE)/f'audit_n{n}';output.mkdir(exist_ok=True)
    def wrapped(self):
        initialize(self);solver=self._GetSolver();solve=solver.SolveSolutionStep
        strategy=solver._GetSolutionStrategy();strategy.SetUpNonconvergedSolutionsFlag(True)
        builder=solver._GetBuilderAndSolver();scheme=solver._GetScheme();part=solver.GetComputingModelPart()
        def context():return strategy.GetSystemMatrix(),strategy.GetSystemVector(),strategy.GetSolutionVector()
        def dofs():return sorted([d for d in builder.GetDofSet() if d.EquationId<len(strategy.GetSystemVector())],key=lambda d:d.EquationId)
        def uvec():return np.array([d.GetSolutionStepValue() for d in dofs()])
        def set_u(u):
            for d,v in zip(dofs(),u):d.SetSolutionStepValue(float(v))
        def build():
            A,R,_=context();KM.UblasSparseSpace().SetToZeroMatrix(A);KM.UblasSparseSpace().SetToZeroVector(R)
            builder.Build(scheme,part,A,R);return to_csr(A).copy(),np.array(R).copy()
        def groups(tangent=False):
            size=len(strategy.GetSystemVector());out={key:np.zeros(size) for key in ('dirichlet','neumann')}
            entries={key:([],[],[]) for key in out}
            for c in part.Conditions:
                if c.IsDefined(KM.ACTIVE) and not c.Is(KM.ACTIVE):continue
                name=c.Info();key='neumann' if 'SbmLoadSolidCondition' in name else 'dirichlet' if 'SbmSolidCondition' in name else None
                if key is None:continue
                ids=np.array([d.EquationId for d in c.GetDofList(part.ProcessInfo)],dtype=int);mask=ids<size;valid=ids[mask]
                localR=KM.Vector()
                if tangent:
                    localK=KM.Matrix();c.CalculateLocalSystem(localK,localR,part.ProcessInfo)
                    rr,cc=np.meshgrid(valid,valid,indexing='ij');entries[key][0].extend(rr.ravel());entries[key][1].extend(cc.ravel());entries[key][2].extend(np.array(localK)[np.ix_(mask,mask)].ravel())
                else:c.CalculateRightHandSide(localR,part.ProcessInfo)
                np.add.at(out[key],valid,np.array(localR)[mask])
            if tangent:return out,{key:coo_matrix((v[2],(v[0],v[1])),shape=(size,size)).tocsr() for key,v in entries.items()}
            return out
        def residual():
            _,R,_=context();KM.UblasSparseSpace().SetToZeroVector(R);builder.BuildRHS(scheme,part,R)
            total=np.array(R).copy();out=groups();out['bulk']=total-out['dirichlet']-out['neumann'];out['total']=total;return out
        def field_error():
            weights=[];errors=[];norms=[]
            for e in self.model['IgaModelPart.StructuralAnalysisDomain'].Elements:
                if e.IsDefined(KM.ACTIVE) and not e.Is(KM.ACTIVE):continue
                if e.GetValue(KM.ACTIVATION_LEVEL)==2:continue
                g=e.GetGeometry();N=np.array(g.ShapeFunctionsValues())[0];xy=N@np.array([[node.X0,node.Y0] for node in g]);v=N@np.array([[node.GetSolutionStepValue(KM.DISPLACEMENT_X),node.GetSolutionStepValue(KM.DISPLACEMENT_Y)] for node in g])
                exact=study.exact.evaluate(xy.reshape(1,2),float(self.time))[0][0];w=float(e.GetValue(KM.INTEGRATION_WEIGHT));assert w>0
                weights.append(w);errors.append((v-exact)@(v-exact));norms.append(exact@exact)
            return float(np.sqrt(np.dot(weights,errors)/np.dot(weights,norms)))
        def projection():
            size=len(strategy.GetSystemVector());ri=[];ci=[];data=[];target=np.zeros(size);cache=[]
            for e in self.model['IgaModelPart.StructuralAnalysisDomain'].Elements:
                if e.IsDefined(KM.ACTIVE) and not e.Is(KM.ACTIVE):continue
                if e.GetValue(KM.ACTIVATION_LEVEL)==2:continue
                g=e.GetGeometry();N=np.array(g.ShapeFunctionsValues())[0];xy=N@np.array([[node.X0,node.Y0] for node in g]);w=float(e.GetValue(KM.INTEGRATION_WEIGHT));assert w>0
                ids=np.array([d.EquationId for d in e.GetDofList(part.ProcessInfo)],dtype=int)
                assert len(ids)==2*len(N) and np.all(ids<size)
                uv=study.exact.evaluate(xy.reshape(1,2),float(self.time))[0][0]
                mass=np.kron(np.outer(N,N),np.eye(2))*w;rr,cc=np.meshgrid(ids,ids,indexing='ij');ri.extend(rr.ravel());ci.extend(cc.ravel());data.extend(mass.ravel())
                np.add.at(target,ids,(N[:,None]*uv[None,:]*w).ravel())
            M=coo_matrix((data,(ri,ci)),shape=(size,size)).tocsr();projected=spsolve(M,target)
            assert np.isfinite(projected).all()
            error=float(np.linalg.norm(M@projected-target)/np.linalg.norm(target));assert error<1e-10,error
            np.savez_compressed(output/'projection.npz',coefficients=projected,target=target,equation_node_ids=[d.Id() for d in dofs()])
            return projected,error
        def directional(values):
            saved=uvec();saved_iteration=int(part.ProcessInfo[KM.NL_ITERATION_NUMBER]);rows=[];checks=[]
            try:
                for index in (0,2,values.shape[1]-2):
                    if index<0 or index+1>=values.shape[1]:continue
                    base=values[:,index];direction=values[:,index+1]-base;set_u(base)
                    part.ProcessInfo[KM.NL_ITERATION_NUMBER]=index+1
                    A,R=build();baseR,bcK=groups(True);baseR['total']=R;baseR['bulk']=R-baseR['dirichlet']-baseR['neumann']
                    K=dict(bcK,total=A,bulk=A-bcK['dirichlet']-bcK['neumann'])
                    # Native sign convention: K = -dR/du.
                    kd={key:mat@direction for key,mat in K.items()}
                    alpha=float((kd['total']@R)/(R@R))
                    checks.append(dict(iteration=index,residual_norm=float(np.linalg.norm(R)),actual_update_norm=float(np.linalg.norm(direction)),effective_step_scale=alpha,
                        relative_K_update_minus_scaled_R=float(np.linalg.norm(kd['total']-alpha*R)/np.linalg.norm(R)),L2_relative_current=field_error()))
                    for eta in (1e-2,1e-4,1e-6,1e-8):
                        set_u(base+eta*direction);plus=residual()
                        set_u(base-eta*direction);minus=residual()
                        for key in ('bulk','dirichlet','neumann','total'):
                            forward=(plus[key]-baseR[key])/eta;backward=(baseR[key]-minus[key])/eta
                            scale=max(np.linalg.norm(kd[key]),1e-14)
                            rows.append(dict(iteration=index,eta=eta,group=key,base_residual_norm=float(np.linalg.norm(baseR[key])),predicted_derivative_norm=float(np.linalg.norm(kd[key])),
                                forward_relative_mismatch=float(np.linalg.norm(forward+kd[key])/scale),backward_relative_mismatch=float(np.linalg.norm(backward+kd[key])/scale),
                                central_relative_mismatch=float(np.linalg.norm((forward+backward)/2+kd[key])/scale),
                                forward_phi_slope=float((plus[key]@plus[key]-baseR[key]@baseR[key])/(2*eta)),predicted_phi_slope=float(-baseR[key]@kd[key])))
                    set_u(base);again=residual()
                    assert np.linalg.norm(again['total']-R)/max(np.linalg.norm(R),1)<1e-10
            finally:
                set_u(saved);part.ProcessInfo[KM.NL_ITERATION_NUMBER]=saved_iteration;build()
            write_csv(output/'directional.csv',rows);(output/'native_directions.json').write_text(json.dumps(checks,indent=2)+'\n')
        def checked():
            if abs(float(self.time)-.00075)>1e-12:return solve()
            start=uvec();initialR=residual()['total'];baseline={}
            if not fresh:
                try:
                    flag=solve();baseline=dict(status='converged',returned=bool(flag))
                except RuntimeError as exc:baseline=dict(status='failed',error=str(exc))
                matrix,stored_dofs=strategy.GetNonconvergedSolutions();values=np.array(matrix).copy()
                np.savez_compressed(output/'native_iterates.npz',coefficients=values,node_ids=[d.Id() for d in dofs()])
                baseline['iterations']=int(part.ProcessInfo[KM.NL_ITERATION_NUMBER]);baseline['time']=float(self.time)
                (output/'baseline.json').write_text(json.dumps(baseline,indent=2)+'\n')
                directional(values)
            set_u(start);repeatR=residual()['total'];history_error=float(np.linalg.norm(repeatR-initialR)/max(1.,np.linalg.norm(initialR)))
            assert history_error<1e-10,history_error
            projected,normal_error=projection();scans=[]
            for a in (0.,.25,.5,.75,1.):
                set_u((1-a)*start+a*projected);r=residual()['total'];scans.append(dict(fraction=a,residual_norm=float(np.linalg.norm(r)),L2_relative_current=field_error()))
            write_csv(output/'projection_scan.csv',scans)
            set_u(projected);warm_initial_R=residual()['total']
            # Re-evaluate the tangent at the projected state before invoking the unchanged native solver.
            build();part.ProcessInfo[KM.NL_ITERATION_NUMBER]=0
            control=dict(history_residual_repeat_relative_error=history_error,mass_normal_equation_relative_error=normal_error,
                original_initial_residual_norm=float(np.linalg.norm(initialR)),projected_initial_residual_norm=float(np.linalg.norm(warm_initial_R)),
                projected_L2_relative=scans[-1]['L2_relative_current'],time=float(self.time))
            (output/'warm_start.json').write_text(json.dumps(control,indent=2)+'\n')
            try:
                result=solve();control['status']='converged';return result
            except RuntimeError as exc:
                control['status']='failed';control['error']=str(exc);raise
            finally:
                control['iterations']=int(part.ProcessInfo[KM.NL_ITERATION_NUMBER]);control['final_residual_norm']=float(np.linalg.norm(residual()['total']));control['final_L2_relative']=field_error()
                (output/'warm_start.json').write_text(json.dumps(control,indent=2)+'\n')
        solver.SolveSolutionStep=checked
    study.base.Analysis.Initialize=wrapped
    load=study.load_parameters
    def limited(*args,**kwargs):
        params,material=load(*args,**kwargs);params['problem_data']['end_time'].SetDouble(.00075);return params,material
    study.load_parameters=limited

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--spans',type=int,required=True);parser.add_argument('--fresh',action='store_true');args=parser.parse_args();install(args.spans,args.fresh)
    sys.argv=[__file__,'--mode','new','--degree','3','--spans',str(args.spans),'--output',str(HERE/'fresh' if args.fresh else HERE)];native.main()
if __name__=='__main__':main()
