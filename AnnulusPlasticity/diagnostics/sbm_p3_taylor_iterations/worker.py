"""Replay native iterations; separate Taylor strain terms and error sources."""
import argparse,importlib.util,json,sys
from math import factorial
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline,PPoly
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(HERE.parent/'sbm_p3_global_causes'))
sys.path.insert(0,str(HERE.parent/'sbm_cyclic_old_vs_new'))
import worker as native
from check_projection_history import strain_matrix
spec=importlib.util.spec_from_file_location('phase_audit',HERE.parent/'sbm_p3_pair_phases/worker.py')
phase=importlib.util.module_from_spec(spec);spec.loader.exec_module(phase)
OriginalPoints=phase.BoundaryPoints
class CurrentPoints(OriginalPoints):
    def __init__(self,study,model):
        super().__init__(study,model)
        nodes=[]
        def collect(p):
            if p.Name=='projection_data':nodes.extend(p.Nodes)
            for child in p.SubModelParts:collect(child)
        collect(model['skin_model_part'].GetSubModelPart('inner'))
        fixed=[]
        for (c,Bs,Bt,s,t),node in zip(self.points,sorted(nodes,key=lambda n:n.Id)):
            g=c.GetGeometry();center=g.Center();d=np.array([node.X0-center.X,node.Y0-center.Y])
            grad=np.array(g.ShapeFunctionDerivatives(1,0))[:,:2]
            degree=round(np.sqrt(len(g)))-1
            for order in range(2,degree+1):
                deriv=np.array(g.ShapeFunctionDerivatives(order,0))
                for k in range(order):
                    a=d[0]**(order-1-k)*d[1]**k/(factorial(order-1-k)*factorial(k))
                    grad[:,0]+=a*deriv[:,k];grad[:,1]+=a*deriv[:,k+1]
            fixed.append((c,Bs,strain_matrix(grad),s,t))
        self.points=fixed
phase.BoundaryPoints=CurrentPoints

class Audit(phase.Audit):
    def __init__(self,study,model,out,n,solver):
        super().__init__(study,model,out)
        self.solver=solver;self.cache=[]
        p=3;breaks=np.linspace(-2.2,2.2,n+1)
        knots=np.r_[[breaks[0]]*p,breaks,[breaks[-1]]*p];count=n+p
        greville=np.array([knots[i+1:i+p+1].mean() for i in range(count)])
        self.global_nodes=list(self.part.Nodes)
        allcoords=np.array([[node.X0,node.Y0] for node in self.global_nodes])
        global_ij=np.argmin(abs(allcoords[:,:,None]-greville),axis=2)
        assert np.max(abs(greville[global_ij]-allcoords))<1e-11
        spl=BSpline(knots,np.eye(count),p)
        polys=[PPoly.from_spline(BSpline(knots,np.eye(count)[i],p)) for i in range(count)]
        for (c,Bs,Bt,s,t),node in zip(self.points.points,self.nodes):
            g=c.GetGeometry();center=np.array([g.Center().X,g.Center().Y]);true=np.array([node.X0,node.Y0]);d=true-center
            mats=[Bs]
            for order in range(2,2*p+1):
                deriv=np.array(g.ShapeFunctionDerivatives(order,0));grad=np.zeros((len(g),2))
                for k in range(order):
                    a=d[0]**(order-1-k)*d[1]**k/(factorial(order-1-k)*factorial(k))
                    grad[:,0]+=a*deriv[:,k];grad[:,1]+=a*deriv[:,k+1]
                mats.append(strain_matrix(grad))
            np.testing.assert_allclose(sum(mats[:3]),Bt,atol=1e-12)
            coords=np.array([[node.X0,node.Y0] for node in g]);ij=np.argmin(abs(coords[:,:,None]-greville),axis=2)
            spans=ij.min(axis=0)+p
            direct=[]
            for i,j in ij:
                ax,ay=polys[i].c[:,spans[0]],polys[j].c[:,spans[1]]
                xx,yy=true-knots[spans]
                direct.append([np.polyval(np.polyder(ax),xx)*np.polyval(ay,yy),np.polyval(ax,xx)*np.polyval(np.polyder(ay),yy)])
            directB=strain_matrix(np.array(direct))
            np.testing.assert_allclose(sum(mats),directB,atol=1e-10)
            vx,vy=spl(true[0]),spl(true[1]);gx,gy=spl.derivative()(true[0]),spl.derivative()(true[1])
            globalgrad=np.column_stack((gx[global_ij[:,0]]*vy[global_ij[:,1]],vx[global_ij[:,0]]*gy[global_ij[:,1]]))
            expected={(i,j) for i in range(count) for j in range(count) if abs(gx[i]*vy[j])+abs(vx[i]*gy[j])>1e-12}
            missing=expected-set(map(tuple,global_ij.tolist()))
            cross=any(np.any((breaks[1:-1]-min(center[i],true[i])>1e-10)&(max(center[i],true[i])-breaks[1:-1]>1e-10)) for i in (0,1))
            self.cache.append((c,center,true,mats,directB,strain_matrix(globalgrad),bool(cross),len(missing),s,t))

    def replay(self,strategy,state):
        matrix,dofs=strategy.GetNonconvergedSolutions();values=np.array(matrix)
        saved=[d.GetSolutionStepValue() for d in dofs]
        KM=self.KM;time=float(self.part.ProcessInfo[KM.TIME]);rows=[];summaries=[]
        truepoints=np.array([v[2] for v in self.cache]);grad=self.study.exact.evaluate(truepoints,time)[1]
        exact=np.column_stack((grad[:,0],grad[:,3],grad[:,1]+grad[:,2]))
        scale=np.maximum(np.linalg.norm(exact,axis=1),1e-14)
        load=abs(float(np.interp(time,self.study.exact.TIMES,self.study.exact.LEVELS)))
        # Only first monotonic loading is used for front indicators; later frames marked unavailable.
        rp=np.sqrt(np.sqrt(3)*self.study.exact.C_FINAL*load/self.study.exact.YIELD) if time<=.001+1e-12 else None
        builder=self.solver._GetBuilderAndSolver();scheme=self.solver._GetScheme();part=self.solver.GetComputingModelPart()
        rhs=strategy.GetSystemVector();rhs_saved=np.array(rhs).copy()
        try:
            for iteration in range(values.shape[1]):
                for dof in dofs:dof.SetSolutionStepValue(float(values[dof.EquationId,iteration]))
                KM.UblasSparseSpace().SetToZeroVector(rhs);builder.BuildRHS(scheme,part,rhs)
                global_residual=float(np.linalg.norm(np.array(rhs)))
                gu=np.array([[n.GetSolutionStepValue(KM.DISPLACEMENT_X),n.GetSolutionStepValue(KM.DISPLACEMENT_Y)] for n in self.global_nodes]).ravel()
                local=[]
                for index,(c,scoord,tcoord,mats,directB,globalB,knot,missing,s,t) in enumerate(self.cache):
                    u=self.points.displacement(c);terms=np.array([m@u for m in mats]);cumulative=terms.cumsum(axis=0)
                    full=directB@u;glob=globalB@gu
                    fs,_,_=self.response(cumulative[0],self.history[c.Id][0]);ft,_,_=self.response(cumulative[2],self.history[c.Id][1])
                    tol=1e-8*max(201.,201.+10000.*self.history[c.Id][0][1],201.+10000.*self.history[c.Id][1][1])
                    robust=(fs < -tol and ft > tol) or (fs > tol and ft < -tol)
                    d=tcoord-scoord
                    roots=np.roots([d@d,2*scoord@d,scoord@scoord-rp*rp]) if rp is not None and d@d else []
                    front=any(abs(z.imag)<1e-12 and 1e-10<float(z.real)<1-1e-10 for z in roots)
                    item=dict(time=time,state=state,iteration=iteration,condition=c.Id,crosses_knot=knot,missing_global_coefficients=missing,
                        exact_front_available=rp is not None,crosses_exact_front=bool(front),robust_phase_mismatch=bool(robust),
                        yield_s=float(fs),yield_t=float(ft),
                        second_term_over_first_order=float(np.linalg.norm(terms[2])/max(np.linalg.norm(cumulative[1]),1e-14)),
                        second_term_over_exact=float(np.linalg.norm(terms[2])/scale[index]),
                        q2_vs_local_continuation=float(np.linalg.norm(cumulative[2]-full)/scale[index]),
                        full_taylor_validation=float(np.linalg.norm(cumulative[-1]-full)/max(np.linalg.norm(full),1e-14)),
                        local_vs_global_spline=float(np.linalg.norm(full-glob)/scale[index]) if not missing else '',
                        q2_vs_global_spline=float(np.linalg.norm(cumulative[2]-glob)/scale[index]) if not missing else '')
                    for k in range(6):
                        item[f'term{k}_over_exact']=float(np.linalg.norm(terms[k])/scale[index])
                        item[f'error_q{k}_over_exact']=float(np.linalg.norm(cumulative[k]-exact[index])/scale[index])
                    rows.append(item);local.append(item)
                sums=dict(time=time,state=state,iteration=iteration,residual_norm=global_residual,pairs=len(local),
                    robust_phase_mismatch=sum(r['robust_phase_mismatch'] for r in local),crosses_exact_front=sum(r['crosses_exact_front'] for r in local),
                    crosses_knot=sum(r['crosses_knot'] for r in local),max_full_taylor_validation=max(r['full_taylor_validation'] for r in local))
                for key in ['second_term_over_first_order','second_term_over_exact','q2_vs_local_continuation','error_q1_over_exact','error_q2_over_exact','error_q5_over_exact']:
                    vals=[r[key] for r in local];sums['median_'+key]=float(np.median(vals));sums['max_'+key]=float(max(vals))
                for cross in (False,True):
                    group=[r for r in local if r['crosses_knot']==cross and not r['missing_global_coefficients']]
                    label='crossing' if cross else 'same_span'
                    sums[label+'_count']=len(group)
                    sums[label+'_median_local_vs_global']=float(np.median([r['local_vs_global_spline'] for r in group])) if group else ''
                    sums[label+'_median_q2_vs_global']=float(np.median([r['q2_vs_global_spline'] for r in group])) if group else ''
                summaries.append(sums)
        finally:
            for dof,value in zip(dofs,saved):dof.SetSolutionStepValue(float(value))
            KM.UblasSparseSpace().SetToZeroVector(rhs);builder.BuildRHS(scheme,part,rhs)
            for i,value in enumerate(rhs_saved):rhs[i]=float(value)
        phase.write_csv(self.out/f'iterations_t{time:.8f}_{state}.csv',summaries)
        phase.write_csv(self.out/f'pairs_t{time:.8f}_{state}_iterations.csv',rows)
        # Validate only absolute/normalized polynomial reconstruction; ratios at zero states are harmless.
        (self.out/f'validation_t{time:.8f}_{state}.json').write_text(json.dumps(dict(max_full_taylor_relative=max(r['full_taylor_validation'] for r in rows)),indent=2))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--spans',type=int,required=True);parser.add_argument('--mode',choices=['baseline','audit'],required=True);args=parser.parse_args()
    if args.mode=='audit':
        study=native.study;initialize=study.base.Analysis.Initialize
        def wrapped(self):
            initialize(self);solver=self._GetSolver();audit=Audit(study,self.model,HERE/'audit'/f'audit_n{args.spans}',args.spans,solver)
            solver._GetSolutionStrategy().SetUpNonconvergedSolutionsFlag(True)
            solve=solver.SolveSolutionStep
            def checked():
                try:converged=solve()
                except RuntimeError:
                    audit.snapshot('failed');audit.replay(solver._GetSolutionStrategy(),'failed');raise
                audit.snapshot('converged');audit.replay(solver._GetSolutionStrategy(),'converged');audit.commit();return converged
            solver.SolveSolutionStep=checked
        study.base.Analysis.Initialize=wrapped
    sys.argv=[__file__,'--mode','new','--degree','3','--spans',str(args.spans),'--output',str(HERE/args.mode)]
    native.main()
if __name__=='__main__':main()
