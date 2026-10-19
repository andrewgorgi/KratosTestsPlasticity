"""Classify native failed trial states using validated independent histories."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'sbm_cyclic_old_vs_new'))
import worker as native_worker
sys.path.insert(0, str(HERE.parent/'sbm_p3_global_causes'))
from branch_points import BoundaryPoints


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def traction(n):
    return np.array([[n[0], 0., n[1]], [0., n[1], n[0]]])


class Audit:
    def __init__(self, study, model, out):
        self.study, self.KM, self.out = study, study.KM, out
        out.mkdir(parents=True, exist_ok=True)
        self.part = model['IgaModelPart']
        self.points = BoundaryPoints(study, model)
        self.history = {c.Id:[[np.zeros(6), 0.], [np.zeros(6), 0.]] for c, *_ in self.points.points}
        self.rows = []
        def collect(part):
            nodes = list(part.Nodes) if part.Name=='projection_data' else []
            for child in part.SubModelParts: nodes += collect(child)
            return nodes
        self.nodes = sorted(collect(model['skin_model_part'].GetSubModelPart('inner')), key=lambda n:n.Id)

    def response(self, strain, old):
        plastic, alpha = old
        eps = np.array([strain[0], strain[1], 0., strain[2], 0., 0.])
        elastic = np.array([67000., 67000., 67000., 33500., 33500., 33500.])
        trial = elastic*(eps-plastic)
        dev = trial.copy(); dev[:3] -= trial[:3].mean()
        norm = np.sqrt(dev[:3]@dev[:3]+2*dev[3:]@dev[3:])
        f = np.sqrt(1.5)*norm-(201.+10000.*alpha)
        dgamma = max(f, 0.)/(np.sqrt(1.5)*(67000.+2*10000./3))
        normal = dev/norm if norm>0 else np.zeros(6)
        stress = trial-67000.*dgamma*normal
        new_plastic = plastic+dgamma*normal*np.array([1., 1., 1., 2., 2., 2.])
        new_alpha = alpha+np.sqrt(2/3)*dgamma
        return f, stress[[0, 1, 3]], [new_plastic, float(new_alpha)]

    def snapshot(self, state):
        KM, time = self.KM, float(self.part.ProcessInfo[self.KM.TIME])
        pairs = []
        norms = np.linalg.norm(self.points.elastic)
        for item, node in zip(self.points.points, self.nodes):
            c, Bs, Bt, s, t = item
            u = self.points.displacement(c)
            ss, Ds = s.response(Bs@u); st, Dt = t.response(Bt@u)
            fs, reference_s, _ = self.response(Bs@u, self.history[c.Id][0])
            ft, reference_t, _ = self.response(Bt@u, self.history[c.Id][1])
            active_s = bool(np.linalg.norm(Ds-self.points.elastic)/norms>1e-8)
            active_t = bool(np.linalg.norm(Dt-self.points.elastic)/norms>1e-8)
            tolerance = 1e-8*max(201., 201.+10000.*self.history[c.Id][0][1],
                                 201.+10000.*self.history[c.Id][1][1])
            near_s, near_t = bool(abs(fs)<=tolerance), bool(abs(ft)<=tolerance)
            robust_s, robust_t = bool(fs>tolerance), bool(ft>tolerance)
            ns = np.array(c.GetValue(KM.NORMAL))[:2]
            nt = np.array(node.GetValue(KM.NORMAL))[:2]
            N = np.array(c.GetGeometry().ShapeFunctionsValues())[0]
            weight = float(c.GetValue(KM.INTEGRATION_WEIGHT))
            force = np.array(node.GetValue(KM.FORCE))[:2]
            expected_R = (N[:,None]*(traction(ns)@ss+(ns@nt)*(force-traction(nt)@st))[None,:]*weight).ravel()
            expected_K = (N[:,None,None]*(-traction(ns)@Ds@Bs+(ns@nt)*traction(nt)@Dt@Bt)[None,:,:]*weight).reshape(len(u),len(u))
            K, R = KM.Matrix(), KM.Vector(); c.CalculateLocalSystem(K,R,self.part.ProcessInfo)
            center = c.GetGeometry().Center()
            row = dict(time=time,state=state,condition=c.Id, surrogate_x=float(center.X),surrogate_y=float(center.Y),
                true_x=float(node.X0),true_y=float(node.Y0),
                active_plastic_s=active_s,active_plastic_t=active_t,
                yield_function_s=float(fs),yield_function_t=float(ft),yield_tolerance=float(tolerance),
                near_yield_s=near_s,near_yield_t=near_t,
                elastic_s_plastic_t_robust=bool(not near_s and not near_t and not robust_s and robust_t),
                plastic_s_elastic_t_robust=bool(not near_s and not near_t and robust_s and not robust_t),
                alpha_committed_s=self.history[c.Id][0][1],alpha_committed_t=self.history[c.Id][1][1],
                residual_validation_error=float(np.linalg.norm(np.array(R)-expected_R)/max(1.,np.linalg.norm(expected_R))),
                tangent_validation_error=float(np.linalg.norm(np.array(K)-expected_K)/max(1.,np.linalg.norm(expected_K))),
                stress_reference_error=float(max(np.linalg.norm(ss-reference_s),np.linalg.norm(st-reference_t))))
            if not near_s: assert active_s==robust_s, row
            if not near_t: assert active_t==robust_t, row
            pairs.append(row)
        summary = dict(time=time,state=state,iteration=int(self.part.ProcessInfo[KM.NL_ITERATION_NUMBER]),points=len(pairs),
            both_elastic=sum(not r['active_plastic_s'] and not r['active_plastic_t'] for r in pairs),
            both_active_plastic=sum(r['active_plastic_s'] and r['active_plastic_t'] for r in pairs),
            elastic_s_plastic_t=sum(not r['active_plastic_s'] and r['active_plastic_t'] for r in pairs),
            plastic_s_elastic_t=sum(r['active_plastic_s'] and not r['active_plastic_t'] for r in pairs),
            robust_elastic_s_plastic_t=sum(r['elastic_s_plastic_t_robust'] for r in pairs),
            robust_plastic_s_elastic_t=sum(r['plastic_s_elastic_t_robust'] for r in pairs),
            near_yield_pairs=sum(r['near_yield_s'] or r['near_yield_t'] for r in pairs),
            accumulated_mismatch=sum((r['alpha_committed_s']>1e-10)!=(r['alpha_committed_t']>1e-10) for r in pairs),
            max_residual_validation_error=max(r['residual_validation_error'] for r in pairs),
            max_tangent_validation_error=max(r['tangent_validation_error'] for r in pairs),
            max_stress_reference_error=max(r['stress_reference_error'] for r in pairs))
        self.rows.append(summary)
        write_csv(self.out/f'pairs_t{time:.8f}_{state}.csv',pairs)
        (self.out/'steps.json').write_text(json.dumps(self.rows,indent=2)+'\n')
        assert summary['max_residual_validation_error']<1e-8, summary
        assert summary['max_stress_reference_error']<1e-7, summary
        if state=='failed': assert summary['max_tangent_validation_error']<1e-7,summary

    def commit(self):
        for c,Bs,Bt,_,_ in self.points.points:
            u=self.points.displacement(c)
            self.history[c.Id]=[self.response(Bs@u,self.history[c.Id][0])[2],
                                self.response(Bt@u,self.history[c.Id][1])[2]]
        self.points.commit()

    def replay(self, strategy, state):
        matrix, dofs = strategy.GetNonconvergedSolutions()
        values = np.array(matrix)
        saved = [d.GetSolutionStepValue() for d in dofs]
        time = float(self.part.ProcessInfo[self.KM.TIME])
        summaries, detailed = [], []
        try:
            for iteration in range(values.shape[1]):
                for dof in dofs:
                    dof.SetSolutionStepValue(float(values[dof.EquationId,iteration]))
                count = dict(time=time,step_status=state,iteration=iteration,points=len(self.points.points),
                    elastic_s_plastic_t=0,plastic_s_elastic_t=0,both_elastic=0,both_active_plastic=0,
                    robust_elastic_s_plastic_t=0,robust_plastic_s_elastic_t=0,near_yield_pairs=0)
                for c,Bs,Bt,s,t in self.points.points:
                    u = self.points.displacement(c)
                    fs,_,_=self.response(Bs@u,self.history[c.Id][0])
                    ft,_,_=self.response(Bt@u,self.history[c.Id][1])
                    _,Ds=s.response(Bs@u);_,Dt=t.response(Bt@u)
                    ps=bool(np.linalg.norm(Ds-self.points.elastic)/np.linalg.norm(self.points.elastic)>1e-8)
                    pt=bool(np.linalg.norm(Dt-self.points.elastic)/np.linalg.norm(self.points.elastic)>1e-8)
                    tol=1e-8*max(201.,201.+10000.*self.history[c.Id][0][1],201.+10000.*self.history[c.Id][1][1])
                    near=bool(abs(fs)<=tol or abs(ft)<=tol)
                    count['elastic_s_plastic_t']+=int(not ps and pt)
                    count['plastic_s_elastic_t']+=int(ps and not pt)
                    count['both_elastic']+=int(not ps and not pt)
                    count['both_active_plastic']+=int(ps and pt)
                    count['robust_elastic_s_plastic_t']+=int(fs < -tol and ft > tol)
                    count['robust_plastic_s_elastic_t']+=int(fs > tol and ft < -tol)
                    count['near_yield_pairs']+=int(near)
                    detailed.append(dict(time=time,step_status=state,iteration=iteration,condition=c.Id,
                        active_plastic_s=ps,active_plastic_t=pt,yield_function_s=float(fs),yield_function_t=float(ft),
                        near_yield=near,robust_elastic_s_plastic_t=bool(fs < -tol and ft > tol),
                        robust_plastic_s_elastic_t=bool(fs > tol and ft < -tol)))
                summaries.append(count)
        finally:
            for dof,value in zip(dofs,saved):dof.SetSolutionStepValue(float(value))
            # Restore trial responses as well; no material state is committed during replay.
            rhs=self.KM.Vector()
            for c,*_ in self.points.points:c.CalculateRightHandSide(rhs,self.part.ProcessInfo)
        write_csv(self.out/f'iterations_t{time:.8f}_{state}.csv',summaries)
        write_csv(self.out/f'iteration_pairs_t{time:.8f}_{state}.csv',detailed)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--spans',type=int,required=True);args=parser.parse_args()
    study=native_worker.study
    initialize=study.base.Analysis.Initialize
    def wrapped(self):
        initialize(self)
        audit=Audit(study,self.model,HERE/f'audit_n{args.spans}')
        solver=self._GetSolver();solve=solver.SolveSolutionStep
        solver._GetSolutionStrategy().SetUpNonconvergedSolutionsFlag(True)
        def checked():
            try: converged=solve()
            except RuntimeError:
                audit.snapshot('failed')
                audit.replay(solver._GetSolutionStrategy(),'failed')
                raise
            audit.snapshot('converged')
            audit.replay(solver._GetSolutionStrategy(),'converged')
            audit.commit()
            return converged
        solver.SolveSolutionStep=checked
    study.base.Analysis.Initialize=wrapped
    sys.argv=[__file__,'--mode','new','--degree','3','--spans',str(args.spans),'--output',str(HERE)]
    native_worker.main()


if __name__=='__main__':main()
