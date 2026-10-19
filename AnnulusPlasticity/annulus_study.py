"""Run/postprocess annular torsion using the installed Kratos SBM/Gap-SBM."""
import argparse
import csv
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from analytic import evaluate, boundary_traction, RI, RO, END_TIME, C_FINAL

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent/'AnalyticalTestRevised'))
import study as shared
KM = shared.KM


class Analysis(shared.Analysis):
    def EvaluateSolution(self, points, time):
        return evaluate(points, time)

    def BoundaryTraction(self, points, time):
        return boundary_traction(points, time)

    def KeepAdvancingSolutionLoop(self):
        # Floating-point accumulation can leave the last time just below end_time
        # and make the base class run an unintended extra load step.
        tolerance = 32*np.finfo(float).eps*max(abs(self.time),abs(self.end_time),np.finfo(float).tiny)
        return self.time < self.end_time - tolerance

    def InitializeSolutionStep(self):
        time = float(self._GetSolver().GetComputingModelPart().ProcessInfo[KM.TIME])
        skin = self.model['skin_model_part']
        skin.ProcessInfo[KM.TIME] = time
        for node in skin.Nodes:
            point = [[node.X0, node.Y0]]
            value = self.EvaluateSolution(point,time)[0][0]
            node.SetValue(KM.DISPLACEMENT,[float(value[0]),float(value[1]),0.])
            # On a linearized inner skin, nodal NORMAL can point into the material.
            # Prescribe traction with the analytic material-outward normal instead.
            force = self.BoundaryTraction(point,time)[0]
            node.SetValue(KM.FORCE,[float(force[0]),float(force[1]),0.])
        # GapSbmSolidCondition reads nonhistorical DISPLACEMENT on the condition.
        part = self.model['IgaModelPart']
        for element in part.Elements:
            element.SetValue(KM.BODY_FORCE,[0.,0.,0.])
        for condition in part.Conditions:
            center = condition.GetGeometry().Center()
            point = [[center.X,center.Y]]
            if np.hypot(center.X,center.Y) < 1e-12:
                continue
            value = self.EvaluateSolution(point,time)[0][0]
            condition.SetValue(KM.DISPLACEMENT,[float(value[0]),float(value[1]),0.])
            force = self.BoundaryTraction(point,time)[0]
            condition.SetValue(KM.FORCE,[float(force[0]),float(force[1]),0.])
        shared.StructuralMechanicsAnalysis.InitializeSolutionStep(self)


def load_parameters(folder,args,spans):
    filename = args.parameters or f'ProjectParameters_{args.boundary}.json'
    parameters = KM.Parameters((folder/filename).read_text())
    geometry = shared.geometry_settings(parameters)
    for i in range(2):
        geometry['polynomial_order'][i].SetInt(args.degree)
        geometry['number_of_knot_spans'][i].SetInt(spans[i])
    if args.dt is not None:
        parameters['solver_settings']['time_stepping']['time_step'].SetDouble(args.dt)
    if getattr(args,'strategy',None) is not None:
        parameters['solver_settings']['solving_strategy_settings']['type'].SetString(args.strategy)
    if getattr(args,'echo_level',None) is not None:
        parameters['solver_settings']['echo_level'].SetInt(args.echo_level)
    return parameters


def solve(folder,args,spans):
    parameters = load_parameters(folder,args,spans)
    if method_has_gap_load(parameters) and not args.gap_load_from_force:
        raise ValueError('GapSbmLoadSolidCondition requires the FORCE patch and --gap-load-from-force.')
    settings = parameters['solver_settings']
    strategy = settings['solving_strategy_settings']['type'].GetString()
    dt = settings['time_stepping']['time_step'].GetDouble()
    print(f'Solver: {strategy}, dt={dt:g}, max_iterations={settings["max_iteration"].GetInt()}',flush=True)
    with shared.working_directory(folder):
        Path('txt_files').mkdir(exist_ok=True)
        model = KM.Model()
        analysis = Analysis(model,parameters)
        try:
            analysis.Run()
        except RuntimeError as error:
            failure = folder/'results'/args.boundary/f'p{args.degree}'/f'failed_spans_{spans[0]}x{spans[1]}'
            failure.mkdir(parents=True,exist_ok=True)
            (failure/'parameters_used.json').write_text(parameters.PrettyPrintJsonString()+'\n')
            info = dict(method=folder.name,spans=list(spans),degree=args.degree,
                        strategy=strategy,dt=dt,error=str(error))
            if model.HasModelPart('IgaModelPart'):
                process_info = model['IgaModelPart'].ProcessInfo
                if process_info.Has(KM.TIME):
                    info['time'] = float(process_info[KM.TIME])
                if process_info.Has(KM.STEP):
                    info['step'] = int(process_info[KM.STEP])
            (failure/'failure.json').write_text(json.dumps(info,indent=2)+'\n')
            raise RuntimeError(f'{folder.name}, p={args.degree}, mesh={spans}, '
                               f'{strategy}, dt={dt:g}: {error}\nDiagnostics: {failure}') from error
        samples = analysis.samples
        # Save actual runtime condition types, not only names in the JSON.
        from collections import Counter
        counts = Counter(condition_name(c) for c in model['IgaModelPart'].Conditions)
    return samples, parameters, dict(counts)


def method_has_gap_load(parameters):
    return 'GapSbmLoadSolidCondition' in parameters.PrettyPrintJsonString()


def condition_name(condition):
    try:
        return KM.CompareElementsAndConditionsUtility.GetRegisteredName(condition)
    except RuntimeError:
        # Modelers also keep auxiliary geometric conditions with no registered name.
        return condition.Info().split('#')[0].strip()


def plot_geometry(folder,args,output):
    """Initialize actual Kratos modelers; do not run any load step."""
    from matplotlib.collections import LineCollection
    parameters = load_parameters(folder,args,args.spans)
    with shared.working_directory(folder):
        Path('txt_files').mkdir(exist_ok=True)
        model = KM.Model()
        analysis = Analysis(model,parameters)
        analysis.Initialize()
    geo = shared.geometry_settings(parameters)
    lower = np.array(geo['lower_point_xyz'].GetVector())[:2]
    upper = np.array(geo['upper_point_xyz'].GetVector())[:2]
    fig,ax = plt.subplots(figsize=(10,8),constrained_layout=True)
    for axis in range(2):
        knots = np.linspace(lower[axis],upper[axis],args.spans[axis]+1)
        for i,knot in enumerate(knots):
            (ax.axvline if axis == 0 else ax.axhline)(
                knot,color='0.82',lw=.5,zorder=0,
                label='Background mesh' if axis == 0 and i == 0 else None)
    angle = np.linspace(0,2*np.pi,721)
    for radius,label in ((RI,'Skin interna (r=1)'),(RO,'Skin esterna (r=2)')):
        ax.plot(radius*np.cos(angle),radius*np.sin(angle),color='black',lw=1.5,label=label,zorder=4)
    saved = {}
    for side,color in [('inner','purple'),('outer','darkgreen')]:
        name = f'IgaModelPart.surrogate_{side}'
        segments = []
        if model.HasModelPart(name):
            for condition in model[name].Conditions:
                geometry = condition.GetGeometry()
                if len(geometry) == 2:
                    segments.append([[geometry[0].X,geometry[0].Y],
                                     [geometry[1].X,geometry[1].Y]])
        if segments:
            ax.add_collection(LineCollection(segments,colors=color,linewidths=1.1,
                                             label=f'Surrogate {side}',zorder=2))
        saved[f'surrogate_{side}'] = np.asarray(segments).reshape(-1,2,2)
    seen = set()
    for suffix,key,color,marker in [('StructuralAnalysisDomain','bulk','tab:red','.'),
                                    ('GapElements','gap','tab:blue','x')]:
        points = []
        name = 'IgaModelPart.'+suffix
        if model.HasModelPart(name):
            for element in model[name].Elements:
                if element.Id in seen or element.GetValue(KM.ACTIVATION_LEVEL) == 2:
                    continue
                seen.add(element.Id)
                center = element.GetGeometry().Center()
                points.append([center.X,center.Y])
        points = np.asarray(points).reshape(-1,2)
        saved[key+'_points'] = points
        if len(points):
            ax.scatter(*points.T,s=8 if key == 'bulk' else 12,color=color,marker=marker,
                       linewidths=.5,zorder=3,label=f'Centri elementi {key} ({len(points)})')
    ax.set(xlabel='x',ylabel='y',aspect='equal',title=f'{folder.name}: p={args.degree}, mesh {args.spans[0]}×{args.spans[1]}',
           xlim=(lower[0]-.05,upper[0]+.05),ylim=(lower[1]-.05,upper[1]+.05))
    ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),fontsize=9)
    output.mkdir(parents=True,exist_ok=True)
    fig.savefig(output/'geometry.png',dpi=200,bbox_inches='tight')
    np.savez(output/'geometry.npz',**saved,lower=lower,upper=upper,spans=args.spans)
    (output/'geometry_parameters.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    print(f'Geometry saved to {output / "geometry.png"}',flush=True)
    if args.show:
        plt.show()
    plt.close(fig)


def errors(samples):
    exact_u, exact_g, _, _ = evaluate(samples['points'],samples['time'])
    w = samples['weights']
    u = samples['displacement']
    bulk = np.isfinite(samples['gradient']).all(axis=1)
    if not np.isfinite(u).all():
        raise RuntimeError('Nonfinite numerical displacement.')
    l2 = np.sum(w[:,None]*(u-exact_u)**2)
    norm = np.sum(w[:,None]*exact_u**2)
    h1 = np.sum(w[bulk,None]*((u[bulk]-exact_u[bulk])**2).sum(axis=1,keepdims=True))
    h1 += np.sum(w[bulk,None]*(samples['gradient'][bulk]-exact_g[bulk])**2)
    h1norm = np.sum(w[bulk,None]*exact_u[bulk]**2)+np.sum(w[bulk,None]*exact_g[bulk]**2)
    if l2 < 0 or h1 < 0:
        raise RuntimeError('Negative squared error with signed gap quadrature; refine gap integration.')
    return dict(L2_relative=float(np.sqrt(l2/norm)),H1_bulk_relative=float(np.sqrt(h1/h1norm)),
                area=float(w.sum()),H1_bulk_area=float(w[bulk].sum()),time=samples['time'])


def plot_solution(samples,output,show):
    points = samples['points']
    exact_u, _, exact_s, exact_alpha = evaluate(points,samples['time'])
    tri = mtri.Triangulation(*points.T)
    # Mask triangles centered in the hole; draw the exact circular boundaries.
    vertices = points[tri.triangles]
    radius = np.linalg.norm(vertices.mean(axis=1),axis=1)
    tri.set_mask(radius < RI)
    fields = [('ux',samples['displacement'][:,0],exact_u[:,0]),
              ('uy',samples['displacement'][:,1],exact_u[:,1])]
    fields += [(label,samples['stress'][:,i],exact_s[:,i])
               for i,label in enumerate(('sigma_xx','sigma_yy','sigma_xy'))]
    fields += [('equivalent_plastic_strain',samples['equivalent_plastic_strain'],exact_alpha)]
    angle = np.linspace(0,2*np.pi,300)
    for label,numeric,exact in fields:
        fig,axes = plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
        low,high = min(numeric.min(),exact.min()),max(numeric.max(),exact.max())
        if high <= low:
            high = low + 1e-12
        for ax,data,title in zip(axes,(numeric,exact,numeric-exact),('Numerico','Analitico','Errore')):
            settings = {'levels':np.linspace(low,high,31),'extend':'both'} if title != 'Errore' else {'levels':30}
            artist = ax.tricontourf(tri,data,cmap='viridis' if 'plastic' in label else 'coolwarm',**settings)
            for radius in (RI,RO):
                ax.plot(radius*np.cos(angle),radius*np.sin(angle),'k-',lw=.6)
            ax.set(title=f'{label}: {title}',xlabel='x',ylabel='y',aspect='equal')
            fig.colorbar(artist,ax=ax)
        fig.savefig(output/f'{label}.png',dpi=180)
    if show:
        plt.show()
    plt.close('all')


def main(method,action):
    parser = argparse.ArgumentParser(description=f'Annular J2 torsion: {method}, {action}')
    parser.add_argument('--boundary',choices=('dirichlet','mixed'),default='dirichlet')
    parser.add_argument('--parameters',help='Alternative JSON relative to the method directory')
    parser.add_argument('--degree',type=int,default=2)
    parser.add_argument('--spans',type=int,nargs=2,default=[20,20])
    # parser.add_argument('--refinements',type=int,nargs='+',default=list(range(9,60)))
    parser.add_argument('--refinements',type=int,nargs='+',default=list([10,20, 40]))#, 80]))
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--dt',type=float)
    parser.add_argument('--strategy',choices=('newton_raphson','line_search'),
                        help='Override the nonlinear strategy specified in the JSON')
    parser.add_argument('--echo-level',type=int,help='Kratos solver verbosity; 1 prints Newton residuals')
    parser.add_argument('--show',action='store_true')
    parser.add_argument('--gap-load-from-force',action='store_true',
                        help='Confirm that the supplied Gap-SBM FORCE patch is compiled in Kratos.')
    args = parser.parse_args()
    if min([args.degree,args.threads,*args.spans,*args.refinements]) < 1:
        parser.error('Degree, threads and spans must be positive.')
    if args.dt is not None and (args.dt <= 0 or args.dt > END_TIME):
        parser.error('--dt must be in (0, END_TIME].')
    if args.echo_level is not None and args.echo_level < 0:
        parser.error('--echo-level must be non-negative.')
    if action != 'plot' and method == 'gap_sbm' and args.boundary == 'mixed' and not args.gap_load_from_force:
        parser.error('Mixed Gap-SBM requires the supplied FORCE patch and --gap-load-from-force; see README.')
    KM.ParallelUtilities.SetNumThreads(args.threads)
    folder = ROOT/method
    output = folder/'results'/args.boundary/f'p{args.degree}'
    output.mkdir(parents=True,exist_ok=True)
    if action == 'plot':
        plot_geometry(folder,args,output)
        return
    meshes = [(n,n) for n in args.refinements] if action == 'convergence' else [args.spans]
    rows = []
    failures = []
    if action == 'convergence':
        (output/'failures.json').write_text('[]\n')
        for filename in ('convergence.csv','convergence.png'):
            (output/filename).unlink(missing_ok=True)
    for spans in meshes:
        print(f'{method}: {args.boundary}, p={args.degree}, spans={spans}',flush=True)
        try:
            samples,parameters,counts = solve(folder,args,spans)
            metrics = errors(samples)
        except RuntimeError as error:
            if action != 'convergence':
                raise
            failures.append(dict(spans_x=spans[0],spans_y=spans[1],error=str(error)))
            (output/'failures.json').write_text(json.dumps(failures,indent=2)+'\n')
            print(f'Simulazione fallita, salto mesh {spans}: {error}',flush=True)
            continue
        row = dict(spans_x=spans[0],spans_y=spans[1],h=4.4/min(spans),**metrics)
        for name in ('L2_relative','H1_bulk_relative'):
            row[name+'_slope'] = (float(np.log(rows[-1][name]/row[name])/np.log(rows[-1]['h']/row['h']))
                                 if rows and row[name]>0 and rows[-1][name]>0 else float('nan'))
        rows.append(row)
        print(json.dumps(row,indent=2),flush=True)
        target = output/f'spans_{spans[0]}x{spans[1]}' if action == 'convergence' else output
        target.mkdir(exist_ok=True)
        (target/'parameters_used.json').write_text(parameters.PrettyPrintJsonString()+'\n')
        (target/'errors.json').write_text(json.dumps(metrics,indent=2)+'\n')
        (target/'conditions.json').write_text(json.dumps(counts,indent=2)+'\n')
        np.savez(target/'solution.npz',**samples)
        if action == 'run_and_post':
            plot_solution(samples,target,args.show)
        else:
            with (output/'convergence.csv').open('w',newline='') as stream:
                writer = csv.DictWriter(stream,fieldnames=list(row))
                writer.writeheader()
                writer.writerows(rows)
    if action == 'convergence':
        print(f'Convergence: {len(rows)} simulazioni riuscite, {len(failures)} fallite.',flush=True)
        if not rows:
            print(f'Nessuna simulazione riuscita; dettagli in {output/"failures.json"}.',flush=True)
            return
        fig,ax = plt.subplots(constrained_layout=True)
        h_values = np.array([r['h'] for r in rows])
        for name in ('L2_relative','H1_bulk_relative'):
            ax.loglog(h_values,[r[name] for r in rows],'o-',label=name)
        reference_h = np.unique(np.sort(h_values))
        if len(reference_h) == 1:
            reference_h = reference_h[0]*np.array([.5,1.])
        for name,order,style in (('H1_bulk_relative',args.degree,'--'),
                                ('L2_relative',args.degree+1,':')):
            error_values = np.array([r[name] for r in rows])
            valid = np.flatnonzero(np.isfinite(error_values) & (error_values > 0))
            if len(valid):
                anchor = valid[np.argmin(h_values[valid])]
                scale = .7*error_values[anchor]
                anchor_h = h_values[anchor]
            else:
                scale,anchor_h = 1.,reference_h[0]
            reference = scale*(reference_h/anchor_h)**order
            ax.loglog(reference_h,reference,style,color='black',
                      label=rf'$h^{{{order}}}$ (ordine {order})')
        ax.set(xlabel='h',ylabel='Errore relativo',title=f'{method}: annulus {args.boundary}')
        ax.grid(True,which='both',alpha=.3)
        ax.legend()
        fig.savefig(output/'convergence.png',dpi=180)
        plt.show()
        plt.close(fig)
