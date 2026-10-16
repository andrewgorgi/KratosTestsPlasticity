"""Runner comune: gli script nelle tre cartelle fissano il metodo geometrico."""
import argparse
import csv
import faulthandler
import json
import math
import os
from contextlib import contextmanager
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from analytic_fields import FIELDS, PHYSICS, evaluator

ROOT = Path(__file__).resolve().parent

@contextmanager
def working_directory(path):
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)

class Analysis(StructuralMechanicsAnalysis):
    def Initialize(self):
        super().Initialize()
        solver = self._GetSolver()
        solve_solution_step = solver.SolveSolutionStep

        def solve_checked():
            converged = solve_solution_step()
            if not converged:
                raise RuntimeError(f'Newton non convergente al tempo {self.time:g}; '
                                   'nessuna norma di errore viene calcolata.')
            return converged

        solver.SolveSolutionStep = solve_checked

    def InitializeSolutionStep(self):
        if self.model.HasModelPart('skin_model_part'):
            skin = self.model['skin_model_part']
            time = self._GetSolver().GetComputingModelPart().ProcessInfo[KM.TIME]
            skin.ProcessInfo[KM.TIME] = time
            name = self.project_parameters['problem_data']['problem_name'].GetString()
            physics = next(p for p in PHYSICS if name.endswith('_'+p))
            field = next(f for f in FIELDS if name.endswith('_'+f+'_'+physics))
            nodes = list(skin.Nodes)
            if nodes:
                values, _, stresses = evaluator(field,physics)([[n.X0,n.Y0] for n in nodes],time)
                for node,value,stress in zip(nodes,values,stresses):
                    # SbmSolidCondition legge il DISPLACEMENT non storico del nodo proiettato.
                    node.SetValue(KM.DISPLACEMENT, [float(value[0]),float(value[1]),0.])
                    # SbmLoadSolidCondition legge la trazione sul bordo vero.
                    # Il medesimo stress manufactured alimenta carico e confronto.
                    normal = node.GetValue(KM.NORMAL)
                    node.SetValue(KM.FORCE, [
                        float(stress[0]*normal[0]+stress[2]*normal[1]),
                        float(stress[2]*normal[0]+stress[1]*normal[1]), 0.])
        super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        # Preserve final samples before Finalize/Clear; no fabricated skin for body-fitted.
        self.samples = sample_model(self.model, self.project_parameters)

def geometry_settings(parameters):
    matches = [m['Parameters'] for m in parameters['modelers'].values() if m['Parameters'].Has('number_of_knot_spans')]
    if len(matches) != 1:
        raise ValueError('Serve esattamente un modeler con number_of_knot_spans.')
    return matches[0]

def sample_model(model, parameters):
    """L2 bulk+gap; H1 sul bulk, la cui derivata e esposta in Python."""
    parts = ['IgaModelPart.StructuralAnalysisDomain']
    if model.HasModelPart('IgaModelPart.GapElements'):
        parts.append('IgaModelPart.GapElements')
    points, disp, grad, weights, stresses, seen = [], [], [], [], [], set()
    is_plastic = parameters['problem_data']['problem_name'].GetString().endswith(('_plastic', '_plastic_hardening'))
    plastic_strains = []
    for name in parts:
        for elem in model[name].Elements:
            if elem.Id in seen or elem.GetValue(KM.ACTIVATION_LEVEL) == 2:
                continue
            seen.add(elem.Id)
            geom = elem.GetGeometry()
            is_gap = name.endswith('GapElements')
            if is_gap:
                # La build locale esporta u e strain al punto vero, ma non la geometria
                # NEIGHBOUR_GEOMETRIES in Python: nessuna interpolazione inventata.
                center = geom.Center()
                point = [center.X, center.Y]
                displacement = np.array(elem.GetValue(KM.DISPLACEMENT))[:2]
                gradient = np.full(4, np.nan)
            else:
                N = np.array(geom.ShapeFunctionsValues())[0]
                dN = np.array(geom.ShapeFunctionDerivatives(1,0))
                nodal = np.array([[n.GetSolutionStepValue(KM.DISPLACEMENT_X), n.GetSolutionStepValue(KM.DISPLACEMENT_Y)] for n in geom])
                point = N @ np.array([[n.X0,n.Y0] for n in geom])
                displacement = N @ nodal
                gradient = (nodal.T @ dN).reshape(-1)
            weight = float(elem.GetValue(KM.INTEGRATION_WEIGHT))
            if weight <= 0 and not is_gap:
                raise RuntimeError(f'Peso non positivo per elemento {elem.Id}: {weight}')
            points.append(point)
            disp.append(displacement)
            grad.append(gradient)
            weights.append(weight)
            stresses.append([elem.GetValue(v) for v in (IGA.CAUCHY_STRESS_XX,IGA.CAUCHY_STRESS_YY,IGA.CAUCHY_STRESS_XY)])
            if is_plastic:
                values = elem.CalculateOnIntegrationPoints(
                    CLA.EQUIVALENT_PLASTIC_STRAIN, model['IgaModelPart'].ProcessInfo)
                if len(values) != 1 or not np.isfinite(values[0]):
                    raise RuntimeError(f'Deformazione plastica equivalente non valida per elemento {elem.Id}: {values}')
                plastic_strains.append(float(values[0]))
    if not points:
        raise RuntimeError('Nessun punto di integrazione attivo.')
    time = model['IgaModelPart'].ProcessInfo[KM.TIME]
    samples = {'points':np.array(points), 'displacement':np.array(disp), 'gradient':np.array(grad),
            'weights':np.array(weights), 'stress':np.array(stresses), 'time':float(time)}
    if is_plastic:
        samples['equivalent_plastic_strain'] = np.array(plastic_strains)
    return samples

def errors(samples, field, physics):
    exact_u, exact_g, _ = evaluator(field,physics)(samples['points'], samples['time'])
    w = samples['weights']
    if not np.all(np.isfinite(samples['displacement'])):
        raise RuntimeError('Spostamenti non finiti: il solve non ha prodotto una soluzione valida.')
    def integral(array):
        return np.sum(w * np.sum(array**2,axis=1))
    norm_u = integral(exact_u)
    error_u = integral(samples['displacement']-exact_u)
    bulk = np.all(np.isfinite(samples['gradient']), axis=1)
    def bulk_integral(array):
        return np.sum(w[bulk]*np.sum(array[bulk]**2,axis=1))
    bulk_error = bulk_integral(samples['displacement']-exact_u) + bulk_integral(samples['gradient']-exact_g)
    bulk_norm = bulk_integral(exact_u) + bulk_integral(exact_g)
    return {'L2_relative':float(np.sqrt(error_u/norm_u)),
            'H1_bulk_relative':float(np.sqrt(bulk_error/bulk_norm)),
            'area':float(w.sum()), 'H1_bulk_area':float(w[bulk].sum()), 'time':samples['time']}

def load_parameters(folder,args):
    path = folder/(args.parameters or f'ProjectParameters_{args.field}_{args.physics}.json')
    parameters = KM.Parameters(path.read_text())
    expected = f'{folder.name}_{args.field}_{args.physics}'
    if parameters['problem_data']['problem_name'].GetString() != expected:
        raise ValueError(f'Il JSON deve avere problem_name={expected}; allineare --field e --physics.')
    if args.spans:
        for i, value in enumerate(args.spans):
            geometry_settings(parameters)['number_of_knot_spans'][i].SetInt(value)
    if args.degree:
        for i in range(2):
            geometry_settings(parameters)['polynomial_order'][i].SetInt(args.degree)
    return parameters

def solve(folder,parameters):
    with working_directory(folder):
        Path('txt_files').mkdir(exist_ok=True)
        model = KM.Model()
        analysis = Analysis(model,parameters)
        analysis.Run()
        return analysis.samples, model

def run_post(folder,args):
    parameters = load_parameters(folder,args)
    samples, _ = solve(folder,parameters)
    output = folder/'results'/f'{args.field}_{args.physics}'
    output.mkdir(parents=True,exist_ok=True)
    metrics = errors(samples,args.field,args.physics)
    (output/'parameters_used.json').write_text(parameters.PrettyPrintJsonString()+'\n')
    (output/'errors.json').write_text(json.dumps(metrics,indent=2)+'\n')
    np.savez(output/'solution.npz', **samples)
    exact_u, _, exact_s = evaluator(args.field,args.physics)(samples['points'],samples['time'])
    # Stessi punti e stessa scala per analitico/numerico, senza deformare le coordinate.
    values = [('ux',samples['displacement'][:,0],exact_u[:,0]),
              ('uy',samples['displacement'][:,1],exact_u[:,1])]
    if samples['stress'].shape[1] >= 3:
        values += [(label,samples['stress'][:,i],exact_s[:,i]) for i,label in enumerate(('sigma_xx','sigma_yy','sigma_xy'))]
    for label,numeric,exact in values:
        fig,axes = plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
        low,high = min(numeric.min(),exact.min()),max(numeric.max(),exact.max())
        for ax,data,title in zip(axes,(numeric,exact,numeric-exact),('Numerico','Analitico','Errore')):
            settings = {'vmin':low,'vmax':high} if title != 'Errore' else {}
            artist = ax.scatter(*samples['points'].T,c=data,s=9,cmap='coolwarm',**settings)
            ax.set(title=f'{label}: {title}',xlabel='x',ylabel='y',aspect='equal')
            fig.colorbar(artist,ax=ax)
        fig.savefig(output/f'{label}.png',dpi=180)
    if args.physics in ('plastic', 'plastic_hardening'):
        fig,ax = plt.subplots(constrained_layout=True)
        plastic_strain = samples['equivalent_plastic_strain']
        # Contour sui punti di integrazione; livelli validi anche se il campo e nullo.
        upper = max(float(plastic_strain.max()), 1e-12)
        artist = ax.tricontourf(*samples['points'].T,plastic_strain,
                               levels=np.linspace(0,upper,31),cmap='viridis',vmin=0)
        ax.set(title=f'Deformazione plastica equivalente: t={samples["time"]:g}',
               xlabel='x',ylabel='y',aspect='equal')
        fig.colorbar(artist,ax=ax,label='Deformazione plastica equivalente')
        fig.savefig(output/'equivalent_plastic_strain.png',dpi=180)
    print(json.dumps(metrics,indent=2))
    plt.show()

def convergence(folder,args):
    rows=[]
    for n in args.refinements:
        parameters=load_parameters(folder,args)
        geo=geometry_settings(parameters)
        for i in range(2):
            geo['number_of_knot_spans'][i].SetInt(n)
        print(f'Convergence: {folder.name}, {args.field}, {args.physics}, '
              f'degree={args.degree}, spans={n}x{n}',flush=True)
        samples,_=solve(folder,parameters)
        widths=[geo['upper_point_xyz'][i].GetDouble()-geo['lower_point_xyz'][i].GetDouble() for i in range(2)]
        row={'spans':n,'h':max(widths)/n,**errors(samples,args.field,args.physics)}
        rows.append(row)
        print(row,flush=True)
        print(f'\nRisultati cumulati ({len(rows)} simulazioni; '
              'H1_relative = H1_bulk_relative):\n'
              f'h = {[r["h"] for r in rows]}\n'
              f'L2_relative = {[r["L2_relative"] for r in rows]}\n'
              f'H1_relative = {[r["H1_bulk_relative"] for r in rows]}\n',
              flush=True)
    for i,row in enumerate(rows):
        for name in ('L2_relative','H1_bulk_relative'):
            row[name+'_slope']=float('nan')
            if i and row[name]>0 and rows[i-1][name]>0:
                row[name+'_slope']=math.log(row[name]/rows[i-1][name])/math.log(row['h']/rows[i-1]['h'])
    output=folder/'results'/f'{args.field}_{args.physics}'
    output.mkdir(parents=True,exist_ok=True)
    with (output/'convergence.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    fig,ax=plt.subplots()
    h_values = np.array([r['h'] for r in rows])
    for name,order in (('L2_relative',args.degree+1),('H1_bulk_relative',args.degree)):
        error_values = np.array([r[name] for r in rows])
        ax.loglog(h_values,error_values,'o-',label=name)
        valid = np.isfinite(error_values) & (error_values > 0)
        if np.any(valid):
            # Ancoraggio vicino alla norma sulla mesh piu fine, con un piccolo offset.
            anchor = np.flatnonzero(valid)[np.argmin(h_values[valid])]
            reference_h = np.sort(h_values)
            reference = 0.7*error_values[anchor]*(reference_h/h_values[anchor])**order
            ax.loglog(reference_h,reference,'--',color='black',label='_nolegend_')
            if reference_h[-1] > reference_h[0]:
                # Cateti con rapporto 1:order nelle coordinate logaritmiche.
                left = reference_h[0]
                right = left*(reference_h[-1]/left)**0.3
                bottom = 0.7*error_values[anchor]*(left/h_values[anchor])**order
                top = bottom*(right/left)**order
                ax.plot([left,right,right,left],[bottom,bottom,top,bottom],
                        color='black',linewidth=1,label='_nolegend_')
                ax.annotate('1',(np.sqrt(left*right),bottom),xytext=(0,-12),
                            textcoords='offset points',ha='center',va='top',color='black')
                ax.annotate(str(order),(right,np.sqrt(bottom*top)),xytext=(5,0),
                            textcoords='offset points',ha='left',va='center',color='black')
    ax.set(xlabel='h = max(hx,hy)',ylabel='Errore relativo',title=f'{folder.name}: {args.field}, {args.physics}')
    ax.grid(True,which='both'); ax.legend()
    fig.savefig(output/'convergence.png',dpi=180)
    plt.show()

def plot_geometry(folder,args):
    parameters=load_parameters(folder,args)
    with working_directory(folder):
        Path('txt_files').mkdir(exist_ok=True)
        model=KM.Model(); analysis=Analysis(model,parameters)
        analysis.Initialize()
        geo=geometry_settings(parameters)
        fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
        for i in range(2):
            knots=np.linspace(geo['lower_point_xyz'][i].GetDouble(),geo['upper_point_xyz'][i].GetDouble(),geo['number_of_knot_spans'][i].GetInt()+1)
            for knot in knots:
                (ax.axvline if i==0 else ax.axhline)(knot,color='0.8',lw=.5)
        used=set()
        for part,color,label in [('IgaModelPart.surrogate_outer','darkred','Surrogate'),('skin_model_part.outer','blue','Contorno vero'),
                                  ('IgaModelPart.SBM_Support_outer','green','Supporto')]:
            if model.HasModelPart(part):
                for condition in model[part].Conditions:
                    geometry=condition.GetGeometry()
                    legend_label=label if label not in used else None
                    if len(geometry)==2:
                        ax.plot([geometry[0].X,geometry[1].X], [geometry[0].Y,geometry[1].Y],
                                color=color,label=legend_label)
                    else:
                        center=geometry.Center()
                        ax.plot(center.X,center.Y,'.',color=color,label=legend_label)
                    used.add(label)
        projection_segments = []
        for condition in model['IgaModelPart'].Conditions:
            if condition.Has(IGA.PROJECTION_NODE_COORDINATES):
                center = condition.GetGeometry().Center()
                projection = condition.GetValue(IGA.PROJECTION_NODE_COORDINATES)
                projection_segments.append(((center.X,center.Y),(projection[0],projection[1])))
        if projection_segments:
            # Coordinate salvate dalle condizioni Kratos dopo Initialize:
            # sono le proiezioni effettivamente usate, inclusi gli angoli.
            segments = np.array(projection_segments)
            ax.add_collection(LineCollection(segments,color='purple',linewidths=.6,
                                             linestyles='--',alpha=.7,zorder=2,
                                             label='Collegamenti bordo–skin'))
            projected_points = np.unique(segments[:,1,:],axis=0)
            ax.scatter(*projected_points.T,color='purple',marker='x',s=24,zorder=4,
                       label='Punti proiettati sulla skin')
        integration_points = []
        seen_elements = set()
        for part in ('IgaModelPart.StructuralAnalysisDomain','IgaModelPart.GapElements'):
            if model.HasModelPart(part):
                for element in model[part].Elements:
                    if element.Id in seen_elements or element.GetValue(KM.ACTIVATION_LEVEL) == 2:
                        continue
                    seen_elements.add(element.Id)
                    center = element.GetGeometry().Center()
                    integration_points.append((center.X,center.Y))
        if integration_points:
            ax.scatter(*np.array(integration_points).T,color='red',s=10,zorder=3,
                       label='Punti di integrazione degli elementi attivi')
        # Contorno fisico comune, incluso body-fitted senza skin_model_part.
        ax.plot([0,4,4,0,0],[0,0,2,2,0],'b-',lw=1,label='Dominio fisico [0,4] x [0,2]')
        ax.set(xlabel='x',ylabel='y',aspect='equal',title=folder.name)
        ax.legend(loc='upper left',bbox_to_anchor=(1.02,1))
        output=folder/'results'/f'{args.field}_{args.physics}'; output.mkdir(parents=True,exist_ok=True)
        fig.savefig(output/'geometry.png',dpi=180)
        plt.show()

def main(method,action):
    faulthandler.enable()
    parser=argparse.ArgumentParser(description=f'{method}: {action}')
    parser.add_argument('--field',choices=FIELDS,default='linear')
    parser.add_argument('--physics',choices=PHYSICS,default='elastic')
    parser.add_argument('--parameters',help='File JSON alternativo, relativo alla cartella del metodo')
    parser.add_argument('--degree',type=int,help='Grado delle funzioni di base; se omesso viene richiesto per le simulazioni')
    parser.add_argument('--threads',type=int,
                        help='Numero di thread Kratos (default dell’ambiente se omesso)')
    parser.add_argument('--spans',type=int,nargs=2)
    # parser.add_argument('--refinements',type=int,nargs='+',default=list(range(10,60)))
    parser.add_argument('--refinements',type=int,nargs='+',default=[10,20,40])
    parser.add_argument('--show',action='store_true')
    args=parser.parse_args()
    if args.threads is not None:
        if args.threads < 1:
            parser.error('--threads deve essere positivo')
        KM.ParallelUtilities.SetNumThreads(args.threads)
    print(f'Thread Kratos: {KM.ParallelUtilities.GetNumThreads()}',flush=True)
    
    if args.degree is not None and args.degree < 1:
        parser.error('--degree deve essere positivo')
    if any(n<1 for n in args.refinements) or (args.spans and min(args.spans)<1):
        parser.error('Il numero di knot span deve essere positivo')
    if action in ('run_and_post','convergence') and args.degree is None:
        while args.degree is None:
            try:
                degree = int(input('Inserisci il grado delle funzioni di base (intero positivo): '))
            except ValueError:
                print('Valore non valido: inserisci un intero positivo.')
                continue
            except EOFError:
                parser.error('Input non disponibile: specifica il grado con --degree.')
            if degree < 1:
                print('Il grado deve essere positivo.')
                continue
            args.degree = degree
    {'run_and_post':run_post,'convergence':convergence,'plot':plot_geometry}[action](ROOT/method,args)
