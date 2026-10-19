"""Exact first-loading strain Taylor audit on previously recorded SBM pairs.
No Kratos solve, no constitutive-state modification.
"""
import csv,json,sys
from pathlib import Path
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
import cyclic_analytic as exact
x,y,dx,dy,A,B=sp.symbols('x y dx dy A B',real=True)
r2=x*x+y*y
gamma=A/r2+B
strain=sp.Matrix([-gamma*x*y/r2,gamma*x*y/r2,gamma*(x*x-y*y)/r2])
terms=[]
v=strain
for k in range(3):
    terms.append(sp.lambdify((x,y,dx,dy,A,B),v/sp.factorial(k),'numpy'))
    v=v.diff(x)*dx+v.diff(y)*dy

def coeff(point,t):
    c=exact.C_FINAL*t/.001
    if np.linalg.norm(point)<np.sqrt(np.sqrt(3)*c/exact.YIELD):
        return c/exact.G+3*c/exact.H,-np.sqrt(3)*exact.YIELD/exact.H
    return c/exact.G,0.

def eps(point,t):
    grad=exact.evaluate(np.asarray(point).reshape(1,2),t)[1][0]
    return np.array([grad[0],grad[3],grad[1]+grad[2]])

summary=[];rows=[]; example=None
for n in (20,40):
    path=HERE.parent/f'sbm_p3_pair_phases/audit_n{n}/pairs_t0.00075000_failed.csv'
    pairs=list(csv.DictReader(path.open()))
    for t in (.00025,.0005,.00075):
        rp=np.sqrt(np.sqrt(3)*exact.C_FINAL*t/.001/exact.YIELD)
        block=[]
        for row in pairs:
            s=np.array([float(row['surrogate_x']),float(row['surrogate_y'])])
            true=np.array([float(row['true_x']),float(row['true_y'])]);d=true-s
            roots=np.roots([d@d,2*s@d,s@s-rp*rp]) if d@d else []
            crossing=any(abs(np.imag(z))<1e-12 and 1e-10<float(np.real(z))<1-1e-10 for z in roots)
            knots=np.linspace(-2.2,2.2,n+1)[1:-1]
            knot_cross=any(np.any((knots-min(s[i],true[i])>1e-10)&(max(s[i],true[i])-knots>1e-10)) for i in (0,1))
            aa,bb=coeff(s,t); target=eps(true,t);approx=np.zeros(3)
            item=dict(mesh=n,time=t,condition=row['condition'],crosses_exact_front=bool(crossing),crosses_knot=bool(knot_cross))
            for k,term in enumerate(terms):
                approx+=np.asarray(term(*s,*d,aa,bb)).ravel()
                item[f'error_q{k}']=float(np.linalg.norm(approx-target)/np.linalg.norm(target))
            rows.append(item);block.append(item)
            if crossing and t==.0005 and example is None:example=(s,true,t)
        for crossing in (False,True):
            group=[r for r in block if r['crosses_exact_front']==crossing]
            summary.append(dict(mesh=n,time=t,front_radius=rp,crosses_front=crossing,count=len(group),
                pairs_crossing_knots=sum(r['crosses_knot'] for r in group),
                median_relative_errors=[float(np.median([r[f'error_q{k}'] for r in group])) if group else None for k in range(3)]))
with (HERE/'pairs.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
if example:
    s,true,t=example;d=true-s;aa,bb=coeff(s,t)
    ss=np.linspace(0,1,401);ref=np.array([eps(s+v*d,t) for v in ss])
    q0=np.array(terms[0](*s,*d,aa,bb)).ravel()
    q1=np.array(terms[1](*s,*d,aa,bb)).ravel()
    q2=np.array(terms[2](*s,*d,aa,bb)).ravel()
    component=np.argmax(np.ptp(ref,axis=0))
    fig,ax=plt.subplots(figsize=(8,4))
    ax.plot(ss,ref[:,component],label='Strain esatto')
    ax.plot(ss,(q0+ss[:,None]*q1)[:,component],'--',label='Taylor ordine 1')
    ax.plot(ss,(q0+ss[:,None]*q1+ss[:,None]**2*q2)[:,component],':',label='Taylor ordine 2')
    rp=np.sqrt(np.sqrt(3)*exact.C_FINAL*t/.001/exact.YIELD)
    for z in np.roots([d@d,2*s@d,s@s-rp*rp]):
        if abs(np.imag(z))<1e-12 and 0<float(np.real(z))<1:ax.axvline(float(np.real(z)),color='grey',label='Fronte plastico')
    ax.set(xlabel='Posizione sul segmento: 0 surrogate, 1 true',ylabel=['epsilon_xx','epsilon_yy','gamma_xy'][component],title=f'Soluzione esatta, primo carico, t={t:g}')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(HERE/'strain_segment.png',dpi=160)
print(json.dumps(summary,indent=2))
