#!/usr/bin/env python3
"""Modelli condizionali di pulizia: nessuna previsione sul campione di Armando.

Python >= 3.10 e NumPy. Eseguire: python model_cleaning.py --out-dir risultati
Le costanti cinetiche sono scenari adimensionali scelti, non dati sperimentali.
L'unico parametro fisico bibliografico usato e' D del PMMA sciolto in acetone.
Il programma verifica conservazione, positivita' e soluzioni analitiche.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np


def exp_generator(generator, duration):
    """Esponenziale di generatore con convenzione colonne, per uniformizzazione.

    La somma di ciascuna colonna deve essere zero. La matrice P=I+G/nu
    e' stocastica. Suddivisione dell'intervallo evita underflow Poisson.
    """
    g = np.asarray(generator, dtype=float)
    if duration < 0:
        raise ValueError('Durata negativa')
    off = g.copy()
    np.fill_diagonal(off, 0)
    if np.min(off) < -1e-15 or np.max(np.abs(g.sum(axis=0))) > 1e-12:
        raise ValueError('Generatore non conservativo o tassi negativi')
    nu = float(np.max(-np.diag(g)))
    eye = np.eye(len(g))
    if nu == 0 or duration == 0:
        return eye
    steps = max(1, math.ceil(nu * duration / 5.0))
    mu = nu * duration / steps
    p = eye + g / nu
    power = eye.copy()
    weight = math.exp(-mu)
    result = weight * eye
    total_weight = weight
    for k in range(1, 1000):
        power = p @ power
        weight *= mu / k
        result += weight * power
        total_weight += weight
        if k > mu and weight < 2e-17:
            break
    else:
        raise RuntimeError('Serie Poisson non convergente')
    if abs(total_weight - 1) > 1e-13:
        raise RuntimeError('Massa Poisson persa')
    return np.linalg.matrix_power(result, steps)


def transport_generator(beta, renewal, release=0.0):
    """Stati x superficie, y bagno, z schermato, w estratto; tutti masse scalate."""
    return np.array([
        [-1.0, beta, release, 0.0],
        [1.0, -beta-renewal, 0.0, 0.0],
        [0.0, 0.0, -release, 0.0],
        [0.0, renewal, 0.0, 0.0],
    ])


def static_closed(beta, tau):
    return (beta + math.exp(-(1+beta)*tau)) / (1+beta)


def fixed_resource_residual(B, tau, n):
    """N bagni: stessi tempo e volume totali; lambda invariata per ipotesi."""
    if n < 1 or B < 0 or tau < 0:
        raise ValueError('Parametri non ammissibili')
    if B == 0:
        return math.exp(-tau)
    loss_per_stage = -math.expm1(-(1+n*B)*tau/n) / (1+n*B)
    if loss_per_stage < 0.5:
        log_fraction = math.log1p(-loss_per_stage)
    else:
        log_fraction = float(np.logaddexp(math.log(n*B), -(1+n*B)*tau/n)) - math.log1p(n*B)
    return math.exp(n * log_fraction)


def fixed_resource_limit(B, tau):
    return math.exp(-tau if B == 0 else math.expm1(-B*tau)/B)


def conformation_generator(a, b, e):
    """L aderente, M mobilizzato, R rimosso: L->M a, M->L b, M->R e.

    R assorbente presuppone asportazione efficace del materiale estratto.
    Non e' un'identificazione chimica dei solventi, ne' una simulazione atomistica.
    """
    return np.array([[-a, b, 0.0], [a, -b-e, 0.0], [0.0, e, 0.0]])


def run_sequence(phases, initial=None):
    state = np.array([1.0, 0.0, 0.0]) if initial is None else np.array(initial, float)
    for g, dt in phases:
        state = exp_generator(g, dt) @ state
    return state


def cycle_residue(d, t, n):
    phase = exp_generator(t, 1/n) @ exp_generator(d, 1/n)
    state = np.linalg.matrix_power(phase, n) @ np.array([1.0, 0.0, 0.0])
    return float(state[0]+state[1])


def plot_latex(curves, cycling, fixed):
    # Tutti i dati sono incorporati: il documento LaTeX non richiede file aggiuntivi.
    colors = ['black', 'blue!75!black', 'red!75!black', 'green!45!black']
    coords = lambda pairs: ' '.join(f'({a:.7g},{b:.7g})' for a,b in pairs)
    lines = [r'\begin{figure}[htbp]', r'\centering',
             r'\begin{tikzpicture}',
             r'\begin{axis}[width=0.94\textwidth,height=6.0cm,xlabel={Tempo adimensionale $\tau=\lambda t$},ylabel={Frazione ancora aderente $x$},xmin=0,xmax=8,ymin=0,ymax=1,grid=major,legend style={font=\footnotesize,at={(0.98,0.98)},anchor=north east},tick label style={font=\small},label style={font=\small}]']
    for color, item in zip(colors, curves):
        lines.append(r'\addplot[thick,'+color+r'] coordinates {'+coords(item['points'])+'};')
        lines.append(r'\addlegendentry{'+item['label']+'}')
    lines += [r'\addplot[black,dashed,domain=0:8,samples=80]{exp(-x)};',
              r'\addlegendentry{Estrazione ideale}', r'\end{axis}',r'\end{tikzpicture}',
              r'\caption{Modello di trasporto: $\beta=1$, con diversi ricambi $r$. Parametri scelti, non stimati sul campione. Il ricambio riduce il riadsorbimento ma non supera la velocit\`a intrinseca di distacco.}',
              r'\label{fig:transport}',r'\end{figure}']
    lines += [r'\begin{figure}[htbp]',r'\centering',r'\begin{tikzpicture}',
              r'\begin{axis}[width=0.94\textwidth,height=6.0cm,xlabel={Numero di alternanze $N$},ylabel={Frazione residua $L+M$},xmin=1,xmax=20,ymin=0,ymax=1,grid=major,legend columns=3,legend style={font=\footnotesize,at={(0.5,1.04)},anchor=south,draw=none},tick label style={font=\small},label style={font=\small}]']
    for color,item in zip(colors,cycling):
        lines.append(r'\addplot[thick,'+color+r',mark=*] coordinates {'+coords(item['points'])+'};')
        label={'complementary':'Complementari','one_way':'Preparazione poi estrazione','same':'Identici'}[item['name']]
        lines.append(r'\addlegendentry{'+label+'}')
    lines += [r'\end{axis}',r'\end{tikzpicture}',
              r'\caption{Modello a stati conformazionali. Ogni curva mantiene invariati il tempo totale nel liquido D e quello nel liquido T. Alternare pi\`u spesso pu\`o aiutare o peggiorare: dipende dai tassi, qui illustrativi. D e T sono stati di trattamento ipotetici, non misure di DMF e THF.}',
              r'\label{fig:cycles}',r'\end{figure}']
    return '\n'.join(lines)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir',type=Path,default=Path('risultati-modello'))
    args = parser.parse_args()
    out = args.out_dir
    out.mkdir(parents=True,exist_ok=True)
    tests = {}
    maximum_mass_error = 0.0
    maximum_closed_error = 0.0
    minimum_entry = 1.0
    # Controlli su ampia griglia di scenari, non ricerca di parametri favorevoli.
    for beta in [0,0.01,0.1,1,10,100]:
        for renewal in [0,0.01,0.1,1,10,100]:
            for release in [0,0.1,1]:
                g = transport_generator(beta,renewal,release)
                for tau in [0,0.01,0.5,2,8]:
                    p = exp_generator(g,tau)
                    maximum_mass_error=max(maximum_mass_error,float(np.max(np.abs(p.sum(axis=0)-1))))
                    minimum_entry=min(minimum_entry,float(p.min()))
                    initial=np.array([1.,0.,0.,0.])
                    state=p@initial
                    assert state[0] >= math.exp(-tau)-1e-11
                    if renewal==0:
                        maximum_closed_error=max(maximum_closed_error,abs(state[0]-static_closed(beta,tau)))
    assert maximum_mass_error < 1e-10
    assert maximum_closed_error < 1e-10
    assert minimum_entry >= -1e-14
    tests.update(mass_error=maximum_mass_error,closed_form_error=maximum_closed_error,minimum_transition_entry=minimum_entry)

    g = conformation_generator(1,0,1)
    t=2.3
    state=run_sequence([(g,t)])
    expected=np.array([math.exp(-t),t*math.exp(-t),1-(1+t)*math.exp(-t)])
    tests['erlang_chain_error']=float(np.max(np.abs(state-expected)))
    assert tests['erlang_chain_error']<1e-12
    limit_error=0.0
    for B in [0,0.01,0.1,1,10]:
        for tau in [0.1,1,5,10]:
            limit_error=max(limit_error,abs(fixed_resource_residual(B,tau,1000000)-fixed_resource_limit(B,tau)))
    tests['finite_stage_limit_error_N_1e6']=limit_error
    assert limit_error<1e-6
    assert fixed_resource_residual(0,100,1)==math.exp(-100)
    assert fixed_resource_residual(0,100,2)==math.exp(-100)
    assert math.isfinite(fixed_resource_residual(1e-20,100,1))

    curves=[]
    for renewal in [0,0.2,2,20]:
        g=transport_generator(1,renewal)
        points=[]
        for tau in np.linspace(0,8,81):
            x=float((exp_generator(g,float(tau))@np.array([1.,0.,0.,0.]))[0])
            points.append([float(tau),x])
        curves.append(dict(r=renewal,label=f'$r={renewal:g}$',points=points))

    # Scenari cinetici dichiarati a priori. Nessuna stima di tassi fisici.
    scenarios=[
        dict(name='complementary',label='Tassi complementari',D=[2,4,0.02],T=[0.02,0.1,4]),
        dict(name='one_way',label='Una preparazione seguita da estrazione',D=[2,0,0],T=[0,0,2]),
        dict(name='same',label='Tassi uguali nei due liquidi',D=[0.4,0.5,0.3],T=[0.4,0.5,0.3]),
    ]
    cycling=[]
    for sc in scenarios:
        d=conformation_generator(*sc['D'])
        t=conformation_generator(*sc['T'])
        points=[[n,cycle_residue(d,t,n)] for n in range(1,21)]
        item=dict(sc,points=points,DT_residual=points[0][1],TD_residual=float(run_sequence([(t,1),(d,1)])[:2].sum()),commutator_norm=float(np.linalg.norm(t@d-d@t)))
        cycling.append(item)
    # Nessuna differenza se i generatori sono identici.
    assert max(p[1] for p in cycling[2]['points'])-min(p[1] for p in cycling[2]['points'])<1e-11
    assert cycling[0]['points'][-1][1]<cycling[0]['points'][0][1]
    assert cycling[1]['points'][-1][1]>cycling[1]['points'][0][1]

    ipa_examples=[]
    d=conformation_generator(2,4,0.02)
    t=conformation_generator(0.02,0.1,4)
    before=run_sequence([(d,1)])
    for name,values in [('mobilizing',[4,0.02,0]),('relocking',[0,4,0])]:
        ipa=conformation_generator(*values)
        state=run_sequence([(ipa,8/30),(t,1)],before)
        baseline=run_sequence([(t,1)],before)
        ipa_examples.append(dict(hypothesis=name,rates=values,with_IPA_residual=float(state[:2].sum()),without_IPA_residual=float(baseline[:2].sum()),warning='Confronto di meccanismi ipotetici; durata totale diversa. Per una prova causale servono controlli temporali abbinati.'))

    fixed=[dict(B=B,tau=5,residual_by_N={str(n):fixed_resource_residual(B,5,n) for n in [1,2,4,8,16,100]},limit=fixed_resource_limit(B,5)) for B in [0.1,1,10]]
    D=2.93e-11
    diffusion=[dict(delta_um=d,t_seconds=(d*1e-6)**2/D) for d in [10,100,1000]]
    redeposit=dict(formula='t_nm <= 0.001 * c_mg_L * h_um / rho_g_cm3',illustrative_c_mg_L=1,illustrative_h_um=10,assumed_rho_g_cm3=1.18,t_nm=0.001*1*10/1.18,c_needed_for_0_5_nm_at_h10=0.5*1.18/(0.001*10),scope='Parametro rho assunto; c e h non misurati. Limite di bilancio sul residuo medio, non previsione locale.')
    report=dict(status='SIMULAZIONI CONDIZIONALI, NON CALIBRATE SUL CAMPIONE',numpy_version=np.__version__,tests=tests,transport=curves,cycling=cycling,IPA_opposing_hypotheses=ipa_examples,fixed_resources=fixed,diffusion=dict(D_m2_s=D,source='Arai et al. Macromolecules 1996, DOI10.1021/ma951274n, tabella3 PMMA atattico Mw9.52e5, acetone25C',warning='Diffusione di catene gia sciolte, non desorbimento o penetrazione nel residuo. Le distanze sono ipotetiche.',values=diffusion),redeposition_bound=redeposit)
    (out/'risultati.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    (out/'figure-incorporate.tex').write_text(plot_latex(curves,cycling,fixed))
    print(json.dumps(dict(tests=tests,cycling=[{k:v for k,v in x.items() if k!='points'} for x in cycling],fixed_resources=fixed,IPA=ipa_examples),indent=2,ensure_ascii=False))


if __name__=='__main__':
    main()
