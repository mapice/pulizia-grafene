"""Tighten all sampled geometries, including initially unconverged starts.

No global-minimum claim. Preserve the initial data for audit.
"""
from clusters_gfn2 import *

def refine(atoms,name,target,steps):
    atoms.calc=calculator()
    opt=BFGS(atoms,logfile=str(OUT/f'{name}.refined.log'),maxstep=.10)
    c=bool(opt.run(fmax=target,steps=steps))
    e=float(atoms.get_potential_energy())
    f=float(np.linalg.norm(atoms.get_forces(),axis=1).max())
    if not np.isfinite(e) or rss_mb()>3000: raise RuntimeError('Resource/numerical limit')
    write(OUT/f'{name}.refined.xyz',atoms)
    return e,c,f

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--pair',default=None);args=parser.parse_args()
    start=time.monotonic()
    if args.pair:
        m=json.loads((OUT/'monomers-refined.json').read_text())
        b=json.loads((OUT/'binding-refined.json').read_text())
        fresh=json.loads((OUT/'binding.json').read_text())['pairs'][args.pair]
        assert args.pair not in b['pairs'],'Preserve existing refined pair'
        b['pairs'][args.pair]=fresh
        names=[]
    else:
        m=json.loads((OUT/'monomers.json').read_text())
        b=json.loads((OUT/'binding.json').read_text());names=list(m)
    for name in names:
        r=m[name]
        a=read(OUT/f'{name}.xyz')
        e,c,f=refine(a,name,.005,250)
        r.update(energy_eV=e,converged=c,max_force_eV_A=f,geometry=name+'.refined.xyz')
        if not c: raise RuntimeError('Unconverged monomer '+name)
    if not args.pair:(OUT/'monomers-refined.json').write_text(json.dumps(m,indent=2)+'\n')
    for key,r in b['pairs'].items():
        if args.pair and key!=args.pair:continue
        trials=[]
        for old in r['trials']:
            if 'energy_eV' not in old: continue
            a=read(OUT/f"{old['name']}.xyz")
            e,c,f=refine(a,old['name'],.008,400)
            new=dict(old,energy_eV=e,converged=c,max_force_eV_A=f)
            # Preserve the intended neutral donor and ensure no bond rupture.
            new['OH_A']=float(np.linalg.norm(a.positions[new['donor_O']]-a.positions[new['donor_H']]))
            new['H_acceptor_A']=float(np.linalg.norm(a.positions[new['acceptor_O']]-a.positions[new['donor_H']]))
            trials.append(new)
        valid=[t for t in trials if t['converged'] and .85<t['OH_A']<1.25]
        if not valid: raise RuntimeError('No converged intact donor '+key)
        best=min(valid,key=lambda t:t['energy_eV'])
        a=read(OUT/f"{best['name']}.refined.xyz")
        n=best['n_acceptor']
        ea=singlepoint(a[:n]); ed=singlepoint(a[n:])
        eb=best['energy_eV']-m[r['acceptor']]['energy_eV']-m[r['donor']]['energy_eV']
        ei=best['energy_eV']-ea-ed
        r.update(trials=trials,best_geometry=best['name']+'.refined.xyz',binding_electronic_kJ_mol=eb*EV_KJMOL,interaction_frozen_fragments_kJ_mol=ei*EV_KJMOL,deformation_kJ_mol=(eb-ei)*EV_KJMOL,OH_A=best['OH_A'],H_acceptor_A=best['H_acceptor_A'],refinement_max_force_eV_A=best['max_force_eV_A'])
        b['refinement']=dict(monomer_fmax_eV_A=.005,complex_fmax_eV_A=.008,global_minimum_proved=False,elapsed_seconds=time.monotonic()-start,peak_rss_mb=rss_mb())
        (OUT/'binding-refined.json').write_text(json.dumps(b,indent=2)+'\n')
        print(key,round(eb*EV_KJMOL,3),[t['converged'] for t in trials],flush=True)
    print(b['refinement'],flush=True)

if __name__=='__main__':main()
