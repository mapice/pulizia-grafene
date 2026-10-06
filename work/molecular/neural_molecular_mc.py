"""Same-Hamiltonian molecular NPT Metropolis + flexible-nuclei HMC.

Translations, rigid rotations, OH torsion rotations and molecular-center
volume moves avoid slow solvent rotations. HMC updates every internal
coordinate with a reversible volume-preserving leapfrog proposal and an
energy acceptance test. No equilibrium or real-time cleaning inference
comes from the move count. Pressure enters via PV, not an approximate
virial. Volume bounds are reported and must be inactive for qualification.
"""
import os
for v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[v]='1'
os.nice(8)
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np
import torch
torch.set_num_threads(1);torch.set_num_interop_threads(1)
torch.set_default_dtype(torch.float32);torch.mps.set_per_process_memory_fraction(.25)
from ase import units
from ase.io import read,write
from ase.md.verlet import VelocityVerlet
from scipy.spatial.transform import Rotation
from mace.calculators import MACECalculator
from neural_structure import SHA
from fused_radial_interaction import wrap_radial_interactions
from chunk_node_products import wrap_node_products

HERE=Path(__file__).resolve().parent
NMOL=16


def whole_groups(a):
    b=a.copy();length=float(b.cell[0,0])
    assert np.allclose(b.cell.array,np.eye(3)*length)
    xyz=b.positions.reshape(NMOL,12,3)
    delta=xyz-xyz[:,:1,:];delta-=length*np.rint(delta/length)
    assert np.linalg.norm(delta,axis=2).max()<length/2,'Grouping outside molecular lifting domain'
    b.positions=(xyz[:,:1,:]+delta).reshape(-1,3)
    return b


def main():
    p=argparse.ArgumentParser();p.add_argument('--moves',type=int,default=1000);p.add_argument('--seed',type=int,default=20261026)
    p.add_argument('--resume-root',default=None)
    p.add_argument('--allow-interrupted-parent',action='store_true')
    p.add_argument('--initial',default=None)
    p.add_argument('--hmc-min-steps',type=int,default=4)
    p.add_argument('--hmc-max-steps',type=int,default=4)
    p.add_argument('--hmc-probability',type=float,default=.025)
    args=p.parse_args();assert 200<=args.moves<=4000
    assert 4<=args.hmc_min_steps<=args.hmc_max_steps<=32
    assert .025<=args.hmc_probability<=.15
    probabilities={'translation':.35/.90*(.925-args.hmc_probability),
                   'rotation':.30/.90*(.925-args.hmc_probability),
                   'OH_torsion':.25/.90*(.925-args.hmc_probability),
                   'volume':.075,'HMC':args.hmc_probability}
    assert abs(sum(probabilities.values())-1)<1e-12
    default_kernel=(args.hmc_min_steps==args.hmc_max_steps==4 and args.hmc_probability==.025)
    kernel_suffix='' if default_kernel else f'_hmc{args.hmc_min_steps}-{args.hmc_max_steps}_p{args.hmc_probability:.3f}'
    control=json.loads((HERE/'results/molecular-mc-controls/control.json').read_text())
    assert control['accepted_numerical_and_analytic_sampler_control']
    if not default_kernel:
        mixture_control=json.loads((HERE/'results/hmc-mixture-controls/control.json').read_text())
        assert mixture_control['accepted_random_length_sampler_control']
        assert mixture_control['harmonic_modes'][1]['all_mode_precision_criterion_met']
    parent=Path(args.resume_root).resolve() if args.resume_root else None
    if parent:
        if not (parent/'complete.json').exists():
            assert args.allow_interrupted_parent,'Explicit interrupted restart requires terminal evidence'
            stopped=json.loads((parent/'interrupted.json').read_text())
            assert stopped['authoritatively_terminal'] and stopped['guard']['termination_reason'] in ['Wall-time cap','RSS cap']
            try:os.kill(int(stopped['guard']['pid']),0)
            except ProcessLookupError:pass
            else:raise RuntimeError('Parent PID still exists; never resume a live parent')
        previous_plan=json.loads((parent/'plan.json').read_text())
        assert previous_plan['model_sha256']==SHA and previous_plan['molecules']==NMOL
        root=parent.parent/(parent.name+f'_continue{args.moves}'+kernel_suffix)
    else:root=HERE/'results/neural-molecular-mc'/(f'HFIP_N16_seed{args.seed}_{args.moves}moves'+kernel_suffix)
    assert not root.exists(),'Preserve completed/interrupted chain and checkpoint'
    root.mkdir(parents=True)
    source=Path(args.initial).resolve() if args.initial else HERE/'results/neural-bulk/HFIP_N16_seed20261016/dynamics-pilot/NVT/last-state.xyz'
    a=whole_groups(read(source));a.info.update(charge=0,spin=1,external_field=[0.,0.,0.])
    assert len(a)==192 and all(a.pbc)
    previous=None
    if parent:
        checkpoint=np.load(parent/'checkpoint.npz')
        a.positions=checkpoint['positions'];a.set_cell(checkpoint['cell'],scale_atoms=False)
        a.set_momenta(checkpoint['momenta'])
        previous=json.loads((parent/'checkpoint.json').read_text())
    path=HERE/'sources/MACE-POLAR-1-M.model';assert hashlib.sha256(path.read_bytes()).hexdigest()==SHA
    torch.serialization.add_safe_globals([slice]);model=torch.load(path,map_location='cpu',weights_only=False).float()
    with torch.no_grad():model.atomic_energies_fn.atomic_energies.zero_()
    wrap_radial_interactions(model,256);wrap_node_products(model,16)
    calc=MACECalculator(models=model,model_type='PolarMACE',device='mps',default_dtype='float32')
    beta=1/(units.kB*298.15);pressure=1e5*units.Pascal
    rng=np.random.default_rng(args.seed);mass=a.get_masses().reshape(NMOL,12)
    if previous:rng.bit_generator.state=previous['rng_state']
    group_mass=mass.sum(axis=1)
    def energy(atoms):
        batch=calc._atoms_to_batch(atoms).to_dict()
        with torch.no_grad():out=model(batch,compute_force=False,compute_stress=False,compute_virials=False,training=False)
        value=float(out['energy'].item())*calc.energy_units_to_eV
        assert np.isfinite(value)
        return value
    current=energy(a)
    if previous:
        assert abs(current-previous['energy_eV'])<.0003
        current=previous['energy_eV']
    # Independent periodic lifting invariant: one complete molecule moves
    # by one cell vector, hence all interactions must be unchanged.
    shifted=a.copy();shifted.positions[:12,0]+=float(a.cell[0,0])
    translation_error=energy(shifted)-current;assert abs(translation_error)<.0003
    plan=dict(model_sha256=SHA,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              T_K=298.15,P_bar=1,molecules=16,atoms=192,moves=args.moves,seed=args.seed,
              volume_proposal='symmetric logV, centers only; Jacobian plus proposal factor exponent17',
              move_probabilities=probabilities,
              HMC_timestep_fs=.25,HMC_leapfrog_steps_uniform_inclusive=[args.hmc_min_steps,args.hmc_max_steps],
              HMC_length_choice_independent_of_configuration=True,volume_box_bounds_A=[12.2,22.],
              periodic_translation_energy_error_eV=translation_error,MPS_memory_fraction=.25,
              chemical_lifting_domain='Each assigned HFIP group is whole within half a box',
              qualified_for_equilibrium_properties=False,qualified_for_cleaning=False)
    plan['parent_chain']=None if parent is None else str(parent)
    plan['completed_parent_chain']=str(parent) if parent and (parent/'complete.json').exists() else None
    plan['parent_completed_or_confirmed_terminal_interrupted']=(None if parent is None else
                      'complete' if (parent/'complete.json').exists() else 'confirmed interrupted')
    plan['preserved_parent_RNG_positions_cell_momenta']=bool(parent)
    if previous:
        plan['initial_checkpoint_move']=int(previous['next_move'])
        plan['initial_attempt_counts']=previous['attempts']
        plan['initial_accept_counts']=previous['accepted']
        plan['parent_checkpoint_npz_sha256']=hashlib.sha256((parent/'checkpoint.npz').read_bytes()).hexdigest()
    (root/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    attempts={k:0 for k in plan['move_probabilities']};accepts=attempts.copy();bounds=0;rows=[];start=time.monotonic()
    offset=0
    if previous:
        attempts=previous['attempts'].copy();accepts=previous['accepted'].copy()
        bounds=previous['volume_bound_rejections'];offset=int(previous['next_move'])
    for step in range(args.moves+1):
        if step%20==0:
            row=dict(attempt=offset+step,potential_eV=current,volume_A3=float(a.get_volume()),
                     instantaneous_density_g_cm3=float(a.get_masses().sum()/(.602214076*a.get_volume())),
                     attempts=attempts.copy(),accepted=accepts.copy(),volume_bound_rejections=bounds,
                     elapsed_seconds=time.monotonic()-start)
            rows.append(row);write(root/'frames.xyz',a,format='extxyz',append=True)
            (root/'observations.json').write_text(json.dumps(rows,indent=2)+'\n')
            np.savez(root/'checkpoint.npz',positions=a.positions,cell=a.cell.array,momenta=a.get_momenta())
            (root/'checkpoint.json').write_text(json.dumps(dict(next_move=offset+step,energy_eV=current,
                       rng_state=rng.bit_generator.state,attempts=attempts,accepted=accepts,
                       volume_bound_rejections=bounds),indent=2)+'\n')
            print(row,flush=True)
        if step==args.moves:break
        chosen=rng.random()
        cumulative=np.cumsum(list(probabilities.values()))
        kind=list(probabilities)[min(int(np.searchsorted(cumulative,chosen,side='right')),4)]
        attempts[kind]+=1;b=whole_groups(a);logjac=0.;kinetic_delta=0.
        mol=int(rng.integers(NMOL));indices=slice(12*mol,12*(mol+1))
        if kind=='translation':b.positions[indices]+=rng.uniform(-.15,.15,size=3)
        elif kind=='rotation':
            axis=rng.normal(size=3);axis/=np.linalg.norm(axis);R=Rotation.from_rotvec(axis*rng.uniform(-.35,.35)).as_matrix()
            center=np.average(b.positions[indices],axis=0,weights=mass[mol])
            b.positions[indices]=center+(b.positions[indices]-center)@R.T
        elif kind=='OH_torsion':
            xyz=b.positions[indices];assert b.numbers[12*mol]==8 and b.numbers[12*mol+10]==1
            axis=xyz[1]-xyz[0];axis/=np.linalg.norm(axis)
            R=Rotation.from_rotvec(axis*rng.uniform(-np.pi,np.pi)).as_matrix()
            b.positions[12*mol+10]=xyz[0]+R@(xyz[10]-xyz[0])
        elif kind=='volume':
            logratio=float(rng.uniform(-.04,.04));scale=np.exp(logratio/3)
            length=float(b.cell[0,0])*scale
            if not 12.2<length<22:
                bounds+=1;continue
            xyz=b.positions.reshape(NMOL,12,3)
            centers=(xyz*mass[:,:,None]).sum(axis=1)/group_mass[:,None]
            b.positions=(xyz+(scale-1)*centers[:,None,:]).reshape(-1,3)
            b.set_cell(b.cell.array*scale,scale_atoms=False);logjac=(NMOL+1)*logratio
        else:
            # Independent Maxwell momenta; leapfrog + momentum flip makes a
            # symmetric Hamiltonian proposal, correcting finite-step error.
            momentum=rng.normal(size=(len(b),3))*np.sqrt(b.get_masses()[:,None]/beta)
            b.set_momenta(momentum);before=float(b.get_kinetic_energy());b.calc=calc
            # A configuration-independent mixture of reversible HMC kernels
            # preserves the same target distribution. Random length avoids
            # a fixed resonance with internal vibrational modes. Default
            # consumes the same RNG sequence as the original 4-step kernel.
            nsteps=(args.hmc_min_steps if args.hmc_min_steps==args.hmc_max_steps
                    else int(rng.integers(args.hmc_min_steps,args.hmc_max_steps+1)))
            dyn=VelocityVerlet(b,.25*units.fs,logfile=None);dyn.run(nsteps)
            kinetic_delta=float(b.get_kinetic_energy())-before
            b.set_momenta(-b.get_momenta());b=whole_groups(b)
        candidate=energy(b)
        exponent=-beta*(candidate-current+pressure*(b.get_volume()-a.get_volume())+kinetic_delta)+logjac
        if np.log(rng.random())<min(0.,exponent):
            a=b;current=candidate;accepts[kind]+=1
    result=dict(moves=args.moves,attempts=attempts,accepted=accepts,volume_bound_rejections=bounds,
                elapsed_seconds=time.monotonic()-start,scope='Single preliminary model chain; convergence/replicas/size checks still required',
                no_physical_time_from_MC_moves=True,liquid_equilibrium_qualified=False,cleaning_solution_established=False)
    (root/'complete.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
