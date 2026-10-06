"""Audit actual biased windows before allowing any free-energy inference.

No numerical bridging of unvisited gaps. Overlap of the distance alone does
not prove equilibration of contact/rotor variables. MBAR estimator is tested
on independent analytic harmonic samples, not calibrated to make a solvent win.
"""
from pathlib import Path
import argparse, json
import numpy as np
from scipy.special import logsumexp
from scipy.optimize import minimize
HERE = Path(__file__).resolve().parent
RT = .00831446261815324 * 298.15


def statistical_inefficiency(x):
    x = np.asarray(x) - np.mean(x)
    n = len(x)
    if np.dot(x, x) < 1e-20:
        return float(n)
    size = 1 << (2*n - 1).bit_length()
    fft = np.fft.rfft(x, n=size)
    corr = np.fft.irfft(fft * fft.conjugate(), n=size)[:n]
    corr /= np.arange(n, 0, -1)
    corr /= corr[0]
    total = 0.
    for i in range(1, n//5 - 1, 2):
        pair = corr[i]+corr[i+1]
        if pair <= 0:
            break
        total += pair
    return float(min(n, max(1., 1+2*total)))


def mbar(samples, centers, stiffness):
    counts = np.array([len(x) for x in samples], dtype=float)
    z = np.concatenate(samples)
    u = .5*np.asarray(stiffness)[:, None]*(z[None, :]-np.asarray(centers)[:, None])**2/RT

    def objective(reduced):
        f = np.r_[0., reduced]
        scores = np.log(counts)[:, None]+f[:, None]-u
        denom = logsumexp(scores, axis=0)
        probabilities = np.exp(scores-denom)
        value = denom.sum()-np.dot(counts, f)
        gradient = probabilities.sum(axis=1)-counts
        return value, gradient[1:]
    fit = minimize(objective, np.zeros(len(counts)-1), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=5000, ftol=1e-15, gtol=1e-7, maxls=50))
    f = np.r_[0., fit.x]
    logw = -logsumexp(np.log(counts)[:, None]+f[:, None]-u, axis=0)
    weights = np.exp(logw-logsumexp(logw))
    _, gradient = objective(fit.x)
    assert np.max(abs(gradient)/counts[1:]) < 2e-6, fit.message
    return z, weights, f


def positive_control():
    rng = np.random.default_rng(20261006)
    a, true_center, k = 30., .82, 750.
    centers = np.arange(.45, 1.451, .1)
    samples = [rng.normal((a*true_center+k*c)/(a+k), np.sqrt(RT/(a+k)), 12000) for c in centers]
    z, weights, f = mbar(samples, centers, np.full(len(centers), k))
    edges = np.arange(.35, 1.551, .015)
    hist, _ = np.histogram(z, bins=edges, weights=weights)
    mid = .5*(edges[1:]+edges[:-1])
    use = (mid>.6)&(mid<1.3)&(hist>0)
    numerical = -RT*np.log(hist[use]); exact = .5*a*(mid[use]-true_center)**2
    numerical -= np.mean(numerical-exact)
    error = float(np.sqrt(np.mean((numerical-exact)**2)))
    assert error < .25, error
    return dict(known_free_energy='30/2*(z-.82)^2 kJ/mol, independent exact Gaussian biased samples',
                profile_RMS_error_kJ_mol=error, windows=len(centers), independent_samples_each=12000,
                scope='Estimator numerical check only, not a graphene or cleaning simulation')


def audit(root):
    windows, series = [], []
    for path in sorted(root.glob('z*/complete.json')):
        r = json.loads(path.read_text())
        folder = path.parent
        x = np.loadtxt(folder/'pullx.xvg', comments=['@', '#'])
        f = np.loadtxt(folder/'pullf.xvg', comments=['@', '#'])
        discard = r.get('pilot_discard_ps', r.get('discard_ps', 0))
        mask = x[:, 0]>=discard
        tail = x[mask, 1]; ftail = f[mask, 1]
        blocks = np.array_split(np.arange(len(tail)), 4)
        g = statistical_inefficiency(tail)
        block_means = [float(tail[b].mean()) for b in blocks]
        windows.append(dict(directory=str(folder.relative_to(HERE)), center_nm=r['center_nm'],
                            duration_ps=r['run_ps'], discarded_ps=discard,
                            z_mean_nm=float(tail.mean()), z_sd_nm=float(tail.std()),
                            block_z_means_nm=block_means,
                            block_force_means_kJ_mol_nm=[float(ftail[b].mean()) for b in blocks],
                            coordinate_statistical_inefficiency=g,
                            nominal_effective_coordinate_samples=len(tail)/g,
                            equilibration_not_implied_by_inefficiency=True))
        series.append(tail)
    overlaps = []
    if len(series)>1:
        edges = np.arange(min(min(x) for x in series)-.005, max(max(x) for x in series)+.011, .005)
        probabilities = [np.histogram(x, edges)[0]/len(x) for x in series]
        for i in range(len(series)-1):
            overlaps.append(dict(left=windows[i]['center_nm'], right=windows[i+1]['center_nm'],
                                 empirical_overlap_mass=float(np.minimum(probabilities[i], probabilities[i+1]).sum()),
                                 Bhattacharyya=float(np.sqrt(probabilities[i]*probabilities[i+1]).sum())))
    result = dict(windows=windows, adjacent_distance_overlaps=overlaps,
                  global_free_energy_profile_accepted=False, cleaning_ranking_accepted=False,
                  reasons=['Distance overlap does not establish stationarity of ester/backbone contacts or hidden torsions.',
                           'Both directions, independent starting conformations, complete coordinate range and model qualification required.'])
    (root/'sampling-audit.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--self-check', action='store_true');args=p.parse_args()
    out=HERE/'results/sampling-audits';out.mkdir(exist_ok=True)
    if args.self_check:
        check=positive_control();(out/'estimator-positive-control.json').write_text(json.dumps(check,indent=2)+'\n')
        print(check)
    for family in ['umbrella', 'umbrella-ladders']:
        root=HERE/'results'/family
        for folder in sorted(root.rglob('windows.json')):
            result=audit(folder.parent)
            print(folder.parent.relative_to(HERE), 'windows',len(result['windows']),
                  'overlaps',[round(r['empirical_overlap_mass'],3) for r in result['adjacent_distance_overlaps']],
                  'acceptedPMF',result['global_free_energy_profile_accepted'])


if __name__=='__main__':
    main()
