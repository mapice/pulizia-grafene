"""Verifiche matematiche del modello; nessun dato sperimentale simulato.

Solo libreria standard. Parametri adimensionali scelti per controllare
identità e limiti, non stime dei residui di AR-P 672.045.
"""

import math


def rates(beta, refresh):
    s = 1.0 + beta + refresh
    disc = math.sqrt(max(0.0, s * s - 4.0 * refresh))
    fast = 0.5 * (s + disc)
    slow = 2.0 * refresh / (s + disc) if refresh else 0.0
    return slow, fast


def exact_x(tau, beta, refresh):
    if beta == 0:
        return math.exp(-tau)
    lo, hi = rates(beta, refresh)
    return ((hi - 1) * math.exp(-lo * tau)
            + (1 - lo) * math.exp(-hi * tau)) / (hi - lo)


def rhs(state, beta, refresh, access):
    x, y, z, exported = state
    return (-x + beta * y + access * z,
            x - (beta + refresh) * y,
            -access * z,
            refresh * y)


def integrate(tau, beta, refresh, access=0.0, initial=(1., 0., 0., 0.)):
    n = max(1000, math.ceil(tau * max(1, beta + refresh, access) * 100))
    dt = tau / n
    state = initial
    max_mass_error = 0.0
    minimum = min(initial)
    mass = sum(initial)
    for _ in range(n):
        a = rhs(state, beta, refresh, access)
        b = rhs(tuple(u + dt * v / 2 for u, v in zip(state, a)),
                beta, refresh, access)
        c = rhs(tuple(u + dt * v / 2 for u, v in zip(state, b)),
                beta, refresh, access)
        d = rhs(tuple(u + dt * v for u, v in zip(state, c)),
                beta, refresh, access)
        state = tuple(u + dt * (v + 2*w + 2*s + t) / 6
                      for u, v, w, s, t in zip(state, a, b, c, d))
        max_mass_error = max(max_mass_error, abs(sum(state) - mass))
        minimum = min(minimum, *state)
    return state, max_mass_error, minimum


def sequential(tau_total, total_partition, n):
    b = n * total_partition
    exponent = -(1 + b) * tau_total / n
    removed = -math.expm1(exponent) / (1 + b)
    return math.exp(n * math.log1p(-removed))


def sequential_limit(tau_total, total_partition):
    if total_partition == 0:
        return math.exp(-tau_total)
    return math.exp(math.expm1(-total_partition * tau_total) / total_partition)


def poor_map(state, a_time, b_time):
    locked, mobile, removed = state
    total_rate = a_time + b_time
    if not total_rate:
        return state
    target = a_time * (locked + mobile) / total_rate
    new_mobile = target + (mobile-target) * math.exp(-total_rate)
    return locked + mobile - new_mobile, new_mobile, removed


def good_map(state, k_time, c_time=0.0):
    locked, mobile, removed = state
    total_rate = k_time + c_time
    if not total_rate:
        return state
    departing = mobile * (-math.expm1(-total_rate))
    return (locked + departing*c_time/total_rate,
            mobile-departing, removed + departing*k_time/total_rate)


def cycles(n, a_total, b_total, k_total, c_total=0.0):
    state = (1., 0., 0.)
    for _ in range(n):
        state = poor_map(state, a_total/n, b_total/n)
        state = good_map(state, k_total/n, c_total/n)
    return state


def run_checks():
    max_error = 0.0
    max_mass_error = 0.0
    for beta in (0., .1, 1., 10.):
        for refresh in (0., .1, 1., 10.):
            for tau in (.01, 1., 5.):
                state, mass_error, minimum = integrate(tau, beta, refresh)
                error = abs(state[0] - exact_x(tau, beta, refresh))
                max_error = max(max_error, error)
                max_mass_error = max(max_mass_error, mass_error)
                assert minimum > -1e-13
                assert error < 5e-10, (beta, refresh, tau, error)
                assert state[0] >= math.exp(-tau) - 1e-12
    for beta in (.1, 1., 10.):
        value = exact_x(2, beta, 0)
        expected = (beta + math.exp(-(1+beta)*2))/(1+beta)
        assert abs(value - expected) < 1e-14
        for tau in (.1, 1., 10.):
            xs = [exact_x(tau, beta, q) for q in (0, .1, 1, 10, 1000)]
            assert all(x1 >= x2 for x1, x2 in zip(xs, xs[1:]))
    state, mass_error, minimum = integrate(10, 2, 3, .1, (1., 0., .6, 0.))
    assert mass_error < 1e-12
    assert minimum >= -1e-13
    assert state[0] >= math.exp(-10)
    assert abs(state[2] - .6*math.exp(-1)) < 1e-12
    assert max_mass_error < 1e-12
    for B in (.01, .1, 1., 10.):
        for tau in (.1, 1., 10.):
            approximate = sequential(tau, B, 1_000_000)
            limit = sequential_limit(tau, B)
            assert abs(approximate - limit) < 2e-7
            assert abs(sequential(tau, B, 1)-exact_x(tau, B, 0)) < 1e-14
    print("Verifiche superate: formula esatta, conservazione, positivita, "
          "limite di distacco, accesso lento, bagni a tempo e volume fissati.")
    print(f"Massimo errore assoluto RK4/formula: {max_error:.3g}")
    print(f"Massimo errore di conservazione: {max_mass_error:.3g}")
    for delta in (1e-5, 1e-4, 1e-3):
        diffusion = 2.93e-11
        print(f"Geometria illustrativa delta={delta:g} m: "
              f"D/delta={diffusion/delta:.3g} m/s; "
              f"delta^2/D={delta*delta/diffusion:.6g} s")
    for a, b, k, c in ((1., 2., 3., 4.), (.1, 0., .8, 0.), (0., 2., .2, .1)):
        for start in ((1., 0., 0.), (0., 1., 0.), (.3, .5, .2)):
            pg = good_map(poor_map(start, a, b), k, c)
            gp = poor_map(good_map(start, k, c), a, b)
            expected = (k/(k+c))*(-math.expm1(-k-c)) \
                * (-math.expm1(-a-b)) * (a*start[0]-b*start[1])/(a+b)
            assert abs((pg[2]-gp[2])-expected) < 1e-14
            assert abs(sum(pg)-sum(start)) < 1e-14
            assert min(pg) >= -1e-14
    for w in (.1, 1., 3.):
        last = 1.
        for n in (1, 2, 5, 10, 100):
            observed = cycles(n, w, 0, w)[2]
            exact = 1-math.exp(-w)*(1+n*(-math.expm1(-w/n)))
            assert abs(observed-exact) < 5e-14
            assert observed < last
            last = observed
    print("Modello a tre stati: formula di ordine esatta, conservazione e "
          "controesempio cicli peggiori verificati.")
    for n in (1, 10):
        print(f"Cicli illustrativi n={n}: regime senza ritorno "
              f"R={cycles(n, 1., 0., 1.)[2]:.9f}; "
              f"regime con ritorno R={cycles(n, 1., 9., 100.)[2]:.9f}")
    example_thickness = 1e-3 * 1e-5 / 1180
    assert abs(example_thickness/1e-9-0.008474576271186441) < 1e-14
    print(f"Esempio limite film uscita: {example_thickness/1e-9:.6g} nm")


if __name__ == "__main__":
    run_checks()
