#!/usr/bin/env python3
"""Independent adversarial checks; does not import either candidate checker.

Exact arithmetic checks audit finite instances, not the universal proofs.
Decimal and floating-point checks are explicitly labelled numerical diagnostics.
"""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
import math
import platform
import sys
import time

import numpy as np
import scipy
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent


def prefix(v):
    return [v[j] - sum(v[:j]) / j for j in range(1, len(v))]


def close_rows(g):
    for i, row in enumerate(g):
        row[i] = -sum(row[j] for j in range(len(g)) if j != i)
    return g


def apply(g, b):
    return [sum(x * y for x, y in zip(row, b)) for row in g]


def transform_by_basis_action(g):
    """Reconstruct each prefix coordinate basis and apply the generator.

    Gauge b[0]=0; for column k, b[k]=1 and b[j]=1/(k+1) for j>k.
    This tests the dynamical coordinate definition directly, without importing
    the supplied Helmert formula or matrix-conjugation implementation.
    """
    n = len(g)
    cols = []
    for k in range(1, n):
        b = [F(0)] * k + [F(1)] + [F(1, k + 1)] * (n - k - 1)
        assert prefix(b) == [F(i == k - 1) for i in range(n - 1)]
        cols.append(prefix(apply(g, b)))
    return [list(row) for row in zip(*cols)]


def hierarchy(n, eps):
    m = n - 3
    z = 1 - 2 * eps - 3 * eps / (m - F(1, 2))
    g = [[F(0) for _ in range(n)] for _ in range(n)]
    g[1][0] = 1
    for j in range(2, n - 1):
        g[j][0], g[j][1] = F(3, 2), F(1)
        for i in range(2, j):
            g[j][i] = eps
    g[-1][0], g[-1][1] = 2 - eps, F(1)
    for i in range(2, n - 1):
        g[-1][i] = z
    return close_rows(g)


def restrict(g, subset):
    h = [[g[i][j] if i != j else F(0) for j in subset] for i in subset]
    return close_rows(h)


def triangular_slow(g):
    rates = [-g[i][i] for i in range(len(g))]
    assert all(rates[i] < rates[i + 1] for i in range(len(g) - 1))
    v = [F(0), F(1)]
    for i in range(2, len(g)):
        v.append(sum(g[i][j] * v[j] for j in range(i)) / (rates[i] - rates[1]))
    assert apply(g, v) == [-rates[1] * x for x in v]
    return v


def exact_hierarchy():
    count = 0
    records = []
    for eps in [F(1, 100), F(1, 15), F(499, 7000)]:
        for n in range(4, 10):
            g = hierarchy(n, eps)
            low = None
            for size in range(2, n):
                for subset in combinations(range(n), size):
                    d = min(prefix(triangular_slow(restrict(g, subset))))
                    assert d > 0
                    low = d if low is None else min(low, d)
                    count += 1
            full = prefix(triangular_slow(g))
            assert all(d > 0 for d in full[:-1]) and full[-1] < 0
            records.append({'N': n, 'epsilon': str(eps), 'minimum_proper_margin': str(low),
                            'full_final_margin': str(full[-1])})
    boundary = min(prefix(triangular_slow(restrict(hierarchy(4, F(1, 14)), [0, 2, 3]))))
    assert boundary == 0
    return {'proper_subsystems_checked': count, 'cases': records,
            'excluded_epsilon_1_over_14_boundary_margin': str(boundary)}


def bosonic_exact(energies, q, kernel):
    n = len(energies)
    g = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                gap = abs(energies[i] - energies[j])
                g[i][j] = kernel(gap) * q ** gap / (1 - q ** gap)
                if i > j:
                    g[i][j] += kernel(gap)
    return close_rows(g)


def exact_concavity():
    records = []
    for n in [3, 4, 6, 8]:
        for energies in [list(range(n)), [i * i for i in range(n)]]:
            span = energies[-1]
            kernels = {
                'strictly_decreasing_affine': lambda d: F(2 * span - d),
                'concave_hump': lambda d: F(2 * span + 2 * span * d - 2 * d * d),
                'zero_at_largest_gap': lambda d: F(span - d),
                'concave_piecewise_linear': lambda d: min(F(1 + d), 1 + F(span, 2)),
                'concave_quadratic_positive': lambda d: F(1 + 2 * span * d - d * d),
            }
            for name, kernel in kernels.items():
                for q in [F(1, 2), F(1, 16)]:
                    g = bosonic_exact(energies, q, kernel)
                    h = transform_by_basis_action(g)
                    low = min(h[i][j] for i in range(n - 1) for j in range(n - 1) if i != j)
                    assert low >= 0
                    if name != 'zero_at_largest_gap':
                        assert low > 0
                    records.append({'N': n, 'energies': energies, 'kernel': name,
                                    'q': str(q), 'minimum_offdiagonal': str(low)})
    # A nonconcave kernel that passes every zero-bath triplet on [0,1].
    # J(a)+J(b)-J(a+b)=1-2ab >= 1/2, but replication defeats the full cone.
    m = 13
    energies = [F(0), F(1, 4)] + [F(1, 2) + F(k - 6, 1000) for k in range(m)] + [F(3, 4)]
    kernel = lambda d: 1 + d * d
    n = len(energies)
    a = [[kernel(energies[i] - energies[j]) if i > j else F(0)
          for j in range(n)] for i in range(n)]
    a = close_rows(a)
    ha = transform_by_basis_action(a)[-1][0]
    assert ha == -F(1, 240)
    # At beta=16 log 2 only the top-to-second and top-to-first gaps enter
    # this thermal matrix element; the displayed value is therefore rational.
    thermal = F(n, 2 * (n - 1)) * (kernel(F(1, 2)) / (2 ** 8 - 1)
                                              - kernel(F(3, 4)) / (2 ** 12 - 1))
    total = ha + thermal
    assert total < 0
    e = np.array([float(x) for x in energies])
    diffs = abs(e[:, None] - e[None, :]); np.fill_diagonal(diffs, 1.)
    weights = 1 + diffs * diffs
    gf = weights / np.expm1(16 * math.log(2) * diffs) + weights * (e[:, None] > e[None, :])
    np.fill_diagonal(gf, 0); gf -= np.diag(gf.sum(axis=1))
    b = np.r_[0., 1., np.full(n - 2, .5)]
    numerical = prefix(gf @ b)[-1]
    assert abs(numerical - float(total)) < 1e-11
    # Smaller distinct-spectrum construction supplied by the physical branch,
    # independently checked here. No minimality claim is made.
    e6 = [F(0), F(1, 2), F(501, 1000), F(502, 1000), F(503, 1000), F(1)]
    a6 = close_rows([[kernel(e6[i] - e6[j]) if i > j else F(0) for j in range(6)] for i in range(6)])
    ha6 = transform_by_basis_action(a6)[-1][0]
    assert ha6 == -F(247, 5000)
    hs6 = F(3, 5) * (F(5, 4) / 255 - F(2) / 65535)
    full6 = ha6 + hs6
    assert full6 < 0
    return {'rational_concave_models': len(records), 'cases': records,
            'six_state_branch_witness': {'origin': 'physical_mpemba_frontier; independently checked here',
                'energies': list(map(str, e6)), 'zero_bath_H_62': str(ha6), 'thermal_H_62': str(hs6),
                'finite_bath_H_62': str(full6), 'bath_inverse_temperature': '16 log 2', 'minimality_claimed': False},
            'nonconcave_cluster_witness': {'kernel': '1+d^2', 'N': n, 'energies': list(map(str, energies)),
                'zero_bath_H_N2': str(ha), 'thermal_H_N2': str(thermal),
                'finite_bath_H_N2': str(total), 'numerical_H_N2': numerical,
                'bath_inverse_temperature': '16 log 2',
                'scope': 'Failure of cone invariance only; no actual thermal crossing claimed.'}}


def finite_time_divergences():
    tests, worst = 0, {'cooling_KL': 0., 'cooling_reverse_KL': 0.,
                       'cooling_chi2': 0., 'heating_chi2': 0., 'stochastic_order': 0.}
    records = []
    for energies in [np.array([0., .3, 1., 2.]), np.array([0., .2, .7, 1.1, 1.8, 3.])]:
        for name, kernel in [('linear', lambda x: x), ('decreasing', lambda x: 7 - x),
                             ('hump', lambda x: 1 + 2 * x - .5 * x * x)]:
            n = len(energies); gap = abs(energies[:, None] - energies[None, :]); np.fill_diagonal(gap, 1)
            j = kernel(gap)
            g = j / np.expm1(gap) + j * (energies[:, None] > energies[None, :])
            np.fill_diagonal(g, 0); g -= np.diag(g.sum(axis=1))
            pi = np.exp(-energies); pi /= pi.sum()
            for t in [0., .001, .01, .1, .5, 2.]:
                p_t = expm(t * g)
                for side, betas in [('cooling', [0., .1, .3, .6, .9, 1.]),
                                    ('heating', [1., 1.1, 1.5, 2., 4.])]:
                    values = []
                    for beta in betas:
                        p = np.exp(-beta * energies); p /= p.sum(); evolved = p @ p_t
                        r = evolved / pi
                        kl = np.sum(evolved * np.log(r)); rev = -np.sum(pi * np.log(r)); chi = np.sum(pi * (r - 1) ** 2)
                        values.append((evolved, kl, rev, chi))
                    for left, right in zip(values, values[1:]):
                        cdfmargin = np.cumsum(right[0] - left[0])[:-1].min()
                        worst['stochastic_order'] = min(worst['stochastic_order'], float(cdfmargin))
                        assert cdfmargin > -2e-13
                        if side == 'cooling':
                            for k, key in [(1, 'cooling_KL'), (2, 'cooling_reverse_KL'), (3, 'cooling_chi2')]:
                                violation = float(right[k] - left[k]); worst[key] = max(worst[key], violation)
                                assert violation < 2e-13
                                tests += 1
                        else:
                            violation = float(left[3] - right[3]); worst['heating_chi2'] = max(worst['heating_chi2'], violation)
                            assert violation < 2e-13
                            tests += 1
                    records.append({'N': n, 'kernel': name, 'time': t, 'side': side})
    # Jensen closure is exact for these concrete non-monotone cone elements.
    # The threshold/sign map is not closed, so the proof does not cover TV.
    b = [F(1, 10), F(7, 5), F(4, 5), F(4)]
    assert min(prefix(b)) >= 0
    sign_b = [F(-1) if x < 1 else F(1) for x in b]
    assert min(prefix(sign_b)) < 0
    return {'numerical_divergence_comparisons': tests, 'worst_roundoff_or_violation': worst,
            'cases': records, 'TV_proof_obstruction': {'b': list(map(str, b)),
                'prefix_margins': list(map(str, prefix(b))),
                'sign_at_one_prefix_margins': list(map(str, prefix(sign_b)))}}


def noncommuting_protocol():
    energies = np.array([0., 1., 2., 4.]); q = F(1, 2)
    exact_g1 = bosonic_exact([0, 1, 2, 4], q, lambda d: F(d))
    exact_g2 = bosonic_exact([0, 1, 2, 4], q, lambda d: F(1))
    g1, g2 = np.array(exact_g1, float), np.array(exact_g2, float)
    p1, p2 = expm(.07 * g1), expm(.13 * g2)
    row_p, k = p1 @ p2, p2 @ p1
    pi = 2. ** -energies; pi /= pi.sum()
    adj = k.T * pi[None, :] / pi[:, None]
    orientation_error = float(abs(adj - row_p).max())
    commutator = float(abs(g1 @ g2 - g2 @ g1).max())
    nonreversibility = float(abs(k - adj).max())
    assert orientation_error < 1e-14
    assert commutator > 1e-3 and nonreversibility > 1e-5
    assert min(x for row in transform_by_basis_action(k) for x in row) > 0
    assert min(x for row in transform_by_basis_action(adj) for x in row) > 0
    def power(r, alpha):
        if alpha == 0:
            return -np.dot(pi, np.log(r))
        if alpha == 1:
            return np.dot(pi * r, np.log(r))
        return (np.dot(pi, r ** alpha) - 1) / (alpha * (alpha - 1))
    largest_violation, largest_density_error, comparisons = 0., 0., 0
    for side, betas, alphas in [('cooling', [0., .1, .2, .4, .6, math.log(2)], [-2, 0, .5, 1, 1.5, 2]),
                               ('heating', [math.log(2), .8, 1., 2., 3.], [2, 3, 4])]:
        vals = []
        for beta in betas:
            p = np.exp(-beta * energies); p /= p.sum()
            r = k @ (p / pi)
            err = float(abs(pi * r - p @ row_p).max())
            largest_density_error = max(largest_density_error, err)
            assert err < 1e-14
            vals.append([power(r, a) for a in alphas])
        for a, b in zip(vals, vals[1:]):
            diff = np.array(b) - np.array(a)
            if side == 'heating':
                diff = -diff
            largest_violation = max(largest_violation, float(diff.max()))
            assert diff.max() < 2e-12
            comparisons += len(alphas)
    return {'energies': energies.tolist(), 'bath_inverse_temperature': 'log 2',
            'kernels': ['d', '1'], 'durations': [.07, .13],
            'generator_commutator_maximum': commutator, 'K_minus_adjoint_maximum': nonreversibility,
            'adjoint_product_orientation_error': orientation_error,
            'row_density_orientation_error': largest_density_error,
            'power_divergence_comparisons': comparisons,
            'largest_violation_or_roundoff': largest_violation,
            'scope': 'Floating-point check of a genuinely noncommuting two-stage protocol; proof remains analytic.'}


def decimal_solve(a, rhs):
    a = [list(row) + [x] for row, x in zip(a, rhs)]
    n = len(a)
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(a[i][k])); a[k], a[pivot] = a[pivot], a[k]
        assert a[k][k] != 0
        s = a[k][k]; a[k] = [x / s for x in a[k]]
        for i in range(n):
            if i != k:
                s = a[i][k]; a[i] = [x - s * y for x, y in zip(a[i], a[k])]
    return [row[-1] for row in a]


def decimal_slow(g):
    n = len(g); q = [[g[i][j] - g[0][j] for j in range(1, n)] for i in range(1, n)]
    # Newton's method for Qx=lambda x with x[0]=1, starting near the downward slow mode.
    lam = D(-1); x = [D(1)] + [D('0.6')] * (n - 2)
    for _ in range(16):
        residual = [sum(q[i][j] * x[j] for j in range(n - 1)) - lam * x[i] for i in range(n - 1)]
        jac = [[-x[i]] + [q[i][j] - (lam if i == j else 0) for j in range(1, n - 1)] for i in range(n - 1)]
        step = decimal_solve(jac, [-r for r in residual])
        lam += step[0]
        x = [D(1)] + [x[j] + step[j] for j in range(1, n - 1)]
        if max(abs(r) for r in residual) < D('1e-65'):
            break
    error = max(abs(sum(q[i][j] * x[j] for j in range(n - 1)) - lam * x[i]) for i in range(n - 1))
    assert error < D('1e-60')
    return lam, [D(0)] + x, error


def decimal_distribution(energies, beta):
    p = [(-beta * e).exp() for e in energies]; z = sum(p)
    return [x / z for x in p]


def decimal_evolve(p, g, t):
    n = len(g); rate = max(-g[i][i] for i in range(n)) + 1
    step = [[D(i == j) + g[i][j] / rate for j in range(n)] for i in range(n)]
    weight = (-rate * t).exp(); total = [weight * x for x in p]; state = p[:]
    for k in range(1, 1000):
        state = [sum(state[i] * step[i][j] for i in range(n)) for j in range(n)]
        weight *= rate * t / k
        total = [x + weight * y for x, y in zip(total, state)]
        if k > rate * t and weight < D('1e-70'):
            assert abs(sum(total) - 1) < D('1e-65')
            return total, k
    raise AssertionError('uniformization did not converge')


def positive_bath_hierarchy():
    with localcontext() as context:
        context.prec = 80
        energies = list(map(D, [0, 1, 2, 1555])); bath = D(10)
        base = hierarchy(4, F(1, 100))
        g = [[D(0)] * 4 for _ in range(4)]
        for i in range(4):
            for j in range(i):
                rate = base[i][j]; g[i][j] = D(rate.numerator) / rate.denominator
                g[j][i] = g[i][j] * (-bath * (energies[i] - energies[j])).exp()
        close_rows(g)
        assert min(g[i][j] for i in range(4) for j in range(4) if i != j) > 0
        pi = decimal_distribution(energies, bath)
        balance_error = max(abs(pi[i] * g[i][j] / (pi[j] * g[j][i]) - 1)
                            for i in range(4) for j in range(4) if i != j)
        assert balance_error < D('1e-65')
        lam, b, residual = decimal_slow(g)
        pi_mean = sum(x * y for x, y in zip(pi, b)); f = [x - pi_mean for x in b]
        proper = []
        for size in [2, 3]:
            for subset in combinations(range(4), size):
                sub = restrict(g, subset); _, sb, err = decimal_slow(sub)
                margin = min(prefix(sb)); assert margin > 0
                proper.append({'subset_one_based': [i + 1 for i in subset], 'minimum_prefix_margin': str(margin), 'residual': str(err)})
        bh, bc = D('0.00001'), D('0.00002')
        ph, pc = decimal_distribution(energies, bh), decimal_distribution(energies, bc)
        ah, ac = sum(x * y for x, y in zip(ph, f)), sum(x * y for x, y in zip(pc, f))
        assert 0 < ah < ac
        tv = lambda p: sum(abs(x - y) for x, y in zip(p, pi)) / 2
        initial = tv(ph) - tv(pc); assert initial > 0
        th, terms_h = decimal_evolve(ph, g, D(10)); tc, terms_c = decimal_evolve(pc, g, D(10))
        late = tv(th) - tv(tc); assert late < 0
        return {'precision_decimal_digits': context.prec, 'energies': list(map(str, energies)), 'beta_b': str(bath),
                'beta_h': str(bh), 'beta_c': str(bc), 'slow_eigenvalue': str(lam),
                'backward_eigenfunction_quotient_residual': str(residual),
                'relative_detailed_balance_residual': str(balance_error),
                'minimum_positive_rate': str(min(g[i][j] for i in range(4) for j in range(4) if i != j)),
                'slow_amplitude_hot': str(ah), 'slow_amplitude_cold': str(ac),
                'initial_TV_hot_minus_cold': str(initial), 'time_10_TV_hot_minus_cold': str(late),
                'uniformization_terms': [terms_h, terms_c], 'proper_subsystems': proper,
                'scope': '80-digit numerical witness, not an interval-arithmetic proof; analytic continuation proves existence.'}


def main():
    started = time.perf_counter()
    result = {'status': 'PASS', 'python': sys.version, 'platform': platform.platform(),
              'numpy': np.__version__, 'scipy': scipy.__version__,
              'hierarchy_exact': exact_hierarchy(), 'concavity_exact': exact_concavity(),
              'finite_time_divergence_diagnostics': finite_time_divergences(),
              'noncommuting_protocol': noncommuting_protocol(),
              'positive_bath_hierarchy': positive_bath_hierarchy()}
    result['elapsed_seconds'] = time.perf_counter() - started
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (ROOT / 'independent_check_results.json').write_text(json.dumps(result, indent=2) + '\n')
    summary = {key: result[key] for key in ['status', 'elapsed_seconds', 'script_sha256']}
    summary['proper_subsystems_checked'] = result['hierarchy_exact']['proper_subsystems_checked']
    summary['rational_concave_models'] = result['concavity_exact']['rational_concave_models']
    summary['nonconcave_cluster_witness'] = result['concavity_exact']['nonconcave_cluster_witness']
    summary['six_state_branch_witness'] = result['concavity_exact']['six_state_branch_witness']
    summary['positive_bath_hierarchy'] = {k: v for k, v in result['positive_bath_hierarchy'].items() if k != 'proper_subsystems'}
    summary['finite_time_divergence_comparisons'] = result['finite_time_divergence_diagnostics']['numerical_divergence_comparisons']
    summary['worst_roundoff_or_violation'] = result['finite_time_divergence_diagnostics']['worst_roundoff_or_violation']
    summary['noncommuting_protocol'] = result['noncommuting_protocol']
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
