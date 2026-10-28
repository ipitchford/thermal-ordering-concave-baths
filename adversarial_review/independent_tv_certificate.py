#!/usr/bin/env python3
"""Independently audit exact heating-inversion certificates.

Unlike the supplied certificate, this script bounds exp(-x) by an alternating
Taylor pair, rather than taking reciprocals of bounds for exp(x). It performs
all sign decisions with Fraction and imports no candidate code.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def prefixes(b):
    return [b[i] - sum(b[:i]) / i for i in range(1, len(b))]


def cone_by_action(g):
    n = len(g); columns = []
    for k in range(1, n):
        basis = [F(0)] * k + [F(1)] + [F(1, k + 1)] * (n - k - 1)
        assert prefixes(basis) == [F(i == k - 1) for i in range(n - 1)]
        derivative = [sum(a * b for a, b in zip(row, basis)) for row in g]
        columns.append(prefixes(derivative))
    return list(map(list, zip(*columns)))


def row_action(p, transition):
    return [sum(p[i] * transition[i][j] for i in range(len(p))) for j in range(len(p))]


def independent_interval(g, p, t, terms=120):
    assert terms % 2 == 0 and all(x >= 0 for x in p) and sum(p) == 1
    rate = F(2); x = rate * t; n = len(g)
    transition = [[F(i == j) + g[i][j] / rate for j in range(n)] for i in range(n)]
    assert min(min(row) for row in transition) >= 0
    assert all(sum(row) == 1 for row in transition)
    assert x < terms + 1
    scalar_term = F(1); exp_negative_even = F(1); state = p[:]; series = p[:]
    for k in range(1, terms + 1):
        scalar_term *= x / k
        exp_negative_even += (-1 if k % 2 else 1) * scalar_term
        state = row_action(state, transition)
        assert all(z >= 0 for z in state) and sum(state) == 1
        series = [a + scalar_term * b for a, b in zip(series, state)]
    first_omitted = scalar_term * x / (terms + 1)
    # After this truncation, the alternating tail terms strictly decrease.
    exp_low, exp_high = exp_negative_even - first_omitted, exp_negative_even
    assert 0 < exp_low < exp_high
    # Positive uniformization tail: all component probabilities lie in [0,1].
    positive_tail = first_omitted / (1 - x / (terms + 2))
    low = [exp_low * v for v in series]
    high = [exp_high * (v + positive_tail) for v in series]
    assert all(0 < a < b for a, b in zip(low, high))
    assert sum(low) < 1 < sum(high)
    return low, high, (exp_low, exp_high)


def outward_interval(a, b, digits=20):
    scale = 10 ** digits
    lo = F((a.numerator * scale) // a.denominator, scale)
    hi = F(-((-b.numerator * scale) // b.denominator), scale)
    assert lo <= a <= b <= hi
    return {'lower_exact': str(lo), 'upper_exact': str(hi),
            'lower_decimal_approximation': float(lo), 'upper_decimal_approximation': float(hi)}


def audit(model):
    weights = list(map(F, model['weights'])); n = len(weights)
    assert all(weights[i] > weights[i + 1] > 0 for i in range(n - 1))
    pi = [w / sum(weights) for w in weights]
    assert pi == list(map(F, model['pi']))
    g = [list(map(F, row)) for row in model['generator']]
    assert all(sum(row) == 0 for row in g)
    assert all(g[i][j] >= 0 for i in range(n) for j in range(n) if i != j)
    assert all(pi[i] * g[i][j] == pi[j] * g[j][i] for i in range(n) for j in range(n))
    h = cone_by_action(g)
    assert h == [list(map(F, row)) for row in model['transformed_quotient']]
    assert all(h[i][j] >= 0 for i in range(n - 1) for j in range(n - 1) if i != j)
    threshold = F(model['threshold']); t = F(model['time'])
    beta_values = [3, 4]; initial = []
    for beta in beta_values:
        p = [w ** beta for w in weights]; z = sum(p); initial.append([v / z for v in p])
    divergence = lambda p: sum(max(v - threshold * pi[i], F(0)) for i, v in enumerate(p))
    d_initial = list(map(divergence, initial)); assert d_initial[1] > d_initial[0]
    intervals = [independent_interval(g, p, t) for p in initial]
    positive, negative = [], []
    for i in range(n):
        cut = threshold * pi[i]
        if all(bounds[0][i] > cut for bounds in intervals):
            positive.append(i)
        elif all(bounds[1][i] < cut for bounds in intervals):
            negative.append(i)
        else:
            raise AssertionError('Threshold sign is not certified by independent bounds')
    assert sorted(positive + negative) == list(range(n))
    lower = sum(intervals[0][0][i] - intervals[1][1][i] for i in positive)
    upper = sum(intervals[0][1][i] - intervals[1][0][i] for i in positive)
    assert lower > 0
    assert [i + 1 for i in positive] == model['positive_set_one_based']
    assert [i + 1 for i in negative] == model['negative_set_one_based']
    if n == 5:
        assert threshold == 1 and min(g[i][j] for i in range(n) for j in range(n) if i != j) > 0
    width = max(b - a for lo, hi, _ in intervals for a, b in zip(lo, hi))
    return {'name': model['name'], 'N': n, 'status': 'PASS', 'all_decisions_exact_rational': True,
            'uniformization_terms': 121, 'exponential_bound_method': 'alternating Taylor even/odd pair for exp(-x)',
            'initial_cold_minus_warm': outward_interval(d_initial[1] - d_initial[0], d_initial[1] - d_initial[0]),
            'final_warm_minus_cold': outward_interval(lower, upper),
            'positive_set_one_based': [i + 1 for i in positive], 'negative_set_one_based': [i + 1 for i in negative],
            'probability_interval_width_upper_bound': outward_interval(F(0), width, 45),
            'scope': 'Heating preparations beta=3,4 toward beta_b=1; abstract reversible cone-preserving generator, not a concave-J physical model.'}


def main():
    source = HERE / 'inspected_sources' / 'metric_sharpness__exact_counterexamples.json'
    models = json.loads(source.read_text())['models']
    results = [audit(model) for model in models]
    output = {'status': 'PASS', 'models': results,
              'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE / 'independent_tv_certificate_results.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
