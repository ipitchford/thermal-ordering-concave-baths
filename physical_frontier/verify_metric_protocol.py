#!/usr/bin/env python3
"""Exact rational audit of the finite-operator divergence/protocol theorem.

This deliberately uses stochastic uniformization *steps*, not floating-point
matrix exponentials. Each step and its weighted adjoint meet the abstract
theorem's hypotheses. The written proof extends the result to all exponential
semigroups and finite schedules. These checks audit orientation, strict signs,
power normalisation and the two distinct adjoints in noncommuting products.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import time

from verify_concave_kernel import (
    eye, gibbs, matrices, mm, quotient_transform, row_step,
)


def mv(a, x):
    return [sum((a[i][j]*x[j] for j in range(len(x))), F(0))
            for i in range(len(a))]


def pi_adjoint(a, pi):
    return [[a[j][i]*pi[j]/pi[i] for j in range(len(pi))]
            for i in range(len(pi))]


def margins(x):
    return [x[j]-sum(x[:j])/j for j in range(1, len(x))]


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), F(0))


def power_f(x, alpha):
    return (x**alpha-alpha*x+alpha-1)/(alpha*(alpha-1))


def power_fp(x, alpha):
    # Includes alpha=0 (reverse KL derivative); alpha=1 would require log.
    return (x**(alpha-1)-1)/(alpha-1)


def petz_fp(x, lam):
    a = 1-lam
    den = a*x+lam
    return (x-1)*(2*den-a*(x-1))/(den*den)


def step(E, qb, J):
    G, _, _ = matrices(E, qb, J)
    H = quotient_transform(G)
    R = 3*max([-G[i][i] for i in range(len(E))]
              + [-H[i][i] for i in range(len(H))])+1
    P = [[F(i == j)+G[i][j]/R for j in range(len(E))]
         for i in range(len(E))]
    # Strict row diagonal dominance certifies invertibility of each step.
    assert all(P[i][i] > F(1, 2) for i in range(len(E)))
    assert all(x >= 0 for row in P for x in row)
    return P


def audit_derivative(E, pi, K, Kstar, q, fp, side):
    """Compare direct differentiation, adjoint transfer and covariance exactly."""
    p = gibbs(E, q)
    mean_E = dot(p, E)
    pprime = [p[i]*(mean_E-E[i]) for i in range(len(E))]
    r = [p[i]/pi[i] for i in range(len(E))]
    rprime = [pprime[i]/pi[i] for i in range(len(E))]
    rt = mv(K, r)
    assert all(x > 0 for x in rt)
    assert [pi[i]*rt[i] for i in range(len(E))] == row_step(p, Kstar)
    phi = [fp(x) for x in rt]
    b = mv(Kstar, phi)
    sign = F(1) if side == "cooling" else F(-1)
    assert all(sign*x >= 0 for x in margins(rt))
    assert all(sign*x >= 0 for x in margins(phi))
    assert all(sign*x >= 0 for x in margins(b))
    direct = dot([pi[i]*phi[i] for i in range(len(E))], mv(K, rprime))
    moved = dot(pprime, b)
    covariance = dot([p[i]*E[i] for i in range(len(E))], b)-mean_E*dot(p, b)
    assert direct == moved == -covariance
    assert sign*direct < 0
    return direct, 7


def main():
    start = time.perf_counter()
    assertions = 0
    derivatives = 0
    nonreversible_products = 0
    cases = []
    example = None
    for n in (3, 4, 5):
        for E in (list(range(n)), [i*(i+1)//2 for i in range(n)]):
            D = E[-1]
            for qb in (F(1, 2), F(2, 3)):
                pi = gibbs(E, qb)
                kernels = [
                    lambda d: F(d),
                    lambda d: F(1),
                    lambda d: F(d)-F(d*d, 2*D),
                ]
                ps = [step(E, qb, J) for J in kernels]
                forward, density = eye(n), eye(n)
                for P in ps:
                    assert pi_adjoint(P, pi) == P
                    assert row_step(pi, P) == pi
                    forward, density = mm(forward, P), mm(P, density)
                    assertions += 4  # Includes invertibility/positivity in step.
                assert pi_adjoint(density, pi) == forward
                assertions += 1
                if density != forward:
                    nonreversible_products += 1
                products = [
                    ("one_step", ps[0], ps[0]),
                    ("three_stages", density, forward),
                    ("reverse_stages", forward, density),
                ]
                for label, K, Kstar in products:
                    assert Kstar == pi_adjoint(K, pi)
                    assert all(sum(row) == 1 for row in K)
                    assert row_step(pi, K) == pi
                    assert all(x >= 0 for row in quotient_transform(K) for x in row)
                    assert all(x >= 0 for row in quotient_transform(Kstar) for x in row)
                    assertions += 5
                    local_derivatives = []
                    for q in (F(1), (1+qb)/2):
                        for alpha in (-3, -1, 0, 2):
                            d, count = audit_derivative(
                                E, pi, K, Kstar, q,
                                lambda x, a=alpha: power_fp(x, a), "cooling")
                            local_derivatives.append(float(d))
                            derivatives += 1
                            assertions += count
                        for lam in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
                            d, count = audit_derivative(
                                E, pi, K, Kstar, q,
                                lambda x, lam=lam: petz_fp(x, lam), "cooling")
                            local_derivatives.append(float(d))
                            derivatives += 1
                            assertions += count
                    for q in (qb/2, qb/3):
                        for alpha in (2, 3, 4):
                            d, count = audit_derivative(
                                E, pi, K, Kstar, q,
                                lambda x, a=alpha: power_fp(x, a), "heating")
                            local_derivatives.append(float(d))
                            derivatives += 1
                            assertions += count
                    # Pairwise distances: rational f_alpha values allow exact
                    # checks in addition to the analytical derivative identity.
                    for side, qs, alphas in [
                        ("cooling", ((1+qb)/2, (1+3*qb)/4), (-3, -1, 2)),
                        ("heating", (qb/2, qb/3), (2, 3, 4)),
                    ]:
                        pair = []
                        for q in qs:
                            p = gibbs(E, q)
                            rt = mv(K, [p[i]/pi[i] for i in range(n)])
                            pair.append([dot(pi, [power_f(x, a) for x in rt])
                                         for a in alphas])
                        if side == "cooling":
                            assert all(a > b for a, b in zip(*pair))
                        else:
                            assert all(a < b for a, b in zip(*pair))
                        assertions += 1
                    # At the bath preparation, r=1 and every tested centred
                    # derivative is exactly zero.
                    assert mv(K, [F(1)]*n) == [F(1)]*n
                    assertions += 1
                    if example is None and E == [0, 1, 2] and qb == F(1, 2) and label == "three_stages":
                        d, _ = audit_derivative(E, pi, K, Kstar, F(3, 4),
                                               lambda x: x-1, "cooling")
                        example = {
                            "E": E, "bath_q": str(qb), "preparation_q": "3/4",
                            "metric": "half Pearson chi-square",
                            "preparation_derivative": str(d),
                            "K": [[str(x) for x in row] for row in K],
                            "pi_adjoint_K": [[str(x) for x in row] for row in Kstar],
                            "K_is_not_pi_self_adjoint": K != Kstar,
                        }
                    cases.append({
                        "N": n, "E": E, "bath_q": str(qb), "protocol": label,
                        "derivative_evaluations": len(local_derivatives),
                        "minimum_absolute_derivative_as_float":
                            min(abs(x) for x in local_derivatives),
                    })
    assert nonreversible_products > 0
    assertions += 1
    result = {
        "status": "PASS",
        "model_cases": len(cases),
        "exact_assertions": assertions,
        "exact_derivative_evaluations": derivatives,
        "nonreversible_three_stage_products": nonreversible_products,
        "normalisations": {
            "power_alpha_2": "one half Pearson chi-square",
            "power_alpha_half": "four times Hellinger squared when H²=half squared root distance",
            "power_alpha_0_limit": "reverse KL, not Renyi order zero",
        },
        "scope": "Exact finite-operator audit; these uniformization steps are not matrix exponentials. The general semigroup and protocol claims rely on the written proof.",
        "example": example,
        "cases": cases,
        "elapsed_seconds": time.perf_counter()-start,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(__file__).with_name("metric_protocol_checks.json").write_text(
        json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("cases", "example")}, indent=2))


if __name__ == "__main__":
    main()
