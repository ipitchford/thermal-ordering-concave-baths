"""E2: exact cone failure and certified relative-entropy cooling crossing for
the complete-graph bosonic generator with super-Ohmic jump weight J(d)=d^3.

Python standard library only. Every decision is made on exact Fractions:
 * the generator, Helmert quotient H, Gibbs states and equilibrium are rational
   (integer energies, q=exp(-beta_b) rational, u=exp(-beta) rational);
 * exp(tG) is enclosed by rational stochastic uniformisation with a rational
   geometric tail bound (as in Appendix A of the manuscript);
 * log is enclosed by the atanh series with a rational geometric tail bound;
 * KL = sum_i pi_i phi(r_i), phi(r)=r log r - r + 1, enclosed exactly using
   convexity/monotonicity of phi on each side of r=1.
Decimal strings in the output are for reading only.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
import json, math

OUT = Path(__file__).resolve().parent

def outward(x, sig=12, up=False):
    """Decimal string with `sig` significant digits rounded outward; exactness checked with Fractions."""
    from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
    with localcontext() as c:
        c.prec = sig; c.rounding = ROUND_CEILING if up else ROUND_FLOOR
        d = +(Decimal(x.numerator) / Decimal(x.denominator))
    fd = F(str(d))
    assert (fd >= x) if up else (fd <= x)
    return str(d)

def dec(x, p=30):
    with localcontext() as c:
        c.prec = p
        return str(Decimal(x.numerator) / Decimal(x.denominator))

def eye(n): return [[F(int(i == j)) for j in range(n)] for i in range(n)]
def mm(a, b): return [[sum((x*y for x, y in zip(r, c)), F(0)) for c in zip(*b)] for r in a]
def rowmul(v, a): return [sum((v[i]*a[i][j] for i in range(len(v))), F(0)) for j in range(len(a[0]))]

def inv(a):
    n = len(a); z = [r.copy() + e for r, e in zip(a, eye(n))]
    for i in range(n):
        j = next(j for j in range(i, n) if z[j][i]); z[i], z[j] = z[j], z[i]
        x = z[i][i]; z[i] = [v/x for v in z[i]]
        for j in range(n):
            if i != j:
                x = z[j][i]; z[j] = [a - x*b for a, b in zip(z[j], z[i])]
    return [r[n:] for r in z]

def cone_matrix(g):
    """H with d(Gb)=H d(b), d_j(b)=b_j - mean(b_1..b_{j-1}) (manuscript Eq. helmert)."""
    n = len(g); t = eye(n-1)
    for i in range(n-1):
        for j in range(i): t[i][j] = -F(1, i+1)
    qm = [[g[i][j] - g[0][j] for j in range(1, n)] for i in range(1, n)]
    return mm(mm(t, qm), inv(t))

def boson_generator(E, q, J):
    """Eq. (boson): G_ij = J(d)/(e^{beta_b d}-1) + 1[E_i>E_j] J(d), q=e^{-beta_b}."""
    n = len(E); g = [[F(0)]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                d = abs(E[i]-E[j]); nb = q**d/(1-q**d)
                g[i][j] = J(d)*(nb + (1 if E[i] > E[j] else 0))
        g[i][i] = -sum(g[i][j] for j in range(n) if j != i)
    return g

def gibbs(E, u):
    w = [u**e for e in E]; z = sum(w); return [x/z for x in w]

def uniformise(g, p0, t, kmax, rate):
    n = len(g); P = eye(n)
    for i in range(n):
        for j in range(n): P[i][j] += g[i][j]/rate
    assert all(v >= 0 for r in P for v in r) and all(sum(r) == 1 for r in P)
    x = rate*t; assert x < kmax + 2
    term = F(1); B = F(1); v = list(p0); S = list(p0)
    for k in range(1, kmax+1):
        term *= x/k; B += term; v = rowmul(v, P); S = [a + term*b for a, b in zip(S, v)]
    T = term*x/(kmax+1)/(1 - x/F(kmax+2))
    return [a/(B+T) for a in S], [(a+T)/B for a in S], T/B

def log_enclosure(y, terms=None):
    """Rational [lo,hi] containing log(y), y>0, via log y = 2 atanh z."""
    z = (y-1)/(y+1); z2 = z*z; s = F(0); pw = z; k = 0
    while True:
        s += pw/(2*k+1); pw *= z2; k += 1
        tail = 2*abs(pw)/((2*k+1)*(1-z2))  # bounds |2 sum_{j>=k} z^{2j+1}/(2j+1)|
        if tail < F(1, 10**70) or (terms and k >= terms): break
    c = 2*s
    return c - tail, c + tail

def phi_enclosure(rlo, rhi):
    """Enclose phi(r)=r log r - r + 1 over r in [rlo,rhi] (all r>0)."""
    def ph(r):
        lo, hi = log_enclosure(r)
        return r*lo - r + 1, r*hi - r + 1
    a, b = ph(rlo), ph(rhi)
    upper = max(a[1], b[1])                      # convex: max at an endpoint
    if rlo <= 1 <= rhi: lower = F(0)
    elif rhi < 1: lower = b[0]                     # decreasing on (0,1)
    else: lower = a[0]                             # increasing on (1,inf)
    return lower, upper

def kl_enclosure(plo, phi_, pi):
    lo = hi = F(0)
    for a, b, w in zip(plo, phi_, pi):
        l, h = phi_enclosure(a/w, b/w); lo += w*l; hi += w*h
    return lo, hi

def main():
    E = [0, 1, 2, 3]
    q = F(2, 3)                     # beta_b = log(3/2)
    u_hot, u_cold = F(9, 10), F(4, 5)  # beta = log(10/9) < log(5/4) < beta_b: both cooling
    assert 1 > u_hot > u_cold > q
    J3 = lambda d: F(d)**3
    g = boson_generator(E, q, J3); n = len(E)
    pi = gibbs(E, q)
    assert all(sum(r) == 0 for r in g)
    assert all(pi[i]*g[i][j] == pi[j]*g[j][i] for i in range(n) for j in range(n))
    assert all(g[i][j] > 0 for i in range(n) for j in range(n) if i != j)
    H = cone_matrix(g)
    neg = [(i+2, j+2, H[i][j]) for i in range(n-1) for j in range(n-1) if i != j and H[i][j] < 0]
    assert neg, 'expected cone failure'
    # contrast: Ohmic J(d)=d on the same spectrum and bath preserves the cone (Theorem 4)
    H1 = cone_matrix(boson_generator(E, q, lambda d: F(d)))
    assert all(H1[i][j] >= 0 for i in range(n-1) for j in range(n-1) if i != j)

    ph, pc = gibbs(E, u_hot), gibbs(E, u_cold)
    rate = F(56); assert rate >= max(-g[i][i] for i in range(n))
    times = [F(0), F(1, 10), F(1, 5), F(3, 10), F(2, 5), F(1, 2)]
    rows = []
    for t in times:
        if t == 0:
            eh, ec = (ph, ph, F(0)), (pc, pc, F(0))
        else:
            eh = uniformise(g, ph, t, 160, rate); ec = uniformise(g, pc, t, 160, rate)
        Kh = kl_enclosure(eh[0], eh[1], pi); Kc = kl_enclosure(ec[0], ec[1], pi)
        dlo, dhi = Kh[0] - Kc[1], Kh[1] - Kc[0]
        sign = '+' if dlo > 0 else ('-' if dhi < 0 else '?')
        assert sign != '?'
        rows.append(dict(t=str(t), KL_hot=[dec(Kh[0]), dec(Kh[1])], KL_cold=[dec(Kc[0]), dec(Kc[1])],
                         hot_minus_cold=[dec(dlo), dec(dhi)], certified_sign=sign,
                         hot_minus_cold_outward=[outward(dlo), outward(dhi, up=True)],
                         ratio_cold_over_hot_lower=dec(Kc[0]/Kh[1], 12),
                         max_prob_width=dec(max(max(b-a for a, b in zip(eh[0], eh[1])),
                                                 max(b-a for a, b in zip(ec[0], ec[1]))), 6)))
        print('t=%-5s KL_hot in [%s, %s]  KL_cold in [%s, %s]  sign(hot-cold)=%s' % (
            t, dec(Kh[0], 12), dec(Kh[1], 12), dec(Kc[0], 12), dec(Kc[1], 12), sign), flush=True)
    signs = [r['certified_sign'] for r in rows]
    assert signs[0] == '+' and signs[4] == '-'
    # bracket of a crossing: last certified '+' before first certified '-'
    first_minus = signs.index('-')
    out = dict(
        model='complete-graph bosonic generator, Eq. (boson), J(d)=d^3 (super-Ohmic)',
        energies=E, q_bath=str(q), beta_b='log(3/2)', u_hot=str(u_hot), u_cold=str(u_cold),
        beta_hot='log(10/9)', beta_cold='log(5/4)',
        generator=[[str(v) for v in r] for r in g],
        equilibrium=[str(v) for v in pi],
        cone_quotient_H=[[str(v) for v in r] for r in H],
        negative_offdiagonal_H_entries=[dict(i=i, j=j, value=str(v), decimal=dec(v, 12)) for i, j, v in neg],
        ohmic_contrast_H=[[str(v) for v in r] for r in H1], ohmic_contrast_metzler=True,
        uniformisation_rate=str(rate), series_terms=161,
        rows=rows,
        crossing_bracket=[str(times[first_minus-1]), str(times[first_minus])],
        certificate='Exact Fractions throughout; exp series and log bounded by rational geometric tails.')
    (OUT/'e2_superohmic_certificate.json').write_text(json.dumps(out, indent=2))
    print('negative H entries:', [(i, j, str(v)) for i, j, v in neg])
    print('crossing certified in t in (%s, %s)' % tuple(out['crossing_bracket']))

if __name__ == '__main__': main()
