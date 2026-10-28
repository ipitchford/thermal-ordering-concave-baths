"""E3: certified finite-bath instance of the hierarchy theorem (Theorem 6) for N=4.

Construction of Sec. 6 with eps=1/100, m=1 (C=199/100, z=23/25), energies
E=(0,1,2,L), L=1555 (> 1+291*16/6 = 777, Eq. (L)), and a rational bath
parameter q=e^{-beta_b}; upward rates by Eq. (finitebath).

Certified (standard library only, exact Fractions; Decimal for display):
 (a) detailed balance, positivity of all rates;
 (b) every proper induced subsystem (6 pairs, 4 triples) has simple spectrum
     and a slow eigenfunction with strictly positive prefix margins, decided
     exactly (pairs: automatic; triples: 2x2 quotient, sign of lambda_+ - H11
     decided without square roots);
 (c) full system: simple spectrum (three disjoint rational sign-change
     brackets of the characteristic polynomial of H), slow-amplitude inversion
     0 < A(beta_h) < A(beta_c) by rational interval evaluation at the
     bracketed slow eigenvalue;
 (d) full system: TV(hot) > TV(cold) at t=0 (exact) and TV(hot) < TV(cold)
     at a rational t_* (uniformisation with rational tail bound and outward
     dyadic rounding), hence a finite-time cooling TV crossing.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from itertools import combinations
from pathlib import Path
import json, math, sys
sys.set_int_max_str_digits(0)

OUT = Path(__file__).resolve().parent
EPS = F(1, 100); M = 1; C = 2-EPS; Z = 1-2*EPS-3*EPS/(M-F(1, 2))
SCALE = 2**360

def dec(x, p=20):
    with localcontext() as c:
        c.prec = p; return str(Decimal(x.numerator)/Decimal(x.denominator))
def floor_d(x): return F(math.floor(x*SCALE), SCALE)
def ceil_d(x): return F(math.ceil(x*SCALE), SCALE)
def eye(n): return [[F(int(i == j)) for j in range(n)] for i in range(n)]
def mm(a, b): return [[sum((x*y for x, y in zip(r, c)), F(0)) for c in zip(*b)] for r in a]
def inv(a):
    n = len(a); z = [r.copy()+e for r, e in zip(a, eye(n))]
    for i in range(n):
        j = next(j for j in range(i, n) if z[j][i]); z[i], z[j] = z[j], z[i]
        x = z[i][i]; z[i] = [v/x for v in z[i]]
        for j in range(n):
            if i != j:
                x = z[j][i]; z[j] = [a-x*b for a, b in zip(z[j], z[i])]
    return [r[n:] for r in z]
def cone_matrix(g):
    n = len(g); t = eye(n-1)
    for i in range(n-1):
        for j in range(i): t[i][j] = -F(1, i+1)
    q = [[g[i][j]-g[0][j] for j in range(1, n)] for i in range(1, n)]
    return mm(mm(t, q), inv(t))

def downward(L):
    E = {1: 0, 2: 1, 3: 2, 4: L}
    a = {(1, 2): F(1), (1, 3): F(3, 2), (2, 3): F(1), (1, 4): C, (2, 4): F(1), (3, 4): Z}
    return E, a
def generator(E, a, q, states):
    n = len(states); G = [[F(0)]*n for _ in range(n)]
    for x, i in enumerate(states):
        for y, j in enumerate(states):
            if i < j:
                G[y][x] = a[(i, j)]
                G[x][y] = a[(i, j)]*q**(E[j]-E[i])
    for x in range(n): G[x][x] = -sum(G[x][y] for y in range(n) if y != x)
    return G

def triple_margin(H):
    """Exact sign of d3/d2 for the slow (largest) eigenvector of a 2x2 H."""
    a, b, c, d = H[0][0], H[0][1], H[1][0], H[1][1]
    D = (a-d)**2+4*b*c
    assert D > 0, 'repeated eigenvalue'
    bc = b*c
    s = 1 if bc > 0 else (0 if bc == 0 and d <= a else (1 if bc == 0 else (1 if d > a else -1)))
    # eigenvector (b, lambda_+ - a); if b==0 use (lambda_+ - d, c)
    if b != 0:
        sb = 1 if b > 0 else -1; return s*sb > 0 and s != 0, None
    s2 = 1  # lambda_+ - d >= 0 always when b==0? (lambda_+ = max(a,d))
    raise AssertionError('degenerate b=0 case not needed here')

def charpoly3(H):
    a = H
    tr = a[0][0]+a[1][1]+a[2][2]
    m2 = (a[0][0]*a[1][1]-a[0][1]*a[1][0])+(a[0][0]*a[2][2]-a[0][2]*a[2][0])+(a[1][1]*a[2][2]-a[1][2]*a[2][1])
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
           +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    return lambda x: x**3-tr*x**2+m2*x-det

def bracket_root(p, lo, hi, iters):
    flo = p(lo); assert flo != 0 and (flo > 0) != (p(hi) > 0)
    for _ in range(iters):
        mid = (lo+hi)/2; fm = p(mid)
        mid = F(math.floor(mid*SCALE), SCALE) if mid.denominator > SCALE else mid
        fm = p(mid)
        if fm == 0: return mid, mid
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    return lo, hi

def iv_mul(x, y):
    c = [x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1]]; return (min(c), max(c))
def iv_add(x, y): return (x[0]+y[0], x[1]+y[1])
def iv_sub(x, y): return (x[0]-y[1], x[1]-y[0])
def iv(c): return (c, c)

def slow_amplitudes(G, pi, preps, lam):
    """Interval amplitudes A(beta)=p_beta^T f - pi^T f for the eigenfunction f
    whose prefix coordinates d solve (H - lam I) d = 0, scaled with d_2 fixed by
    a cross product; lam is an interval."""
    H = cone_matrix(G); n = len(G)
    R = [[iv_sub(iv(H[i][j]), lam) if i == j else iv(H[i][j]) for j in range(3)] for i in range(3)]
    best = None
    for (i0, i1) in ((0, 1), (0, 2), (1, 2)):
        r0, r1 = R[i0], R[i1]
        dd = [iv_sub(iv_mul(r0[1], r1[2]), iv_mul(r0[2], r1[1])),
              iv_sub(iv_mul(r0[2], r1[0]), iv_mul(r0[0], r1[2])),
              iv_sub(iv_mul(r0[0], r1[1]), iv_mul(r0[1], r1[0]))]
        qual = min((abs(x[0]+x[1]) - (x[1]-x[0])) / (abs(x[0]+x[1]) + (x[1]-x[0]) + F(1, 10**4000)) for x in dd)
        if best is None or qual > best[0]: best = (qual, dd, (i0, i1))
    d = best[1]
    # lift: f_1 = 0, f_j = d_j + mean(f_1..f_{j-1})
    f = [iv(F(0))]
    for j in range(2, n+1):
        s = (F(0), F(0))
        for v in f: s = iv_add(s, v)
        f.append(iv_add(d[j-2], (s[0]/(j-1), s[1]/(j-1))))
    def amp(p):
        acc = (F(0), F(0))
        for pi_, pv, fv in zip(pi, p, f):
            acc = iv_add(acc, iv_mul(iv(pv-pi_), fv))
        return acc
    return d, [amp(p) for p in preps]

def uniformise_iv(G, p0, t, rate, kmax):
    n = len(G); P = eye(n)
    for i in range(n):
        for j in range(n): P[i][j] += G[i][j]/rate
    assert all(v >= 0 for r in P for v in r) and all(sum(r) == 1 for r in P)
    Pd = [[(floor_d(v), ceil_d(v)) for v in r] for r in P]
    x = rate*t; assert x < kmax+2
    lo = [floor_d(v) for v in p0]; hi = [ceil_d(v) for v in p0]
    Slo, Shi = list(lo), list(hi); term = F(1); B = F(1)
    for k in range(1, kmax+1):
        term *= x/k; B += term
        lo = [floor_d(sum(lo[i]*Pd[i][j][0] for i in range(n))) for j in range(n)]
        hi = [ceil_d(min(F(1), sum(hi[i]*Pd[i][j][1] for i in range(n)))) for j in range(n)]
        Slo = [floor_d(s+term*v) for s, v in zip(Slo, lo)]
        Shi = [ceil_d(s+term*v) for s, v in zip(Shi, hi)]
    T = term*x/(kmax+1)/(1-x/F(kmax+2))
    return [s/(B+T) for s in Slo], [(s+T)/B for s in Shi]

def tv_iv(lo, hi, pi):
    pos = []
    for i in range(len(pi)):
        if lo[i] > pi[i]: pos.append(i)
        elif hi[i] < pi[i]: pass
        else: raise AssertionError('component sign not certified')
    return sum(lo[i]-pi[i] for i in pos), sum(hi[i]-pi[i] for i in pos), pos

def run(L, q, uh, uc, tstar):
    E, a = downward(L); full = [1, 2, 3, 4]
    G = generator(E, a, q, full)
    w = [q**E[s] for s in full]; pi = [x/sum(w) for x in w]
    assert all(pi[i]*G[i][j] == pi[j]*G[j][i] for i in range(4) for j in range(4))
    assert all(G[i][j] > 0 for i in range(4) for j in range(4) if i != j)
    rec = dict(L=L, q=str(q), u_hot=str(uh), u_cold=str(uc))
    # exact closed form for the ground-slow-top triple (binding constraint on q_0)
    Ht = cone_matrix(generator(E, a, q, [1, 2, 4]))
    assert Ht[1][0] == (EPS-q)/2 + q**(L-1)*(1-C*q)/4 and Ht[0][1] == q**(L-1)*(1-C*q)
    rec['triple_124_H32_formula'] = '(eps-q)/2 + q^(L-1)(1-Cq)/4 verified exactly'
    # (b) proper subsystems
    subs = []
    for k in (2, 3):
        for S in combinations(full, k):
            g = generator(E, a, q, list(S))
            if k == 2:
                subs.append(dict(states=list(S), margins_positive=True, reason='two-state: eigenfunction (-pi2,pi1) increasing'))
                continue
            H = cone_matrix(g); ok, info = triple_margin(H)
            assert ok, ('margin fails', S)
            subs.append(dict(states=list(S), margins_positive=True, H_decimal=[[dec(v, 8) for v in r] for r in H]))
    rec['proper_subsystems'] = subs
    # (c) full slow amplitude inversion
    H = cone_matrix(G); p = charpoly3(H)
    # brackets from float estimates
    import cmath
    est = sorted(_roots_float(H))
    brs = []
    for r in est:
        lo, hi = F(r)-F(1, 10**6), F(r)+F(1, 10**6)
        brs.append(bracket_root(p, lo, hi, 220))
    for lo, hi in brs:
        assert lo == hi and p(lo) == 0 or (p(lo) > 0) != (p(hi) > 0), 'bracket sign'
    assert brs[0][1] < brs[1][0] and brs[1][1] < brs[2][0], 'brackets overlap'
    # three disjoint sign-change brackets of a cubic: exactly one simple root in each
    lam = brs[2]  # slow = largest (closest to zero)
    ph = [x/sum(uh**E[s] for s in full) for x in (uh**E[s] for s in full)]
    pc = [x/sum(uc**E[s] for s in full) for x in (uc**E[s] for s in full)]
    d, (Ah, Ac) = slow_amplitudes(G, pi, [ph, pc], lam)
    sign = 1 if d[0][0] > 0 else -1
    assert d[0][0] > 0 or d[0][1] < 0
    Ah = (Ah[0], Ah[1]) if sign > 0 else (-Ah[1], -Ah[0]); Ac = (Ac[0], Ac[1]) if sign > 0 else (-Ac[1], -Ac[0])
    assert Ah[0] > 0 and Ac[0] > Ah[1], 'slow amplitude inversion not certified'
    rec['full_eigenvalue_brackets'] = [[dec(lo, 25), dec(hi, 25)] for lo, hi in brs]
    rec['slow_amplitude_hot'] = [dec(Ah[0], 15), dec(Ah[1], 15)]
    rec['slow_amplitude_cold'] = [dec(Ac[0], 15), dec(Ac[1], 15)]
    rec['slow_amplitude_ratio_cold_over_hot_lower'] = dec(Ac[0]/Ah[1], 15)
    # (d) TV crossing
    tv0 = [sum(max(x-y, F(0)) for x, y in zip(pp, pi)) for pp in (ph, pc)]
    assert tv0[0] > tv0[1]
    rate = F(4); assert rate >= max(-G[i][i] for i in range(4))
    eh = uniformise_iv(G, ph, tstar, rate, 200); ec = uniformise_iv(G, pc, tstar, rate, 200)
    a_, b_ = tv_iv(*eh, pi), tv_iv(*ec, pi)
    lo, hi = a_[0]-b_[1], a_[1]-b_[0]
    assert hi < 0, 'reversal at t* not certified'
    rec.update(TV0_hot=dec(tv0[0]), TV0_cold=dec(tv0[1]), TV0_gap=dec(tv0[0]-tv0[1]), t_star=str(tstar),
               TVt_hot=[dec(a_[0]), dec(a_[1])], TVt_cold=[dec(b_[0]), dec(b_[1])],
               TVt_hot_minus_cold=[dec(lo), dec(hi)], positive_set=[i+1 for i in a_[2]],
               min_upward_rate_log10=float(-(L-2)*math.log10(1/q)), beta_b_decimal=math.log(1/q))
    return rec

def _roots_float(H):
    # float eigenvalues via characteristic polynomial (estimates only; certified by brackets)
    import cmath
    a = [[float(v) for v in r] for r in H]
    tr = a[0][0]+a[1][1]+a[2][2]
    m2 = (a[0][0]*a[1][1]-a[0][1]*a[1][0])+(a[0][0]*a[2][2]-a[0][2]*a[2][0])+(a[1][1]*a[2][2]-a[1][2]*a[2][1])
    det = (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
    # Newton from Gershgorin-ish starts
    roots = []
    f = lambda x: x**3-tr*x**2+m2*x-det; df = lambda x: 3*x*x-2*tr*x+m2
    for s in (-10.0, -3.0, -1.5, -0.5, 0.5):
        x = s
        for _ in range(200):
            dx = df(x)
            if dx == 0: break
            x -= f(x)/dx
        if abs(f(x)) < 1e-9 and all(abs(x-r) > 1e-6 for r in roots): roots.append(x)
    assert len(roots) == 3, roots
    return roots

QS = (F(1, 10), F(1, 20), F(1, 50), F(1, 100), F(1, 1000))
def main():
    L = 1555; uh, uc = 1-F(1, 10**5), 1-F(2, 10**5); ts = F(10)
    results = []
    for q in QS:
        Ht = cone_matrix(generator(*downward(L), q, [1, 2, 4]))
        assert Ht[1][0] == (EPS-q)/2 + q**(L-1)*(1-C*q)/4
        try:
            r = run(L, q, uh, uc, ts); r['certified'] = True
        except AssertionError as e:
            r = dict(L=L, q=str(q), certified=False, failure=str(e))
        print(json.dumps({k: v for k, v in r.items() if k != 'proper_subsystems'}), flush=True)
        results.append(r)
    (OUT/'e3_hierarchy_n4_certificate.json').write_text(json.dumps(dict(
        construction='Sec. 6, eps=1/100, m=1, C=199/100, z=23/25, E=(0,1,2,1555); upward rates a_ij q^(E_j-E_i)',
        preparations='u=e^{-beta}: u_hot=1-1e-5, u_cold=1-2e-5 (beta ~ 1e-5, 2e-5), both hotter than the bath (cooling side), beta_h<beta_c<beta_b',
        results=results), indent=2))

if __name__ == '__main__': main()
