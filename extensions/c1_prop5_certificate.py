"""C1: exact component intervals for Proposition 5 (five-state TV heating
reversal), rational bounds for the sign quantity, rational certificates that
bracket BOTH crossings, a 60-digit Decimal scan of the reversal window, and a
robustness check under small rational perturbations of w_5 and eta.

Adapted from metric_sharpness/certify_counterexamples.py (copied, not imported).
Standard library only. Exact decisions use Fractions; Decimal is used only for
the high-precision scan (numerical, labelled as such).
"""
from fractions import Fraction as F
from decimal import Decimal, getcontext, localcontext
from pathlib import Path
import json

OUT = Path(__file__).resolve().parent

def eye(n): return [[F(int(i == j)) for j in range(n)] for i in range(n)]
def rowmul(v, a): return [sum((v[i]*a[i][j] for i in range(len(v))), F(0)) for j in range(len(v))]
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
def dec(x, p=40):
    with localcontext() as c:
        c.prec = p; return str(Decimal(x.numerator)/Decimal(x.denominator))
def thermal(w, beta):
    raw = [x**beta for x in w]; z = sum(raw); return [x/z for x in raw]

G4 = [[-F(113,552), F(203,1104), F(1,48), F(0)],
      [F(145,184), -F(313,368), F(0), F(1,16)],
      [F(5,8), F(0), -F(11,16), F(1,16)],
      [F(0), F(7,8), F(1,8), -F(1)]]

def model(w5, eta):
    w = list(map(F, [60, 14, 2, 1])) + [w5]; pi = [x/sum(w) for x in w]
    g = [r+[F(0)] for r in G4] + [[F(0)]*5]
    for i in range(5):
        for j in range(5): g[i][j] += eta*(pi[j]-int(i == j))
    n = 5
    assert all(sum(r) == 0 for r in g)
    assert all(g[i][j] > 0 for i in range(n) for j in range(n) if i != j)
    assert all(pi[i]*g[i][j] == pi[j]*g[j][i] for i in range(n) for j in range(n))
    h = cone_matrix(g)
    assert all(h[i][j] >= 0 for i in range(n-1) for j in range(n-1) if i != j)
    return w, pi, g, h

def uniformise(g, p0, t, kmax=100, rate=F(2)):
    n = len(g); P = eye(n)
    for i in range(n):
        for j in range(n): P[i][j] += g[i][j]/rate
    assert all(v >= 0 for r in P for v in r) and all(sum(r) == 1 for r in P)
    x = rate*t; assert x < kmax+2
    term = F(1); B = F(1); v = p0.copy(); s = p0.copy()
    for k in range(1, kmax+1):
        term *= x/k; B += term; v = rowmul(v, P); s = [a+term*b for a, b in zip(s, v)]
    T = term*x/(kmax+1)/(1-x/F(kmax+2))
    return [a/(B+T) for a in s], [(a+T)/B for a in s], T

def outward(x, sig=12, up=False):
    """Decimal string with `sig` significant digits rounded outward; exactness checked with Fractions."""
    from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
    with localcontext() as c:
        c.prec = sig; c.rounding = ROUND_CEILING if up else ROUND_FLOOR
        d = +(Decimal(x.numerator) / Decimal(x.denominator))
    fd = F(str(d))
    assert (fd >= x) if up else (fd <= x)
    return str(d)

def tv_interval(lo, hi, pi):
    """Rigorous TV enclosure; requires certified sign of every p_i - pi_i."""
    pos = []
    for i in range(len(pi)):
        if lo[i] > pi[i]: pos.append(i)
        elif hi[i] < pi[i]: pass
        else: raise AssertionError('component sign not certified')
    return sum(lo[i]-pi[i] for i in pos), sum(hi[i]-pi[i] for i in pos), pos

def certify_at(w5, eta, t):
    w, pi, g, h = model(w5, eta)
    p3, p4 = thermal(w, 3), thermal(w, 4)
    e3 = uniformise(g, p3, t); e4 = uniformise(g, p4, t)
    a = tv_interval(e3[0], e3[1], pi); b = tv_interval(e4[0], e4[1], pi)
    lo, hi = a[0]-b[1], a[1]-b[0]   # TV(p3(t)) - TV(p4(t))
    sign = '+' if lo > 0 else ('-' if hi < 0 else '?')
    return dict(lo=lo, hi=hi, sign=sign, pos3=a[2], pos4=b[2], e3=e3, e4=e4, pi=pi)

# ---------- 60-digit Decimal scan (numerical corroboration) ----------
getcontext().prec = 60
def dmodel(w5, eta, kmax):
    w, pi, g, h = model(w5, eta)
    D = lambda x: Decimal(x.numerator)/Decimal(x.denominator)
    P = [[D(int(i == j)+g[i][j]/2) for j in range(5)] for i in range(5)]
    pid = [D(x) for x in pi]
    vs = []
    for beta in (3, 4):
        v = [D(x) for x in thermal(w, beta)]; seq = [v]
        for k in range(kmax):
            v = [sum(v[i]*P[i][j] for i in range(5)) for j in range(5)]; seq.append(v)
        vs.append(seq)
    return pid, vs
def dTVdiff(dm, t):
    pid, vs = dm; x = 2*Decimal(t); wk = (-x).exp(); acc = [[Decimal(0)]*5 for _ in range(2)]
    K = len(vs[0])-1
    for k in range(K+1):
        if k: wk = wk*x/k
        for s in range(2):
            vk = vs[s][k]
            for i in range(5): acc[s][i] += wk*vk[i]
    tv = [sum(max(a-p, Decimal(0)) for a, p in zip(acc[s], pid)) for s in range(2)]
    return tv[0]-tv[1]
def bisect(dm, a, b, tol=Decimal('1e-45')):
    fa = dTVdiff(dm, a)
    while b-a > tol:
        c = (a+b)/2; fc = dTVdiff(dm, c)
        if (fc > 0) == (fa > 0): a, fa = c, fc
        else: b = c
    return (a+b)/2
def scan(dm, t0, t1, step):
    roots = []; t = Decimal(t0); prev = dTVdiff(dm, t); best = (prev, t)
    while t < Decimal(t1):
        t2 = t+Decimal(step); cur = dTVdiff(dm, t2)
        if cur > best[0]: best = (cur, t2)
        if (cur > 0) != (prev > 0): roots.append(bisect(dm, t, t2))
        t, prev = t2, cur
    return roots, best

def main():
    W5, ETA, TS = F(5763, 100000), F(1, 100), F(1047, 100)
    res = {}
    # (1) exact component intervals at t*
    c = certify_at(W5, ETA, TS); assert c['sign'] == '+'
    print('Proposition 5, t* = 1047/100: rational component intervals (decimal rendering of exact endpoints)')
    comp = []
    for name, e in (('beta=3', c['e3']), ('beta=4', c['e4'])):
        for i in range(5):
            row = dict(prep=name, i=i+1, lo=dec(e[0][i]), hi=dec(e[1][i]), pi=dec(c['pi'][i]),
                       width=dec(e[1][i]-e[0][i], 6), sign='+' if e[0][i] > c['pi'][i] else '-',
                       lo_exact=str(e[0][i]), hi_exact=str(e[1][i]))
            comp.append(row)
            print('  %s p_%d in [%s, %s]  pi=%s  sign %s  width %s' % (name, i+1, row['lo'], row['hi'], row['pi'][:20], row['sign'], row['width']))
    print('  sign quantity TV(p3)-TV(p4) in [%s, %s]' % (dec(c['lo']), dec(c['hi'])))
    res['t_star'] = dict(t='1047/100', components=comp, positive_set_beta3=[i+1 for i in c['pos3']],
        positive_set_beta4=[i+1 for i in c['pos4']], sign_quantity_lower_exact=str(c['lo']),
        sign_quantity_upper_exact=str(c['hi']), sign_quantity_lower=dec(c['lo']), sign_quantity_upper=dec(c['hi']),
        tail_bound_T=dec(c['e3'][2], 10))
    # (2) rational brackets for both crossings
    brackets = []
    for t in (F(1046, 100), F(1047, 100), F(1048, 100)):
        r = certify_at(W5, ETA, t)
        brackets.append(dict(t=str(t), lower=dec(r['lo'], 20), upper=dec(r['hi'], 20), sign=r['sign'],
                             lower_outward=outward(r['lo']), upper_outward=outward(r['hi'], up=True),
                             positive_set_beta3=[i+1 for i in r['pos3']], positive_set_beta4=[i+1 for i in r['pos4']]))
        print('  certified TV(p3)-TV(p4) at t=%s in [%s, %s] sign %s' % (t, dec(r['lo'], 15), dec(r['hi'], 15), r['sign']))
    assert [b['sign'] for b in brackets] == ['-', '+', '-']
    res['certified_brackets'] = brackets
    # (3) 60-digit scan
    dm = dmodel(W5, ETA, 260)
    roots, best = scan(dm, '0', '60', '0.002')
    print('  60-digit scan on [0,60], step 0.002: sign changes at', [str(r)[:22] for r in roots])
    print('  max of TV(p3)-TV(p4):', str(best[0])[:14], 'at t ~', best[1])
    assert len(roots) == 2
    res['decimal_scan'] = dict(precision_digits=60, interval='[0,60]', step='0.002', series_terms=261,
        crossings=[str(r)[:50] for r in roots], window_length=str(roots[1]-roots[0])[:20],
        max_reversal=str(best[0])[:20], max_reversal_t=str(best[1]),
        initial_gap_TV4_minus_TV3=str(-dTVdiff(dm, Decimal(0)))[:20])
    # (4) robustness: perturb w5 and eta by small rationals
    rob = []
    for dw in (F(-1, 10**5), F(-1, 10**6), F(0), F(1, 10**6), F(1, 10**5)):
        for de in (F(-1, 10**5), F(0), F(1, 10**5)):
            if dw == 0 and de == 0: continue
            w5, eta = W5+dw, ETA+de
            dmp = dmodel(w5, eta, 140)
            rts, bst = scan(dmp, '10.3', '10.7', '0.001')
            entry = dict(w5=str(w5), eta=str(eta), window=[str(r)[:12] for r in rts], max_reversal=str(bst[0])[:12])
            if len(rts) == 2:
                tm = F(str(bst[1]))  # rational time at the scan maximum
                r = certify_at(w5, eta, tm); entry['certified_at_t'] = str(tm)
                entry['certified_lower'] = dec(r['lo'], 12); entry['certified_sign'] = r['sign']
            rob.append(entry); print('  robustness', entry, flush=True)
    res['robustness'] = rob
    (OUT/'c1_prop5_certificate.json').write_text(json.dumps(res, indent=2))

if __name__ == '__main__': main()
