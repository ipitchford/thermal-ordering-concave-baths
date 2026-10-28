"""E1 audit: exact sign checks of the curvature factor of standard jump weights
(stdlib Fractions). J'' = (positive factor) * Q, with Q polynomial/rational:
  d^s e^{-d/c}:  Q(y)=y^2-2sy+s(s-1), y=d/c          threshold y=s+sqrt(s) (0<s<=1)
  Drude c d/(d^2+c^2): Q = -(3c^2-d^2)                threshold d=sqrt(3)c
  d e^{-d^2/c^2}:     Q = -(3-2y^2), y=d/c            threshold y=sqrt(3/2)
  d^s:                Q = s(s-1)                       concave iff 0<=s<=1
Checks signs at rational points just inside/outside each threshold, plus a
second-difference (midpoint) test of J itself using rational exp bounds.
"""
from fractions import Fraction as F
import json
from pathlib import Path

def exp_bounds(x, n=60):
    # x rational, |x|<=10: partial sum bounds of e^{x} for x>=0, reciprocal for x<0
    neg = x < 0; x = abs(x); s = F(1); t = F(1)
    for k in range(1, n+1): t *= x/k; s += t
    tail = t*x/(n+1)/(1-x/F(n+2))
    lo, hi = s, s+tail
    return (1/hi, 1/lo) if neg else (lo, hi)

def Q_exp(s, y): return y*y-2*s*y+s*(s-1)
res = {}
# exponential cutoff, s in {1/4, 1/2, 1}: threshold s+sqrt(s) is rational for these s
for s, root in ((F(1, 4), F(3, 4)), (F(1, 2), None), (F(1), F(2))):
    if root is not None:
        assert Q_exp(s, root) == 0
        assert Q_exp(s, root-F(1, 10**6)) < 0 < Q_exp(s, root+F(1, 10**6))
    # s=1/2: threshold 1/2+1/sqrt2 irrational; check rational bracket
    if s == F(1, 2):
        a, b = F(12071, 10000), F(12072, 10000)  # 1/2+sqrt(1/2)=1.20710678...
        assert Q_exp(s, a) < 0 < Q_exp(s, b)
    assert all(Q_exp(s, F(k, 100)) < 0 for k in range(1, int((root or F(12071, 10000))*100)))
    res['exp_cutoff_s=%s' % s] = 'threshold %s*w_c confirmed' % (root if root is not None else '1/2+sqrt(1/2)')
for s in (F(3, 2), F(3)):
    assert Q_exp(s, F(1, 1000)) > 0  # convex near zero for s>1
    res['exp_cutoff_s=%s' % s] = 'convex near 0: not concave on any [0,D]'
# Drude: sign(J'') = -sign(3c^2-d^2), c=1
a, b = F(17320, 10000), F(17321, 10000)
assert 3-a*a > 0 > 3-b*b
res['drude'] = 'threshold sqrt(3) w_c confirmed (1.7320 < sqrt3 < 1.7321)'
# Gaussian Ohmic
a, b = F(12247, 10000), F(12248, 10000)
assert 3-2*a*a > 0 > 3-2*b*b
res['ohmic_gaussian'] = 'threshold sqrt(3/2) w_c confirmed'
# direct midpoint test on J itself for Ohmic J=d e^{-d} (c=1): concavity violated across (2.4, 3.6)
def J_bounds(d): lo, hi = exp_bounds(-d); return d*lo, d*hi
x, h = F(3), F(3, 5)
l1, h1 = J_bounds(x-h); l2, h2 = J_bounds(x+h); lm, hm = J_bounds(x)
assert (l1+l2)/2 > hm  # midpoint below chord: not concave on [0,3.6]
x, h = F(1), F(1, 2)
l1, h1 = J_bounds(x-h); l2, h2 = J_bounds(x+h); lm, hm = J_bounds(x)
assert (h1+h2)/2 < lm  # concave sample inside [0,2]
res['ohmic_midpoint'] = 'J(3)<[J(2.4)+J(3.6)]/2 certified; J(1)>[J(.5)+J(1.5)]/2 certified'
(Path(__file__).resolve().parent/'e1_threshold_check.json').write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
