#!/usr/bin/env python3
"""Exact finite audits for the concave-kernel theorem; standard library only.

The general result is proved analytically in concave_kernel_theorem.md. These
checks audit formulas, orientations, scaling, boundary cases and small instances.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import random
import time


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def mm(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def close_rows(a):
    for i in range(len(a)):
        a[i][i] = -sum(a[i][j] for j in range(len(a)) if j != i)
    return a


def matrices(E, qb, J):
    n = len(E)
    S, A = [[F(0) for _ in E] for _ in E], [[F(0) for _ in E] for _ in E]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = abs(E[i]-E[j])
            S[i][j] = J(d)*qb**d/(1-qb**d)
            A[i][j] = J(d) if i > j else F(0)
    close_rows(S); close_rows(A)
    G = [[S[i][j]+A[i][j] for j in range(n)] for i in range(n)]
    return G, S, A


def quotient_transform(G):
    n = len(G)
    Q = [[G[i][j]-G[0][j] for j in range(1,n)] for i in range(1,n)]
    T, Ti = eye(n-1), eye(n-1)
    for i in range(n-1):
        for j in range(i):
            T[i][j] = -F(1,i+1)
            Ti[i][j] = F(1,j+2)
    assert mm(T,Ti) == eye(n-1)
    return mm(mm(T,Q),Ti)


def helmert_transform(G):
    n = len(G)
    U = [[F(-1) if k<j-1 else F(j-1) if k==j-1 else F(0)
          for k in range(n)] for j in range(2,n+1)]
    out = [[F(0) for _ in U] for _ in U]
    for i in range(n-1):
        for j in range(n-1):
            numerator = sum((U[i][a]*G[a][b]*U[j][b]
                             for a in range(n) for b in range(n)), F(0))
            out[i][j] = numerator/((i+1)*(j+2))
    return out


def gibbs(E, q):
    p = [q**e for e in E]
    z = sum(p)
    return [x/z for x in p]


def row_step(p, P):
    return [sum((p[i]*P[i][j] for i in range(len(p))), F(0))
            for j in range(len(p))]


def main():
    start = time.perf_counter()
    rng = random.Random(20261026)
    checks = 0
    records = []
    for n in [3,4,5,7]:
        random_E = [0]
        for _ in range(n-1):
            random_E.append(random_E[-1]+rng.randint(1,4))
        for E in [list(range(n)), [i*i for i in range(n)], random_E]:
            D = E[-1]
            kernels = {
                'linear': lambda d: F(d),
                'constant': lambda d: F(1),
                'concave_quadratic': lambda d: F(d)-F(d*d,2*D),
                'concave_exponential': lambda d: F(1)-F(1,2)**d,
            }
            for name,J in kernels.items():
                for qb in [F(1,2), F(2,3)]:
                    G,S,A = matrices(E,qb,J)
                    H,HS,HA = map(quotient_transform,(G,S,A))
                    assert H == helmert_transform(G)
                    checks += 1
                    pi = gibbs(E,qb)
                    assert all(sum(G[i]) == 0 for i in range(n))
                    assert all(pi[i]*G[i][j] == pi[j]*G[j][i]
                               for i in range(n) for j in range(n))
                    checks += 2
                    h = lambda d: J(d)*qb**d/(1-qb**d)
                    for ii in range(n-1):
                        for jj in range(n-1):
                            if ii == jj:
                                continue
                            i,j = ii+2,jj+2
                            assert H[ii][jj] > 0
                            assert HS[ii][jj] >= 0 and HA[ii][jj] >= 0
                            checks += 2
                            if i < j:
                                upper = h(E[j-1]-E[i-1])-sum(h(E[j-1]-E[k]) for k in range(i-1))/(i-1)
                                upper += sum((h(E[l]-E[i-1])-sum(h(E[l]-E[k]) for k in range(i-1))/(i-1)
                                              for l in range(j,n)),F(0))/j
                                assert HS[ii][jj] == upper
                                assert HS[jj][ii] == F((i-1)*j,(j-1)*i)*upper
                                assert HA[ii][jj] == 0
                                checks += 3
                            else:
                                sj = sum(J(E[j-1]-E[l]) for l in range(j-1))
                                fj = lambda x: (j-1)*J(x-E[j-1])-sum(J(x-E[l]) for l in range(j-1))
                                numerator = j*(fj(E[i-1])+sj)+sum(fj(E[i-1])-fj(E[k-1]) for k in range(j+1,i))
                                assert HA[ii][jj] == numerator/((i-1)*j)
                                checks += 1
                    # Exact uniformization-step audit. R controls both G and H,
                    # so P is stochastic and its quotient maps the orthant to itself.
                    R = max([-G[i][i] for i in range(n)]+[-H[i][i] for i in range(n-1)])+1
                    P = [[F(i==j)+G[i][j]/R for j in range(n)] for i in range(n)]
                    assert all(x>=0 for row in P for x in row)
                    ph,pc = gibbs(E,F(3,4)),gibbs(E,F(1,4))
                    for _ in range(4):
                        ph,pc = row_step(ph,P),row_step(pc,P)
                        assert sum(ph)==sum(pc)==1
                        assert all(sum(pc[:k])>sum(ph[:k]) for k in range(1,n))
                        checks += 2
                    records.append({'N':n,'E':E,'kernel':name,'bath_q':str(qb),
                                    'minimum_offdiagonal_H':str(min(H[i][j] for i in range(n-1) for j in range(n-1) if i!=j))})
    G,_,_ = matrices([0,1,2],F(1,2),lambda d:F(d))
    assert G[1][0]==2 and G[2][0]==F(8,3)
    checks += 1
    _,_,A = matrices([0,1,2],F(1,2),lambda d:F(d*d))
    assert quotient_transform(A)[1][0] == -1
    checks += 1
    G,_,_ = matrices([0,1,2],F(1,16),lambda d:F(d*d))
    nonconcave_finite = quotient_transform(G)[1][0]
    assert nonconcave_finite < 0
    checks += 1
    output = {'status':'PASS','matrix_cases':len(records),'exact_assertions':checks,
              'finite_nonconcave_H_32':str(nonconcave_finite),
              'scope':'Exact finite formula and orientation audit; universal statements rely on the written proofs.',
              'cases':records,'elapsed_seconds':time.perf_counter()-start,
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('concave_kernel_checks.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='cases'},indent=2))


if __name__ == '__main__':
    main()
