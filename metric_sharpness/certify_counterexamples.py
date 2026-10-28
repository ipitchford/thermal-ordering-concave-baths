"""Exact rational certificates for heating divergence inversions.

No floating-point decisions are used. The only transcendental factor in the
uniformisation formula is bounded by rational partial sums and a geometric
remainder bound. Decimal output is for readability only.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
import json

OUT=Path(__file__).resolve().parent


def mm(a,b):
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*b)] for row in a]


def rowmul(v,a):
    return [sum((v[i]*a[i][j] for i in range(len(v))),F(0)) for j in range(len(v))]


def eye(n):
    return [[F(int(i==j)) for j in range(n)] for i in range(n)]


def inv(a):
    n=len(a);z=[r.copy()+e for r,e in zip(a,eye(n))]
    for i in range(n):
        j=next(j for j in range(i,n) if z[j][i])
        z[i],z[j]=z[j],z[i]
        x=z[i][i];z[i]=[v/x for v in z[i]]
        for j in range(n):
            if i!=j:
                x=z[j][i];z[j]=[a-x*b for a,b in zip(z[j],z[i])]
    return [r[n:] for r in z]


def cone_matrix(g):
    n=len(g);t=eye(n-1)
    for i in range(n-1):
        for j in range(i): t[i][j]=-F(1,i+1)
    q=[[g[i][j]-g[0][j] for j in range(1,n)] for i in range(1,n)]
    return mm(mm(t,q),inv(t))


def decimal(x):
    with localcontext() as ctx:
        ctx.prec=50
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def thermal(pi,beta):
    raw=[x**beta for x in pi];z=sum(raw)
    return [x/z for x in raw]


def uniformise(g,p0,t,kmax=100,rate=F(2)):
    """Return rational unnormalised series and exp-series common bounds."""
    n=len(g);p=eye(n)
    for i in range(n):
        for j in range(n):p[i][j]+=g[i][j]/rate
    assert all(v>=0 for r in p for v in r)
    assert all(sum(r)==1 for r in p)
    x=rate*t
    assert x<kmax+2
    term=F(1);expo=F(1);v=p0.copy();s=p0.copy()
    for k in range(1,kmax+1):
        term*=x/k;expo+=term;v=rowmul(v,p)
        s=[a+term*b for a,b in zip(s,v)]
    tail=term*x/(kmax+1)/(1-x/F(kmax+2))
    # exp(x) lies in [expo, expo+tail]. Each component's unnormalised
    # omitted series lies in [0,tail] since P^k is stochastic.
    lo=[a/(expo+tail) for a in s]
    hi=[(a+tail)/expo for a in s]
    return dict(series=s,expo=expo,tail=tail,lo=lo,hi=hi)


def check_model(name,weights,g,threshold,t=F(1047,100)):
    n=len(weights);pi=[x/sum(weights) for x in weights]
    assert all(sum(r)==0 for r in g)
    assert all(g[i][j]>=0 for i in range(n) for j in range(n) if i!=j)
    assert all(pi[i]*g[i][j]==pi[j]*g[j][i] for i in range(n) for j in range(n))
    h=cone_matrix(g)
    assert all(h[i][j]>=0 for i in range(n-1) for j in range(n-1) if i!=j)
    initial=[thermal(pi,b) for b in [3,4]]
    evolved=[uniformise(g,p,t) for p in initial]
    above=[];below=[]
    for i in range(n):
        cut=threshold*pi[i]
        if all(e['lo'][i]>cut for e in evolved):above.append(i)
        elif all(e['hi'][i]<cut for e in evolved):below.append(i)
        else:raise AssertionError('Common threshold sign not certified')
    assert above and below
    # With a common threshold set, the subtraction of k*pi cancels
    # between preparations. Bound the two probability sums directly.
    lower=sum(evolved[0]['lo'][i]-evolved[1]['hi'][i] for i in above)
    upper=sum(evolved[0]['hi'][i]-evolved[1]['lo'][i] for i in above)
    assert lower>0
    def divergence(p):
        return sum(max(p[i]-threshold*pi[i],F(0)) for i in range(n))-max(1-threshold,F(0))
    d0=[divergence(p) for p in initial]
    assert d0[1]>d0[0]
    dt=[[sum(e[key][i]-threshold*pi[i] for i in above)-max(1-threshold,F(0))
         for key in ['lo','hi']] for e in evolved]
    minrate=min(g[i][j] for i in range(n) for j in range(n) if i!=j)
    return dict(name=name,beta_b=1,beta_warm=3,beta_cold=4,
        weights=[str(v) for v in weights],pi=[str(v) for v in pi],
        energy='E_i = log(w_1/w_i), so bath beta_b=1 and all energies are strictly ordered.',
        generator=[[str(v) for v in r] for r in g],
        transformed_quotient=[[str(v) for v in r] for r in h],
        threshold=str(threshold),time=str(t),uniformisation_rate='2',terms=101,
        positive_set_one_based=[i+1 for i in above],negative_set_one_based=[i+1 for i in below],
        all_offdiagonal_rates_strictly_positive=bool(minrate>0),minimum_offdiagonal_rate=str(minrate),
        initial_divergences=[decimal(v) for v in d0],
        initial_cold_minus_warm_gap=decimal(d0[1]-d0[0]),
        final_divergence_intervals=[[decimal(v) for v in pair] for pair in dt],
        final_warm_minus_cold_gap_lower=decimal(lower),final_warm_minus_cold_gap_upper=decimal(upper),
        final_likelihood_ratio_intervals=[[[decimal(e['lo'][i]/pi[i]),decimal(e['hi'][i]/pi[i])] for i in range(n)] for e in evolved],
        maximum_probability_interval_width=decimal(max(e['hi'][i]-e['lo'][i] for e in evolved for i in range(n))),
        unnormalised_exponential_series_tail_bound=decimal(evolved[0]['tail']),
        certificate='All assertions and interval endpoints computed as exact fractions; decimal rendering only.')


def main():
    g4=[[-F(113,552),F(203,1104),F(1,48),F(0)],
        [F(145,184),-F(313,368),F(0),F(1,16)],
        [F(5,8),F(0),-F(11,16),F(1,16)],
        [F(0),F(7,8),F(1,8),-F(1)]]
    out4=check_model('four-state hinge-divergence heating inversion',list(map(F,[60,14,2,1])),g4,F(9992521,10000000))
    print('four-state certificate passed; final gap lower',out4['final_warm_minus_cold_gap_lower'],flush=True)
    w5=list(map(F,[60,14,2,1]))+[F(5763,100000)]
    pi5=[x/sum(w5) for x in w5]
    g5=[r+[F(0)] for r in g4]+[[F(0)]*5]
    eta=F(1,100)
    for i in range(5):
        for j in range(5):g5[i][j]+=eta*(pi5[j]-int(i==j))
    out5=check_model('five-state total-variation heating inversion, all rates positive',w5,g5,F(1))
    out5['reset_rate']=str(eta)
    print('five-state certificate passed; final gap lower',out5['final_warm_minus_cold_gap_lower'],flush=True)
    (OUT/'exact_counterexamples.json').write_text(json.dumps(dict(models=[out4,out5]),indent=2))


if __name__=='__main__':main()
