#!/usr/bin/env python3
"""Exact finite checks of the hierarchy construction; not a substitute for proof."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def rates(n,epsilon=F(1,100)):
    m=n-3
    z=1-2*epsilon-3*epsilon/(m-F(1,2))
    a=[[F(0) for _ in range(n)] for _ in range(n)]
    a[0][1]=F(1)
    for j in range(2,n-1):
        a[0][j]=F(3,2);a[1][j]=F(1)
        for i in range(2,j):a[i][j]=epsilon
    a[0][-1]=2-epsilon;a[1][-1]=F(1)
    for i in range(2,n-1):a[i][-1]=z
    return a,z

def restriction(a,S):
    return [[a[i][j] for j in S] for i in S]

def slow_response(a):
    n=len(a)
    escape=[sum(a[i][j] for i in range(j)) for j in range(n)]
    assert all(escape[j]>escape[j-1] for j in range(1,n))
    g=[F(0),F(1)]
    for j in range(2,n):
        g.append(sum(a[i][j]*g[i] for i in range(j))/(escape[j]-escape[1]))
    for j in range(n):
        assert sum(a[i][j]*(g[i]-g[j]) for i in range(j))==-escape[1]*g[j]
    return g,escape

def prefix_margins(g):
    return [g[j]-sum(g[:j])/j for j in range(1,len(g))]

def all_time_certificate(n,S,epsilon,z):
    """Check positive coefficients of the exact all-time proof, not time samples."""
    m=n-3; C=2-epsilon
    r=sum(2<=j<n-1 for j in S)
    if 0 in S and 1 in S and n-1 in S:
        obstruction=epsilon*(1-F(r)/(m-F(1,2)))
        kappa=F(r)*(2*z-F(1,2)+epsilon)/(r+2)
        assert obstruction>0 and kappa>=0
        return 'both_lower_states'
    if 0 in S and 1 not in S and n-1 in S and r>=1:
        b=r*z/(C+r*z-F(3,2))
        assert b-F(r,r+1)>0 and F(1,r+1)>0
        return 'ground_without_slow'
    if 0 not in S and 1 not in S and n-1 in S and r>=2:
        b=(r-1)*z/(r*z-epsilon)
        assert b-F(r-1,r)>0 and F(1,r)>0
        return 'middle_ground'
    return 'common_or_two_exponentials'

def main():
    rows=[];total=0;classes={}
    for n in range(4,13):
        m=n-3;a,z=rates(n);g,escape=slow_response(a)
        full=prefix_margins(g)
        P=196*m*m+294*m-199
        assert all(v>0 for v in full[:-1]) and full[-1]<0
        assert full[-1]==-F(1,(m+2)*P)
        threshold=F(m+1,2)+F(P*(m+3)**2,6)
        L=(2*threshold).numerator//(2*threshold).denominator+1
        E=list(map(F,range(n-1)))+[F(L)]
        derivative=(sum(g)*sum(E)-n*sum(x*y for x,y in zip(g,E)))/n**2
        expected=(F(L)-F(m+1,2)-F(P*n*n,6))/(P*n*n)
        assert derivative==expected and derivative>0
        count=0; min_margin=None
        for size in range(2,n):
            for S in combinations(range(n),size):
                sub=restriction(a,S);gs,_=slow_response(sub);margins=prefix_margins(gs)
                assert min(margins)>0,(n,S,margins)
                min_margin=min(margins) if min_margin is None else min(min_margin,*margins)
                label=all_time_certificate(n,S,F(1,100),z)
                classes[label]=classes.get(label,0)+1
                count+=1
        expected_count=2**n-n-2
        assert count==expected_count
        total+=count
        rows.append({'N':n,'proper_subsystems_checked':count,'min_proper_prefix_margin':str(min_margin),
                     'full_final_prefix_margin':str(full[-1]),'escape_rates':list(map(str,escape)),
                     'energy_top':L,'uniform_derivative':str(derivative)})
    # Distinct rational choices test the open parameter formulas without
    # pretending finite coverage proves the continuous interval.
    parameter_rows=[]
    for e in [F(1,50),F(1,20),F(1,15)]:
        for n in range(4,9):
            a,z=rates(n,e);g,_=slow_response(a)
            assert prefix_margins(g)[-1]<0
            count=0
            for size in range(2,n):
                for S in combinations(range(n),size):
                    gs,_=slow_response(restriction(a,S))
                    assert min(prefix_margins(gs))>0
                    all_time_certificate(n,S,e,z);count+=1
            parameter_rows.append({'N':n,'epsilon':str(e),'proper_subsystems_checked':count})
    result={'status':'PASS','method':'fractions.Fraction, triangular back-substitution, exact survival proof coefficients',
            'primary_range':'N=4,...,12; epsilon=1/100','primary_proper_subsystems':total,
            'all_time_certificate_classes':classes,'instances':rows,'other_parameter_checks':parameter_rows,
            'limits':['Finite arithmetic checks do not prove the unrestricted theorem.',
                      'Positive-bath continuation is proved analytically in hierarchy_theorem.md, not by this script.',
                      'No strong Mpemba cancellation or finite-bath all-time local TV claim is tested.']}
    (ROOT/'hierarchy_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
