#!/usr/bin/env python3
"""Decimal illustrations of the positive-bath theorem; not interval certificates."""
import json
from decimal import Decimal as D, getcontext
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from verify_hierarchy import rates
ROOT=Path(__file__).resolve().parent
getcontext().prec=80
getcontext().Emin=-999999999
getcontext().Emax=999999999

def number(q):return D(q.numerator)/q.denominator
def mv(A,x):return [sum(a*b for a,b in zip(row,x)) for row in A]
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def solve(A,b):
 n=len(b);a=[row[:]+[v] for row,v in zip(A,b)]
 for k in range(n):
  pivot=max(range(k,n),key=lambda i:abs(a[i][k]));a[k],a[pivot]=a[pivot],a[k]
  assert a[k][k]
  p=a[k][k];a[k]=[x/p for x in a[k]]
  for i in range(n):
   if i==k:continue
   c=a[i][k];a[i]=[x-c*y for x,y in zip(a[i],a[k])]
 return [row[-1] for row in a]
def build(n,E,beta):
 a,_=rates(n);M=[[D(0)]*n for _ in range(n)]
 for i in range(n):
  for j in range(i+1,n):
   M[i][j]=number(a[i][j]);M[j][i]=M[i][j]*(-beta*(E[j]-E[i])).exp()
 for j in range(n):M[j][j]=-sum(M[i][j] for i in range(n) if i!=j)
 return M
def restrict(M,S):
 A=[[M[i][j] if i!=j else D(0) for j in S] for i in S]
 for j in range(len(S)):A[j][j]=-sum(A[i][j] for i in range(len(S)) if i!=j)
 return A
def slow(M):
 n=len(M);G=list(map(list,zip(*M)))
 Q=[[G[i][j]-G[0][j] for j in range(1,n)] for i in range(1,n)]
 lam=Q[0][0]
 for step in range(100):
  B=[[Q[i][j]-(lam if i==j else 0) for j in range(1,n-1)] for i in range(1,n-1)]
  h=[D(1)]+solve(B,[-Q[i][0] for i in range(1,n-1)])
  new=dot(Q[0],h)
  if abs(new-lam)<D('1e-65'):
   lam=new;break
  lam=new
 else:raise AssertionError('Slow root fixed-point iteration failed')
 g=[D(0)]+h
 residual=max(abs(dot(Q[i],h)-lam*h[i]) for i in range(n-1))
 assert residual<D('1e-60')
 margins=[g[j]-sum(g[:j])/j for j in range(1,n)]
 return lam,g,margins,residual
def gibbs(E,beta):
 w=[(-beta*x).exp() for x in E];Z=sum(w);return [x/Z for x in w]
def evolve(M,preparations,t):
 # Positive Poisson uniformization, generalised from the fresh review audit.
 n=len(M);rate=max(-M[j][j] for j in range(n))
 P=[[D(i==j)+M[i][j]/rate for j in range(n)] for i in range(n)]
 assert min(x for row in P for x in row)>=0
 x=rate*t;weight=(-x).exp();terms=[p[:] for p in preparations]
 states=[[weight*y for y in p] for p in preparations];k=0
 while True:
  nxt=weight*x/(k+1)
  if k+2>x:
   tail=nxt/(1-x/(k+2))
   if tail<D('1e-65'):return states,tail,k
  k+=1;terms=[mv(P,p) for p in terms];weight=nxt
  states=[[v+weight*u for v,u in zip(state,term)] for state,term in zip(states,terms)]
  assert k<5000
def tv(p,q):return sum(abs(x-y) for x,y in zip(p,q))/2
def mode_isolation_bound(M,a):
 # Exact eigenvector conditioning for the triangular zero-bath quotient,
 # followed by a numerical Bauer-Fike perturbation bound. This records why
 # the root near the first escape rate is still the slow root.
 n=len(a);k=n-1;esc=[sum(a[i][j] for i in range(j)) for j in range(n)]
 V=[[F(0)]*k for _ in range(k)]
 for c in range(k):
  V[c][c]=F(1)
  for r in range(c+1,k):
   V[r][c]=sum(a[i+1][r+1]*V[i][c] for i in range(r))/(esc[r+1]-esc[c+1])
 invcols=[solve(V,[F(i==j) for i in range(k)]) for j in range(k)]
 Vinv=list(map(list,zip(*invcols)))
 cond=max(sum(abs(x) for x in row) for row in V)*max(sum(abs(x) for x in row) for row in Vinv)
 G=list(map(list,zip(*M)))
 Q=[[G[i][j]-G[0][j] for j in range(1,n)] for i in range(1,n)]
 Q0=[[a[j][i] if i!=j else -esc[i] for j in range(1,n)] for i in range(1,n)]
 perturb=max(sum(abs(Q[i][j]-number(Q0[i][j])) for j in range(k)) for i in range(k))
 radius=number(cond)*perturb
 gap=min([esc[1]]+[esc[j]-esc[j-1] for j in range(2,n)])
 assert radius<number(gap)/3
 return radius/number(gap)
def main():
 out=[]
 for n in range(4,9):
  m=n-3;P=196*m*m+294*m-199;L=(P*n*n+3*(m+1))//3+1
  E=list(map(D,range(n-1)))+[D(L)]
  beta_b=D(20);M=build(n,E,beta_b);lam,g,margins,residual=slow(M)
  assert margins[-1]<0
  a,_=rates(n);isolation_ratio=mode_isolation_bound(M,a)
  local_min=D('Infinity');count=0
  for size in range(2,n):
   for S in combinations(range(n),size):
    sub=restrict(M,S)
    _,_,local,res=slow(sub);assert min(local)>0
    isolation_ratio=max(isolation_ratio,mode_isolation_bound(sub,[[a[i][j] for j in S] for i in S]))
    local_min=min(local_min,*local);residual=max(residual,res);count+=1
  pi=gibbs(E,beta_b);centre=dot(g,pi);f=[x-centre for x in g]
  beta_h=D(1)/(100*L);beta_c=D(2)/(100*L)
  ph=gibbs(E,beta_h);pc=gibbs(E,beta_c);ah=dot(f,ph);ac=dot(f,pc)
  assert 0<ah<ac
  initial=tv(ph,pi)-tv(pc,pi)
  states,tail,terms=evolve(M,[ph,pc],D(20));later=tv(states[1],pi)-tv(states[0],pi)
  assert initial>0 and later>0
  derivative=(sum(g)*sum(E)-n*dot(g,E))/n**2
  assert derivative>0 and min(M[i][j] for i in range(n) for j in range(n) if i!=j)>0
  out.append({'N':n,'L':L,'beta_b':'20','beta_h':str(beta_h),'beta_c':str(beta_c),
   'proper_subsystems_checked':count,'min_proper_prefix_margin':str(local_min),
   'full_final_prefix_margin':str(margins[-1]),'full_slow_eigenvalue':str(lam),
   'uniform_amplitude_derivative':str(derivative),'cold_minus_hot_slow_amplitude':str(ac-ah),
   'initial_TV_hot_minus_cold':str(initial),'time20_TV_cold_minus_hot':str(later),
   'min_off_diagonal_rate':str(min(M[i][j] for i in range(n) for j in range(n) if i!=j)),
   'min_rate_log10':str(min(M[i][j].log10() for i in range(n) for j in range(n) if i!=j)),
   'quotient_eigenfunction_residual_max':str(residual),'bauer_fike_radius_over_min_gap_max':str(isolation_ratio),'time20_poisson_tail_bound':str(tail),
   'time20_last_poisson_index':terms})
 result={'status':'PASS','precision_decimal_digits':getcontext().prec,'decimal_Emin':getcontext().Emin,'decimal_Emax':getcontext().Emax,
  'method':'Decimal quotient-eigenvector iteration and positive Poisson uniformization; not interval arithmetic',
  'instances':out,
  'limits':['Finite computations illustrate the analytic existence theorem.',
   'No all-time proper-subsystem TV assertion at positive bath temperature.',
   'The engineered energy scales and tiny reverse rates are not claimed experimentally realistic.']}
 (ROOT/'positive_bath_numerics.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
