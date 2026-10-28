from fractions import Fraction as F
from pathlib import Path
import json
import hashlib
from itertools import permutations

def matmul(a,b):
 return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0)) for j in range(len(b[0]))] for i in range(len(a))]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def characteristic(a):
 n=len(a);b=eye(n);coeff=[F(1)]
 for k in range(1,n+1):
  ab=matmul(a,b);c=-sum(ab[i][i] for i in range(n))/k;coeff.append(c)
  b=[[ab[i][j]+(c if i==j else 0) for j in range(n)] for i in range(n)]
 return coeff

def model(E,q=F(1,2),family='linear'):
 n=len(E);G=[[F(0) for _ in E] for _ in E]
 for i in range(n):
  for j in range(n):
   if i==j:continue
   d=abs(E[i]-E[j]);J=F(d) if family=='linear' else 1-q**d
   G[i][j]=J*q**d/(1-q**d)+(J if i>j else 0)
 for i in range(n):G[i][i]=-sum(G[i])
 Q=[[G[i][j]-G[0][j] for j in range(1,n)] for i in range(1,n)]
 T=eye(n-1);Ti=eye(n-1)
 for i in range(n-1):
  for j in range(i):T[i][j]=-F(1,i+1);Ti[i][j]=F(1,j+2)
 H=matmul(matmul(T,Q),Ti)
 assert matmul(T,Ti)==eye(n-1)
 assert all(H[i][j]>=0 for i in range(n-1) for j in range(n-1) if i!=j)
 return G,H

if __name__=='__main__':
 out=[]
 for family in ['linear','matched_exponential']:
  for E in [[0,1,2],[0,1,2,3],[0,1,2,3,4]]:
   G,H=model(E,family=family)
   out.append({'family':family,'energies':E,'q':'1/2','G':[[str(x) for x in r] for r in G], 'H':[[str(x) for x in r] for r in H], 'characteristic_coefficients':[str(x) for x in characteristic(G)], 'trace':str(sum(G[i][i] for i in range(len(E)))),'metzler':True})
  
 enc=lambda a:json.dumps(a,separators=(',',':')).encode()
 for r in out:
  a=r['G'];n=len(a);trans=list(map(list,zip(*a)));rev=[list(reversed(x)) for x in reversed(a)]
  r['dimensions']=[n,n];r['rank']=n-1;r['determinant']='0'
  r['row_sums']=[str(sum(map(F,x))) for x in a]
  r['column_sums']=[str(sum(F(a[i][j]) for i in range(n))) for j in range(n)]
  r['sorted_entry_multiset']=[str(x) for x in sorted(F(x) for row in a for x in row)]
  r['hashes']={'matrix':hashlib.sha256(enc(a)).hexdigest(),'transpose':hashlib.sha256(enc(trans)).hexdigest(),'reversal':hashlib.sha256(enc(rev)).hexdigest(),'minimum_over_simultaneous_permutations':min(hashlib.sha256(enc([[a[i][j] for j in perm] for i in perm])).hexdigest() for perm in permutations(range(n)))}
  r['zero_pattern']=[[int(F(x)==0) for x in row] for row in a]
  r['positivity']='Strict positive off-diagonal rates, negative diagonal entries'
  r['nonzero_spectrum']=['(-14+sqrt(7))/3','(-14-sqrt(7))/3'] if r['family']=='linear' and n==3 else 'Roots of the recorded characteristic polynomial after removing the factor lambda'
 path=Path(__file__).with_name('small_fingerprints.json');path.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
