import numpy as np,json
from scipy.linalg import eigh,expm
from pathlib import Path
rng=np.random.default_rng(20261026)
records=[];crossings=[];min_order=0
for trial in range(400):
 n=int(rng.integers(3,11));E=np.r_[0,np.cumsum(np.exp(rng.uniform(-2,1,n-1)))];E*=rng.uniform(.5,12)/E[-1]
 alpha=rng.choice([0,.2,.5,.8,1]);D=abs(E[:,None]-E[None,:]);np.fill_diagonal(D,1)
 J=D**alpha;h=J/np.expm1(D);G=h+J*(E[:,None]>E[None,:]);np.fill_diagonal(G,0);G-=np.diag(G.sum(axis=1))
 pi=np.exp(-E);pi/=pi.sum();R=np.sqrt(pi);S=R[:,None]*G/R[None,:];w,U=eigh(S);gap=-w[-2]
 bh,bc=sorted(rng.uniform(0,1,2));ph=np.exp(-bh*E);ph/=ph.sum();pc=np.exp(-bc*E);pc/=pc.sum()
 for t in np.r_[0,np.logspace(-4,1.5,45)/gap]:
  P=(U*np.exp(w*t))@U.T;qh=(ph/R)@P*R;qc=(pc/R)@P*R
  if t==0:continue
  dh=.5*abs(qh-pi).sum();dc=.5*abs(qc-pi).sum()
  order=np.cumsum(qc-qh)[:-1].min();min_order=min(min_order,float(order))
  if dc-dh>1e-10:
   crossings.append({'trial':trial,'N':n,'E':E.tolist(),'alpha':float(alpha),'beta_h':float(bh),'beta_c':float(bc),'time':float(t),'hot_tv':float(dh),'cold_tv':float(dc),'min_cdf_order':float(order)})
   break
 if crossings:break
out={'trials':trial+1,'TV_crossings':crossings,'min_stochastic_order_gap':min_order,'scope':'numerical probe, not theorem'}
Path(__file__).with_name('transient_probe.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
