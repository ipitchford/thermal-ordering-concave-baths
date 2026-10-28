from pathlib import Path
import numpy as np, scipy.linalg as la, json, time
rng=np.random.default_rng(261027)
root=Path(__file__).resolve().parent
start=time.monotonic(); best=None; checked=0; bad=0
for case in range(1800):
 n=int(rng.integers(3,13));gaps=10**rng.uniform(-3,0,n-1); E=np.r_[0,np.cumsum(gaps)];E/=E[-1]
 beta=10**rng.uniform(-1,np.log10(20)); pi=np.exp(-beta*E);pi/=pi.sum()
 knots=np.sort(rng.uniform(.01,.99,5)); weights=10**rng.uniform(-2,2,5); ends=10**rng.uniform(-4,1,2)
 def J(d):return ends[0]*(1-d)+ends[1]*d+np.sum(weights*(np.minimum(d,knots)-knots*d))
 G=np.zeros((n,n))
 for i in range(n):
  for j in range(i+1,n):
   d=E[j]-E[i]; h=J(d)/np.expm1(beta*d);G[i,j]=h;G[j,i]=h+J(d)
 G[np.diag_indices(n)]=-G.sum(axis=1);G/=np.max(-np.diag(G))
 S=np.sqrt(pi)[:,None]*G/np.sqrt(pi)[None,:]
 assert np.max(np.abs(S-S.T))<2e-10
 vals,U=la.eigh((S+S.T)/2);vals[-1]=0
 gamma=-vals[-2]; times=np.r_[0,np.geomspace(.001,35/gamma,70)]
 for direction in ['cooling','heating']:
  if direction=='cooling':
   near=beta*rng.uniform(.03,.99);far=near*rng.uniform(0,.95)
  else:
   near=beta+10**rng.uniform(-2,1);far=near+10**rng.uniform(-2,2)
  p=np.exp(-E[:,None]*np.array([far,near])[None,:]);p/=p.sum(axis=0)
  coeff=U.T@(p/np.sqrt(pi)[:,None]);result=np.einsum('ij,tj,jk->tik',np.sqrt(pi)[:,None]*U,np.exp(times[:,None]*vals[None,:]),coeff)
  err=max(abs(result.sum(axis=1)-1).max(),max(0,-result.min()))
  if err>1e-7:bad+=1;continue
  tv=.5*np.abs(result-pi[None,:,None]).sum(axis=1);delta=tv[:,1]-tv[:,0];idx=np.argmax(delta)
  assert delta[0]<1e-10
  checked+=1
  if best is None or delta[idx]>best['gap']:
   best={'gap':float(delta[idx]),'time':float(times[idx]),'case':case,'direction':direction,'n':n,'E':E.tolist(),'beta_b':float(beta),'beta_far':float(far),'beta_near':float(near),'knots':knots.tolist(),'weights':weights.tolist(),'ends':ends.tolist(),'G':G.tolist(),'tv_far':float(tv[idx,0]),'tv_near':float(tv[idx,1]),'normalisation_error':float(err)}
  if delta[idx]>1e-7 and tv[idx,1]>1e-6:
   # Independently check promising point by a direct matrix exponential.
   exact=la.expm(G.T*times[idx])@p
   dt=.5*np.abs(exact[:,1]-pi).sum()-.5*np.abs(exact[:,0]-pi).sum()
   if dt>1e-7:
    best['direct_expm_gap']=float(dt);print('COUNTEREXAMPLE',json.dumps(best),flush=True)
    (root/'result.json').write_text(json.dumps({'status':'counterexample_candidate','checked':checked,'bad':bad,'best':best},indent=2));raise SystemExit
 if case%200==0:print('completed',case,'checked',checked,'largest gap',best['gap'],flush=True)
result={'status':'no_crossing_found','checked':checked,'discarded_numerically_unstable':bad,'best':best,'seconds':time.monotonic()-start,'scope':'1800 sampled concave piecewise-linear kernels, 3-12 levels, one cooling and one heating pair each, 71 sampled times; not a proof or interval certificate.'}
(root/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
