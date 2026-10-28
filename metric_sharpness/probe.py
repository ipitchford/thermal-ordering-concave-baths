"""Exploratory LP search; numerical witnesses are not interval certificates."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import eigh, expm
from scipy.optimize import linprog

OUT = Path(__file__).resolve().parent
RNG = np.random.default_rng(20261026)


def transforms(n):
    t = np.eye(n-1)
    for j in range(1, n-1):
        t[j, :j] = -1/(j+1)
    return t, np.linalg.inv(t)


def build(pi, edges, c):
    n = len(pi)
    g = np.zeros((n,n))
    for (i,j),v in zip(edges,c):
        g[i,j] = v/pi[i]
        g[j,i] = v/pi[j]
    g[np.diag_indices(n)] = -g.sum(axis=1)
    return g


def geometry(pi):
    n = len(pi)
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    t,ti = transforms(n)
    off = ~np.eye(n-1,dtype=bool)
    amap = []
    for k in range(len(edges)):
        c = np.zeros(len(edges)); c[k] = 1
        g = build(pi,edges,c)
        h = t @ (g[1:,1:]-g[0,1:]) @ ti
        amap.append(h[off])
    return edges,np.array(amap).T,t,ti


def thermal(pi,beta):
    energy = -np.log(pi)
    logp = np.log(pi)[None,:]*beta[:,None]
    logp -= logp.max(axis=1)[:,None]
    p = np.exp(logp); p /= p.sum(axis=1)[:,None]
    dp = -p*(energy[None,:]-(p@energy)[:,None])
    return p,dp


def evaluate(pi,g,case):
    root = np.sqrt(pi)
    g = g/max(-g.diagonal())
    ev,u = eigh(root[:,None]*g/root[None,:])
    gap = -ev[-2]
    if gap < 1e-9:
        return None
    times = np.unique(np.r_[0,np.geomspace(.0003, min(30/gap, 30000),100)])
    betac = np.unique(np.r_[np.linspace(0,1,26)[:-1],.99,.999])
    betah = 1+np.geomspace(.001,100,45)
    beta = np.r_[betac,betah]
    p,dp = thermal(pi,beta)
    pe = (p/root) @ u
    dpe = (dp/root) @ u
    # T x B x N arrays
    decay = np.exp(times[:,None]*ev[None,:])
    pt = np.einsum('bk,tk,ik->tbi',pe,decay,u)*root
    dpt = np.einsum('bk,tk,ik->tbi',dpe,decay,u)*root
    rt = np.maximum(pt/pi,1e-250)
    # Remove round-off in the equilibrium mass component.
    dpt -= dpt.sum(axis=2)[:,:,None]*pi
    # A positive badscore indicates the wrong temperature-derivative sign.
    specs=[('cooling',3),('cooling',4),('cooling',10),('cooling',100),
           ('heating',.1),('heating',.5),('heating',1),('heating',1.5),
           ('cooling','TV'),('heating','TV'),('cooling','hinge'),('heating','hinge')]
    best=None
    for side,alpha in specs:
        sl = slice(0,len(betac)) if side=='cooling' else slice(len(betac),len(beta))
        r = rt[:,sl,:]; v = dpt[:,sl,:]
        threshold=None
        if alpha=='TV':
            deriv = .5*np.sum(np.sign(r-1)*v,axis=2)
        elif alpha=='hinge':
            # Any positive part f_k(x)=(x-k)_+ is a convex divergence generator,
            # up to the irrelevant subtraction of f_k(1). Check every rank cut.
            order = np.argsort(r,axis=2)
            sortedv = np.take_along_axis(v,order,axis=2)
            tail = np.cumsum(sortedv[:,:,::-1],axis=2)[:,:,::-1]
            sr = np.take_along_axis(r,order,axis=2)
            valid = (sr[:,:,1:]-sr[:,:,:-1]) > 1e-8
            scoreall = tail[:,:,1:]*(1 if side=='cooling' else -1)
            scoreall = np.where(valid,scoreall,-np.inf)
            idx3=np.unravel_index(np.argmax(scoreall),scoreall.shape)
            deriv = np.zeros(scoreall.shape[:2]); deriv[idx3[:2]]=scoreall[idx3]*(1 if side=='cooling' else -1)
            threshold = float((sr[idx3[0],idx3[1],idx3[2]]+sr[idx3[0],idx3[1],idx3[2]+1])/2)
        elif alpha==1:
            deriv=np.sum(np.log(r)*v,axis=2)
        else:
            # Positive rescaling preserves sign and prevents overflow for large alpha.
            logphi=(alpha-1)*np.log(r)
            shift=logphi.max(axis=2,keepdims=True)
            phi=np.exp(logphi-shift)/(alpha-1)
            deriv=np.sum(phi*v,axis=2)
        score=deriv*(1 if side=='cooling' else -1)
        # Suppress numerical equilibrium floor and negligible initial derivatives.
        dist=np.sum(abs(pt[:,sl,:]-pi),axis=2)
        score=np.where(dist>1e-7,score,-np.inf)
        idx=np.unravel_index(np.argmax(score),score.shape)
        val=float(score[idx])
        if val>1e-8 and (best is None or val>best['score']):
            bi=idx[1] if side=='cooling' else idx[1]+len(betac)
            best=dict(case=case,side=side,alpha=alpha,score=val,time=float(times[idx[0]]),
                      beta=float(beta[bi]),pi=pi.tolist(),G=g.tolist(),
                      eigenvalues=ev.tolist(),r=rt[idx[0],bi].tolist(),
                      p_derivative=dpt[idx[0],bi].tolist(),threshold=threshold)
    return best


def run(count=1200):
    start=time.time()
    hits=[]; tested=0; lp_failed=0
    for case in range(count):
        n=3+case%6
        if case%4==0:
            weights=2.**(-np.arange(n))
        elif case%4==1:
            weights=np.exp(-np.cumsum(np.r_[0,RNG.uniform(.01,2,n-1)]))
        elif case%4==2:
            weights=np.exp(-np.cumsum(np.r_[0,RNG.uniform(.001,.2,n-1)]))
        else:
            weights=np.exp(-np.cumsum(np.r_[0,RNG.uniform(.1,4,n-1)]))
        pi=weights/weights.sum()
        edges,amap,t,ti=geometry(pi)
        objective=RNG.normal(size=len(edges))
        objective*=np.exp(RNG.uniform(-2,2,len(edges)))
        opt=linprog(objective,A_ub=-amap,b_ub=np.zeros(len(amap)),
                    A_eq=np.ones((1,len(edges))),b_eq=[1],bounds=(0,None),method='highs')
        if not opt.success:
            lp_failed+=1;continue
        c=opt.x
        reset=np.array([pi[i]*pi[j] for i,j in edges]);reset/=reset.sum()
        mix=10**RNG.uniform(-5,-.3)
        c=(1-mix)*c+mix*reset
        g=build(pi,edges,c)
        h=t@(g[1:,1:]-g[0,1:])@ti
        if min(h[~np.eye(n-1,dtype=bool)]) < -1e-5:
            raise RuntimeError('LP infeasible numerical result')
        tested+=1
        hit=evaluate(pi,g,case)
        if hit:
            hit['conductances']=c.tolist();hit['edges']=edges;hit['H']=h.tolist()
            hits.append(hit)
            print('HIT',case,n,hit['side'],hit['alpha'],hit['score'],flush=True)
            (OUT/'probe_hits.json').write_text(json.dumps(hits,indent=2))
            if hit['alpha']!='hinge':
                break
        if case%50==0:
            print('progress',case,'hits',len(hits),'seconds',round(time.time()-start,1),flush=True)
    out=dict(seed=20261026,cases_requested=count,tested=tested,lp_failed=lp_failed,
             elapsed_seconds=time.time()-start,hits=hits,
             limitation='Exploratory floating-point scan; negative result is not an exclusion theorem.')
    (OUT/'probe_result.json').write_text(json.dumps(out,indent=2))
    print('done',tested,'hits',len(hits),flush=True)


if __name__=='__main__':
    run()
