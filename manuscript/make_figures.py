#!/usr/bin/env python3
"""Numerical illustrations of proved thermal ordering; no data are experimental."""
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/thermal_order_mpl")
import json
import numpy as np
import scipy.linalg as la
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":9,
                    "axes.spines.top":False,"axes.spines.right":False,
                    "savefig.facecolor":"white"})
PLUM, TEAL, OCHRE = "#753F68", "#22766E", "#B06B2E"
E = np.array([0.,1.,2.,4.,7.])
BB = np.log(2.)
def gibbs(beta):
    p = np.exp(-beta*E)
    return p/p.sum()
PI = gibbs(BB)
def build(kind):
    G=np.zeros((len(E),len(E)))
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]
            w=d if kind=="linear" else 1-np.exp(-BB*d)
            up=w/np.expm1(BB*d)
            G[i,j]=up;G[j,i]=up+w
    G[np.diag_indices(len(E))]=-G.sum(axis=1)
    G/=max(-np.diag(G))
    assert np.max(np.abs(PI[:,None]*G-(PI[:,None]*G).T))<1e-13
    return G
G1,G2=build("linear"),build("metropolis")
def kl(p):
    r=p/PI
    return float(np.sum(PI*(r*np.log(r)-r+1)))
ratios=np.linspace(0,1,161)
betas=BB*ratios
static={t:np.array([kl(la.expm(G1.T*t)@gibbs(b)) for b in betas]) for t in [0.,1.,3.]}
for values in static.values():assert np.max(np.diff(values))<2e-13

stage_duration=.8
switches=np.arange(0,6.4+1e-10,stage_duration)
times=np.linspace(0,6.4,257)
prep_ratios=[.1,.4,.75]
def density_at(t,beta):
    p=gibbs(beta)
    for k,start in enumerate(switches[:-1]):
        duration=min(stage_duration,max(0,t-start))
        if duration<=0:break
        p=la.expm((G1 if k%2==0 else G2).T*duration)@p
    assert abs(p.sum()-1)<3e-14 and p.min()>0
    return p
scheduled=np.array([[kl(density_at(t,BB*r)) for r in prep_ratios] for t in times])
assert np.min(scheduled[:,0]-scheduled[:,1])>0
assert np.min(scheduled[:,1]-scheduled[:,2])>0
fig,ax=plt.subplots(1,2,figsize=(7.2,3.1),layout="constrained")
for t,c in zip(static,[PLUM,TEAL,OCHRE]):ax[0].plot(ratios,static[t],color=c,label=f"t = {t:g}")
ax[0].set(xlabel=r"Preparation inverse temperature, $\beta/\beta_b$",ylabel="Relative entropy to equilibrium (nats)",title="(a) Fixed linear jump weight")
ax[0].legend(frameon=False)
for j,(r,c) in enumerate(zip(prep_ratios,[PLUM,TEAL,OCHRE])):ax[1].semilogy(times,scheduled[:,j],color=c,label=rf"$\beta/\beta_b={r:g}$")
for s in switches[1:-1]:ax[1].axvline(s,color="#C7C4BD",lw=.6,zorder=0)
ax[1].set(xlabel="Time (common rate units)",ylabel="Relative entropy to equilibrium (nats)",title="(b) Alternating admissible generators")
ax[1].legend(frameon=False)
fig.savefig(OUT/"thermal_order.png",dpi=450)
fig.savefig(OUT/"thermal_order.svg")
plt.close(fig)

m=5;r=np.arange(m+1);F=1-r/(m-.5)
fig,ax=plt.subplots(figsize=(5.6,2.7),layout="constrained")
ax.axhline(0,color="#5A5A56",lw=.8)
ax.plot(r,F,color="#A8AAA5",lw=1)
ax.scatter(r[:-1],F[:-1],color=TEAL,s=45,zorder=3,label="Proper retained systems")
ax.scatter(r[-1:],F[-1:],color=PLUM,s=65,marker="D",zorder=3,label="Full eight-level system")
ax.set(xticks=r,xlabel="Retained middle levels, r",ylabel=r"Exact forcing margin $F(r)/\varepsilon$",ylim=(-.25,1.15))
ax.legend(frameon=False,loc="upper right")
fig.savefig(OUT/"hierarchy_margin.png",dpi=450)
fig.savefig(OUT/"hierarchy_margin.svg")
plt.close(fig)

np.savetxt(OUT/"thermal_order_curves.csv",np.column_stack([times,scheduled]),delimiter=",",header="time,KL_beta_ratio_0.1,KL_beta_ratio_0.4,KL_beta_ratio_0.75",comments="")
metadata={"energies":E.tolist(),"beta_b":float(BB),"generator_normalisation":"Each generator divided by its maximum escape rate; the resulting two matrices define fixed common time units.","G_linear":G1.tolist(),"G_metropolis":G2.tolist(),"commutator_Frobenius_norm":float(la.norm(G1@G2-G2@G1)),"stage_duration":stage_duration,"preparation_beta_ratios":prep_ratios,"minimum_scheduled_ordering_gap":float(np.min(scheduled[:,:-1]-scheduled[:,1:])),"interpretation":"Numerical illustration of theorem predictions; not an empirical dataset, interval proof, or claim of scheduling optimality."}
(OUT/"figure_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
print(json.dumps(metadata,indent=2))
