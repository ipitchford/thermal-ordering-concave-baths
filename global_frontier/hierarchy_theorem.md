# Arbitrarily large minimal Mpemba systems

Research note, 26 October 2026. This is a new candidate theorem within the present session. Its proof and scope are set out for independent audit. The novelty judgement remains bounded by `novelty_report.md`; no publication or REF rating is asserted.

## 1. Statement and meaning

**Theorem.** For every integer N >= 4 there is a strictly ordered energy spectrum and a complete, detailed-balance N-state generator at a finite positive bath temperature, with a simple spectrum, for which:

1. The full system has a strict slow-mode cooling Mpemba effect and an actual finite-time total-variation crossing between two finite Gibbs preparation temperatures.
2. Every proper induced subsystem of at least two states has a strictly decreasing signed slow Gibbs amplitude on all nonnegative preparation inverse temperatures. Consequently, every proper subsystem excludes a slow-mode Mpemba effect between any two preparations on the same side of its bath temperature, and has no nontrivial strong-mode zero there.
3. These statements survive sufficiently small detailed-balance-preserving rate perturbations.

At zero bath temperature the same construction has a stronger local property: every proper subsystem has conventional total-variation ordering at every time for every strictly ordered energy spectrum. The full system violates that ordering for some spectrum.

An induced subsystem retains internal off-diagonal rates and recomputes its diagonal escape rates. The statement does not concern marginal dynamics or a Schur-complement effective generator.

**Consequence.** There is no dimension bound k such that every slow-mode thermal Mpemba system has an induced anomalous subsystem of at most k states. This is a failure of hereditary anomaly witnessing. It is not an information-theoretic impossibility: all pairwise rate restrictions jointly reveal the whole generator.

## 2. The rates

Put m=N-3. The ordered states are ground 1, slow state 2, m middle states 3,...,N-1 and top state N. Choose

\[
0<\varepsilon<1/14,\qquad
\delta=\varepsilon,\qquad
C=2-\varepsilon,\qquad
z=1-2\varepsilon-\frac{3\varepsilon}{m-1/2}.
\]

Every downward rate is strictly positive. In the convention a_ij is j -> i, set

\[
a_{12}=1,
\quad a_{1j}=3/2,
\quad a_{2j}=1\quad(3\le j<N),
\quad a_{ij}=\delta\quad(3\le i<j<N),
\]

and

\[
a_{1N}=C,\qquad a_{2N}=1,\qquad a_{jN}=z\quad(3\le j<N).
\]

At zero bath temperature upward rates vanish. Three inequalities used below hold uniformly for m>=1:

\[
z\ge1-8\varepsilon>\varepsilon=\delta,
\qquad
z-C+3/2\ge 1/2-7\varepsilon>0,
\qquad
2z-1/2+\varepsilon\ge3/2-15\varepsilon>0.
\]

The full escape rates are 0,1,5/2,5/2+delta,...,5/2+(m-1)delta,3-epsilon+m z. They are strictly ordered. Every subsystem also has a simple spectrum, by the four cases below: its middle escape rates are spaced by delta, its first excited escape rate is 1, 3/2 or delta, and its top escape rate exceeds all retained middle escape rates because z>delta and C>3/2.

## 3. The Gibbs comparison lemma used

If a finite observable b has nondecreasing arithmetic prefix means in increasing energy order, its Gibbs expectation

\[
F(\beta)=\frac{\sum_i b_i e^{-\beta E_i}}{\sum_i e^{-\beta E_i}}
\]

decreases on beta>=0 for every ordered energy vector. It decreases strictly if b is nonconstant. This is the classical increasing-in-mean/Chebyshev ordering argument already proved and attributed in Lemma 2 of the reviewed manuscript. A negative final prefix margin can be converted into a positive derivative by choosing the final energy sufficiently large.

At zero bath temperature, total variation from equilibrium is survival probability. If h_i(t) is the probability of not yet reaching the subsystem's ground state at time t, then

\[
D_{\rm TV}(e^{Mt}p_E(\beta),e_{\rm ground})=\sum_i p_{E,i}(\beta)h_i(t).
\]

Thus increasing prefix means of the survival vector h(t) imply conventional total-variation ordering at that time for every ordered energy spectrum.

## 4. Every proper subsystem: all four classes

Only the retained count of middle states matters for the survival responses; which middle labels were retained does not change the reduced rates after reindexing.

### A. Subsystems retaining both 1 and 2

Let r be the number of retained middle states. If the top is retained, properness means r<=m-1. Without the top there is no such restriction.

The ground response is zero, and

\[
h_2(t)=e^{-t},\qquad
h_M(t)=\frac23e^{-t}+\frac13e^{-5t/2}
\]

for every middle state. Middle-to-middle transitions cancel in this equation because all middle survival responses agree. This assertion follows by induction in increasing middle order, or by substitution into the backward equations with the shared initial value 1.

The middle prefix comparisons are positive:

\[
2h_M-h_2=\frac13e^{-t}+\frac23e^{-5t/2}>0.
\]

If the top is present, put Gamma=C+1+r z and

\[
d(t)=h_N(t)-\frac{h_2(t)+r h_M(t)}{r+2},
\quad
\kappa=\frac{r(2z-1/2+\varepsilon)}{r+2}\ge0.
\]

The top backward equation yields

\[
d'=-\Gamma d+
\frac{F(r)}{r+2}e^{-t}+
\frac{\kappa}{3}e^{-5t/2},
\qquad d(0)=\frac1{r+2},
\]

where the exact linear obstruction is

\[
F(r)=\varepsilon+r\left(\frac{z-1}{3}+\frac{2\varepsilon}{3}\right)
=\varepsilon\left(1-\frac{r}{m-1/2}\right).
\]

For every proper subsystem F(r)>0. Variation of constants therefore gives d(t)>0 for every finite t>=0. All survival prefix means increase strictly.

The normalised slow response, with ground value 0 and state-2 value 1, is

\[
g=(0,1,\underbrace{2/3,\ldots,2/3}_{r},\eta_r),
\quad
\eta_r=\frac{1+(2/3)r z}{C+r z}.
\]

Its last prefix margin is

\[
\eta_r-\frac{1+(2/3)r}{r+2}
=\frac{F(r)}{(r+2)(C+r z)}>0.
\]

The same response without the final coordinate covers subsystems without the top.

### B. Subsystems retaining 1 but not 2

Let r be the number of middle states. If r=0, the subsystem containing the top is two-state and conventional. For r>=1, every middle survival response is e^{-3t/2}. The top escape rate is Gamma=C+r z>3/2, and

\[
h_N(t)=b e^{-3t/2}+(1-b)e^{-\Gamma t},
\quad b=\frac{r z}{C+r z-3/2}.
\]

The final prefix comparison is

\[
d(t)=e^{-3t/2}\left[b-\frac{r}{r+1}+(1-b)e^{-(\Gamma-3/2)t}\right].
\]

The bracket is affine in a variable ranging over [0,1]. Its endpoint values are 1/(r+1)>0 and

\[
b-\frac{r}{r+1}
=\frac{r(z-C+3/2)}{(r+1)(C+r z-3/2)}>0.
\]

All prefix means increase. The slow response is (0,1,...,1,b), with the same strict comparisons. Without the top, the conclusion is immediate.

### C. Subsystems retaining 2 but not 1

State 2 becomes the ground. Every other retained state has direct absorption rate 1 into it. Hence every excited survival response is exactly e^{-t}, whatever the internal transitions. The slow response is (0,1,...,1), with strict prefix comparisons. Two-state cases are included.

### D. Subsystems retaining neither 1 nor 2

If two or more middle states remain, the least retained middle state becomes ground and every other middle state has direct absorption rate delta. If r is the number of retained middle states, their non-ground survival responses are e^{-delta t}. When the top is retained its escape rate is Gamma=r z>delta, and

\[
h_N(t)=b e^{-\delta t}+(1-b)e^{-r zt},
\quad b=\frac{(r-1)z}{r z-\delta}.
\]

The final prefix comparison has positive endpoint values 1/r and

\[
b-\frac{r-1}{r}
=\frac{(r-1)\delta}{r(r z-\delta)}>0.
\]

This proves the claim as in case B. With one middle state and the top the subsystem is two-state; without the top all excited survival responses agree. Singletons have no nonzero relaxation mode and are trivially conventional.

These four cases exhaust all induced subsystems. They prove the all-time zero-bath statement, and their slow responses have strictly increasing prefix means.

## 5. The full system fails precisely at the last comparison

For r=m,

\[
F(m)=-\frac{\varepsilon}{2m-1}<0.
\]

Thus the top slow response lies strictly below the preceding prefix mean. All earlier slow prefix comparisons remain positive.

For a fully explicit choice take epsilon=delta=1/100 and put

\[
P_m=196m^2+294m-199,
\quad
z=\frac{49m-26}{25(2m-1)},
\quad
\eta_m=\frac{2m+3}{3(m+2)}-\frac1{(m+2)P_m}.
\]

Choose the shared spectrum

\[
E=(0,1,2,\ldots,m+1,L),
\quad
L>\frac{m+1}{2}+\frac{P_m(m+3)^2}{6}.
\]

The full Gibbs slow amplitude A(beta)=p_E(beta)^T g has derivative

\[
A'(0)=
\frac{L-(m+1)/2-P_m(m+3)^2/6}{P_m(m+3)^2}>0.
\]

For clarity, this follows by writing n=m+3,

\[
\sum_{i<N}E_i=(m+1)(m+2)/2,
\qquad
\sum_{i<N}g_iE_i=((m+1)(m+2)+1)/3,
\]

and differentiating the normalised Gibbs expectation at beta=0. Continuity gives two finite 0<beta_h<beta_c with 0<A(beta_h)<A(beta_c). Initially the hotter Gibbs state has less ground probability and larger total variation from the absorbing ground. At late times its smaller slow amplitude gives a smaller total variation. A finite-time crossing follows.

This is ordinary amplitude inversion. Both amplitudes are nonzero, so no exponential improvement of the asymptotic decay rate is claimed.

## 6. Simultaneous positive-bath continuation

Fix one N, the above energy spectrum and two finite inverse preparation temperatures inside the increasing interval. Keep every downward rate fixed. For q=e^{-beta_b}>0 set

\[
a_{ji}(q)=a_{ij}q^{E_j-E_i}\quad(i<j).
\]

Every full-system rate is strictly positive and the common bath Gibbs law satisfies detailed balance. Each induced subsystem satisfies detailed balance with the restricted, renormalised Gibbs law.

At q=0 every subsystem has a simple zero eigenvalue and a simple slow eigenvalue, with a positive gap to the next mode. All eigenvalues are in fact distinct by the rate comparisons in Section 2. Standard continuity of simple eigenpairs gives a continuous real slow eigenfunction for each subsystem. Normalise it by its lowest two components, whose difference is nonzero at q=0. Subtracting a constant does not change a prefix margin.

Every proper subsystem has finitely many strictly positive prefix margins at q=0. There are finitely many proper subsystems for fixed N. Therefore a common q_0>0 exists such that all these margins remain positive for 0<q<q_0, and all full/subsystem spectra remain simple. The Gibbs comparison lemma now makes each centred proper-subsystem slow amplitude strictly decreasing for every beta>=0. It vanishes exactly at beta_b, so its absolute value has conventional ordering separately on either side of the bath.

The full-system derivative and chosen positive amplitude ordering also persist. Taking q still smaller makes beta_b>beta_c. The initial total-variation inequality persists from q=0. The simple slow mode then gives a late-time reversal, hence a finite-time crossing at positive bath temperature.

The same finitely many strict margins, spectral gaps and full amplitude comparisons persist in a sufficiently small open neighbourhood of the fixed downward rates, with upward rates assigned by detailed balance. This proves the claimed robustness.

The argument does not give a uniform q_0 valid simultaneously for all N. It also does not claim that all-time proper-subsystem total-variation ordering survives at finite bath temperature: only their global signed slow-amplitude monotonicity has been continued.

## 7. Limits and physical interpretation

The construction proves that anomaly witnessing by induced state subsets has unbounded minimal order. It strengthens a single triplet counterexample into an arbitrary-order structural obstruction. It does not establish a physically constrained Hamiltonian/bath realisation or an experimentally substantial speedup.

With epsilon fixed at 1/100, the full prefix deficit is 1/[(m+2)P_m]=O(m^-3); the displayed sufficient top-energy scale grows as O(m^4). These are properties of this construction, not proved optimal lower bounds. They limit an immediate practical interpretation and should be stated rather than hidden.

Every induced triplet also strictly fails the original necessary-rate screen at zero bath temperature: its B2 margin is one of 1/2, epsilon, delta, z, 3/2-C+z, or another positive constant from the listed rates. Strictness continues to sufficiently small positive bath temperatures. This is consistent with, and extends, the original four-level screen failure.

The theorem does not preclude a valid global criterion computed from local rates. In particular, a list of all pairwise restrictions already determines the generator.

### Lumpability and the size of the dynamical core

At zero bath temperature the partition {1}, {2}, {all retained middles}, {top} is strongly lumpable. Its column generator is

\[
\overline M_r=
\begin{pmatrix}
0&1&3/2&C\\
0&-1&1&1\\
0&0&-5/2&r z\\
0&0&0&-(C+1+r z)
\end{pmatrix}.
\]

The middle block's Gibbs preparation weight is the sum of its individual Boltzmann weights. Removing middle levels changes both its multiplicity and the aggregate top-to-middle rate. The hierarchy therefore has an arbitrarily large minimal *induced-level* witness but a four-block survival dynamics. It does not exhibit an increasing number of independently necessary relaxation modes. Nondegenerate middle energies generally break this exact lumpability at positive bath temperature. The novelty gate was rerun under lumpability, multiplicity and entropic/degeneracy aliases after this simplification.

## 8. Finite computations

`verify_hierarchy.py` uses only exact rational arithmetic. It reconstructs slow eigenvectors by triangular back-substitution for all 8,086 proper subsystems of N=4,...,12 at epsilon=1/100. Every proper prefix margin is positive and every full final prefix margin is negative. It checks the positive coefficients/endpoints in the all-time survival argument, rather than relying on time samples. Another 1,368 proper subsystem cases cover epsilon=1/50, 1/20 and 1/15 for N=4,...,8. The general proof remains the source of the unrestricted theorem.

`check_positive_bath.py` independently illustrates the finite-temperature continuation with 80-digit Decimal arithmetic, beta_b=20, fixed downward rates and all strictly positive upward rates. For N=4,...,8 it checks every proper subsystem's slow prefix margins, confirms the full amplitude inversion between beta_h=1/(100L) and beta_c=2/(100L), and finds a positive time-20 TV reversal through Poisson uniformization. The chosen integer L is slightly above twice the analytic threshold. The N=4 and N=8 amplitude gaps are approximately 1.05e-6 and 1.23e-8 respectively. These are numerical illustrations, not interval certificates. The reverse rates in these examples are extremely small, with the N=8 minimum rate of order 10^(-1143540). Their positivity is retained by arbitrary-exponent Decimal arithmetic, and no experimental feasibility is inferred.
