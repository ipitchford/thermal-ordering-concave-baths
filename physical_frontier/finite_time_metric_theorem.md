# All-time divergence ordering under a cone-preserving thermal protocol

Research proof note, 26 October 2026. The coordinating agent proposed this cone/Jensen use of divergence differentiation and the protocol extension; this branch verified and expanded them. This note extends *concave_kernel_theorem.md*. The KL preparation-derivative/covariance identity itself is already present in Ohga, Hayakawa and Ito, Appendix B, Result S2 and Lemma S3, Eqs. (S14)–(S20), and is not claimed as new. Novelty and external significance require separate assessment. The statement below is not a theorem for every distance to equilibrium.

## 1. General statement

Fix \(N\ge2\) strictly ordered energies \(E_1<\cdots<E_N\), a finite positive bath inverse temperature \(\beta_b\), and its full-support Gibbs distribution
\[
\pi_i=p_{\beta_b,i},\qquad
p_{\beta,i}=\frac{e^{-\beta E_i}}{\sum_j e^{-\beta E_j}}.
\]
Preparation inverse temperatures are finite and nonnegative. Let
\[
\mathcal C=\left\{b\in\mathbb R^N:
b_j\ge\frac1{j-1}\sum_{i<j}b_i,\quad j=2,\ldots,N\right\}.
\]
Constants belong to both \(\mathcal C\) and \(-\mathcal C\), and the intersection consists exactly of constants.

Write the weighted inner product as
\(\langle u,v\rangle_\pi=\sum_i\pi_i u_i v_i\).
Let \(K\) be a Markov operator on likelihood-ratio columns: \(K\) has nonnegative entries, \(K\mathbf1=\mathbf1\), and \(\pi^T K=\pi^T\). Its \(\pi\)-adjoint is
\[
K^*_\pi=\operatorname{diag}(\pi)^{-1}K^T\operatorname{diag}(\pi).
\]
Assume the two cone conditions
\[
K\mathcal C\subseteq\mathcal C,\qquad
K^*_\pi\mathcal C\subseteq\mathcal C. \tag{1}
\]
For \(r_\beta=p_\beta/\pi\), define the evolved distribution \(p_\beta^K\) by \(p_{\beta,i}^K=\pi_i(Kr_\beta)_i\). It is positive and normalised.

For a convex \(C^3\) function \(f:(0,\infty)\to\mathbb R\) with \(f(1)=0\), use the convention
\[
D_f(p\Vert\pi)=\sum_i\pi_i f(p_i/\pi_i). \tag{2}
\]
Adding an affine multiple of \(x-1\) to \(f\) changes neither this divergence nor the result.

**Theorem C (finite-protocol divergence ordering).**

1. If \(f''\ge0\) and \(f'''\le0\), then
   \[
   \frac{\partial}{\partial\beta}D_f(p_\beta^K\Vert\pi)\le0,
   \qquad 0\le\beta\le\beta_b. \tag{3}
   \]
   Thus every pair of cooling Gibbs preparations retains its conventional distance ordering after this same protocol.
2. If \(f''\ge0\) and \(f'''\ge0\), then
   \[
   \frac{\partial}{\partial\beta}D_f(p_\beta^K\Vert\pi)\ge0,
   \qquad \beta\ge\beta_b. \tag{4}
   \]
   Thus every pair of heating Gibbs preparations retains its conventional distance ordering after the same protocol.
3. If \(K\) is invertible and \(f'\) is strictly increasing, then these inequalities are strict at every \(\beta\ne\beta_b\) in the stated ranges. They vanish at \(\beta=\beta_b\).

The hypotheses can equivalently be stated using an increasing concave derivative \(f'\) in part 1 or an increasing convex derivative in part 2; the third derivative is a convenient sufficient smoothness formulation. The two signs of \(f'''\) distinguish the cooling and heating statements. They cannot be silently interchanged.

### 1.1 Homogeneous dynamics

If \(G\) is a row-generator on observables, reversible with respect to \(\pi\), and \(e^{tG}\mathcal C\subseteq\mathcal C\), take \(K=e^{tG}\). Reversibility gives \(K^*_\pi=K\), and this exponential is invertible at every finite \(t\). Consequently Theorem C applies at **every finite time**, including \(t=0\). At equilibrium in the infinite-time limit all differences may vanish.

### 1.2 Arbitrary finite switching schedules

Suppose a protocol uses any finite sequence of row-generators \(G_1,\ldots,G_m\), each reversible with respect to the **same** \(\pi\), each preserving \(\mathcal C\), for durations \(t_k\ge0\). Put \(P_k=e^{t_kG_k}\).

For the row-distribution convention the evolution is
\[
p_\beta^T\longmapsto p_\beta^T P_1P_2\cdots P_m.
\]
Detailed balance at each stage makes its likelihood-ratio propagator and adjoint
\[
K=P_mP_{m-1}\cdots P_1,\qquad
K^*_\pi=P_1P_2\cdots P_m. \tag{5}
\]
Both products preserve \(\mathcal C\). All factors are invertible, so strictness holds at every finite endpoint. The total product need not be reversible and the stages need not commute. Applying the same argument to any partial final stage gives the ordering throughout the protocol.

This statement compares two initial thermal states under **the same prescribed schedule**, which is independent of their preparation temperature. It neither compares two differently chosen protocols nor excludes a speed advantage from control. It does not allow changing the bath temperature, changing the energy order, or adding a Hamiltonian drive within the stated proof.

## 2. Proof

### 2.1 Gibbs comparison for the prefix-mean cone

For completeness, let \(b\in\mathcal C\) and put \(A_k=k^{-1}\sum_{i\le k}b_i\). The defining condition is equivalent to \(A_1\le\cdots\le A_N\). For positive decreasing weights \(w_i\), Abel summation gives
\[
\sum_{i=1}^k w_i(b_i-A_k)
=\sum_{\ell=1}^{k-1}(w_\ell-w_{\ell+1})
  \ell(A_\ell-A_k)\le0. \tag{6}
\]
Thus the weighted mean of \(b\) on a prefix is at most its arithmetic mean. Every component \(b_j\) with \(j>k\) is at least \(A_{j-1}\ge A_k\), so its weighted tail mean is at least \(A_k\).

For Gibbs weights at \(\beta\ge0\), write \(P_k=\sum_{i\le k}p_{\beta,i}\), and let \(L_k,R_k\) be the corresponding weighted prefix and tail means of \(b\). The cut decomposition of covariance is
\[
\operatorname{Cov}_{p_\beta}(E,b)
=\sum_{k=1}^{N-1}(E_{k+1}-E_k)P_k(1-P_k)(R_k-L_k)\ge0. \tag{7}
\]
All factors except possibly the last are strictly positive. If \(b\) is nonconstant, let \(j\) be the first component differing from \(b_1\). The cone inequality makes \(b_j>b_1\), and all later components are at least \(b_1\). Therefore the cut after \(j-1\) has \(R_{j-1}>L_{j-1}=b_1\); the covariance is strictly positive. For \(b\in-\mathcal C\), the signs reverse by applying the same argument to \(-b\).

Equivalently, the Gibbs expectation of every nonconstant \(b\in\mathcal C\) strictly decreases with \(\beta\) at every finite nonnegative inverse temperature. This is the existing Gibbs comparison lemma used in the original manuscript, included here to make the new inference self-contained; the ordering and summation-by-parts argument are not claimed as new.

### 2.2 Composition closure

If \(b\in\mathcal C\) and \(\phi\) is increasing and concave on an interval containing its components, then for every \(j\),
\[
\phi(b_j)\ge
\phi\!\left(\frac1{j-1}\sum_{i<j}b_i\right)
\ge\frac1{j-1}\sum_{i<j}\phi(b_i).
\tag{8}
\]
The first inequality uses monotonicity, the second Jensen's inequality. Hence \(\phi(b)\in\mathcal C\). If \(b\in-\mathcal C\) and \(\phi\) is increasing and convex, both inequalities reverse, giving \(\phi(b)\in-\mathcal C\).

The normalised initial density ratio has the form
\[
r_{\beta,i}=\frac{Z(\beta_b)}{Z(\beta)}
 e^{(\beta_b-\beta)E_i}. \tag{9}
\]
It is increasing in the ordered energies for cooling and decreasing for heating. Therefore \(r_\beta\in\mathcal C\) when \(\beta\le\beta_b\), and \(r_\beta\in-\mathcal C\) when \(\beta\ge\beta_b\).

### 2.3 Moving the propagator across the inner product

Differentiate (2), holding \(K,\pi\) and the energies fixed. This is the general \(f\)-divergence version of the known KL preparation-derivative/covariance calculation in Ohga–Hayakawa–Ito, Appendix B. The substantive sufficient condition here is the two uses of cone preservation around the Jensen composition step:
\[
\begin{aligned}
\frac{\partial}{\partial\beta}D_f(p_\beta^K\Vert\pi)
&=\langle f'(Kr_\beta),K\partial_\beta r_\beta\rangle_\pi\\
&=\langle K^*_\pi f'(Kr_\beta),\partial_\beta r_\beta\rangle_\pi\\
&=\sum_i p'_{\beta,i}
       [K^*_\pi f'(Kr_\beta)]_i\\
&=-\operatorname{Cov}_{p_\beta}
  \!\left(E,K^*_\pi f'(Kr_\beta)\right).
\end{aligned} \tag{10}
\]

For cooling, the first cone condition puts \(Kr_\beta\) in \(\mathcal C\). An increasing concave \(f'\) and (8) put \(f'(Kr_\beta)\) in \(\mathcal C\). The second cone condition puts \(b=K^*_\pi f'(Kr_\beta)\) in \(\mathcal C\). Equation (7) proves (3). The heating proof uses \(-\mathcal C\) throughout, which is invariant because the operators are linear, and the convex version of (8), proving (4).

For strictness, \(\beta\ne\beta_b\) makes \(r_\beta\) nonconstant. An invertible \(K\) fixing constants cannot send a nonconstant vector to a constant. A strictly increasing \(f'\) preserves nonconstancy. An invertible \(K^*_\pi\) fixing constants again preserves nonconstancy. The strict part of (7) now gives the claimed derivative sign. At \(\beta=\beta_b\), \(r_\beta=\mathbf1\) and \(b=f'(1)\mathbf1\), so (10) is zero. This proves the theorem.

## 3. Standard divergences and exact normalisations

The finite-time cooling theorem includes the following functions.

| Divergence | \(f(x)\) in (2) | \(f'''(x)\) |
|---|---|---|
| \(D_{\mathrm{KL}}(p\Vert\pi)\) | \(x\log x-x+1\) | \(-x^{-2}\) |
| \(D_{\mathrm{KL}}(\pi\Vert p)\) | \(x-1-\log x\) | \(-2x^{-3}\) |
| Squared Hellinger \(H^2(p,\pi)=\tfrac12\sum_i(\sqrt{p_i}-\sqrt{\pi_i})^2\) | \(\tfrac12(\sqrt x-1)^2\) | \(-\tfrac38x^{-5/2}\) |
| Half Pearson chi-square \(\tfrac12\chi^2(p\Vert\pi)\) | \(\tfrac12(x-1)^2\) | \(0\) |

The square root of Hellinger squared is an increasing transform, so the same ordering applies to the Hellinger distance. Both orientations of relative entropy here are cooling results. This argument supplies a heating theorem for Pearson chi-square but does not generally supply one for either orientation of relative entropy.

### 3.1 Power divergences

Use the established Cressie–Read power-divergence family, with our index \(\alpha\) equal to their index \(\lambda+1\). For \(\alpha\notin\{0,1\}\), define its normalised generator
\[
f_\alpha(x)=
\frac{x^\alpha-\alpha x+\alpha-1}{\alpha(\alpha-1)}.
\tag{11}
\]
Its continuous limits at \(\alpha=1\) and \(\alpha=0\) are respectively the forward and reverse KL generators in the table. Direct differentiation gives
\[
f_\alpha''(x)=x^{\alpha-2}>0,\qquad
f_\alpha'''(x)=(\alpha-2)x^{\alpha-3}. \tag{12}
\]
Thus the power divergence
\[
\mathcal P_\alpha(p\Vert\pi)
=\frac{\sum_i\pi_i(p_i/\pi_i)^\alpha-1}
       {\alpha(\alpha-1)}
\tag{13}
\]
is strictly ordered throughout every finite cone-preserving protocol:

- **Cooling:** every real \(\alpha\le2\), with the continuous limits at 0 and 1 included.
- **Heating:** every \(\alpha\ge2\).

Full support and finite preparation temperatures make all these finite for each finite protocol. At \(\alpha=2\), \(\mathcal P_2=\tfrac12\chi^2(p\Vert\pi)\), the common boundary valid on both sides of the bath. At \(\alpha=\tfrac12\), \(\mathcal P_{1/2}=4H^2\) under the table's Hellinger convention.

### 3.2 Rényi divergences

For \(\alpha>0\), \(\alpha\ne1\), use
\[
D_\alpha^{\mathrm R}(p\Vert\pi)
=\frac{\log \sum_i \pi_i(p_i/\pi_i)^\alpha}{\alpha-1}
=\frac{\log[1+\alpha(\alpha-1)\mathcal P_\alpha]}
       {\alpha-1}. \tag{14}
\]
The derivative of this expression with respect to \(\mathcal P_\alpha\) is
\[
\frac{\alpha}{1+\alpha(\alpha-1)\mathcal P_\alpha}>0. \tag{15}
\]
Hence Rényi divergence is strictly conventionally ordered for cooling when \(0<\alpha\le2\), including forward KL at \(\alpha=1\), and for heating when \(\alpha\ge2\).

Do not identify the reverse KL limit of (11) at \(\alpha=0\) with Rényi divergence of order zero. For full-support distributions \(D_0^{\mathrm R}=0\); it is not a strict distance ordering. The real-negative power-divergence extension also does not justify transferring the same ordering to a negative-order Rényi convention because the derivative in (15) changes sign.

### 3.3 Optional operator-convex-family corollary, restricted to populations

The standard operator-convex representation on the positive half-line expresses a normalised non-affine generator as a positive mixture of
\[
g_\lambda(x)=\frac{(x-1)^2}{(1-\lambda)x+\lambda},
\qquad 0\le\lambda\le1.
\tag{16}
\]
This representation is established mathematics; see Hiai, Mosonyi, Petz and Bény, and the normalised form in Salazar, Eqs. (17)–(20). Direct scalar differentiation gives
\[
g_\lambda''(x)=\frac{2}{[(1-\lambda)x+\lambda]^3},
\qquad
g_\lambda'''(x)=\frac{-6(1-\lambda)}
 {[(1-\lambda)x+\lambda]^4}\le0.
\tag{17}
\]
Therefore every member of this **classical population reduction** of the normalised operator-convex Petz family obeys the all-time cooling ordering. One may either apply Theorem C to the generators or integrate the strictly ordered kernel divergences in (16). This does not establish a theorem for noncommuting density matrices. It does give an explicit physical no-crossing class for the same population divergence family whose general ordering Salazar tests through a kernel profile. The family is still smaller than all scalar convex functions.

## 4. Physical application

Under the complete-graph, equal-pair-coupling bosonic model of *concave_kernel_theorem.md*, every continuous nonnegative concave jump weight \(J\) on the entire energy-gap interval gives a reversible, \(\mathcal C\)-preserving generator at each finite positive bath temperature. Theorem C therefore upgrades that physical class from a slow-mode exclusion to **all-time cooling order for KL, reverse KL, Hellinger, and the stated power/Rényi families**, as well as all-time chi-square order on either side of the bath.

One may also switch arbitrarily often through any prescribed **finite** schedule of such concave \(J_k\), preserving the same energies, bath and common pairwise coupling assumption at every stage. Multiplying a stage by a nonnegative rate scale is allowed. The resulting product obeys the same inequalities even when its factors do not commute.

The earlier energy-threshold result also extends throughout such schedules: the backward observable product \(P_1\cdots P_m\) preserves \(\mathcal C\), so every increasing energy observable has conventionally ordered Gibbs expectations. This energy-distribution order is distinct from a divergence to equilibrium.

## 5. Scope and exclusions

The theorem is a sufficient criterion for these divergence orderings. It is not an iff classification of actual Mpemba behaviour. The separate necessity construction for concavity classifies preservation of the entire prefix-mean cone, not thermal crossings.

The result does **not** prove finite-bath total-variation ordering or thermomajorization ordering. Total variation corresponds to \(f(x)=\tfrac12|x-1|\), whose derivative is not globally concave or globally convex; the composition step does not apply. Thermomajorization would imply ordering for all convex \(f\), a stronger statement than the one proved here. No result outside the derivative class follows from the argument alone.

The theorem concerns a finite classical population process. Its abstract operator form allows any graph satisfying (1), but the proved bosonic application assumes a complete graph and common pairwise coupling. No general conclusion for state-dependent matrix elements or sparse graphs follows from that application. Coherent quantum evolution, negative preparation temperatures, varying energy levels, varying bath temperatures, or a control schedule selected separately for each initial temperature also require additional arguments.

The proof is elementary once the physical cone-preservation result is available. Its potential significance is the physical class and its robust finite-time/protocol conclusion, not a claim that Jensen's inequality, adjoint differentiation, power divergences, or positive-cone methods are new.

## 6. Verification and prior-art status

The analytical proof above has passed an independent internal adversarial derivation of the cone composition, strictness, product orientation and power/Hellinger/Rényi normalisations. The exact checker *verify_metric_protocol.py* passes 6,493 tracked rational assertions in 36 finite-operator cases, including 864 derivative evaluations and 12 nonreversible three-stage products. It verifies the direct derivative, adjoint and covariance formulas, both cone signs, and pairwise power-divergence ordering.

The exact operators in that checker are stochastic uniformization steps, not exponentials; they audit the general \(K\) theorem. The analytical proof establishes continuous-time semigroups and arbitrary finite switching schedules. The bounded source comparison is summarised in `PRIOR_ART.md` at the package root. No external specialist review is claimed.

## Primary sources inspected for this extension

- Ohga, N., Hayakawa, H., and Ito, S. *Microscopic theory of Mpemba effects and a no-Mpemba theorem for monotone many-body systems*, arXiv:2410.06623v2, Appendix B, Result S2, Lemma S3 and Eqs. (S14)–(S20). https://arxiv.org/html/2410.06623v2
- Vu, T. V., and Hayakawa, H. *Thermomajorization Mpemba Effect*, Physical Review Letters **134**, 107101 (2025). The all-convex-\(f\) criterion is stronger than this note's conclusion. https://arxiv.org/html/2410.06686v2
- Cressie, N., and Read, T. R. C. *Multinomial Goodness-of-Fit Tests*, Journal of the Royal Statistical Society B **46**, 440–464 (1984). Power-divergence family and parameter convention. https://doi.org/10.1111/j.2517-6161.1984.tb01318.x
- Hiai, F., Mosonyi, M., Petz, D., and Bény, C. *Quantum \(f\)-divergences and error correction*, arXiv:1008.2529v6, Section 8. The integral representation is a method antecedent. https://arxiv.org/pdf/1008.2529
- Salazar, D. S. P. *Divergence Geometry of Quantum Multi-Mpemba Effects*, arXiv:2609.20320, Eqs. (2), (17)–(27). It studies family-wide order through the known operator-convex kernel representation, not a concave-bath Gibbs no-crossing theorem. https://arxiv.org/pdf/2609.20320
- Summer, A., Moroder, M., Bettmann, L. P., Turkeshi, X., Marvian, I., and Goold, J. *Resource-Theoretical Unification of Mpemba Effects: Classical and Quantum*, arXiv:2507.16976v2. The inspected framework uses resource monotones and slow modes. https://arxiv.org/pdf/2507.16976
- Wang, Y., and Lu, Z. *Slow is fast: raising barriers to accelerate thermal relaxation*, arXiv:2607.11877v1. Its protocol optimisation compares controls for a given preparation; our theorem compares preparations under a shared restricted control. https://arxiv.org/html/2607.11877v1
