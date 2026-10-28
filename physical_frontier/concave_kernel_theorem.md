# Concave bosonic jump weights preserve thermal ordering

Research note, 26 October 2026. This is an analytical extension of the current manuscript. Its bounded prior-art assessment is summarised in `PRIOR_ART.md` at the package root (the original development report is not shipped). It is not an external review.

## 1. Result and physical assumptions

Let `E_1 < ... < E_N`, let `D = E_N-E_1`, and fix a bath inverse temperature `beta_b > 0`. A row-generator `G` acts on observables. Assume a complete transition graph and the same coupling prefactor for every pair. Its off-diagonal rates are

\[
G_{ij}=\frac{J(|E_i-E_j|)}{e^{\beta_b|E_i-E_j|}-1}
       +\mathbf 1_{E_i>E_j}J(|E_i-E_j|),\qquad i\ne j.
\]

Rows sum to zero. Here `J` is the effective spontaneous-emission jump weight, incorporating any frequency dependence of the matrix element. It must not be equated with a named bath spectral density before matching that convention. The rates satisfy detailed balance with `pi_i proportional exp(-beta_b E_i)`.

Define the observable cone

\[
\mathcal C=\left\{b:\ b_j-\frac1{j-1}\sum_{i<j}b_i\ge0,
\quad j=2,\ldots,N\right\}.
\]

**Theorem A (concave jump weights).** If `J` is continuous, nonnegative and concave on `[0,D]`, then `exp(tG) C` is contained in `C` for all `t>=0`. Consequently, for every pair `0<=beta_h<beta_c`, every increasing energy observable `b`, and every `t>=0`,

\[
p_E(\beta_h)^T e^{tG}b\ \ge\ p_E(\beta_c)^T e^{tG}b.
\]

Equivalently, the two evolving energy distributions remain in first-order stochastic order at every time. This includes every energy threshold, the mean energy and every increasing function of energy. It needs neither a simple-spectrum assumption nor a restriction placing both preparations on the same side of the bath temperature.

If `J(d)>0` at every nonzero gap, then the slow nonzero eigenvalue is simple automatically, and its eigenfunction lies strictly inside `C` up to an overall sign and an additive constant. Its Gibbs amplitude is strictly monotone and vanishes only at the bath inverse temperature. Thus neither weak slow-mode inversion nor strong slow-mode cancellation occurs between two nonnegative preparation inverse temperatures on the same side of the bath.

**Zero-bath corollary.** The cone result also holds in the downward limit. If every excited state can reach the ground state, its equilibrium is the ground-state point mass. Because total variation then equals the survival probability `1-p_1(t)`, total-variation ordering is conventional at every time.

**Finite-time divergence extension.** The separate proof in *finite_time_metric_theorem.md* uses reversibility, this cone invariance and Jensen's inequality to prove all-time cooling ordering for KL, reverse KL, Hellinger, and power divergences with index at most 2; Rényi indices are restricted to \(0<\alpha\le2\). Chi-square is ordered on both sides of the bath. The same result survives finite switching schedules of possibly noncommuting concave-kernel generators at a common bath and fixed energies. The known KL preparation-derivative identity is credited to Ohga, Hayakawa and Ito.

**Scope.** At positive bath temperature neither proof establishes total-variation or all-convex-\(f\) thermomajorization ordering. The physical application assumes a complete graph with common pairwise coupling; it does not cover arbitrary sparse graphs, state-dependent coupling strengths, or quantum coherences.

## 2. Coordinates and exact transformation

For `j=2,...,N`, let the Helmert contrast `u^(j)` have entries `-1` before position `j`, entry `j-1` at `j`, and zero afterwards. These vectors are mutually orthogonal and orthogonal to constants. Set

\[
d_j(b)=\frac{(u^{(j)})^Tb}{j-1}
      =b_j-\frac1{j-1}\sum_{i<j}b_i.
\]

Modulo constants, `b = sum_{j=2}^N d_j u^(j)/j`. Therefore the exact generator in these coordinates is

\[
H_{ij}=\frac{(u^{(i)})^TG u^{(j)}}{(i-1)j},\qquad i,j=2,\ldots,N.
\]

Cone preservation is equivalent to the off-diagonal entries of `H` being nonnegative. This is the standard equivalence between preservation of an orthant by a continuous semigroup and the Metzler property, transported through the chosen coordinates.

## 3. Proof of Theorem A

Split `G=S+A`. The symmetric thermal generator `S` has rates

\[
h(d)=\frac{J(d)}{e^{\beta_b d}-1},
\]

and `A` is the downward generator with rates `J(E_i-E_j)` for `i>j`.

### 3.1 Symmetric thermal part

Concavity and `J(0)>=0` imply that `J(d)/d` is nonincreasing for `d>0`. The function `d/(e^(beta_b d)-1)` is positive and strictly decreasing. Their product `h` is nonincreasing, and is strictly decreasing wherever `J` is positive.

For `2<=i<j<=N`, direct expansion gives

\[
H^S_{ij}=
h(E_j-E_i)-\frac1{i-1}\sum_{k<i}h(E_j-E_k)
 +\frac1j\sum_{\ell>j}\left[
h(E_\ell-E_i)-\frac1{i-1}\sum_{k<i}h(E_\ell-E_k)\right].
\]

Every bracket is nonnegative because its first gap is smaller than each gap in the average. Symmetry of `S` gives

\[
H^S_{ji}=\frac{(i-1)j}{(j-1)i}H^S_{ij}\ge0.
\]

Thus every off-diagonal thermal entry is nonnegative. When all relevant jump weights are positive, the first bracket is strictly positive, so all these entries are strictly positive.

### 3.2 Downward part

The upper entries `H^A_ij`, `i<j`, vanish. For a fixed column `j`, write

\[
S_j=\sum_{\ell<j}J(E_j-E_\ell),\qquad
F_j(x)=(j-1)J(x-E_j)-\sum_{\ell<j}J(x-E_\ell),\quad x\ge E_j.
\]

Concavity makes every increment `J(z+a)-J(z)` nonincreasing in `z`. Hence `F_j` is nondecreasing. Also `F_j(E_j)+S_j=(j-1)J(0)>=0`.

For `i>j`, direct multiplication gives

\[
\begin{aligned}
(i-1)jH^A_{ij}
&=(i-1)F_j(E_i)-\sum_{k=j+1}^{i-1}F_j(E_k)+jS_j\\
&=j\{F_j(E_i)+S_j\}
  +\sum_{k=j+1}^{i-1}\{F_j(E_i)-F_j(E_k)\}\ge0.
\end{aligned}
\]

Therefore `H=H^S+H^A` is Metzler. This proves invariance of `C`.

### 3.3 Thermal preparation ordering

The Gibbs monotonicity lemma in the original manuscript states that, for `b in C`, the expectation `p_E(beta)^T b` is nonincreasing on `beta>=0`; it is strictly decreasing if `b` is nonconstant. Its proof uses Abel summation and a covariance identity. Apply it to `e^(tG)b`, which belongs to `C` whenever `b` does. Every increasing observable belongs to `C`. In particular, every indicator `b_i=1_{i>=k}` belongs to `C`, giving the stated energy-threshold and stochastic ordering.

### 3.4 Slow-mode conclusion

If `J` is positive at the gaps, the finite-temperature `H` has strictly positive off-diagonal entries. Its spectral bound is a simple real eigenvalue with a strictly positive eigenvector. The quotient construction identifies its spectrum with the nonzero spectrum of `G`. Detailed balance makes that spectrum real, and irreducibility places it below zero. The reconstructed eigenfunction has strictly increasing prefix means. Subtracting its stationary mean produces a genuine slow eigenfunction `f` with `pi^T f=0` without changing those prefix differences. Its amplitude `p_E(beta)^T f` strictly decreases and is zero at `beta_b`. Magnitudes therefore have conventional ordering on each side of the bath. For `N=2`, the same conclusion follows from its single nonzero mode.

## 4. Sharpness for universal preservation of this cone

The next statement was developed with the coordinating agent. It concerns the cone property, not an iff criterion for actual Mpemba behaviour.

**Theorem B (universal cone classification).** Fix `D>0` and continuous nonnegative `J` on `[0,D]`. The following are equivalent:

1. `J` is concave.
2. For every finite ordered spectrum in an interval of width at most `D` and every finite `beta_b>0`, the bosonic generator above preserves `C`.

Theorem A proves sufficiency. For necessity, nonconcavity supplies `h>0` and `0<=x<y`, with `h+y<=D`, such that

\[
A(y)>A(x),\qquad A(z)=J(z+h)-J(z).
\]

Choose `N=M+3` levels: `0`, `h`, `M` distinct levels clustered near `h+x`, and the top level `h+y`. For the zero-bath downward generator, its transformed entry in row `N`, column `2` tends, as the cluster shrinks, to

\[
H^A_{N2}\longrightarrow
\frac{2J(h)+M A(x)-(M+2)A(y)}{2(M+2)}.
\]

The numerator becomes negative for sufficiently large `M`. The inequality is strict, so a sufficiently small cluster of distinct levels retains it. At those fixed distinct energies, the thermal part tends to zero as `beta_b` tends to infinity. Thus the same transformed entry is negative for all sufficiently cold finite baths. A non-Metzler entry prevents cone preservation, proving necessity.

This argument gives no implication from nonconcavity to an actual thermal crossing. The observable that witnesses failure of cone invariance need not be an increasing energy observable, and the energies selecting Gibbs preparation are tied to the generator.

## 5. Exact separators and examples

For `E=(0,1,2)`, `beta_b=log 2`, `J(d)=d`,

\[
G=\begin{pmatrix}-5/3&1&2/3\\2&-3&1\\8/3&2&-14/3\end{pmatrix},\qquad
H=\begin{pmatrix}-23/6&1/3\\1/4&-11/2\end{pmatrix}.
\]

The characteristic polynomial is `lambda(lambda^2+28lambda/3+21)`. Its nonzero roots are `(-14 +/- sqrt(7))/3`.

The usual monotone coupling on the energy order fails: a chain starting at energy 2 reaches energy 0 initially at rate `8/3`, while one starting at energy 1 does so at rate `2`. For small time, the initially higher point mass therefore puts more probability below both starting states. Theorem A nonetheless orders all Gibbs preparations. On a complete graph with strictly ordered energies, an energy-compatible partial order satisfying comparability of every edge must be this total order, so changing the partial order does not repair the standard coupling hypothesis.

For `J(d)=d^2`, `E=(0,1,2)` and zero bath, `H_32=-1`. This is an exact elementary failure of the cone criterion outside concavity.

For a power-law jump weight `J(d)=d^alpha`, `0<=alpha<=1`, Theorem A applies at every energy range and every positive bath temperature. If an exponential cutoff is included, `J(d)=d^alpha exp(-d/d_c)` with `0<alpha<=1`, its second derivative is nonpositive for

\[
0\le d\le d_c(\alpha+\sqrt\alpha).
\]

Thus Theorem A holds when the entire Bohr-frequency range lies in that interval. It does not justify discarding the cutoff for a larger energy range.

At a fixed bath, `J(d)=1-exp(-beta_b d)` reduces the rates exactly to complete-graph Metropolis dynamics: all downward rates are 1 and upward rates are `exp(-beta_b d)`. Its no-Mpemba property is already proved by Klich and Vucelja. This is a verified collision for that special case and must be cited, not claimed as new.

### 5.1 A nonconcave kernel invisible to every three-level cone test

The independent reviewer suggested \(J(d)=1+d^2\) on \([0,1]\) as a useful distinction between triplet and global cone tests. The following six-level refinement is exact. For any three levels with adjacent gaps \(a,b>0\), \(a+b\le1\), the downward transformed entry is
\[
H^A_{32}=\tfrac12\{J(a)+J(b)-J(a+b)\}
=\tfrac12(1-2ab)\ge\tfrac14.
\]
Also \(J(d)/d=d^{-1}+d\) is nonincreasing on \((0,1]\), so the symmetric thermal part has nonnegative transformed off-diagonals at every positive bath temperature. Thus every three-level spectrum in this range preserves the cone at every bath temperature despite this kernel's convexity.

For the six-level spectrum
\[
E=(0,\tfrac12,\tfrac{501}{1000},
       \tfrac{502}{1000},\tfrac{503}{1000},1),
\]
the downward entry is \(H^A_{62}=-247/5000\). At the finite bath \(\beta_b=10\log2\), its thermal part is
\[
H^S_{62}=\tfrac35\left\{\frac{5/4}{2^5-1}
-\frac{2}{2^{10}-1}\right\}=\frac{157}{6820}.
\]
Consequently
\[
H_{62}=-\frac{44977}{1705000}<0.
\]
The closely spaced middle-level thermal rates cancel from this particular transformed entry. This example proves a failure of the global cone criterion despite universal three-level cone preservation for the same physical kernel; it does not establish an actual thermal Mpemba crossing. No minimality claim for six levels is made.

## 6. Numerical and exact verification scope

`physical_fingerprints.py` produces exact rational generators and transformed matrices at 3, 4 and 5 levels for the linear and matched-exponential jump weights. `small_fingerprints.json` retains their characteristic polynomials and entries.

`verify_concave_kernel.py` passes 5,859 exact finite assertions across 96 rational models, with 3, 4, 5 and 7 states, three energy designs, four concave kernels and two bath parameters. It compares the quotient and Helmert transforms, checks detailed balance, verifies the separate thermal and downward formulas, and checks four exact uniformization steps of stochastic ordering for Gibbs preparations. The nonconcave quadratic example has `H_32=-327/340` at the finite bath `beta_b=log 16`. The initial script run exposed an empty Python sum that coerced a rational expression to floating point; adding the rational zero corrected the checker before the successful run. No proof formula changed.

`transient_probe.py` is a separate numerical stress probe over 400 random spectra, sizes 3-10 and concave power weights. It found no transient total-variation reversal; the smallest sampled stochastic-order gap was about `-2.1e-15`, consistent with roundoff. This failed search is not a theorem about total variation and is not used in the proof.

## 7. Contribution judgement

Theorem A converts a bath-shape restriction into an all-time ordering result without diagonalisation. Theorem B gives a sharp classification for the chosen cone across all sizes, spectra and bath temperatures. Both are broader than a single counterexample. Their physical applicability depends critically on the complete transition graph and uniform matrix elements. A paper combining this exact physical class with the global no-local-screening result would still need full prior-art and adversarial proof review.

## Primary references inspected

- Avitan, I., Factor, R., & Gelbwaser-Klimovsky, D. (2026). Necessary conditions for the Markovian Mpemba effect (v1), Eq. (5) and discussion of spectral exclusions. https://arxiv.org/html/2603.04567v1
- Klich, I., & Vucelja, M. (2019). Solution of the Metropolis dynamics on a complete graph with application to the Markov chain Mpemba effect (v2), Secs. II and IV.1. https://arxiv.org/html/1812.11962v2
- Ohga, N., Hayakawa, H., & Ito, S. (2024). Microscopic theory of Mpemba effects and a no-Mpemba theorem for monotone many-body systems (v2), Result 2 and proof. https://arxiv.org/html/2410.06623v2
- Atkins, J. E., Boman, E. G., & Hendrickson, B. (1998). A spectral algorithm for seriation and the consecutive ones problem. SIAM Journal on Computing, 28(1), 297-310. Theorem 3.2 uses adjacent differences for symmetric Robinson Laplacians, a related positive-transform method. https://doi.org/10.1137/S0097539795285771
- Lorek, P., & Szekli, R. (2012). Strong stationary duality for Möbius monotone Markov chains. Queueing Systems, 71, 79-95. Observable-cone preservation and transformed positivity are established general methods. https://doi.org/10.1007/s11134-012-9284-z
- Vu, T. V., & Hayakawa, H. (2025). Thermomajorization Mpemba effect. Physical Review Letters, 134, 107101. Its stronger metric-independent notion is not claimed here. https://doi.org/10.1103/PhysRevLett.134.107101
