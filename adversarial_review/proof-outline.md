# Independent proof outline

This outline records the reasoning independently checked during the review. It reorganises the candidate arguments without editing their source files. All claims use finite, strictly ordered energy spectra and finite nonnegative preparation inverse temperatures unless an explicit zero-bath limit is stated.

## Common Gibbs comparison lemma

For b in the arithmetic-prefix cone, the arithmetic means A_k of the first k components increase with k. Abel summation shows that a decreasing positive weighting of a prefix has mean at most A_k. Every later component is at least A_k, hence every weighted tail mean is at least every corresponding prefix mean. Decompose E into its positive adjacent increments to obtain

\[
\operatorname{Cov}_{p_\beta}(E,b)
=\sum_{k=1}^{N-1}(E_{k+1}-E_k)P_k(1-P_k)(R_k-L_k)\ge0.
\]

For nonconstant b the first departure from its initial component is upward. The cut immediately before that departure is strict, making the covariance positive. Therefore the Gibbs expectation decreases strictly with β. The negative cone reverses the signs.

An independent equivalent check uses the contrast basis u^(j). Its pairwise covariance has numerator

\[
jq_j\sum_{a<j}q_a(E_j-E_a)
+\sum_{b>j}q_b\left[
\sum_{a<j}q_a(E_b-E_a)-(j-1)q_j(E_b-E_j)\right],
\]

where q_i=exp(−βE_i). Each bracket is positive because q_a≥q_j and E_a<E_j. Every nonconstant cone vector has at least one positive contrast coefficient. This independently confirms the covariance sign, including β=0 and arbitrary upper energy gaps.

## Hierarchy: six steps

1. **Triangular process.** Choose the downward rates with m=N−3 middle states and 0<ε<1/14. Direct comparison gives increasing escape rates in every retained-state restriction, so every spectrum is simple at zero bath.
2. **Proper-subset split.** Divide all subsets into four classes according to the presence of the two distinguished low states. The middle states have a common survival response within each class because their internal transition terms cancel.
3. **Survival-prefix positivity.** In the class retaining both distinguished states, the top-minus-prefix-mean comparison obeys a first-order scalar equation with positive initial condition and positive forcing F(r), κ for every proper r≤m−1. In the other classes, after factoring the slow exponential, the last comparison is affine in an exponential parameter ranging over [0,1], and both endpoints are positive. The remaining prefix comparisons are immediate. Apply the common Gibbs lemma to survival probabilities to obtain all-time zero-bath TV order.
4. **Full-system sign reversal.** Set r=m. Exactly one slow-prefix comparison becomes negative, through F(m)<0. At ε=1/100 its deficit is −1/[(m+2)P_m]. Put all earlier energies at 0,1,…,m+1 and choose the last energy beyond the displayed bound. Direct covariance calculation gives A′(0)>0, hence a finite pair of positive inverse temperatures has inverted slow amplitudes.
5. **A real crossing.** The hotter preparation initially has larger TV from the ground. Its positive slow amplitude is smaller. The simple dominant decay term therefore reverses TV order at sufficiently late finite time. Continuity in time gives a crossing.
6. **Positive-bath continuation.** Add Gibbs-compatible upward rates. For fixed N and energies, simple spectra and the finite list of strict proper-prefix inequalities persist near the downward limit. These finite eigenfunction inequalities imply signed-amplitude monotonicity on the entire β≥0 half-line. Preserve the full pair's strict amplitude and initial-distance inequalities as well. Choose the bath colder than both preparations. This proves the finite-bath theorem and its open perturbation neighbourhood.

The proof has one main theorem, one common ordering lemma and four nonredundant subsystem cases. The important direct algebra was re-derived independently, including the forcing coefficient, rational deficit and A′(0). No all-N conclusion rests on enumeration.

## Concave-kernel theorem: five steps

1. **Choose coordinates.** Contrast coordinates identify the prefix cone modulo constants with the positive orthant. The generator is represented by H_ij=(u^(i))ᵀGu^(j)/[(i−1)j]. Cone preservation is exactly the Metzler condition.
2. **Thermal part.** Concavity and nonnegativity of J give J(d)/d nonincreasing. Its product with d/(exp(β_b d)−1) is decreasing. Direct transformation of the symmetric thermal rates yields nonnegative off-diagonal entries. Positivity at the gaps makes them strict.
3. **Downward part.** The transformed upper entries vanish. Concavity makes each negative sum of fixed-length increments nondecreasing. The transformed lower numerator is a sum of its nonnegative increments plus the nonnegative boundary term (j−1)J(0).
4. **Consequences.** Sum the two parts. The common Gibbs lemma orders all evolved increasing energy observables. If the gap weights are positive, the transformed matrix has positive off-diagonal entries; Perron–Frobenius on the quotient gives a simple slow mode with positive prefix coordinates. Lift and centre it, then apply the Gibbs lemma to obtain the same-side slow-amplitude exclusion. At zero bath, ground survival is TV.
5. **Universal necessity.** A continuous nonconcave J has an increasing translated increment somewhere. A sufficiently large cluster of retained levels amplifies that violation into one negative transformed downward entry. Keep the levels distinct, fix them, and make the finite bath sufficiently cold. This violates cone invariance. The proof does not produce a Gibbs crossing, nor does it classify a fixed N or fixed bath.

The exact linear-J three-level model separates this cone criterion from ordinary stochastic monotonicity: the higher state initially reaches the ground faster than the middle state. That prevents a total-energy-order monotone coupling while leaving the Gibbs-family ordering intact.

## Divergence and protocol theorem: five steps

1. **Work with densities.** For π=p_βb and r_β=p_β/π, write the evolved density as Kr_β. Assume K and K* both preserve the prefix cone. K is Markov and π-stationary, so the evolved distribution remains positive and normalised.
2. **Use composition closure.** Jensen proves that an increasing concave scalar map sends the prefix cone to itself. An increasing convex map sends the negative cone to itself. The initial r_β lies in the first cone for cooling and in the second for heating.
3. **Differentiate correctly.** Holding K, π and the energies fixed, move K across the π-weighted inner product. The divergence derivative becomes −Cov_pβ(E,K*f′(Kr_β)). The two cone assumptions and the closure step determine its sign.
4. **Prove strictness separately.** Invertible Markov K and K* fix constants and cannot send a nonconstant vector to a constant. Strictly increasing f′ also preserves nonconstancy. The strict Gibbs covariance lemma applies for every finite β≠β_b. At β=β_b the ratio is one and the derivative vanishes.
5. **Pass to a schedule.** Each reversible stage is its own π-adjoint. A row-distribution schedule P₁…P_m has density propagator P_m…P₁ and adjoint P₁…P_m. Both preserve the cone and both are invertible at finite durations, regardless of commutativity. Apply the theorem to every partial endpoint.

The standard divergence formulas are corollaries of the curvature signs of f′, not extra dynamical assumptions. Cooling power parameters extend to all real α≤2, but Rényi monotone transfer requires α>0. No step gives TV or an all-convex-f order.

## Logical dependency summary

The Gibbs-prefix lemma supports the hierarchy, energy-ordering and divergence conclusions. The concave-kernel result supplies one explicit physical family satisfying the abstract two-operator criterion. Its necessity classifies that cone property only. The hierarchy supplies a separate negative result for unrestricted detailed-balance rates. Combining the theorems is scientifically coherent but does not turn either result's hypotheses into the other's.
