# Review summary (rating-free public record)

All review was internal and AI-assisted. The original reports are kept privately. This summary keeps the substance of the objections and how each was resolved, without scores or ratings.

## Round 1

Five reviewer roles (editor, methodology, domain, perspective, devil's advocate), an outside-model correctness audit and a clean replay of every shipped script. No critical issue was raised, and every exact checker passed on a fresh extraction. The merged fix list asked for:

- **Proof completeness:** a step-by-step Lemma 1; the normalisation of Theorem 2; Corollary 3 stated as a corollary of closure; a precise operator-convex integral representation; a proved increment lemma and an explicit constant in the necessity proof of Theorem 4, with the joint all-N/cold-bath quantifier stated; fixed eigenfunction normalisation, existential finite-bath threshold and stated topology in the hierarchy theorem; an appendix of matrix derivations.
- **Certificate disclosure:** printed rational intervals and sign bounds for the total-variation example, effect size, robustness and terminology; float tolerance only for the floating-point audit jobs.
- **Prior art:** state plainly that Avitan, Factor and Gelbwaser-Klimovsky (2026) already exclude the slow-mode effect for Ohmic and sub-Ohmic spectra, and narrow the claimed increment; cite stochastic-order classics; add missing Mpemba and heating/cooling-asymmetry literature; update stale references.
- **Physical relevance:** classify standard bath spectral densities; certify an actual super-Ohmic crossing; if tractable, certify a finite-bath N=4 hierarchy instance; state the sparse-graph scope.
- **Readability:** plain-first abstract; a picture of crossing curves; definitions before use; one scope-and-limits section.

Each item and its disposition is in [REVISION_RESPONSE.md](REVISION_RESPONSE.md).

## Verification of the revision

- **Re-review.** Every item resolved, except the finite-bath N=4 instance, which is partial as disclosed (one certified q, no explicit threshold; acceptable because the theorem stays existential). Spot checks reproduced the Proposition 8 generator and KL values and the Proposition 6 second derivatives.
- **Proof audit,** with its own exact code: Lemma 5, the explicit necessity constant, Proposition 6, Corollary 7, Appendix C, Lemma 1 and the hierarchy normalisation passed. It asked for one clarity fix in Proposition 8, wording in Proposition 9, and the removal of one false phrase in Remark 11 (a bracketed eigenvalue had been called an exact root).
- **Final outside-model audit.** It found Theorems 2 and 4 essentially correct. It asked for an exact check of the exponent convention mapping to Avitan et al., self-contained certificate details for Proposition 8 and Remark 11, narrower paraphrases of what Theorem 2 excludes, the cutoff wording in Proposition 6(a), the spectral-width condition in Corollary 7, and an explicit finite minimum in the hierarchy continuation.

All of these were fixed in [FINAL_FIXES.md](FINAL_FIXES.md) (F1–F9), without a further review round.

## Open after review

Total-variation order for the concave physical class; maximal Rényi ranges; quantum coherences; sparse graphs and arbitrary couplings; an explicit q_0(N). The prior-art search was bounded.
