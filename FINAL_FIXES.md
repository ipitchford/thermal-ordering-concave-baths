# Final fixes (from reviews/proof-audit.md, reviews/re-review.md, sol2/sol-audit.md)
F1 Remark 11: -391/100 is NOT an exact root (sympy). Say eigenvalue lies in a certified rational bracket around -3.91; spectrum simple since brackets disjoint.
F2 Prop 9 proof + App A: at t=10.46 and 10.48 the preparations have different positive sets; say each TV is bounded on its own certified set; print the exact rational sign bounds at t=10.46 and 10.48.
F3 Prop 8: state slow mode normalisation used for amplitudes; make certificate self-contained: print rational KL-difference bounds at the certified times (supplementary table ok).
F4 Prop 6(a): replace "any cutoff lying above D" by "with no cutoff factor on [0,D], e.g. a hard cutoff above D" (smooth cutoffs change curvature).
F5 Corollary 7 takeaway: condition uses full spectral width (top minus bottom), not largest spacing; narrow to cooling Gibbs preparations in the stated model.
F6 Avitan mapping: quote their Eq. (5) Delta^{s-1} explicitly (arXiv version) when stating their s = our s+1... check direction carefully and state it exactly.
F7 Theorem 10 finite-q continuation: explicitly take minimum over the finite family of subsystem margins and spectral gaps.
F8 Temper "physical dividing line": concavity characterises the universal cone criterion, not existence of crossings in general.
F9 Optional: small schematic of two crossing relaxation curves in Sec 1.1 if cheap (matplotlib, illustrative, labelled).
