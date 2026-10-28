# Assurance boundary

## What was done

1. **Written proofs** in [manuscript/thermal_order.tex](manuscript/thermal_order.tex), with derivations of the main matrix formulas in Appendix C.
2. **Exact certificates.** Eight standard-library checkers in rational or outward-rounded interval arithmetic ([README.md](README.md) table), compared exactly with retained outputs by `run_checks.py`. A separately implemented floating-point audit (`adversarial_review/independent_checks.py`) and a separate interval TV certificate (`adversarial_review/independent_tv_certificate.py`) use different code paths.
3. **Internal AI-assisted review.** Round 1: an editor, methodology, domain, perspective and devil's-advocate reviewer and an outside-model correctness audit, plus a clean replay of every script from a fresh extraction. One revision round, including four extensions. Verification: a re-review of all changes, a targeted proof audit with its own exact code, and a final outside-model audit of the new proofs; their findings were fixed as [FINAL_FIXES.md](FINAL_FIXES.md). A rating-free summary is in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
4. **Replay.** From a fresh extraction of the frozen archive, the default checks and the negative controls were rerun in normal mode and with `-O -B`; see [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). CI replays the default checks on Linux.
5. **Semantic negative controls.** [negative_controls.py](negative_controls.py) applies five claim-relevant corruptions (concave weight in place of d³, swapped preparations, convex test weight, a far-perturbed Proposition 9 model, a certificate time outside the reversal window); each must be rejected.
6. **Citations.** Every DOI in the bibliography resolves and matches its Crossref record ([CITATION_AUDIT.md](CITATION_AUDIT.md)); arXiv entries were checked against arXiv records during review.

## What was not done

- No human specialist or journal peer review, no independent reproduction by an unaffiliated party, no proof-assistant formalisation.
- **Priority is bounded** ([PRIOR_ART.md](PRIOR_ART.md)). The closest antecedent, Avitan, Factor and Gelbwaser-Klimovsky (2026), already excludes the slow-mode Mpemba effect for Ohmic and sub-Ohmic spectra; the physical message is not new here. The increment is the all-time, all-temperature, divergence-family and necessity statements for the stated model class.
- The model class is narrow: complete transition graph, common pair coupling, classical populations.
- Final editor fixes were made without a further review round.
- Agreement between AI reviewers is not independent evidence.
