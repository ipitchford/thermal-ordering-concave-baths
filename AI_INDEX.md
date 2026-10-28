# Agent navigation

## Identity

- Release slug: `thermal-ordering-concave-baths`. Version `1.0.0-candidate`, tag `v1.0.0-candidate`. DOI 10.5281/zenodo.23721564; concept DOI 10.5281/zenodo.23721563.
- Repository: https://github.com/ipitchford/thermal-ordering-concave-baths. Release page: https://evidencepress.org/releases/thermal-ordering-concave-baths/
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md)).
- Manuscript: [PDF](manuscript/thermal_order.pdf) (27 pages) and its self-contained [LaTeX source](manuscript/thermal_order.tex).
- Machine-readable claims: [CLAIMS.json](CLAIMS.json). Human overview: [README.md](README.md).
- Assurance limits: [ASSURANCE.md](ASSURANCE.md). Review record: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md), [REVISION_RESPONSE.md](REVISION_RESPONSE.md), [FINAL_FIXES.md](FINAL_FIXES.md).
- Prior art and sources: [PRIOR_ART.md](PRIOR_ART.md), [RESEARCH_GATES.json](RESEARCH_GATES.json), [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md). Provenance: [PROVENANCE.md](PROVENANCE.md).
- Citation: [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json). Licences: [LICENSES.md](LICENSES.md).

## Scientific map

Object: a finite classical Markov jump process on energy levels with one bath; Gibbs preparations p_beta; the cone C of observables whose arithmetic prefix means are nondecreasing in energy order. Question: when do two Gibbs preparations on the same side of the bath stay ordered in a divergence from equilibrium at every time (no Markovian Mpemba crossing)?

- Lemma 1: Gibbs likelihood-ratio comparison extends from increasing observables to C.
- Theorem 2 and Corollary 3: if K and its pi-adjoint preserve C, cooling order holds at every time for f-divergences with increasing concave f' (heating: convex); also for common schedules of noncommuting reversible stages. Table 2 lists covered divergences; total variation is not covered.
- Theorem 4 and Lemma 5: for complete-graph bosonic rates with common coupling, cone preservation for every spectrum and bath temperature is equivalent to concavity of the jump weight J.
- Proposition 6 and Corollary 7: concavity thresholds of standard spectral densities (Ohmic width <= 2 w_c, Drude-Lorentz <= sqrt3 w_c; super-Ohmic never). Section 4.5 states the rate convention and the mapping to Avitan et al.
- Proposition 8: certified super-Ohmic (J=d^3) relative-entropy crossing in (1/5,3/10).
- Proposition 9, Table 4, Appendix A: certified TV reversal in a cone-preserving five-state generator.
- Theorem 10, Remark 11, Appendix B: for every N>=4 a crossing invisible to every proper induced subsystem; certified N=4 instance.
- Section 7.1: scope and limits.

## Evidence and replay

- Runner: [run_checks.py](run_checks.py) (standard library; exact comparison with retained outputs; use `--output`). Retained runner output: [reproduction/run_summary.json](reproduction/run_summary.json).
- Core checkers: [physical_frontier/verify_concave_kernel.py](physical_frontier/verify_concave_kernel.py) -> [concave_kernel_checks.json](physical_frontier/concave_kernel_checks.json); [physical_frontier/verify_metric_protocol.py](physical_frontier/verify_metric_protocol.py) -> [metric_protocol_checks.json](physical_frontier/metric_protocol_checks.json); [global_frontier/verify_hierarchy.py](global_frontier/verify_hierarchy.py) -> [hierarchy_verification.json](global_frontier/hierarchy_verification.json); [metric_sharpness/certify_counterexamples.py](metric_sharpness/certify_counterexamples.py) -> [exact_counterexamples.json](metric_sharpness/exact_counterexamples.json).
- Extension certificates: [extensions/e1_threshold_check.py](extensions/e1_threshold_check.py), [extensions/e2_superohmic_certificate.py](extensions/e2_superohmic_certificate.py), [extensions/e3_hierarchy_n4_certificate.py](extensions/e3_hierarchy_n4_certificate.py), [extensions/c1_prop5_certificate.py](extensions/c1_prop5_certificate.py) with their JSON outputs; [extensions/c1_make_table.py](extensions/c1_make_table.py) builds Table 4.
- Separate implementations: [adversarial_review/independent_checks.py](adversarial_review/independent_checks.py) (floating point, `--independent`) and [adversarial_review/independent_tv_certificate.py](adversarial_review/independent_tv_certificate.py).
- Proof notes: [physical_frontier/concave_kernel_theorem.md](physical_frontier/concave_kernel_theorem.md), [physical_frontier/finite_time_metric_theorem.md](physical_frontier/finite_time_metric_theorem.md), [global_frontier/hierarchy_theorem.md](global_frontier/hierarchy_theorem.md), [metric_sharpness/counterexamples.md](metric_sharpness/counterexamples.md) (development records; the manuscript is authoritative).
- Semantic negative controls: [negative_controls.py](negative_controls.py) (5 corruptions; refuses -O). Environment: [pyproject.toml](pyproject.toml), [uv.lock](uv.lock), [ENVIRONMENT.txt](ENVIRONMENT.txt). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml). Records: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). Integrity: [SHA256SUMS](SHA256SUMS).

## Trust boundary

All review was internal and AI-assisted. There was no peer review, no independent reproduction and no formal verification. Finite checks audit formulas and signs; they do not prove the all-N theorems. Priority holds only within a bounded search; Avitan et al. (2026) already exclude the slow-mode effect for Ohmic and sub-Ohmic spectra. Crossing locations from scans, slow-mode amplitudes and positive-bath Decimal illustrations are numerical only.
