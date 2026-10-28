# Thermal ordering under Markovian relaxation: a finite-time criterion and its limits

**Unrefereed research candidate**, version 1.0.0-candidate (28 October 2026). DOI [10.5281/zenodo.23721564](https://doi.org/10.5281/zenodo.23721564); concept DOI 10.5281/zenodo.23721563. Release page: https://evidencepress.org/releases/thermal-ordering-concave-baths/

Creator: Anonymous. Text and data CC0-1.0; code MIT ([LICENSES.md](LICENSES.md)).

## The question in plain words

Put two copies of the same system, one hotter and one cooler, in contact with the same cold bath. Normally the hotter copy stays farther from equilibrium at every later time. The *Markovian Mpemba effect* is the surprising case in which the two relaxation curves cross, so the initially hotter copy ends up closer to equilibrium first. This package asks when such crossings can be ruled out **for all times at once**, without solving the dynamics, and where that guarantee stops.

## Results and essential limits

1. **Finite-time divergence criterion (Theorem 2).** If a likelihood-ratio Markov operator and its equilibrium adjoint preserve the arithmetic-prefix-mean cone, Gibbs cooling order is retained at every time for convex divergence generators with increasing concave derivative (both relative entropies, squared Hellinger, Rényi 0<α≤2). Heating uses a convex derivative. The conclusion includes common prescribed finite schedules of reversible stages, even when they do not commute (Corollary 3).
2. **Physical cone characterisation (Theorem 4).** For complete-graph bosonic population rates with common pair coupling, a continuous nonnegative spontaneous jump weight J preserves the cone for every finite spectrum and every positive bath exactly when it is concave. Necessity needs all N and arbitrarily cold baths jointly; at fixed N or fixed bath temperature concavity is not necessary. This is a statement about the cone, not about actual Mpemba crossings.
3. **Standard bath spectral densities (Proposition 6, Corollary 7).** With J the spontaneous emission weight: J = d^s is concave iff 0≤s≤1; d^s e^(-d/ω_c) iff 0<s≤1 and spectral width ≤ (s+√s)ω_c (Ohmic: ≤ 2ω_c); Drude–Lorentz iff width ≤ √3 ω_c; super-Ohmic never. The mapping to other papers' spectral-density conventions (for example Avitan et al.) is convention-dependent and not claimed beyond the stated observation.
4. **Certified super-Ohmic crossing (Proposition 8).** J = d³, E = (0,1,2,3), e^(-β_b) = 2/3: exact cone failure and a certified forward-KL crossing of two cooling preparations in t ∈ (1/5, 3/10). The crossing time, slow-mode amplitudes and strong-effect indication are numerical only.
5. **Total-variation boundary (Proposition 9).** A five-state positive-rate reversible cone-preserving generator has certified heating total-variation crossings in (10.46, 10.47) and (10.47, 10.48), with exact rational sign bounds (Table 4). The effect is about 1e-4 relative and lies outside the concave physical class.
6. **Induced-subsystem hierarchy (Theorem 10, Remark 11).** For every N≥4, a finite-bath N-state generator has a cooling slow-mode inversion and TV crossing although every proper induced subsystem excludes slow-mode inversion. The finite-bath threshold q_0(N) is existential; a certified N=4 instance exists at q = 1/1000 (and 1/100), but no explicit q_0(4). The construction needs β_b ≳ log(1/ε) and extreme energy/rate scales.

Not proved: total-variation order for the finite-bath concave physical class, maximality of the Rényi ranges, quantum coherences, sparse graphs (harmonic ladder, birth–death) and arbitrary pair couplings.

## Where to start

- Manuscript: [manuscript/thermal_order.pdf](manuscript/thermal_order.pdf) (27 pages) and its [LaTeX source](manuscript/thermal_order.tex).
- Claims: [CLAIMS.json](CLAIMS.json). Status: [STATUS.md](STATUS.md). Assurance limits: [ASSURANCE.md](ASSURANCE.md). Agent guide: [AI_INDEX.md](AI_INDEX.md).
- Review record without scores: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md), [REVISION_RESPONSE.md](REVISION_RESPONSE.md), [FINAL_FIXES.md](FINAL_FIXES.md). Prior art: [PRIOR_ART.md](PRIOR_ART.md), [CITATION_AUDIT.md](CITATION_AUDIT.md).
- Replay: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [negative_controls.py](negative_controls.py).

## Reproduce the exact checks

With Python 3.11 or later (standard library only), from any directory, in a copy of the package:

```bash
python3 run_checks.py --output /tmp/thermal-replay
python3 negative_controls.py          # do not use -O: the certificates use assert
sha256sum -c SHA256SUMS
```

Do not run the checks with `python -O` or `PYTHONOPTIMIZE`: the certificates use `assert`, and `extensions/e3_hierarchy_n4_certificate.py` records the rejected parameter values by catching assertion failures, so under `-O` it would report them as certified (the comparison with the retained output then fails, which is how the replay detects it).

`run_checks.py` copies the scripts into an isolated temporary tree, runs four core exact checks and four extension certificates, and compares each JSON result exactly with the retained output (elapsed-time and platform metadata excluded). Without `--output` it writes its fresh results to `reproduction/`, overwriting the shipped reproduction record.

| Check | Evidence |
|---|---|
| `physical_frontier/verify_concave_kernel.py` | 5,859 exact assertions across 96 rational models |
| `physical_frontier/verify_metric_protocol.py` | 6,493 exact assertions, including 12 nonreversible three-stage products |
| `global_frontier/verify_hierarchy.py` | All 8,086 proper subsets for N=4..12 at ε=1/100, plus 1,368 parameter variations |
| `metric_sharpness/certify_counterexamples.py` | Rational admissibility, cone signs and uniformisation bounds for the boundary examples |
| `extensions/e1_threshold_check.py` | Concavity thresholds of Proposition 6 |
| `extensions/e2_superohmic_certificate.py` | Proposition 8: exact generator, cone failure, rational KL enclosures (about 20 s) |
| `extensions/e3_hierarchy_n4_certificate.py` | Remark 11 at q=1/1000 and 1/100 (about 15 s) |
| `extensions/c1_prop5_certificate.py` | Table 4 intervals and sign bounds, 14 perturbed models (about 30 s) |

Optional floating-point audit and figures (NumPy, SciPy, Matplotlib, pinned in [uv.lock](uv.lock)):

```bash
uv run --frozen python3 run_checks.py --independent --output /tmp/thermal-replay-ind
uv run --frozen python3 adversarial_review/independent_tv_certificate.py   # writes its result beside the script
uv run --frozen python3 run_checks.py --numerical --output /tmp/thermal-replay-num   # Decimal illustrations, not certificates
```

These jobs compare with relative tolerance 1e-9 and absolute tolerance 1e-12. Figures: `manuscript/make_figures.py` and `manuscript/make_schematic.py` (illustrative, not data). Compile the manuscript from `manuscript/` with two runs of `pdflatex`.

These are finite checks, not an enumeration proof for unbounded dimension; the analytical arguments in the paper establish the unrestricted theorems.

## Package map

- `manuscript/`: paper source, PDF, figure scripts and figure data.
- `physical_frontier/`, `global_frontier/`, `metric_sharpness/`: development proof notes, exact checkers and retained outputs. They predate the final revision and some refer to development records that are not shipped (prior-art reports, target files); their substance is in [PRIOR_ART.md](PRIOR_ART.md).
- `extensions/`: certificates added in the 28 October revision.
- `adversarial_review/`: the independently implemented audit scripts and their results (the audit report itself is summarised in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)).
- `metric_probe/`: an unsuccessful numerical search for TV reversals in concave physical kernels. A negative search is not a no-crossing theorem.
- `reproduction/`: the retained output of `run_checks.py`.

## Cite

Anonymous (2026). *Thermal ordering under Markovian relaxation: a finite-time criterion and its limits*. Evidence Press research candidate, version 1.0.0-candidate. https://doi.org/10.5281/zenodo.23721564
