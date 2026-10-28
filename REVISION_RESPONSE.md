# Revision response: thermal_order (round 1)

Revised bundle: `revised/thermal_order_research/` (copied from `in/thermal_order_research_bundle/thermal_order_research_2026-10-26/`; the input is untouched).
Manuscript: `manuscript/thermal_order.tex` / `.pdf`, 26 pages (was 16).

Numbering in the revised PDF:
- Lemma 1: Gibbs comparison.
- Theorem 2: two-adjoint ordering.
- Corollary 3: schedules.
- Theorem 4: concavity characterisation.
- Lemma 5: increments.
- Proposition 6: bath spectral densities.
- Corollary 7: physical reading.
- Proposition 8: super-Ohmic crossing.
- Proposition 9: TV reversal (was Prop 5).
- Theorem 10: hierarchy (was Thm 6).
- Remark 11: certified N=4 instance.
- Tables: Table 1 positioning, Table 2 divergences, Table 3 super-Ohmic certificate, Table 4 Prop 9 certificate.
- Appendices: A certificate, B proper subsystems, C derivations (new).

## Proof completeness

| Item | What changed and where |
|---|---|
| P1 | Proof of Lemma 1 (Sec. 2) rewritten in labelled steps. It first gives the derivative-as-covariance. It then derives Eq. (3) from the indicator expansion `E_i = E_1 + Σ_k ΔE_k 1{i>k}` and the event covariance identity. It shows `R_k ≥ μ_k ≥ L_k` (Eq. (4) and the tail bound). Strictness comes from the first index where b differs from b_1, at the cut `k=j−1`. |
| P2 | Sec. 3, new Eq. (8) after Eq. (7). It states `p^K ≥ 0`, `1ᵀp^K = πᵀKr = 1`, and `r=1 ⇒ p^K = π`. It also gives the population map `diag(π) K diag(π)^{-1}` and warns against confusing K with the population transition matrix. |
| P3 | Corollary 3 is now introduced as "a corollary of the closure of the cone-preserving class under products and adjoints", not a separate result. The proof gives the partial-stage operator `K(t) = e^{τG_j} P_{j−1}⋯P_1` and its adjoint (Eq. (14)), and notes π-self-adjointness of each stage. |
| P4 | Sec. 3.1 now states the operator-convex representation precisely (Eq. (18)): domain (0,∞), real analyticity (so f''(1) is finite), affine normalisation `f(1)=f'(1)=0`, endpoint atoms λ=0,1 identified, and `g_λ''(1)=2`, so `f''(1)=2` ⇔ probability measure. It cites Hiai et al. "and references therein", with no unverified theorem number. A fallback sentence says the conclusion holds for the listed families by direct differentiation alone. |
| P5 | Theorem 4 is renamed "Uniform-in-(N,β_b) cone characterisation". A new paragraph right after it says (ii) needs all N and arbitrarily cold baths jointly, and that concavity is not necessary at fixed N or β_b. It cites the 1+d² example: every triplet is fine; the exact six-level failure is at β_b = 10 log 2; a floating-point search found no failure for ≤10 levels at β_b = 1/4, 1/2, 1, with failures only for ≥12 clustered levels (labelled numerical). New Lemma 5 proves "continuous J concave ⇔ every increment x ↦ J(x+Δ)−J(x) is nonincreasing", via the chord-slope argument and midpoint concavity. The explicit choice `M = 1 + max{0, ⌊2(J(Δ)−I_Δ(y))/(I_Δ(y)−I_Δ(x))⌋}` (Eq. (30)) replaces "sufficiently large M". |
| P6 | Theorem 10(ii) now fixes the slow eigenfunction by `π_Sᵀf = 0`, `f_2 − f_1 = 1` (states reindexed by energy) and states `f ∈ int C`. A new sentence after the theorem says that q_0(N) exists but no explicit value is proved. Sec. 6.3 says "by finite-dimensional perturbation theory there exists q_0(N)" and "the argument is existential". Theorem 10(iii) now specifies the topology: an open set of independent downward rates in R_{>0}^{N(N−1)/2}, with upward rates set by detailed balance. Appendix B defines r at the start of each case; case B.4 says r includes the new ground. |
| P7 | New Appendix C derives: the quotient formula Eq. (21); H^S_ij, Eq. (23), by explicit pair expansion; H^W, Eq. (26), via `W u^{(j)}`; the necessity limit, Eq. (29); and the amplitude slope a'(0), Eq. (58), via `Cov(E,g⁰)=1/6` and `Cov(E,e_N)`. The a'(0) identity was checked symbolically (sympy) during revision. |

## Certificate disclosure

| Item | What changed and where |
|---|---|
| C1 | Proposition 9 statement: "two colder-than-bath Gibbs preparations undergoing heating"; both crossings are in rational brackets (10.46,10.47) and (10.47,10.48). New Table 4 (generated from `extensions/c1_prop5_certificate.json`) gives: all 10 component intervals at t* = 1047/100 (common digits, widths ≤ 5.02e−36); the exact rational lower and upper bounds for TV(p3) − TV(p4); and the certified signs −,+,− at t = 10.46, 10.47, 10.48. The text after Table 4 adds: the 60-digit scan (two crossings only, window 0.0144, labelled numerical); the effect size (6.0e−8, about 1e−4 relative); the initial gap; exact η-independence of the window; robustness under 14 perturbed models, each certified at t*; and the boundary-Metzler note (exact zero entries). The mechanism (the TV indicator leaves C) is added. Appendix A points to the extension script. |
| C2 | `run_checks.py`: a new `close()` compares numbers with rel 1e−9 / abs 1e−12, used only for the `independent_audit` and `positive_bath_numerics` jobs. Exact jobs still use `==`. Each record now has a `comparison` field. The default set adds four stdlib extension certificates (`extensions/`). The timeout is raised to 900 s. |

## Prior art and positioning

| Item | What changed and where |
|---|---|
| L1 | Intro §1.2, second paragraph, adapted from `refs/prior_art_paragraphs.tex` (2) and made more cautious. Avitan et al. already exclude the **slow-mode** effect for s ≤ 2 (their parametrisation). Their N-level result is proved only at zero bath with strong triplet violation, is conjectural (numerics N ≤ 7) at positive bath, and is analytic only in the high/low-temperature limits. The text says "We do not claim that their class coincides with our concave class". It lists the increments: all β_b/N/spectra, all-time order, the divergence family and schedules, and exact necessity, which concerns the cone. A new Table 1 row covers Avitan. Convention check (Sec. 4.5): the rates recorded from their v1 are Δ^{s−1}e^{−Δ/Δ_c}(n+1), so their s corresponds to the exponent s−1 of our J. Their Ohmic s=1 corresponds to J ∝ e^{−d/ω_c}, which is not concave, so the classes differ. E2 (J = d³) corresponds to their s = 4, outside their s ≤ 2, so there is **no tension**. Beyond this, no mapping is claimed. The Scope subsection repeats the disclaimer. |
| L2 | Sec. 2 has a new stochastic-order paragraph: the likelihood-ratio order of the Gibbs family, increasing-in-mean (Chebyshev/Biernacki) sequences, and Daley 1968, Keilson–Kester 1977, Massey 1987 and Shaked–Shanthikumar 2007. Daley is also cited at the stochastic-monotonicity example (Eq. (32)), and Keilson–Kester and Massey at the Metzler criterion (Sec. 4.1). |
| L3 | Intro §1.2 has a chronological paragraph covering Lu–Raz, Klich–Raz–Hirschberg–Vucelja (weak/strong, index; also cited in the weak/strong definitions), Kumar–Bechhoefer, Walker–Vucelja, Klich–Vucelja, Ohga, Vu–Hayakawa, Summer, Avitan and Gelbwaser. A new paragraph after Table 2 covers the heating/cooling asymmetry (Lapolla–Godec with erratum, Van Vu–Hasegawa, Ibáñez et al.). Klich2019 is replaced by the J. Phys. A version of record. Hiai2017 is updated to Rev. Math. Phys. (2011), and "T. Van Vu" now matches the record. All new bibitems were copied verbatim from `refs/new_bibitems.tex`; only a line-break hint (`\allowbreak`) was added in the Klich title. |
| L4 | Novelty statements now use exact combinations and quantifiers (intro increments list, Table 1). The schedule result is demoted to a corollary of closure. The Scope subsection separates three things: the absence of an anomalous induced subsystem (proved); the three-level necessary conditions of Avitan (not refuted); and global inference from pair data (not prevented). The four-block lumpability remains in §6.4 and is repeated in Scope. |

## Physical relevance and extensions

| Item | What changed and where |
|---|---|
| E1 | New §4.5 "Standard bath spectral densities". It states the convention (J = 2π\|g\|²𝒥, Bose factor explicit, ω^k caveat). Proposition 6 gives the pure power law, exponential cutoff with threshold (s+√s)ω_c (Ohmic 2ω_c; s=0 and s>1 never), Drude–Lorentz √3 ω_c, and Gaussian √(3/2) ω_c, each with a full proof. Corollary 7 states that sufficiency is universal and necessity only existential. The old one-paragraph d^α remark is replaced. The exact audit is `extensions/e1_threshold_check.py`. |
| E2 | New §4.6, Proposition 8: J = d³, E = (0,1,2,3), e^{−β_b} = 2/3. It gives the exact G and H (all off-diagonal entries of H negative, so the cone fails) and a certified forward-KL crossing in (1/5, 3/10), with Table 3. A paragraph "What is certified and what is numerical" separates the exact items from the float items (crossing time ≈ 0.2423, slow amplitudes, strong-effect indication). It also notes that the Ohmic contrast on the same spectrum is Metzler. Reproduced by `extensions/e2_superohmic_certificate.py`. |
| E3 | New Remark 11 in §6.3. It gives a certified N = 4 instance at q = 1/1000, also passing at q = 1/100 on a knife edge. It states that **no explicit q_0(4)** is proved and no interval of q is certified. It gives the binding triple formula (Eq. (60)) and q* = ε + O(ε^{L−1}), and says the construction needs β_b ≳ log(1/ε). The q = 1/50 local anomaly is labelled numerical. The general finite-bath statement in Theorem 10 stays existential. Reproduced by `extensions/e3_hierarchy_n4_certificate.py`. |
| E4 | The "Model class" bullet in §7.1 Scope says the complete graph with common coupling excludes the harmonic ladder and birth–death chains, and that this is open. |

## Readability

| Item | What changed and where |
|---|---|
| R1 | The abstract is rewritten plain-first: the two-copies picture, the Mpemba effect in one sentence, then the results and limits in words. |
| R2 | §1.1 is a concrete picture of two crossing relaxation curves, using the certified numbers from Prop 8. It defines weak (slow-mode), strong and all-time order, with KlichRaz2019. §1.2 is the chronological literature and the Avitan positioning. §1.3 lists the "Main results in plain words" (5 items). §1.4 is a roadmap. |
| R3 | §2 has a new paragraph "Distances, modes and subsystems" that defines TV, slow mode/amplitude, weak/strong and induced subsystem before first use. Double notation is resolved in four places: prefix means are now μ_k; the emission part is now G = S + W (H^W) instead of A; the increment is now I_Δ with step Δ instead of A, h; and in §6 the slow-eigenfunction top coordinate is now θ_r, θ_m (η stays the reset rate in Prop 9) and the amplitude is a(β) instead of A(β). |
| R4 | One "Scope and limits" subsection (§7.1) collects the disclaimers. The self-instruction sentence is deleted. The checker counts move to "Reproducibility and provenance". The AI-assistance statement is on the title page (author footnote, "Anonymous") and in the reproducibility section. The Table 2 row for total variation now reads "Not guaranteed (Proposition 9)". Figure 1 axes now read "Relative entropy to equilibrium (nats)" and "Preparation inverse temperature, β/β_b". The normalisation caveat is added to the caption and to the text. A one-sentence "In plain words" takeaway follows Lemma 1, Theorems 2, 4, 10, Corollary 7, and Propositions 8, 9. |

## Bundle changes
- New `extensions/` contains `e1_threshold_check.py`, `e2_superohmic_certificate.py`, `e3_hierarchy_n4_certificate.py`, `c1_prop5_certificate.py`, `c1_make_table.py`, the retained JSON outputs and `c1_prop5_table.tex`. These are copied from `p3-thermal/ext/`.
- `run_checks.py`: four new default exact jobs and tolerant comparison for the float jobs.
- `reproduction/` is regenerated by the default run with 8 checks.
- `README.md`: the results list (6 items with limits), the check table (8 rows plus the tolerance note) and the "Not proved" line are updated.
- `manuscript/make_figures.py`: axis labels only. The figures are regenerated and the curve CSV is kept, because it was identical up to round-off.
- `manuscript/PDF_QA.md` and `pdf_qa_checks.json` are updated.
- `MANIFEST.sha256` and `PACKAGE_CONTENTS.json` are regenerated with the same scheme: every file except the two indexes, plus PACKAGE_CONTENTS in the manifest, sorted by path.
- `academic_review/` and `adversarial_review/` are unchanged.

## Checks (run from outside the bundle, absolute path)
- `uv run python3 -I -B .../run_checks.py`: PASS. All 8 checks pass: 4 core plus 4 extensions, all exact.
- `uv run --with numpy --with scipy python3 -I -B .../run_checks.py --independent`: PASS. All 9 pass; the independent audit uses the tolerant comparison.
- `--numerical` was also run: PASS (positive_bath_numerics, tolerant comparison).
- The PDF builds with pdfLaTeX with 0 warnings, 0 overfull boxes and 0 undefined references. All 26 pages were rendered and inspected.

## Not done / caveats
- No explicit q_0(N), including q_0(4), was proved. The Remark certifies only two q values.
- The fixed-(N,β_b) threshold N*(β_b) for 1+d² is reported only as a numerical search result from the devil's-advocate review; it is not certified.
- The Teza et al. Phys. Rep. review was not added: its record was not fetched, per REFS_AUDIT §D.
- Theorem numbers for the operator-convex representation in Hiai et al. were not cited, because they could not be verified.

# Final fixes (FINAL_FIXES.md, F1-F9)

| Item | What changed and where |
|---|---|
| F1 | Remark 11: "-3.91 (an exact root)" is removed. The text now says each eigenvalue lies in a certified rational bracket, that the first bracket ends at -3.91, which is not an exact root, and that the spectrum is simple because the brackets are disjoint. |
| F2 | Proposition 9 proof: new Eq. (45) gives the exact rational sign bounds at t=10.46 (about -5.845e-9) and t=10.48 (about -7.871e-10). It states that the positive sets of the two preparations differ at 10.46 ({1,2,4} for p3 vs {1,2} for p4); at 10.48 the script certifies {1,2,3,4} for both. Each TV is bounded on its own certified set. Appendix A explains the per-preparation version. `extensions/c1_prop5_certificate.py` now writes outward-rounded bounds, checked against exact Fractions, and the per-preparation positive sets for all three times. Note: the proof audit listed p4's set at 10.48 as {1,2,4}; the certified script output is {1,2,3,4}, and the printed sign bound is unaffected. |
| F3 | Table 3 has a new column with the certified rational bounds on KL(hot)-KL(cold) at all six times. These come from the new `hot_minus_cold_outward` field in `extensions/e2_superohmic_certificate.py`. The slow-mode normalisation used for the numerical amplitudes (pi^T f = 0, ||f||_2 = 1, f_2 - f_1 > 0) is now stated in §4.6. |
| F4 | Proposition 6(a): "with no cutoff factor on [0,D] (for example a hard cutoff above D; a smooth cutoff changes the curvature...)". |
| F5 | The Corollary 7 takeaway is narrowed to cooling Gibbs preparations in the complete-graph, common-coupling model. It now uses the full spectral width (top minus bottom) for both the Ohmic and the Drude-Lorentz thresholds. |
| F6 | §4.5 quotes Avitan et al. Eq. (5) (arXiv v1) exactly: a_kj = f(Delta) ∝ Delta^{s-1} e^{-Delta/Delta_c}[1/(e^{beta_b Delta}-1)+1]. Matching this to the downward rate J(n_b+1) gives J ∝ Delta^{s-1}e^{-Delta/Delta_c}, so their s equals our exponent s' plus one. Their Ohmic case and the s=4 consequences are kept. |
| F7 | §6.3 continuation: the text now defines mu_N as the minimum of all proper-subsystem slow margins and eigenvalue gaps at q=0, gives a per-subsystem q_S from continuity, and sets q_0(N) = min_S q_S. |
| F8 | "Physical dividing line" / "physically meaningful" are removed from the Prop 8 takeaway and from Sec. 7. The text now says concavity characterises the universal cone criterion, not the existence of crossings in general. |
| F9 | New Figure 1 in §1.1 is an illustrative schematic of two crossing relaxation curves (sums of exponentials, labelled "Schematic only"). It is generated by the new `manuscript/make_schematic.py`. Later figures are renumbered. |

The PDF has 27 pages and builds with zero warnings or overfull boxes. Pages 2, 12, 14, 16 and 20 (the changed pages) were inspected. README, MANIFEST.sha256 and PACKAGE_CONTENTS.json were regenerated. The extension JSONs were regenerated, and run_checks.py compares them exactly.
