# Verification and scope

The result is a candidate theorem requiring independent review. Its unrestricted claims rest on the proof in `hierarchy_theorem.md`; the finite checks below corroborate formulas and boundary cases.

## Reproduce

Both programs use Python's standard library only:

```bash
python3 verify_hierarchy.py
python3 check_positive_bath.py
```

They write the corresponding JSON and print matching results. The retained `.log` files are the recorded outputs from the current run.

## Exact arithmetic

`verify_hierarchy.py` reconstructs every downward generator and every induced restriction, checks strictly ordered escape rates, solves the slow backward eigenfunction by rational back-substitution and checks all prefix margins. It checks exact positive forcing coefficients or positive endpoint values in the all-time survival argument, rather than sampling times.

- N=4,...,12, epsilon=delta=1/100: all **8,086** proper subsystems of at least two states passed.
- N=4,...,8, epsilon=delta in {1/50,1/20,1/15}: all **1,368** additional proper subsystems passed.
- Every full system has a negative final prefix margin and positive earlier margins.
- The explicit uniform-temperature derivative formula matches direct rational covariance computation.

These are coverage counts, not confidence scores or an enumeration proof for all N.

## Positive bath illustrations

`check_positive_bath.py` uses 80-digit Decimal arithmetic with an expanded exponent range. This retains strictly positive upward rates even when their exponents are far beyond ordinary binary floating-point range. It checks N=4,...,8, bath inverse temperature 20, fixed downward rates, the shared displayed energy spectrum, and every proper induced restriction.

The script computes the slow branch through a quotient eigenvector iteration, checks its residual, and uses the exact zero-bath eigenvector condition number with a numerical Bauer–Fike perturbation estimate to separate its spectral neighbourhood from all other modes. Every radius is less than one third of the minimum zero-bath spectral spacing. This is a numerical stability check, not directed-rounding interval certification.

The program evaluates the two full-system preparations beta_h=1/(100L), beta_c=2/(100L). Their initial TV ordering is conventional, their positive slow amplitudes are inverted, and a positive time-20 TV reversal appears in a separate positive Poisson uniformization calculation. The Poisson tail estimate is recorded, while Decimal rounding itself is not interval-enclosed.

All **456** proper subsystems in the five displayed positive-bath cases retain positive slow prefix margins. The general theorem's positive-bath existence and simultaneous proper-subsystem exclusion follow analytically from strict inequalities and continuity, rather than extrapolation from these cases.

## Limits

The zero-bath process has a four-block lumpable survival dynamics. Its arbitrary minimal induced-level order does not demonstrate an unbounded number of dynamically essential modes. The numerical energy ratios and reverse rates are extreme and have no claimed experimental realisation. There is no strong-mode cancellation or large practical speedup in this family. At positive bath temperature only slow-amplitude exclusion has been proved for proper subsystems, not absence of every transient crossing.

The initial attempt to use mpmath found that dependency unavailable. The completed implementation instead uses Decimal, Gaussian elimination and uniformization from the standard library. No package was installed and no dependency fallback is hidden.

## Research artefacts

- `target.yaml`: frozen theorem, scope and normalisations.
- `fingerprints.json`: exact small matrices and invariants.
- `novelty_report.md`: thirty-query bounded prior-art gate and structural rerun after lumpability.
- `hierarchy_theorem.md`: complete proof and significance limits.
- `verify_hierarchy.py`, `hierarchy_verification.json`, `hierarchy_verification.log`: exact finite checks.
- `check_positive_bath.py`, `positive_bath_numerics.json`, `positive_bath_numerics.log`: finite-bath numerical illustrations.
