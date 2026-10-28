# Replay receipt

**Date:** 28 October 2026. **Input:** a fresh extraction of the frozen 1.0.0-candidate archive. **Runner:** the shared local research queue, macOS arm64, CPython 3.11 via uv.

| Run | Mode | Result |
|---|---|---|
| `run_checks.py --output <tmp>` (8 exact jobs) | normal, `-B` | all PASS, exact match with retained outputs |
| same, children with `PYTHONOPTIMIZE=1` | `-O -B` | 7 of 8 match exactly; `e3_hierarchy_n4_certificate` mismatches as expected, because its rejected q values are recorded by catching assertion failures. `-O` is unsupported; the exact comparison detects it. |
| `run_checks.py --independent` | normal, numpy/scipy from uv.lock | PASS within rel 1e-9 / abs 1e-12 |
| `negative_controls.py` | normal | 4 baselines pass; 5 corruptions rejected |
| `negative_controls.py` | `-O` | refused, as designed |

The publication gate repeats this replay from the frozen extraction; its log is retained privately with the publication receipts. Earlier, during review, every shipped script (default, independent, numerical, separate TV certificate) was replayed from a fresh extraction of the pre-revision bundle.
