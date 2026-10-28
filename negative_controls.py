#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Semantic negative controls for the exact certificates.

Each control copies one certificate script into a temporary directory, applies a
single claim-relevant corruption (for example, replacing the super-Ohmic cube law
by a concave weight), runs it, and requires the corrupted run to be REJECTED:
either an assertion fails (non-zero exit) or its output no longer matches the
retained certificate. The unmodified baseline must pass first. Assertions are
load-bearing, so this script refuses to run under python -O.
"""
import json, os, shutil, subprocess, sys, tempfile
from pathlib import Path
if not __debug__:
    sys.exit('negative_controls.py must not run under -O (the certificates use assert)')
ROOT = Path(__file__).resolve().parent
SKIP = {'elapsed_seconds', 'python', 'platform', 'numpy', 'scipy'}
def norm(v):
    if isinstance(v, dict): return {k: norm(x) for k, x in v.items() if k not in SKIP}
    if isinstance(v, list): return [norm(x) for x in v]
    return v
CONTROLS = [
 ('super-Ohmic weight replaced by concave J(d)=d (cone failure must disappear)', 'extensions/e2_superohmic_certificate.py', 'extensions/e2_superohmic_certificate.json',
  "J3 = lambda d: F(d)**3", "J3 = lambda d: F(d)"),
 ('hot and cold preparations swapped in Proposition 8', 'extensions/e2_superohmic_certificate.py', 'extensions/e2_superohmic_certificate.json',
  "u_hot, u_cold = F(9, 10), F(4, 5)", "u_hot, u_cold = F(4, 5), F(9, 10)"),
 ('concave test weight made convex in the cone checker', 'physical_frontier/verify_concave_kernel.py', 'physical_frontier/concave_kernel_checks.json',
  "'concave_quadratic': lambda d: F(d)-F(d*d,2*D)", "'concave_quadratic': lambda d: F(d)+F(d*d,2*D)"),
 ('Proposition 9 weight w5 perturbed far outside the robust range', 'metric_sharpness/certify_counterexamples.py', 'metric_sharpness/exact_counterexamples.json',
  "w5=list(map(F,[60,14,2,1]))+[F(5763,100000)]", "w5=list(map(F,[60,14,2,1]))+[F(1,2)]"),
 ('Table 4 certificate time moved outside the reversal window', 'extensions/c1_prop5_certificate.py', 'extensions/c1_prop5_certificate.json',
  "W5, ETA, TS = F(5763, 100000), F(1, 100), F(1047, 100)", "W5, ETA, TS = F(5763, 100000), F(1, 100), F(1000, 100)"),
]
def run(script, expected, old=None, new=None):
    with tempfile.TemporaryDirectory(prefix='neg-') as td:
        td = Path(td); folder = td / Path(script).parent; folder.mkdir(parents=True)
        for p in (ROOT / Path(script).parent).glob('*.py'): shutil.copy2(p, folder / p.name)
        target = td / script
        if old is not None:
            text = target.read_text()
            if text.count(old) != 1: sys.exit(f'control anchor not found exactly once in {script}: {old}')
            target.write_text(text.replace(old, new))
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        r = subprocess.run([sys.executable, '-B', str(target)], cwd=td, capture_output=True, text=True, env=env, timeout=1800)
        out = td / expected
        same = r.returncode == 0 and out.is_file() and norm(json.loads(out.read_text())) == norm(json.loads((ROOT / expected).read_text()))
        return r.returncode, same
failed = False
for script, expected in sorted({(c[1], c[2]) for c in CONTROLS}):
    code, same = run(script, expected)
    print(f'baseline {script}: exit {code}, matches retained output: {same}')
    if not same: failed = True
for name, script, expected, old, new in CONTROLS:
    code, same = run(script, expected, old, new)
    verdict = 'REJECTED' if not same else 'NOT REJECTED'
    how = 'assertion failure' if code else ('output differs from retained certificate' if not same else 'accepted')
    print(f'{verdict}: {name} ({how})')
    if same: failed = True
print('FAIL' if failed else f'PASS: baselines pass; all {len(CONTROLS)} corruptions rejected')
sys.exit(1 if failed else 0)
