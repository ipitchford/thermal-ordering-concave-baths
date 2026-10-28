#!/usr/bin/env python3
"""Reproduce the manuscript's finite checks in an isolated temporary tree.

This runner does not overwrite the retained research outputs. The default
checks (four core checks and four extension certificates) use only Python's
standard library and are compared exactly with the retained JSON. The
floating-point jobs (--independent, --numerical) are compared with tolerance
rel 1e-9 / abs 1e-12 on numbers, because round-off differs across platforms. --independent adds the separately
implemented audit, requiring NumPy and SciPy. --numerical adds the optional
high-precision positive-bath hierarchy illustrations. These finite checks do
not constitute a computer proof of the unrestricted analytical theorems.
"""
import argparse
import hashlib
import math
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent


def normalise(value):
    if isinstance(value, dict):
        return {k: normalise(v) for k, v in value.items()
                if k not in {'elapsed_seconds', 'python', 'platform', 'numpy', 'scipy'}}
    if isinstance(value, list):
        return [normalise(v) for v in value]
    return value


REL_TOL = 1e-9
ABS_TOL = 1e-12


def close(got, expected):
    """Tolerant structural comparison used only for floating-point audit jobs."""
    if isinstance(got, bool) or isinstance(expected, bool):
        return got == expected
    if isinstance(got, (int, float)) and isinstance(expected, (int, float)):
        if isinstance(got, int) and isinstance(expected, int):
            return got == expected
        return math.isclose(got, expected, rel_tol=REL_TOL, abs_tol=ABS_TOL)
    if isinstance(got, dict) and isinstance(expected, dict):
        return got.keys() == expected.keys() and all(close(got[k], expected[k]) for k in got)
    if isinstance(got, list) and isinstance(expected, list):
        return len(got) == len(expected) and all(close(a, b) for a, b in zip(got, expected))
    return got == expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'reproduction')
    parser.add_argument('--independent', action='store_true')
    parser.add_argument('--numerical', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    jobs = [
        ('concave_kernel', 'physical_frontier/verify_concave_kernel.py',
         'physical_frontier/concave_kernel_checks.json'),
        ('metric_protocol', 'physical_frontier/verify_metric_protocol.py',
         'physical_frontier/metric_protocol_checks.json'),
        ('hierarchy', 'global_frontier/verify_hierarchy.py',
         'global_frontier/hierarchy_verification.json'),
        ('metric_boundary', 'metric_sharpness/certify_counterexamples.py',
         'metric_sharpness/exact_counterexamples.json'),
        ('ext_bath_thresholds', 'extensions/e1_threshold_check.py',
         'extensions/e1_threshold_check.json'),
        ('ext_superohmic_kl_crossing', 'extensions/e2_superohmic_certificate.py',
         'extensions/e2_superohmic_certificate.json'),
        ('ext_hierarchy_n4_finite_bath', 'extensions/e3_hierarchy_n4_certificate.py',
         'extensions/e3_hierarchy_n4_certificate.json'),
        ('ext_tv_certificate_table', 'extensions/c1_prop5_certificate.py',
         'extensions/c1_prop5_certificate.json'),
    ]
    tolerant = set()
    if args.independent:
        jobs.append(('independent_audit', 'adversarial_review/independent_checks.py',
                     'adversarial_review/independent_check_results.json'))
        tolerant.add('independent_audit')
    if args.numerical:
        jobs.append(('positive_bath_numerics', 'global_frontier/check_positive_bath.py',
                     'global_frontier/positive_bath_numerics.json'))
        tolerant.add('positive_bath_numerics')
    records = []
    with tempfile.TemporaryDirectory(prefix='thermal-order-checks-') as tmp:
        work = Path(tmp)
        for folder in sorted({Path(s).parent for _, s, _ in jobs}):
            target = work / folder
            target.mkdir(parents=True)
            for source in (ROOT / folder).glob('*.py'):
                shutil.copy2(source, target / source.name)
        for name, script, result_path in jobs:
            started = time.perf_counter()
            env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                       PYTHONDONTWRITEBYTECODE='1')
            log = out / f'{name}.log'
            with log.open('w', encoding='utf-8') as stream:
                completed = subprocess.run(
                    [sys.executable, str(work / script)], cwd=work,
                    stdout=stream, stderr=subprocess.STDOUT, env=env,
                    timeout=900, check=False)
            record = {
                'check': name, 'script': script,
                'script_sha256': hashlib.sha256((ROOT / script).read_bytes()).hexdigest(),
                'exit_code': completed.returncode,
                'elapsed_seconds': time.perf_counter() - started,
            }
            generated = work / result_path
            if completed.returncode == 0 and generated.is_file():
                got = json.loads(generated.read_text())
                expected = json.loads((ROOT / result_path).read_text())
                if name in tolerant:
                    record['comparison'] = f'tolerant (rel {REL_TOL:g}, abs {ABS_TOL:g})'
                    match = close(normalise(got), normalise(expected))
                else:
                    record['comparison'] = 'exact'
                    match = normalise(got) == normalise(expected)
                record['matches_retained_output_except_runtime_metadata'] = match
                shutil.copy2(generated, out / f'{name}.json')
            else:
                record['matches_retained_output_except_runtime_metadata'] = False
            record['status'] = ('PASS' if completed.returncode == 0 and
                record['matches_retained_output_except_runtime_metadata'] else 'FAIL')
            records.append(record)
            print(f"{name}: {record['status']}", flush=True)
    summary = {'status': 'PASS' if all(r['status'] == 'PASS' for r in records) else 'FAIL',
               'scope': 'Finite formula and certificate reproduction, not an automated proof of all-N theorems.',
               'checks': records}
    (out / 'run_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f"Reports: {out}")
    return 0 if summary['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
