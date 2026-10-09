"""Ensure the documented verification CLI executes outside the repository."""
import json
import os
from pathlib import Path
import subprocess
import sys


def test_verification_cli_records_independent_benchmarks(tmp_path):
    script = Path(__file__).resolve().parents[1] / 'subprojects/ou-calibration/verify.py'
    output = tmp_path / 'results'
    subprocess.run(
        [sys.executable, str(script), '--paths', '5000', '--output', str(output)],
        cwd=tmp_path, env={**os.environ, 'MPLCONFIGDIR': str(tmp_path / 'mpl')},
        check=True, capture_output=True, text=True, timeout=120,
    )
    result = json.loads((output / 'results.json').read_text())
    assert result['all_checks_passed']
    assert len(result['checks']) == 57
    assert len(result['variance_quadrature']) == 9
    assert len(result['likelihood']) == 8
    assert all(len(row['optimizer_runs']) == 2 for row in result['likelihood'])
    assert all(0 <= row['coverage_hits'] <= 5000 for row in result['moments_and_forecasts'])
    assert (output / 'forecast_mse_verification.png').read_bytes().startswith(b'\x89PNG')
    assert '57/57' in (output / 'report.md').read_text()
