"""Audit experiment denominators and replay from the public CLI."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


def test_study_cli_preserves_attempts_and_reproducible_records(tmp_path):
    root = Path(__file__).resolve().parents[1]
    script = root/'subprojects/ou-inference-study/run.py'
    config = json.loads((root/'subprojects/ou-inference-study/ci-config.json').read_text())
    config.update(replications=8, bootstrap_replications=19)
    config_path = tmp_path/'config.json'
    config_path.write_text(json.dumps(config))
    outputs = [tmp_path/'first', tmp_path/'second']
    for output in outputs:
        subprocess.run([sys.executable, str(script), '--config', str(config_path), '--output', str(output)],
            cwd=tmp_path, env={**os.environ, 'MPLCONFIGDIR': str(tmp_path/'mpl')},
            check=True, capture_output=True, text=True, timeout=120)
    archive = outputs[0]/'replications.csv.gz'
    assert archive.read_bytes() == (outputs[1]/archive.name).read_bytes()
    summary = json.loads((outputs[0]/'summary.json').read_text())
    assert summary['record_archive_sha256'] == hashlib.sha256(archive.read_bytes()).hexdigest()
    with gzip.open(archive, 'rt') as file:
        rows = list(csv.DictReader(file))
    assert len(rows) == 40
    for scenario in summary['scenarios']:
        relevant = [r for r in rows if r['scenario'] == scenario['id']]
        assert len(relevant) == scenario['attempted'] == 8
        for method, metrics in scenario['methods'].items():
            available = [r for r in relevant if r[method+'_available'] == 'True']
            hits = sum(r[method+'_hit'] == 'True' for r in relevant)
            assert len(available) == metrics['available']
            assert metrics['available']+metrics['unavailable'] == 8
            assert metrics['success_and_coverage']['coverage'] == pytest.approx(hits/8)
            if available:
                assert metrics['conditional']['coverage'] == pytest.approx(hits/len(available))
    assert (outputs[0]/'coverage_and_noise.pdf').read_bytes().startswith(b'%PDF')
    assert (outputs[0]/'report.md').exists()
