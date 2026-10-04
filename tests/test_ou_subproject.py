"""Verify the subproject's public entry point and its output accounting."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys


def test_cli_runs_outside_repo_and_records_every_attempt(tmp_path):
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root/'subprojects/ou-calibration/config.json').read_text())
    config.update(replications=20, training_durations=[10, 30])
    config_path = tmp_path/'config.json'
    config_path.write_text(json.dumps(config))
    output = tmp_path/'results'
    env = {**os.environ, 'MPLCONFIGDIR': str(tmp_path/'mpl')}
    subprocess.run([sys.executable, str(root/'subprojects/ou-calibration/run.py'),
                    '--config', str(config_path), '--output', str(output)],
                   cwd=tmp_path, env=env, check=True, capture_output=True, text=True)
    summary = json.loads((output/'summary.json').read_text())
    assert len(summary['scenarios']) == 4
    with (output/'replications.csv').open() as file:
        rows = list(csv.DictReader(file))
    assert len(rows) == 80
    for scenario in summary['scenarios']:
        relevant = [r for r in rows if r['scenario'] == scenario['name']]
        assert len(relevant) == scenario['attempted'] == 20
        assert sum(r['success'] == 'True' for r in relevant) == scenario['successful']
        assert scenario['successful'] + scenario['failures'] == 20
        for forecast in scenario['forecasts']:
            assert all(forecast[key]['n'] == scenario['successful'] for key in
                       ('plugin', 'oracle', 'shift_plugin', 'shift_oracle'))
            assert forecast['oracle']['coverage'] == forecast['shift_oracle']['coverage']
    assert (output/'report.md').exists()
    assert len(list((output/'figures').glob('*.png'))) >= 3
