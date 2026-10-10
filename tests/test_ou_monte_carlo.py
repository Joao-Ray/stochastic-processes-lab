"""Independent references for stationary rank calibration and cell inversion."""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys

import numpy as np
import pytest

from src.ou_monte_carlo import (
    PhiCell, cell_at, cell_tv_bound, invert_cells, phi_cells, profile_statistic,
    rank_pvalue, stationary_kl, stationary_null_paths, stationary_rank_test,
    test_cell as evaluate_cell,
)


def test_conservative_rank_enumeration_with_ties():
    # Exchangeability: every member is equally likely to be observed.
    for values in [np.arange(20.), np.array([0.]*7+[1.]*6+[2.]*7)]:
        ps = [rank_pvalue(v, np.delete(values, i))[0] for i, v in enumerate(values)]
        for t in np.linspace(0, 1, 101):
            assert np.mean(np.array(ps) <= t) <= np.floor(20*t)/20+1e-15
    assert rank_pvalue(2., [2., 2., 1.]) == (.75, 2)


def test_stationary_paths_against_direct_recursion():
    phi, n, count, seed = .72, 12, 5, 907
    rng = np.random.default_rng(seed)
    expected = np.empty((count, n+1))
    expected[:, 0] = rng.standard_normal(count)
    eps = rng.standard_normal((count, n))
    for k in range(n):
        expected[:, k+1] = phi*expected[:, k]+np.sqrt(1-phi**2)*eps[:, k]
    actual = stationary_null_paths(phi, n, count, np.random.default_rng(seed))
    np.testing.assert_allclose(actual, expected, atol=1e-14)


def test_profile_against_independent_least_squares_and_affine_invariance():
    x = stationary_null_paths(.8, 40, 1, np.random.default_rng(901))[0]
    candidate = .6
    reg = np.column_stack([np.ones(len(x)-1), x[:-1]])
    coef = np.linalg.lstsq(reg, x[1:], rcond=None)[0]
    best = np.clip(coef[1], 0, 1)
    def residual_sse(slope):
        intercept = np.mean(x[1:]-slope*x[:-1])
        return np.sum((x[1:]-intercept-slope*x[:-1])**2)
    expected = (len(x)-1)*np.log(residual_sse(candidate)/residual_sse(best))
    assert profile_statistic(x, candidate) == pytest.approx(expected, abs=1e-12)
    for shifted in [4+3*x, -9-2*x]:
        assert profile_statistic(shifted, candidate) == pytest.approx(expected, abs=1e-11)
        original = stationary_rank_test(x, candidate, np.random.default_rng(4), 99)
        transformed = stationary_rank_test(shifted, candidate, np.random.default_rng(4), 99)
        assert original['pvalue'] == transformed['pvalue']


@pytest.mark.parametrize('phi,reference,n', [(0., .3, 3), (.8, .2, 10), (.2, .8, 10), (.999, .99, 20), (.9+1e-9, .9, 10)])
def test_path_kl_against_full_gaussian_covariance(phi, reference, n):
    distance = np.abs(np.arange(n+1)[:, None]-np.arange(n+1))
    p, q = phi**distance, reference**distance
    expected = .5*(np.trace(np.linalg.solve(q, p))-(n+1)+np.linalg.slogdet(q)[1]-np.linalg.slogdet(p)[1])
    assert stationary_kl(phi, reference, n) == pytest.approx(expected, rel=2e-9, abs=2e-10)
    assert stationary_kl(reference, reference, n) == 0.


def test_kl_endpoint_supremum_and_numerical_stability():
    low, high, center, n = .81, .85, .83, 18
    bound = cell_tv_bound(low, high, center, n)
    for phi in np.linspace(low, high, 51):
        assert np.sqrt(stationary_kl(float(phi), center, n)/2) <= bound
    # The near-equality quadratic limit is n*(1+psi^2)/(2*(1-psi^2)^2).
    psi, h = .95, 1e-9
    delta = (psi+h)-psi
    approximation = n*(1+psi**2)*delta**2/(2*(1-psi**2)**2)
    assert stationary_kl(psi+h, psi, n) == pytest.approx(approximation, rel=1e-6)
    assert cell_tv_bound(.99, 1., .995, n) == 1.


def test_grid_full_coverage_and_untested_tail():
    cells = phi_cells(4, resolution=.2, u_max=2.)
    assert cells[0].lower == cells[0].center == 0.
    assert cells[-1].upper == cells[-1].tv_bound == 1.
    assert all(a.upper == b.lower for a, b in zip(cells, cells[1:]))
    for phi in np.linspace(0, np.nextafter(1., 0.), 100):
        _, cell = cell_at(float(phi), cells)
        assert cell.lower <= phi < cell.upper
    series = [.2, -.4, .7, .1, -.2]
    result = invert_cells(series, .5, 42, 19, resolution=.2, u_max=2.)
    assert result['decisions'][-1]['accepted'] is True
    assert result['decisions'][-1]['simulated'] is False
    assert result['hull']['lower'] == 0.
    assert result['components'][0]['lower'] == 0.
    # Components must represent all accepted cells without filling rejected gaps.
    for decision in result['decisions']:
        middle = (decision['lower']+decision['upper'])/2
        theta = -np.log(middle)/.5
        membership = any(c['lower'] <= theta and (c['upper'] is None or theta <= c['upper']) for c in result['components'])
        assert membership == decision['accepted']


def test_api_rejects_invalid_or_underbounded_cells():
    series = [1., -.4, .3, .7, -.1]
    with pytest.raises(ValueError): stationary_rank_test(series, 1., np.random.default_rng(1), 19)
    with pytest.raises(ValueError): stationary_null_paths(.5, 2, 1, np.random.default_rng(1))
    with pytest.raises(ValueError): rank_pvalue(np.nan, [1.])
    with pytest.raises(ValueError): phi_cells(3, u_max=100.)
    with pytest.raises(ValueError): evaluate_cell(series, PhiCell(.2, .3, .25, 0.), np.random.default_rng(1), 19)
    with pytest.raises(ValueError): evaluate_cell(series, PhiCell(.9, 1., .9, 1.), np.random.default_rng(1), 1)
    with pytest.raises(ValueError): invert_cells(series, 0., 1)
    with pytest.raises(ValueError): stationary_rank_test([1.]*5, .5, np.random.default_rng(1), 19)


def test_finite_sample_cli_replay_and_archive_denominators(tmp_path):
    root = Path(__file__).resolve().parents[1]
    folder = root/'subprojects/ou-inference-study/finite-sample'
    config = json.loads((folder/'ci-config.json').read_text())
    config.update(replications=8, simulations=19)
    cfg = tmp_path/'config.json'
    cfg.write_text(json.dumps(config))
    destinations = [tmp_path/'one', tmp_path/'two']
    for output in destinations:
        subprocess.run([sys.executable, str(folder/'run.py'), '--config', str(cfg), '--output', str(output)],
            cwd=tmp_path, env={**os.environ, 'MPLCONFIGDIR': str(tmp_path/'mpl')},
            check=True, capture_output=True, text=True, timeout=120)
    first = destinations[0]
    archive = first/'replications.csv.gz'
    assert archive.read_bytes() == (destinations[1]/archive.name).read_bytes()
    assert (first/'examples.json').read_bytes() == (destinations[1]/'examples.json').read_bytes()
    result = json.loads((first/'summary.json').read_text())
    assert result['record_archive_sha256'] == hashlib.sha256(archive.read_bytes()).hexdigest()
    for name, digest in result['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest() == digest
    with gzip.open(archive, 'rt') as file:
        rows = list(csv.DictReader(file))
    assert len(rows) == len(config['scenarios'])*8
    for scenario in result['scenarios']:
        selected = [row for row in rows if row['scenario'] == scenario['id']]
        for method in ('rank', 'cell', 'profile'):
            metric = scenario[method+'_coverage']
            assert metric['n'] == len(selected) == 8
            assert metric['value'] == sum(row[method+'_hit'] == 'True' for row in selected)/8
        for row in selected:
            assert float(row['rank_pvalue'])*20 == pytest.approx(1+int(row['rank_exceedances']))
    assert (first/'coverage.pdf').read_bytes().startswith(b'%PDF')
    examples = json.loads((first/'examples.json').read_text())
    assert len(examples) == 2
    assert all(e['hull']['lower'] == 0. for e in examples)
    assert (first/'accepted_sets.png').exists()
