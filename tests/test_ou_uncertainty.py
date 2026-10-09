"""Independent profile optimization, recursion, and boundary benchmarks."""
import numpy as np
import pytest
from scipy.optimize import minimize
from scipy.stats import chi2

from src.ou_inference import OUFitError, conditional_ou_nll, fit_ou
from src.ou_process import simulate_ou_exact
from src.ou_uncertainty import (bootstrap_ou, profile_theta_interval,
    theta_interval_from_phi, wald_theta_interval)


def test_profile_endpoints_match_independent_nuisance_optimization():
    dt = .2
    series = simulate_ou_exact(.7, 1.5, .8, 1.5, dt, 500, 1, np.random.default_rng(517))[0]
    fitted = fit_ou(series, dt)
    interval = profile_theta_interval(series, dt)
    assert 0 < interval.lower < fitted.theta < interval.upper < np.inf
    for theta in [interval.lower, interval.upper]:
        reference = minimize(lambda p: conditional_ou_nll(series, dt, theta, p[0], np.exp(p[1])),
                             [1., np.log(.6)], method='BFGS')
        assert 2*(reference.fun-fitted.negative_log_likelihood) == pytest.approx(chi2.ppf(.95, 1), abs=1e-6)


@pytest.mark.parametrize('series', [[0, 1, 3, 6, 11, 20], [1, -.3, .7, -.8, .9, -.7]])
def test_profile_boundary_supremum_does_not_create_a_point_fit(series):
    with pytest.raises(OUFitError):
        fit_ou(series, .1)
    interval = profile_theta_interval(series, .1)
    assert not interval.empty
    assert interval.lower == 0 or np.isinf(interval.upper)


def test_bootstrap_matches_manual_recursion_and_keeps_non_ou_slopes():
    series = np.array([1., .8, .7, .9, .4, .5])
    dt, count, seed = .5, 99, 654
    fitted = fit_ou(series, dt)
    result = bootstrap_ou(series, dt, np.random.default_rng(seed), count)
    innovations = np.random.default_rng(seed).standard_normal((count, len(series)-1))*np.sqrt(fitted.innovation_variance)
    slopes = []
    for innovation in innovations:
        path = [series[0]]
        for noise in innovation:
            path.append(fitted.intercept+fitted.phi*path[-1]+noise)
        design = np.column_stack([np.ones(len(path)-1), path[:-1]])
        slopes.append(np.linalg.lstsq(design, path[1:], rcond=None)[0][1])
    np.testing.assert_allclose(result['slopes'], slopes, atol=1e-13)
    bad = sum(not 0 < phi < 1 for phi in slopes)
    assert bad > 0
    assert result['inadmissible_slopes'] == bad
    assert result['bootstrap_mean_phi'] == pytest.approx(np.mean(slopes))
    assert len(result['likelihood_ratios']) == count


def test_time_unit_changes_scale_intervals_and_correction():
    series = simulate_ou_exact(1., 0., 1., 0., .1, 200, 1, np.random.default_rng(777))[0]
    for method in (profile_theta_interval, wald_theta_interval):
        first, second = method(series, .1), method(series, .4)
        np.testing.assert_allclose([second.lower, second.upper], np.array([first.lower, first.upper])/4)
    first = bootstrap_ou(series, .1, np.random.default_rng(19), 99)
    second = bootstrap_ou(series, .4, np.random.default_rng(19), 99)
    np.testing.assert_allclose(first['slopes'], second['slopes'], atol=1e-13)
    assert second['corrected_theta'] == pytest.approx(first['corrected_theta']/4)


def test_drift_inference_is_invariant_under_affine_state_normalization():
    series = simulate_ou_exact(1., 0., 1., 0., .1, 200, 1, np.random.default_rng(719))[0]
    transformed = 3*series-7
    for method in (profile_theta_interval, wald_theta_interval):
        original, rescaled = method(series, .1), method(transformed, .1)
        np.testing.assert_allclose([original.lower, original.upper], [rescaled.lower, rescaled.upper], atol=1e-11)
    original = bootstrap_ou(series, .1, np.random.default_rng(19), 99)
    rescaled = bootstrap_ou(transformed, .1, np.random.default_rng(19), 99)
    np.testing.assert_allclose(original['slopes'], rescaled['slopes'], atol=1e-11)
    assert original['critical_value'] == pytest.approx(rescaled['critical_value'], abs=1e-11)


@pytest.mark.parametrize('bounds, empty, unbounded', [((-1, -.1), True, False),
    ((1.1, 2), True, False), ((-.1, .9), False, True), ((.3, 1.1), False, False)])
def test_slope_confidence_sets_keep_empty_and_unbounded_cases(bounds, empty, unbounded):
    interval = theta_interval_from_phi(*bounds, .1)
    assert interval.empty == empty
    assert np.isinf(interval.upper) == unbounded
    if empty:
        assert not interval.contains(1.)


@pytest.mark.parametrize('dt, level, cutoff', [(0, .95, 1), (.1, 1, 1), (.1, .95, -1)])
def test_profile_rejects_invalid_arguments(dt, level, cutoff):
    with pytest.raises(ValueError):
        profile_theta_interval([1, .8, .7, .9, .4, .5], dt, level, cutoff)
