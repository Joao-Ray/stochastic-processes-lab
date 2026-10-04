"""Check inference against a numerical likelihood and independent benchmarks."""
import numpy as np
import pytest
from scipy.optimize import minimize

from src.ou_inference import OUFitError, fit_ou, forecast_ou
from src.ou_process import simulate_ou_exact


def test_fit_matches_numerically_optimized_conditional_likelihood():
    series = simulate_ou_exact(.7, 1.5, .8, 1.5, .2, 1500, 1,
                               np.random.default_rng(123))[0]
    fit = fit_ou(series, .2)

    def objective(parameters):
        theta, mu, sigma = np.exp(parameters[0]), parameters[1], np.exp(parameters[2])
        phi = np.exp(-theta * .2)
        q = sigma**2 * -np.expm1(-2 * theta * .2) / (2 * theta)
        residual = series[1:] - mu - phi * (series[:-1] - mu)
        return .5 * (residual.size * np.log(2 * np.pi * q) + np.dot(residual, residual) / q)

    optimized = minimize(objective, [np.log(.5), 1., np.log(.6)], method='BFGS',
                         options={'gtol': 1e-5})
    assert optimized.fun == pytest.approx(fit.negative_log_likelihood, abs=1e-6)
    np.testing.assert_allclose([np.exp(optimized.x[0]), optimized.x[1], np.exp(optimized.x[2])],
                               [fit.theta, fit.mu, fit.sigma], rtol=1e-4)


def test_fitted_variance_uses_likelihood_denominator_and_time_units():
    series = simulate_ou_exact(1., 0., .5, 0., .5, 500, 1, np.random.default_rng(7))[0]
    fit = fit_ou(series, .5)
    residual = series[1:] - fit.intercept - fit.phi * series[:-1]
    assert fit.innovation_variance == pytest.approx(np.dot(residual, residual) / 500)
    doubled = fit_ou(series, 1.)
    assert doubled.theta == pytest.approx(fit.theta / 2)
    assert doubled.sigma == pytest.approx(fit.sigma / np.sqrt(2))
    assert doubled.mu == pytest.approx(fit.mu)
    assert doubled.half_life == pytest.approx(2 * fit.half_life)


def test_forecast_limits_match_one_step_transition_and_zero_horizon():
    theta, mu, sigma, x, h = .7, 1.5, .8, -1., .5
    mean, lower, upper = forecast_ou(theta, mu, sigma, x, np.array([0., h]))
    assert mean[0] == lower[0] == upper[0] == x
    assert mean[1] == pytest.approx(mu + np.exp(-theta*h)*(x-mu))
    assert upper[1] - mean[1] == pytest.approx(
        1.959963984540054 * sigma * np.sqrt((1-np.exp(-2*theta*h))/(2*theta)))


@pytest.mark.parametrize('series', [np.ones(10), 2.**np.arange(10), (-.5)**np.arange(10)])
def test_degenerate_or_non_ou_fits_are_not_clipped(series):
    with pytest.raises(OUFitError):
        fit_ou(series, .1)


@pytest.mark.parametrize('series, dt', [([1,2,3], .1), ([1,2,np.nan,3], .1),
                                       ([1,2,3,4], 0), ([[1,2,3,4]], .1)])
def test_invalid_data_rejected(series, dt):
    with pytest.raises(ValueError):
        fit_ou(np.asarray(series), dt)


@pytest.mark.parametrize('horizons, level', [([-1], .95), ([1], 1.), ([], .95)])
def test_invalid_forecasts_rejected(horizons, level):
    with pytest.raises(ValueError):
        forecast_ou(.7, 1., .8, 0., np.asarray(horizons), level)
