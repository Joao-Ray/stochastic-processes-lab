"""Independent benchmarks for the methods used in supplementary experiments."""

import numpy as np
import pytest

from src.markov_chain import n_step_transition, total_variation_distance
from src.monte_carlo import replicate_control_variate_integrals
from src.ou_process import (
    ou_euler_moments, ou_mean, ou_variance, simulate_ou_exact,
)


def test_exact_ou_deterministic_solution_even_at_a_coarse_step():
    times = np.arange(5) * 3.0
    paths = simulate_ou_exact(1, 2, 0, -1, 3, 4, 2, np.random.default_rng(12))
    np.testing.assert_allclose(paths, np.tile(ou_mean(times, 1, 2, -1), (2, 1)))


def test_exact_ou_endpoint_matches_continuous_moments():
    arguments = dict(theta=0.7, mu=1.5, sigma=0.8, x0=-1, dt=0.5,
                     n_steps=8, n_paths=40_000)
    first = simulate_ou_exact(**arguments, rng=np.random.default_rng(93))
    second = simulate_ou_exact(**arguments, rng=np.random.default_rng(93))
    assert first.shape == (40_000, 9)
    assert np.array_equal(first, second)
    variance = ou_variance(4, 0.7, 0.8)
    assert abs(first[:, -1].mean() - ou_mean(4, 0.7, 1.5, -1)) < 6 * np.sqrt(variance / 40_000)
    assert abs(first[:, -1].var(ddof=1) - variance) < 6 * variance * np.sqrt(2 / 39_999)


def test_euler_moments_match_geometric_series_and_refine_to_sde():
    theta, mu, sigma, x0, dt, steps = 0.7, 1.5, 0.8, -1, 0.5, 8
    decay = 1 - theta * dt
    mean, variance = ou_euler_moments(theta, mu, sigma, x0, dt, steps)
    assert mean == pytest.approx(mu + (x0 - mu) * decay**steps)
    assert variance == pytest.approx(sigma**2 * dt * (1 - decay**(2 * steps)) / (1 - decay**2))
    fine_mean, fine_variance = ou_euler_moments(theta, mu, sigma, x0, 0.0005, 8000)
    assert fine_mean == pytest.approx(ou_mean(4, theta, mu, x0), abs=1e-4)
    assert fine_variance == pytest.approx(ou_variance(4, theta, sigma), abs=1e-4)
    assert ou_euler_moments(theta, mu, sigma, x0, dt, 0) == (x0, 0)


@pytest.mark.parametrize("alpha", [0.2, 0.02, 1.0])
def test_markov_mixing_distance_matches_eigenvalue_formula(alpha):
    matrix = np.array([[1 - alpha, alpha], [alpha, 1 - alpha]])
    for n in (0, 1, 2, 10, 50):
        distribution = np.array([1., 0.]) @ n_step_transition(matrix, n)
        assert total_variation_distance(distribution, np.array([0.5, 0.5])) == pytest.approx(
            0.5 * abs(1 - 2 * alpha)**n, abs=1e-13)


def test_control_variate_cancels_a_known_integrand_exactly():
    plain, adjusted, plain_se, adjusted_se = replicate_control_variate_integrals(
        lambda u: u**2, lambda u: u**2, 1 / 3, 1., 200, 20, np.random.default_rng(7))
    np.testing.assert_allclose(adjusted, 1 / 3, atol=1e-14)
    assert np.all(adjusted_se < 1e-14)
    assert np.all(plain_se > 0)
    assert plain.std() > 0


def test_zero_control_coefficient_equals_plain_and_is_reproducible():
    args = (lambda u: np.exp(-u**2), lambda u: u**2, 1 / 3, 0., 100, 10)
    first = replicate_control_variate_integrals(*args, np.random.default_rng(19))
    second = replicate_control_variate_integrals(*args, np.random.default_rng(19))
    np.testing.assert_array_equal(first, second)
    np.testing.assert_array_equal(first[0], first[1])
    np.testing.assert_array_equal(first[2], first[3])


@pytest.mark.parametrize("action", [
    lambda: simulate_ou_exact(0, 0, 1, 0, .1, 10, 1, np.random.default_rng(1)),
    lambda: simulate_ou_exact(1, 0, -1, 0, .1, 10, 1, np.random.default_rng(1)),
    lambda: simulate_ou_exact(1, 0, 1, 0, 0, 10, 1, np.random.default_rng(1)),
    lambda: simulate_ou_exact(1, 0, 1, 0, .1, 0, 1, np.random.default_rng(1)),
    lambda: ou_euler_moments(1, 0, 1, 0, 2., 10),
    lambda: total_variation_distance(np.array([.5, .5]), np.array([1.])),
    lambda: total_variation_distance(np.array([-.1, 1.1]), np.array([.5, .5])),
    lambda: replicate_control_variate_integrals(
        lambda u: u, lambda u: u, .5, 1, 1, 5, np.random.default_rng(1)),
    lambda: replicate_control_variate_integrals(
        lambda u: u, lambda u: np.nan * u, .5, 1, 10, 5, np.random.default_rng(1)),
])
def test_invalid_supplementary_inputs(action):
    with pytest.raises(ValueError):
        action()
