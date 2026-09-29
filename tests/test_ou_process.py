"""Checks for Ornstein--Uhlenbeck theory and Euler simulation."""

import numpy as np
import pytest

from src.ou_process import (
    ou_autocorrelation,
    ou_half_life,
    ou_mean,
    ou_variance,
    sample_autocorrelation,
    simulate_ou_euler,
    stationary_variance,
)


def test_theoretical_moments_and_half_life() -> None:
    theta, mu, sigma, x0 = 0.7, 1.5, 0.8, -1.0
    assert ou_mean(0, theta, mu, x0) == pytest.approx(x0)
    assert ou_variance(0, theta, sigma) == pytest.approx(0)
    assert ou_mean(100, theta, mu, x0) == pytest.approx(mu)
    assert ou_variance(100, theta, sigma) == pytest.approx(
        stationary_variance(theta, sigma)
    )
    half_life = ou_half_life(theta)
    displacement = ou_mean(half_life, theta, mu, x0) - mu
    assert displacement == pytest.approx((x0 - mu) / 2)


def test_euler_deterministic_recursion_when_sigma_is_zero() -> None:
    paths = simulate_ou_euler(
        theta=1.0,
        mu=2.0,
        sigma=0.0,
        x0=0.0,
        dt=0.1,
        n_steps=4,
        n_paths=2,
        rng=np.random.default_rng(1),
    )
    expected = np.array([0.0, 0.2, 0.38, 0.542, 0.6878])
    assert np.allclose(paths, np.tile(expected, (2, 1)))


def test_simulation_is_reproducible_and_has_expected_shape() -> None:
    arguments = dict(theta=0.5, mu=1.0, sigma=0.3, x0=0.0, dt=0.02,
                     n_steps=100, n_paths=50)
    first = simulate_ou_euler(**arguments, rng=np.random.default_rng(4))
    second = simulate_ou_euler(**arguments, rng=np.random.default_rng(4))
    assert first.shape == (50, 101)
    assert np.array_equal(first, second)
    assert np.all(first[:, 0] == 0)


def test_long_run_autocorrelation_matches_ou_theory() -> None:
    theta, dt = 0.8, 0.01
    path = simulate_ou_euler(theta, 0.0, 0.7, 0.0, dt, 100_000, 1,
                             np.random.default_rng(5))[0, 10_000:]
    lag = 50
    empirical = sample_autocorrelation(path, lag)[lag]
    assert empirical == pytest.approx(ou_autocorrelation(lag * dt, theta), abs=0.04)


@pytest.mark.parametrize(
    "action",
    [
        lambda: simulate_ou_euler(0, 0, 1, 0, 0.1, 10, 1, np.random.default_rng()),
        lambda: simulate_ou_euler(1, 0, -1, 0, 0.1, 10, 1, np.random.default_rng()),
        lambda: simulate_ou_euler(1, 0, 1, 0, 2.0, 10, 1, np.random.default_rng()),
        lambda: ou_variance(-1, 1, 1),
        lambda: sample_autocorrelation(np.ones(10), 2),
    ],
)
def test_invalid_inputs_raise_value_error(action) -> None:
    with pytest.raises(ValueError):
        action()
