"""Checks for Monte Carlo estimators and convergence utilities."""

import numpy as np
import pytest

from src.monte_carlo import (
    estimate_integral,
    estimate_pi,
    replicate_pi_estimates,
    root_mean_squared_error,
    running_mean,
)


def test_pi_estimate_is_reproducible_and_close_to_pi() -> None:
    first = estimate_pi(300_000, np.random.default_rng(40))
    second = estimate_pi(300_000, np.random.default_rng(40))
    assert first == second
    assert first == pytest.approx(np.pi, abs=0.01)


def test_replicated_pi_estimates_have_expected_shape() -> None:
    first = replicate_pi_estimates(500, 80, np.random.default_rng(41), batch_size=13)
    second = replicate_pi_estimates(500, 80, np.random.default_rng(41), batch_size=13)
    assert first.shape == (80,)
    assert np.array_equal(first, second)
    assert np.all((first >= 0) & (first <= 4))


def test_integral_estimator_for_x_squared() -> None:
    estimate, standard_error = estimate_integral(
        lambda x: x**2, 0.0, 1.0, 300_000, np.random.default_rng(42)
    )
    assert estimate == pytest.approx(1 / 3, abs=0.002)
    assert standard_error > 0


def test_running_mean_and_rmse_have_known_values() -> None:
    assert np.allclose(running_mean(np.array([1.0, 3.0, 2.0])), [1.0, 2.0, 2.0])
    assert root_mean_squared_error(np.array([1.0, 3.0]), 2.0) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "action",
    [
        lambda: estimate_pi(0, np.random.default_rng()),
        lambda: replicate_pi_estimates(10, 0, np.random.default_rng()),
        lambda: estimate_integral(lambda x: x, 1, 0, 10, np.random.default_rng()),
        lambda: estimate_integral(lambda x: 1.0, 0, 1, 10, np.random.default_rng()),
        lambda: running_mean(np.array([])),
        lambda: root_mean_squared_error(np.array([np.nan]), 0),
    ],
)
def test_invalid_inputs_raise_value_error(action) -> None:
    with pytest.raises(ValueError):
        action()
