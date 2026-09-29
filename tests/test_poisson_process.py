"""Checks for Poisson-process simulation invariants."""

import numpy as np
import pytest

from src.poisson_process import (
    simulate_arrival_times,
    simulate_count_paths,
    simulate_interarrival_times,
)


def test_count_paths_start_at_zero_and_are_nondecreasing() -> None:
    grid = np.linspace(0, 3, 13)
    paths = simulate_count_paths(2.0, grid, 100, np.random.default_rng(8))
    assert paths.shape == (100, 13)
    assert np.all(paths[:, 0] == 0)
    assert np.all(np.diff(paths, axis=1) >= 0)
    assert np.array_equal(
        paths, simulate_count_paths(2.0, grid, 100, np.random.default_rng(8))
    )


def test_interarrival_sample_matches_exponential_scale() -> None:
    waits = simulate_interarrival_times(4.0, 100_000, np.random.default_rng(12))
    assert np.all(waits > 0)
    assert np.mean(waits) == pytest.approx(0.25, abs=0.003)


def test_arrival_times_respect_horizon_and_are_strictly_increasing() -> None:
    arrivals = simulate_arrival_times(3.0, 5.0, np.random.default_rng(15))
    assert arrivals.size > 0
    assert np.all(np.diff(arrivals) > 0)
    assert np.all((arrivals > 0) & (arrivals <= 5.0))
    assert np.array_equal(
        arrivals, simulate_arrival_times(3.0, 5.0, np.random.default_rng(15))
    )


@pytest.mark.parametrize(
    "action",
    [
        lambda: simulate_count_paths(0, np.array([0.0, 1.0]), 2, np.random.default_rng()),
        lambda: simulate_count_paths(1, np.array([0.0, 0.0]), 2, np.random.default_rng()),
        lambda: simulate_interarrival_times(1, 0, np.random.default_rng()),
        lambda: simulate_arrival_times(1, -1, np.random.default_rng()),
    ],
)
def test_invalid_inputs_raise_value_error(action) -> None:
    with pytest.raises(ValueError):
        action()
