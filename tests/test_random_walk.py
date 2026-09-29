"""Checks of random-walk invariants and first-hitting-time conventions."""

from itertools import product

import numpy as np
import pytest

from src.random_walk import (
    endpoint_statistics,
    simulate_first_hitting_times,
    simulate_random_walk,
)


def test_paths_start_at_zero_and_move_one_unit() -> None:
    paths = simulate_random_walk(40, 200, np.random.default_rng(7))
    assert paths.shape == (200, 41)
    assert np.all(paths[:, 0] == 0)
    assert np.all(np.isin(np.diff(paths, axis=1), (-1, 1)))
    assert np.array_equal(paths, simulate_random_walk(40, 200, np.random.default_rng(7)))


def test_exact_four_step_moments_by_enumeration() -> None:
    endpoints = np.array([sum(steps) for steps in product((-1, 1), repeat=4)])
    assert np.mean(endpoints) == 0
    assert np.var(endpoints) == 4
    mean, sample_variance = endpoint_statistics(endpoints)
    assert mean == 0
    assert sample_variance == pytest.approx(64 / 15)


def test_hitting_times_are_first_possible_lattice_times_or_censored() -> None:
    times = simulate_first_hitting_times(3, 500, 15, np.random.default_rng(21))
    observed = times[times >= 0]
    assert observed.size > 0
    assert np.all((observed >= 3) & (observed <= 15))
    assert np.all(observed % 2 == 1)
    assert np.array_equal(
        times, simulate_first_hitting_times(3, 500, 15, np.random.default_rng(21))
    )
    assert np.all(simulate_first_hitting_times(3, 20, 2, np.random.default_rng(1)) == -1)


@pytest.mark.parametrize(
    "action",
    [
        lambda: simulate_random_walk(0, 2, np.random.default_rng()),
        lambda: endpoint_statistics(np.array([0])),
        lambda: simulate_first_hitting_times(0, 2, 10, np.random.default_rng()),
    ],
)
def test_invalid_inputs_raise_value_error(action) -> None:
    with pytest.raises(ValueError):
        action()
