"""Checks for finite-state Markov-chain calculations and simulations."""

import numpy as np
import pytest

from src.markov_chain import (
    empirical_transition_matrix,
    n_step_transition,
    simulate_chains,
    stationary_distribution,
    validate_transition_matrix,
)


WEATHER = np.array([[0.8, 0.2], [0.4, 0.6]])


def test_two_step_transition_matches_manual_total_probability() -> None:
    expected = np.array([[0.72, 0.28], [0.56, 0.44]])
    assert np.allclose(n_step_transition(WEATHER, 2), expected)
    assert np.array_equal(n_step_transition(WEATHER, 0), np.eye(2))


def test_stationary_distribution_solves_balance_equation() -> None:
    stationary = stationary_distribution(WEATHER)
    assert np.allclose(stationary, np.array([2 / 3, 1 / 3]))
    assert np.allclose(stationary @ WEATHER, stationary)
    assert np.isclose(stationary.sum(), 1)


def test_simulation_is_reproducible_and_approaches_stationarity() -> None:
    paths = simulate_chains(WEATHER, 0, 30, 50_000, np.random.default_rng(31))
    assert paths.shape == (50_000, 31)
    assert np.all(paths[:, 0] == 0)
    assert np.mean(paths[:, -1] == 0) == pytest.approx(2 / 3, abs=0.01)
    assert np.array_equal(
        paths, simulate_chains(WEATHER, 0, 30, 50_000, np.random.default_rng(31))
    )


def test_empirical_transition_matrix_counts_observed_moves() -> None:
    paths = np.array([[0, 0, 1, 1], [0, 1, 0, 0]])
    expected = np.array([[0.5, 0.5], [0.5, 0.5]])
    assert np.allclose(empirical_transition_matrix(paths, 2), expected)


@pytest.mark.parametrize(
    "action",
    [
        lambda: validate_transition_matrix(np.array([[0.5, 0.6], [0.4, 0.6]])),
        lambda: n_step_transition(WEATHER, -1),
        lambda: stationary_distribution(np.eye(2)),
        lambda: simulate_chains(WEATHER, 2, 5, 10, np.random.default_rng()),
        lambda: empirical_transition_matrix(np.array([0, 0, 0]), 2),
    ],
)
def test_invalid_or_underdetermined_inputs_raise_value_error(action) -> None:
    with pytest.raises(ValueError):
        action()
