"""Utilities for finite-state, time-homogeneous Markov chains."""

from __future__ import annotations

import numpy as np


def validate_transition_matrix(matrix: np.ndarray, atol: float = 1e-12) -> np.ndarray:
    """Return a validated floating-point row-stochastic transition matrix."""
    transition = np.asarray(matrix, dtype=float)
    if transition.ndim != 2 or transition.shape[0] != transition.shape[1]:
        raise ValueError("transition matrix must be square")
    if transition.shape[0] < 1 or np.any(~np.isfinite(transition)):
        raise ValueError("transition matrix must contain finite values")
    if np.any(transition < -atol) or not np.allclose(
        transition.sum(axis=1), 1.0, atol=atol, rtol=0
    ):
        raise ValueError("transition matrix must have nonnegative rows summing to one")
    return transition


def n_step_transition(matrix: np.ndarray, n_steps: int) -> np.ndarray:
    """Return the matrix of ``n_steps``-step transition probabilities."""
    transition = validate_transition_matrix(matrix)
    if not isinstance(n_steps, (int, np.integer)) or n_steps < 0:
        raise ValueError("n_steps must be a nonnegative integer")
    return np.linalg.matrix_power(transition, int(n_steps))


def stationary_distribution(matrix: np.ndarray) -> np.ndarray:
    """Solve ``pi @ P = pi`` with entries of ``pi`` summing to one.

    A unique solution is required. Reducible chains can have many stationary
    distributions, in which case this function raises ``ValueError``.
    """
    transition = validate_transition_matrix(matrix)
    n_states = transition.shape[0]
    system = np.vstack((transition.T - np.eye(n_states), np.ones(n_states)))
    target = np.r_[np.zeros(n_states), 1.0]
    if np.linalg.matrix_rank(system) < n_states:
        raise ValueError("stationary distribution is not unique")
    distribution, *_ = np.linalg.lstsq(system, target, rcond=None)
    distribution[np.abs(distribution) < 1e-14] = 0.0
    if np.any(distribution < -1e-10):
        raise ValueError("failed to obtain a probability distribution")
    distribution = np.clip(distribution, 0.0, None)
    return distribution / distribution.sum()


def total_variation_distance(first: np.ndarray, second: np.ndarray) -> float:
    """Return half the L1 distance between two finite distributions."""
    left, right = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    if left.ndim != 1 or right.shape != left.shape or left.size == 0:
        raise ValueError("distributions must be nonempty vectors of equal size")
    for values in (left, right):
        if (np.any(~np.isfinite(values)) or np.any(values < 0)
                or not np.isclose(values.sum(), 1.0, atol=1e-12, rtol=0)):
            raise ValueError("distributions must be finite, nonnegative and sum to one")
    return float(0.5 * np.sum(np.abs(left - right)))


def simulate_chains(
    matrix: np.ndarray,
    initial_state: int,
    n_steps: int,
    n_paths: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Simulate paths and include the initial state in column zero."""
    transition = validate_transition_matrix(matrix)
    n_states = transition.shape[0]
    if not isinstance(initial_state, (int, np.integer)) or not 0 <= initial_state < n_states:
        raise ValueError("initial_state is outside the state space")
    if not isinstance(n_steps, (int, np.integer)) or n_steps < 0:
        raise ValueError("n_steps must be a nonnegative integer")
    if not isinstance(n_paths, (int, np.integer)) or n_paths < 1:
        raise ValueError("n_paths must be a positive integer")

    paths = np.empty((int(n_paths), int(n_steps) + 1), dtype=np.int64)
    paths[:, 0] = int(initial_state)
    cumulative = np.cumsum(transition, axis=1)
    cumulative[:, -1] = 1.0
    for step in range(1, int(n_steps) + 1):
        uniforms = rng.random(int(n_paths))
        paths[:, step] = np.sum(
            uniforms[:, None] > cumulative[paths[:, step - 1]], axis=1
        )
    return paths


def empirical_transition_matrix(paths: np.ndarray, n_states: int) -> np.ndarray:
    """Estimate one-step transition probabilities from one or more paths."""
    states = np.asarray(paths)
    if states.ndim == 1:
        states = states[None, :]
    if states.ndim != 2 or states.shape[1] < 2:
        raise ValueError("paths must contain at least two time points")
    if not isinstance(n_states, (int, np.integer)) or n_states < 1:
        raise ValueError("n_states must be positive")
    if np.any(states < 0) or np.any(states >= n_states):
        raise ValueError("paths contain a state outside the state space")

    counts = np.zeros((n_states, n_states), dtype=np.int64)
    np.add.at(counts, (states[:, :-1].ravel(), states[:, 1:].ravel()), 1)
    row_totals = counts.sum(axis=1, keepdims=True)
    if np.any(row_totals == 0):
        raise ValueError("every state must have at least one observed departure")
    return counts / row_totals
