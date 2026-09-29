"""Simulation helpers for the one-dimensional simple symmetric random walk."""

from __future__ import annotations

import numpy as np


def simulate_random_walk(
    n_steps: int, n_paths: int, rng: np.random.Generator
) -> np.ndarray:
    """Return integer positions with shape ``(n_paths, n_steps + 1)``.

    Column zero is the common starting position, S_0 = 0. Each later
    increment is independently +1 or -1 with probability one half.
    """
    if n_steps < 1 or n_paths < 1:
        raise ValueError("n_steps and n_paths must be positive")

    increments = 2 * rng.integers(0, 2, size=(n_paths, n_steps), dtype=np.int8) - 1
    paths = np.empty((n_paths, n_steps + 1), dtype=np.int64)
    paths[:, 0] = 0
    paths[:, 1:] = np.cumsum(increments, axis=1, dtype=np.int64)
    return paths


def endpoint_statistics(endpoints: np.ndarray) -> tuple[float, float]:
    """Return the empirical mean and unbiased sample variance of endpoints."""
    values = np.asarray(endpoints)
    if values.ndim != 1 or values.size < 2:
        raise ValueError("endpoints must be a one-dimensional array of size >= 2")
    return float(np.mean(values)), float(np.var(values, ddof=1))


def simulate_first_hitting_times(
    level: int, n_paths: int, max_steps: int, rng: np.random.Generator
) -> np.ndarray:
    """Simulate first visits to positive ``level`` from S_0 = 0.

    Returns one time per path. A value of -1 means the path did not hit
    the level by ``max_steps``; it is right-censored, not an infinite time.
    Only unfinished paths are advanced at each step.
    """
    if level < 1 or n_paths < 1 or max_steps < 1:
        raise ValueError("level, n_paths, and max_steps must be positive")

    hitting_times = np.full(n_paths, -1, dtype=np.int64)
    active_ids = np.arange(n_paths)
    active_positions = np.zeros(n_paths, dtype=np.int64)

    for step in range(1, max_steps + 1):
        if active_ids.size == 0:
            break
        increments = 2 * rng.integers(0, 2, size=active_ids.size, dtype=np.int8) - 1
        active_positions += increments
        hit = active_positions == level
        hitting_times[active_ids[hit]] = step
        active_ids = active_ids[~hit]
        active_positions = active_positions[~hit]

    return hitting_times
