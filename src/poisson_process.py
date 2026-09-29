"""Simulation helpers for a homogeneous Poisson process."""

from __future__ import annotations

import numpy as np


def _validate_rate(rate: float) -> None:
    if not np.isfinite(rate) or rate <= 0:
        raise ValueError("rate must be a positive finite number")


def simulate_count_paths(
    rate: float,
    time_grid: np.ndarray,
    n_paths: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Simulate Poisson-process counts on an increasing grid.

    The returned array has shape ``(n_paths, len(time_grid))``. The first
    grid point must be zero, so every path starts from N(0) = 0.
    """
    _validate_rate(rate)
    times = np.asarray(time_grid, dtype=float)
    if times.ndim != 1 or times.size < 2 or times[0] != 0:
        raise ValueError("time_grid must be one-dimensional and start at zero")
    intervals = np.diff(times)
    if np.any(~np.isfinite(times)) or np.any(intervals <= 0):
        raise ValueError("time_grid must contain finite, strictly increasing values")
    if n_paths < 1:
        raise ValueError("n_paths must be positive")

    increments = rng.poisson(rate * intervals, size=(n_paths, intervals.size))
    paths = np.empty((n_paths, times.size), dtype=np.int64)
    paths[:, 0] = 0
    paths[:, 1:] = np.cumsum(increments, axis=1)
    return paths


def simulate_interarrival_times(
    rate: float, size: int, rng: np.random.Generator
) -> np.ndarray:
    """Draw independent exponential waiting times with mean ``1 / rate``."""
    _validate_rate(rate)
    if size < 1:
        raise ValueError("size must be positive")
    return rng.exponential(scale=1 / rate, size=size)


def simulate_arrival_times(
    rate: float, horizon: float, rng: np.random.Generator
) -> np.ndarray:
    """Return all arrival times up to and including ``horizon``."""
    _validate_rate(rate)
    if not np.isfinite(horizon) or horizon <= 0:
        raise ValueError("horizon must be a positive finite number")

    arrivals: list[float] = []
    current_time = 0.0
    while True:
        current_time += float(rng.exponential(scale=1 / rate))
        if current_time > horizon:
            break
        arrivals.append(current_time)
    return np.asarray(arrivals, dtype=float)
