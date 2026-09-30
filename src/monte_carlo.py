"""Reusable Monte Carlo estimators and error summaries."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


def _validate_positive_integer(value: int, name: str) -> int:
    if not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def pi_indicators(n_samples: int, rng: np.random.Generator) -> np.ndarray:
    """Return indicators for uniform points falling inside a quarter circle."""
    size = _validate_positive_integer(n_samples, "n_samples")
    x = rng.random(size)
    y = rng.random(size)
    return x * x + y * y <= 1.0


def estimate_pi(n_samples: int, rng: np.random.Generator) -> float:
    """Estimate pi by quarter-circle area sampling."""
    return float(4 * np.mean(pi_indicators(n_samples, rng)))


def replicate_pi_estimates(
    n_samples: int,
    n_replications: int,
    rng: np.random.Generator,
    batch_size: int = 64,
) -> np.ndarray:
    """Generate independent pi estimates in memory-bounded batches."""
    samples = _validate_positive_integer(n_samples, "n_samples")
    replications = _validate_positive_integer(n_replications, "n_replications")
    batch = _validate_positive_integer(batch_size, "batch_size")
    estimates = np.empty(replications, dtype=float)
    for start in range(0, replications, batch):
        stop = min(start + batch, replications)
        current_batch = stop - start
        x = rng.random((current_batch, samples))
        y = rng.random((current_batch, samples))
        estimates[start:stop] = 4 * np.mean(x * x + y * y <= 1.0, axis=1)
    return estimates


def estimate_integral(
    function: Callable[[np.ndarray], np.ndarray],
    lower: float,
    upper: float,
    n_samples: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Estimate an integral by uniform sampling and return estimate and SE."""
    samples = _validate_positive_integer(n_samples, "n_samples")
    if not np.isfinite(lower + upper) or upper <= lower:
        raise ValueError("integration bounds must be finite with lower < upper")
    points = rng.uniform(lower, upper, size=samples)
    values = np.asarray(function(points), dtype=float)
    if values.shape != points.shape or np.any(~np.isfinite(values)):
        raise ValueError("function must return one finite value per sample")
    width = upper - lower
    estimate = width * np.mean(values)
    standard_error = width * np.std(values, ddof=1) / np.sqrt(samples)
    return float(estimate), float(standard_error)


def running_mean(values: np.ndarray) -> np.ndarray:
    """Return the cumulative sample mean after every observation."""
    samples = np.asarray(values, dtype=float)
    if samples.ndim != 1 or samples.size < 1 or np.any(~np.isfinite(samples)):
        raise ValueError("values must be a nonempty finite one-dimensional array")
    return np.cumsum(samples) / np.arange(1, samples.size + 1)


def root_mean_squared_error(estimates: np.ndarray, truth: float) -> float:
    """Return the root mean squared error around a known target."""
    values = np.asarray(estimates, dtype=float)
    if values.ndim != 1 or values.size < 1 or np.any(~np.isfinite(values)):
        raise ValueError("estimates must be a nonempty finite one-dimensional array")
    if not np.isfinite(truth):
        raise ValueError("truth must be finite")
    return float(np.sqrt(np.mean((values - truth) ** 2)))
