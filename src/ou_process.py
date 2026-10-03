"""Theory and Euler--Maruyama simulation for the OU process."""

from __future__ import annotations

import numpy as np


def _validate_parameters(theta: float, sigma: float) -> None:
    if not np.isfinite(theta) or theta <= 0:
        raise ValueError("theta must be a positive finite number")
    if not np.isfinite(sigma) or sigma < 0:
        raise ValueError("sigma must be a nonnegative finite number")


def ou_mean(times: np.ndarray | float, theta: float, mu: float, x0: float) -> np.ndarray:
    """Return ``E[X_t | X_0=x0]`` for an Ornstein--Uhlenbeck process."""
    _validate_parameters(theta, 0.0)
    time = np.asarray(times, dtype=float)
    if np.any(~np.isfinite(time)) or np.any(time < 0) or not np.isfinite(mu + x0):
        raise ValueError("times must be nonnegative and parameters must be finite")
    return mu + (x0 - mu) * np.exp(-theta * time)


def ou_variance(times: np.ndarray | float, theta: float, sigma: float) -> np.ndarray:
    """Return ``Var[X_t | X_0=x0]`` for deterministic ``X_0``."""
    _validate_parameters(theta, sigma)
    time = np.asarray(times, dtype=float)
    if np.any(~np.isfinite(time)) or np.any(time < 0):
        raise ValueError("times must be finite and nonnegative")
    return sigma**2 / (2 * theta) * (1 - np.exp(-2 * theta * time))


def stationary_variance(theta: float, sigma: float) -> float:
    """Return the variance of the stationary normal distribution."""
    _validate_parameters(theta, sigma)
    return float(sigma**2 / (2 * theta))


def ou_autocorrelation(lags: np.ndarray | float, theta: float) -> np.ndarray:
    """Return the stationary autocorrelation ``exp(-theta * lag)``."""
    _validate_parameters(theta, 0.0)
    lag = np.asarray(lags, dtype=float)
    if np.any(~np.isfinite(lag)) or np.any(lag < 0):
        raise ValueError("lags must be finite and nonnegative")
    return np.exp(-theta * lag)


def ou_half_life(theta: float) -> float:
    """Return the time required for a mean displacement to halve."""
    _validate_parameters(theta, 0.0)
    return float(np.log(2) / theta)


def simulate_ou_euler(
    theta: float,
    mu: float,
    sigma: float,
    x0: float,
    dt: float,
    n_steps: int,
    n_paths: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Simulate OU paths using the Euler--Maruyama method.

    The initial value is included in column zero. The restriction
    ``theta * dt < 2`` keeps the autoregressive Euler recursion stable.
    """
    _validate_parameters(theta, sigma)
    if not np.isfinite(mu + x0) or not np.isfinite(dt) or dt <= 0:
        raise ValueError("mu and x0 must be finite and dt must be positive")
    if theta * dt >= 2:
        raise ValueError("Euler step is unstable: require theta * dt < 2")
    if not isinstance(n_steps, (int, np.integer)) or n_steps < 1:
        raise ValueError("n_steps must be a positive integer")
    if not isinstance(n_paths, (int, np.integer)) or n_paths < 1:
        raise ValueError("n_paths must be a positive integer")

    paths = np.empty((int(n_paths), int(n_steps) + 1), dtype=float)
    paths[:, 0] = x0
    noise_scale = sigma * np.sqrt(dt)
    for step in range(int(n_steps)):
        current = paths[:, step]
        paths[:, step + 1] = (
            current
            + theta * (mu - current) * dt
            + noise_scale * rng.standard_normal(int(n_paths))
        )
    return paths


def simulate_ou_exact(
    theta: float,
    mu: float,
    sigma: float,
    x0: float,
    dt: float,
    n_steps: int,
    n_paths: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample the exact OU transition law at equally spaced grid points.

    This is exact in distribution on the grid, up to floating-point and
    pseudorandom sampling error. It does not interpolate a continuous path.
    Unlike Euler, it has no step-size stability restriction.
    """
    _validate_parameters(theta, sigma)
    if not np.isfinite(mu + x0) or not np.isfinite(dt) or dt <= 0:
        raise ValueError("mu and x0 must be finite and dt must be positive")
    if not isinstance(n_steps, (int, np.integer)) or n_steps < 1:
        raise ValueError("n_steps must be a positive integer")
    if not isinstance(n_paths, (int, np.integer)) or n_paths < 1:
        raise ValueError("n_paths must be a positive integer")
    decay = np.exp(-theta * dt)
    noise_scale = sigma * np.sqrt(-np.expm1(-2 * theta * dt) / (2 * theta))
    paths = np.empty((int(n_paths), int(n_steps) + 1), dtype=float)
    paths[:, 0] = x0
    for step in range(int(n_steps)):
        paths[:, step + 1] = (
            mu + decay * (paths[:, step] - mu)
            + noise_scale * rng.standard_normal(int(n_paths))
        )
    return paths


def ou_euler_moments(
    theta: float, mu: float, sigma: float, x0: float, dt: float, n_steps: int
) -> tuple[float, float]:
    """Return the finite-step Euler mean and variance for deterministic X_0.

    These are moments of the discretized recursion, not the continuous SDE.
    The recurrence avoids cancellation in the geometric-series formula.
    """
    _validate_parameters(theta, sigma)
    if not np.isfinite(mu + x0) or not np.isfinite(dt) or dt <= 0:
        raise ValueError("mu and x0 must be finite and dt must be positive")
    if theta * dt >= 2:
        raise ValueError("Euler step is unstable: require theta * dt < 2")
    if not isinstance(n_steps, (int, np.integer)) or n_steps < 0:
        raise ValueError("n_steps must be a nonnegative integer")
    decay = 1 - theta * dt
    mean, variance = float(x0), 0.0
    for _ in range(int(n_steps)):
        mean = mu + decay * (mean - mu)
        variance = decay**2 * variance + sigma**2 * dt
    return float(mean), float(variance)


def sample_autocorrelation(series: np.ndarray, max_lag: int) -> np.ndarray:
    """Estimate autocorrelation at integer lags from one series."""
    values = np.asarray(series, dtype=float)
    if values.ndim != 1 or values.size < 2 or np.any(~np.isfinite(values)):
        raise ValueError("series must be a finite one-dimensional array")
    if not isinstance(max_lag, (int, np.integer)) or not 0 <= max_lag < values.size:
        raise ValueError("max_lag must be between zero and len(series) - 1")
    centered = values - np.mean(values)
    denominator = np.dot(centered, centered)
    if denominator == 0:
        raise ValueError("autocorrelation is undefined for a constant series")
    result = np.empty(int(max_lag) + 1)
    result[0] = 1.0
    for lag in range(1, int(max_lag) + 1):
        result[lag] = np.dot(centered[:-lag], centered[lag:]) / denominator
    return result
