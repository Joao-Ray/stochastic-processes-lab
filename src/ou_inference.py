"""Conditional Gaussian OU estimation and forecasts on a regular time grid."""
from dataclasses import dataclass

import numpy as np
from scipy.stats import norm

from src.ou_process import ou_mean, ou_variance


class OUFitError(ValueError):
    """The unconstrained AR(1) fit has no admissible interior OU mapping."""


@dataclass(frozen=True)
class OUFit:
    theta: float
    mu: float
    sigma: float
    phi: float
    intercept: float
    innovation_variance: float
    dt: float
    n_transitions: int
    negative_log_likelihood: float

    @property
    def half_life(self) -> float:
        return float(np.log(2) / self.theta)


def fit_ou(series: np.ndarray, dt: float) -> OUFit:
    """Fit exact OU transitions conditional on the first observation.

    OLS gives the conditional Gaussian likelihood optimum for the AR(1)
    intercept and slope. The residual variance uses SSE/n (the MLE), not
    SSE/(n-2). Mapping to OU is valid only if 0 < phi < 1 and variance > 0.
    Boundary/non-OU fits are reported explicitly, never clipped or discarded.
    No stationary likelihood term for the first observation is included.
    """
    values = np.asarray(series, dtype=float)
    if values.ndim != 1 or values.size < 4 or np.any(~np.isfinite(values)):
        raise ValueError("series must have at least four finite observations")
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be positive and finite")
    previous, following = values[:-1], values[1:]
    centered = previous - previous.mean()
    denominator = np.dot(centered, centered)
    if denominator == 0:
        raise OUFitError("constant regressors do not identify OU parameters")
    phi = float(np.dot(centered, following - following.mean()) / denominator)
    if not 0 < phi < 1:
        raise OUFitError(f"unconstrained AR(1) slope {phi:.6g} is outside (0, 1)")
    intercept = float(following.mean() - phi * previous.mean())
    residuals = following - intercept - phi * previous
    variance = float(np.mean(residuals**2))
    if np.linalg.norm(residuals) <= 32 * np.finfo(float).eps * max(np.linalg.norm(following), 1.):
        raise OUFitError("zero residual variance has no finite Gaussian likelihood maximum")
    theta = float(-np.log(phi) / dt)
    mu = float(intercept / (1 - phi))
    sigma = float(np.sqrt(variance * 2 * theta / -np.expm1(-2 * theta * dt)))
    n = following.size
    nll = float(n / 2 * (np.log(2 * np.pi * variance) + 1))
    return OUFit(theta, mu, sigma, phi, intercept, variance, float(dt), n, nll)


def forecast_ou(
    theta: float, mu: float, sigma: float, last_value: float,
    horizons: np.ndarray, level: float = .95,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return conditional means and central Gaussian prediction limits.

    With estimated parameters, these are plug-in limits: parameter uncertainty
    is excluded. With true parameters, they have exact marginal coverage.
    Horizons are elapsed times, not numbers of observations.
    """
    times = np.asarray(horizons, dtype=float)
    if times.ndim != 1 or times.size == 0:
        raise ValueError("horizons must be a nonempty one-dimensional array")
    if not np.isfinite(level) or not 0 < level < 1:
        raise ValueError("level must lie strictly between zero and one")
    mean = ou_mean(times, theta, mu, last_value)
    variance = ou_variance(times, theta, sigma)
    radius = norm.ppf((1 + level) / 2) * np.sqrt(variance)
    return mean, mean - radius, mean + radius
