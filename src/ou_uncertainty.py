"""Conditional OU drift confidence sets and reproducible parametric bootstrap.

Likelihood ratios are exact calculations; their chi-square calibration is
asymptotic. Bootstrap calibration is a plug-in approximation, not exact inference.
"""
from dataclasses import dataclass

import numpy as np
from scipy.signal import lfilter
from scipy.stats import chi2, norm

from src.ou_inference import fit_ou


@dataclass(frozen=True)
class ARStatistics:
    n: int
    phi: float
    intercept: float
    xx: float
    sse: float

    def sse_at(self, phi):
        return self.sse + self.xx * (np.asarray(phi) - self.phi)**2


@dataclass(frozen=True)
class ThetaInterval:
    lower: float
    upper: float
    empty: bool = False

    def contains(self, theta):
        return bool(not self.empty and self.lower <= theta <= self.upper)

    @property
    def width(self):
        return 0. if self.empty else self.upper - self.lower


def _dt_level(dt, level):
    if not np.isfinite(dt) or dt <= 0 or not np.isfinite(level) or not 0 < level < 1:
        raise ValueError('dt must be positive and level must be in (0, 1)')


def ar_statistics(series):
    values = np.asarray(series, dtype=float)
    if values.ndim != 1 or values.size < 4 or not np.isfinite(values).all():
        raise ValueError('at least four finite one-dimensional observations required')
    u, y = values[:-1], values[1:]
    centered = u - u.mean()
    xx = float(centered @ centered)
    if xx <= 0:
        raise ValueError('constant regressors do not identify a slope')
    phi = float(centered @ (y-y.mean()) / xx)
    intercept = float(y.mean()-phi*u.mean())
    residual = y-intercept-phi*u
    sse = float(residual @ residual)
    if sse <= 0 or np.linalg.norm(residual) <= 32*np.finfo(float).eps*max(np.linalg.norm(y), 1.):
        raise ValueError('positive residual variance required')
    return ARStatistics(y.size, phi, intercept, xx, sse)


def theta_interval_from_phi(lower, upper, dt):
    """Intersect a slope interval with (0,1); preserve empty/unbounded sets.

    Boundary endpoints denote the closure of the set for reporting. No slope
    estimate is clipped into the physical OU domain.
    """
    _dt_level(dt, .95)
    if not np.isfinite([lower, upper]).all() or lower > upper:
        raise ValueError('slope bounds must be finite and ordered')
    if upper <= 0 or lower >= 1:
        return ThetaInterval(0., 0., True)
    lower, upper = max(lower, 0.), min(upper, 1.)
    return ThetaInterval(float(-np.log(upper)/dt),
                         float(-np.log(lower)/dt) if lower > 0 else np.inf)


def profile_theta_interval(series, dt, level=.95, critical_value=None):
    """Invert the profiled conditional likelihood over the physical OU domain."""
    _dt_level(dt, level)
    cutoff = float(chi2.ppf(level, 1)) if critical_value is None else critical_value
    if not np.isfinite(cutoff) or cutoff < 0:
        raise ValueError('critical value must be finite and nonnegative')
    stats = ar_statistics(series)
    # This is a boundary supremum for the likelihood ratio, not a point estimate.
    best_phi = min(max(stats.phi, 0.), 1.)
    best_sse = float(stats.sse_at(best_phi))
    radius_sq = (best_sse-stats.sse + best_sse*np.expm1(cutoff/stats.n))/stats.xx
    radius = np.sqrt(max(radius_sq, 0.))
    return theta_interval_from_phi(stats.phi-radius, stats.phi+radius, dt)


def wald_theta_interval(series, dt, level=.95):
    """Delta-method Wald interval using the conditional regression Hessian."""
    _dt_level(dt, level)
    fit, stats = fit_ou(series, dt), ar_statistics(series)
    se = np.sqrt(fit.innovation_variance/stats.xx)/(dt*fit.phi)
    radius = norm.ppf((1+level)/2)*se
    return ThetaInterval(float(max(0., fit.theta-radius)), float(fit.theta+radius))


def _batch_statistics(paths):
    u, y = paths[:, :-1], paths[:, 1:]
    u_mean, y_mean = u.mean(axis=1), y.mean(axis=1)
    centered = u-u_mean[:, None]
    xx = np.sum(centered**2, axis=1)
    if np.any(xx <= 0):
        raise ValueError('degenerate bootstrap regressor')
    phi = np.sum(centered*(y-y_mean[:, None]), axis=1)/xx
    intercept = y_mean-phi*u_mean
    residual = y-intercept[:, None]-phi[:, None]*u
    sse = np.sum(residual**2, axis=1)
    if np.any(sse <= 0):
        raise ValueError('degenerate bootstrap residuals')
    return phi, xx, sse


def bootstrap_ou(series, dt, rng, replications=399, level=.95):
    """Bootstrap conditional on observed X_0, retaining every AR slope.

    Basic slope intervals are intersected with (0,1) only after construction.
    Bootstrap LR minima can be boundary suprema. No inadmissible replicate is
    resampled or omitted. The generating fit must be an interior OU fit.
    """
    _dt_level(dt, level)
    if type(replications) is not int or replications < 19:
        raise ValueError('at least 19 bootstrap replications required')
    fit = fit_ou(series, dt)
    n = fit.n_transitions
    innovation = np.sqrt(fit.innovation_variance)*rng.standard_normal((replications, n))
    paths = np.empty((replications, n+1))
    paths[:, 0] = np.asarray(series)[0]
    paths[:, 1:] = (fit.mu + lfilter([1.], [1., -fit.phi], innovation, axis=1)
                     + (paths[:, :1]-fit.mu)*fit.phi**np.arange(1, n+1))
    slopes, xx, sse = _batch_statistics(paths)
    constrained_slopes = np.clip(slopes, 0., 1.)
    best_sse = sse+xx*(constrained_slopes-slopes)**2
    null_sse = sse+xx*(fit.phi-slopes)**2
    lr = np.maximum(n*np.log(null_sse/best_sse), 0.)
    critical_value = float(np.quantile(lr, level, method='linear'))
    alpha = (1-level)/2
    low, high = np.quantile(slopes, [alpha, 1-alpha], method='linear')
    basic = theta_interval_from_phi(2*fit.phi-high, 2*fit.phi-low, dt)
    calibrated = profile_theta_interval(series, dt, level, critical_value)
    corrected_phi = 2*fit.phi-float(slopes.mean())
    corrected_theta = float(-np.log(corrected_phi)/dt) if 0 < corrected_phi < 1 else None
    return dict(basic=basic, calibrated_profile=calibrated,
        corrected_theta=corrected_theta, corrected_phi=corrected_phi,
        bootstrap_mean_phi=float(slopes.mean()), critical_value=critical_value,
        inadmissible_slopes=int(np.count_nonzero((slopes <= 0) | (slopes >= 1))),
        replications=replications, slopes=slopes, likelihood_ratios=lr, fit=fit)
