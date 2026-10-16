"""Stationary Gaussian OU + independent measurement error: ML/REML and limits.

These are established likelihood methods, not finite-sample calibrated CIs.
Equal-time replicates share a latent state and have independent sensor errors.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import cholesky_banded, cho_solve_banded, cho_factor, cho_solve, eigvalsh
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2

MIN_STEP = 1e-5
MAX_STEP = 20.


@dataclass(frozen=True)
class NoiseSample:
    means: np.ndarray
    within_sse: float
    offset: float
    sites: int
    repeats: int

    @property
    def size(self):
        return self.sites*self.repeats


def prepare_sample(observations):
    y = np.asarray(observations, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    if y.ndim != 2 or y.shape[0] < 4 or y.shape[1] < 1 or not np.isfinite(y).all():
        raise ValueError('at least four finite sites, one or more replicates per site required')
    offset = float(y.mean())
    centered = y-offset
    means = centered.mean(axis=1)
    within = float(np.sum((centered-means[:, None])**2))
    if np.sum(centered**2) <= np.finfo(float).tiny:
        raise ValueError('nonconstant observations required')
    if y.shape[1] > 1 and within <= 0:
        raise ValueError('zero within-site variance is a singular replicated-data likelihood; use site means')
    return NoiseSample(means, within, offset, *y.shape)


def _profile(sample, phi, gamma, restricted):
    """Profile location/total variance using a tridiagonal precision identity."""
    if not np.isfinite([phi, gamma]).all() or not 0 <= phi < 1 or not 0 <= gamma <= 1:
        raise ValueError('phi in [0,1), gamma in [0,1] required')
    m, k = sample.sites, sample.repeats
    noise = 1-gamma
    if k > 1 and noise == 0:
        return dict(nll=np.inf, mu=None, scale=None)
    q = (1-phi)*(1+phi)
    # Q=q*R_phi^{-1}; H=gamma*q*I+(noise/k)*Q.
    diagonal = np.full(m, 1+phi*phi)
    diagonal[[0, -1]] = 1.
    band = np.zeros((2, m))
    band[0] = gamma*q+(noise/k)*diagonal
    band[1, :-1] = -(noise/k)*phi
    chol = cholesky_banded(band, lower=True, check_finite=False)
    rhs = np.column_stack([sample.means, np.ones(m)])
    product = diagonal[:, None]*rhs
    product[1:] -= phi*rhs[:-1]
    product[:-1] -= phi*rhs[1:]
    precision_rhs = cho_solve_banded((chol, True), product, check_finite=False)
    information_mean = float(np.sum(precision_rhs[:, 1]))
    mu_centered = float(np.sum(precision_rhs[:, 0])/information_mean)
    residual = sample.means-mu_centered
    quadratic = float(residual @ (precision_rhs[:, 0]-mu_centered*precision_rhs[:, 1]))
    if k > 1:
        quadratic += sample.within_sse/noise
    if quadratic <= 0 or information_mean <= 0:
        raise FloatingPointError('nonpositive profiled quadratic or mean information')
    logdet = float(2*np.log(chol[0]).sum()-np.log(q)+m*np.log(k))
    if k > 1:
        logdet += m*(k-1)*np.log(noise)
    degrees = sample.size-int(restricted)
    scale = quadratic/degrees
    objective = .5*(logdet+degrees*(np.log(2*np.pi*scale)+1))
    if restricted:
        objective += .5*np.log(information_mean)
    return dict(nll=float(objective), mu=float(sample.offset+mu_centered), scale=float(scale))


def profile_shape(observations, phi, gamma, restricted=False):
    return _profile(prepare_sample(observations), phi, gamma, restricted)


def _gamma_profile(sample, phi, restricted):
    if phi == 0 and sample.repeats == 1:
        return dict(**_profile(sample, 0., 0., restricted), gamma=0., converged=True)
    upper = 1. if sample.repeats == 1 else 1-1e-9
    grid = np.linspace(0, upper, 13)
    values = [_profile(sample, phi, float(g), restricted)['nll'] for g in grid]
    candidates = [(v, float(g), True) for v, g in zip(values, grid)]
    for i in range(1, len(grid)-1):
        if values[i] <= values[i-1] and values[i] <= values[i+1]:
            result = minimize_scalar(lambda g: _profile(sample, phi, float(g), restricted)['nll'],
                bounds=(grid[i-1], grid[i+1]), method='bounded', options=dict(xatol=1e-9))
            candidates.append((float(result.fun), float(result.x), bool(result.success)))
    # Search endpoint neighborhoods even when the best grid point is a boundary.
    for low, high in [(grid[0], grid[1]), (grid[-2], grid[-1])]:
        result = minimize_scalar(lambda g: _profile(sample, phi, float(g), restricted)['nll'],
            bounds=(low, high), method='bounded', options=dict(xatol=1e-9))
        candidates.append((float(result.fun), float(result.x), bool(result.success)))
    value, gamma, converged = min(candidates)
    return dict(**_profile(sample, phi, gamma, restricted), gamma=gamma, converged=converged)


def profile_phi(observations, phi, restricted=False):
    """Numerical nuisance profiling, without claiming global certification."""
    return _gamma_profile(prepare_sample(observations), phi, restricted)


def fit_noisy_ou(observations, dt, restricted=False, fixed_gamma=None):
    """Multi-start bounded shape fit; preserve white/cap/nonconvergence flags.

    Drift is unavailable on white-noise limits or the search caps. An optimizer
    convergence flag is numerical, not a proof of global optimality.
    """
    sample = prepare_sample(observations)
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError('positive finite dt required')
    if fixed_gamma is not None and (not np.isfinite(fixed_gamma) or not 0 <= fixed_gamma <= 1 or (fixed_gamma == 1 and sample.repeats > 1)):
        raise ValueError('valid fixed gamma required; replicated data need positive noise')
    upper = 1. if sample.repeats == 1 else 1-1e-9
    def objective(z):
        gamma = float(z[1]) if fixed_gamma is None else fixed_gamma
        return _profile(sample, float(np.exp(-z[0])), gamma, restricted)['nll']
    steps = [.003, .03, .15, .5, 1.5, 5.]
    if fixed_gamma is None:
        # Profile the nuisance at each drift seed: a joint grid can miss a
        # narrow low-signal mode that even differential evolution may miss.
        seeds = [np.array([a, _gamma_profile(sample, float(np.exp(-a)), restricted)['gamma']]) for a in steps]
    else:
        seeds = [np.array([a]) for a in steps]
    seeds.sort(key=objective)
    candidates = []
    bounds = [(MIN_STEP, MAX_STEP)]+([(0., upper)] if fixed_gamma is None else [])
    for seed in seeds[:3]:
        result = minimize(objective, seed, method='L-BFGS-B', bounds=bounds,
            options=dict(ftol=1e-11, gtol=1e-6, maxiter=150))
        if np.isfinite(result.fun):
            candidates.append((float(result.fun), float(np.exp(-result.x[0])),
                float(result.x[1]) if fixed_gamma is None else fixed_gamma, bool(result.success)))
    if not candidates:
        raise FloatingPointError('no finite optimization result')
    nll, phi, gamma, converged = min(candidates, key=lambda row: row[0])
    if fixed_gamma is None:
        # Near the Brownian/REML ridge, noise/innovation coordinates avoid the
        # almost collinear derivatives of a and gamma. No truth enters this step.
        a = float(-np.log(phi))
        ratio = (1-gamma)/max(gamma*(1-phi*phi), np.finfo(float).tiny)
        z0 = [np.log(a), np.clip(np.log(max(ratio, np.exp(-25.))), -25., 30.)]
        def transformed(z):
            candidate_phi = float(np.exp(-np.exp(z[0])))
            candidate_gamma = float(min(upper, 1/(1+np.exp(z[1])*(1-candidate_phi*candidate_phi))))
            return _profile(sample, candidate_phi, candidate_gamma, restricted)['nll']
        alternate = minimize(transformed, z0, method='L-BFGS-B',
            bounds=[(np.log(MIN_STEP), np.log(MAX_STEP)), (-25.,30.)],
            options=dict(ftol=1e-12, gtol=1e-6, maxiter=150))
        if alternate.fun < nll-1e-8:
            nll = float(alternate.fun)
            phi = float(np.exp(-np.exp(alternate.x[0])))
            gamma = float(min(upper, 1/(1+np.exp(alternate.x[1])*(1-phi*phi))))
            converged = bool(alternate.success)
        # Explicitly check the lower-drift search cap; preserve it as a limit,
        # never manufacture a finite drift estimate from that cap.
        cap_phi = float(np.exp(-MIN_STEP))
        cap_fit = _gamma_profile(sample, cap_phi, restricted)
        if cap_fit['nll'] <= nll+1e-7:
            nll, phi, gamma, converged = cap_fit['nll'], cap_phi, cap_fit['gamma'], cap_fit['converged']
    # Explicit phi=0 limiting model. For one replicate it is white for all gamma.
    boundary = _gamma_profile(sample, 0., restricted) if fixed_gamma is None else dict(
        **_profile(sample, 0., fixed_gamma, restricted), gamma=fixed_gamma, converged=True)
    if boundary['nll'] <= nll+1e-8:
        nll, phi, gamma, converged = boundary['nll'], 0., boundary['gamma'], boundary['converged']
    raw_a = float(-np.log(phi)) if phi > 0 else np.inf
    white = phi == 0 or gamma <= 1e-8
    cap = phi > 0 and (raw_a <= MIN_STEP*1.001 or raw_a >= MAX_STEP*.9999)
    grid_better = objective(seeds[0]) < nll-1e-6
    converged = bool(converged and not grid_better)
    identified = bool(converged and not white and not cap)
    details = _profile(sample, phi, gamma, restricted)
    theta = raw_a/dt if identified else None
    return dict(**details, phi=phi, gamma=gamma, theta=theta, identified=identified,
        optimizer_converged=converged, white_limit=bool(white), search_cap=bool(cap),
        noise_boundary=bool(gamma >= upper-1e-7), restricted=bool(restricted),
        latent_variance=None if white and sample.repeats == 1 else gamma*details['scale'],
        noise_variance=None if white and sample.repeats == 1 else (1-gamma)*details['scale'],
        sites=sample.sites, repeats=sample.repeats, dt=float(dt))


def profile_membership(observations, dt, theta, fit, level=.95):
    """Asymptotic profile-LR membership, not a finite-sample guarantee."""
    if not np.isfinite([dt, theta, level]).all() or dt <= 0 or theta <= 0 or not 0 < level < 1:
        raise ValueError('positive dt/theta and level in (0,1) required')
    if fit['dt'] != dt:
        raise ValueError('fit and membership must use the same dt')
    null = profile_phi(observations, float(np.exp(-theta*dt)), fit['restricted'])
    lr = 2*(null['nll']-fit['nll'])
    available = bool(fit['optimizer_converged'] and null['converged'] and lr >= -1e-6)
    return dict(available=available, hit=bool(available and lr <= chi2.ppf(level, 1)),
        lr=float(max(0., lr)) if available else None, null_gamma=null['gamma'],
        optimizer_gap=bool(lr < -1e-6))


def information_matrix(theta, variance, noise_variance, dt, sites, repeats=1):
    """Expected Gaussian Fisher information for (theta,V,eta^2), interior noise.

    Its Schur complement is local information; the CR bound is for regular
    unbiased estimators, not a universal lower bound on empirical RMSE.
    """
    if not np.isfinite([theta, variance, noise_variance, dt]).all() or min(theta, variance, noise_variance, dt) <= 0 or type(sites) is not int or sites < 4 or type(repeats) is not int or repeats < 1:
        raise ValueError('positive parameters, >=4 sites and positive integer repeats required')
    distance = np.abs(np.arange(sites)[:, None]-np.arange(sites))*dt
    r = np.exp(-theta*distance)
    covariance = variance*r+noise_variance/repeats*np.eye(sites)
    derivatives = [-variance*distance*r, r, np.eye(sites)/repeats]
    factor = cho_factor(covariance, lower=True, check_finite=False)
    products = [cho_solve(factor, d, check_finite=False) for d in derivatives]
    info = np.array([[.5*np.einsum('ij,ji->', a, b) for b in products] for a in products])
    info[2, 2] += sites*(repeats-1)/(2*noise_variance**2)
    efficient = float(info[0, 0]-info[0, 1:] @ np.linalg.solve(info[1:, 1:], info[1:, 0]))
    if efficient <= 0:
        raise FloatingPointError('singular efficient information')
    return dict(matrix=info.tolist(), efficient_theta=efficient, local_se_bound=float(1/np.sqrt(efficient)))


def information_limit(phi0, phi1, gamma, sites, alpha=.05):
    """Two-point bounds in a fixed-total-variance submodel, one reading per site."""
    if not np.isfinite([phi0, phi1, gamma, alpha]).all() or not 0 <= phi0 < 1 or not 0 <= phi1 < 1 or not 0 <= gamma <= 1 or not 0 < alpha < .5 or type(sites) is not int or sites < 2:
        raise ValueError('valid stationary slopes, gamma, alpha and integer sites required')
    distance = np.abs(np.arange(sites)[:, None]-np.arange(sites))
    r0, r1 = phi0**distance, phi1**distance
    p, q = gamma*r0+(1-gamma)*np.eye(sites), gamma*r1+(1-gamma)*np.eye(sites)
    def divergence(a, b):
        delta = eigvalsh(a, b, check_finite=False)-1
        small = np.abs(delta) < 1e-4
        terms = np.empty_like(delta)
        x = delta[small]
        terms[small] = x*x/2-x**3/3+x**4/4-x**5/5+x**6/6
        terms[~small] = delta[~small]-np.log1p(delta[~small])
        return float(max(0., terms.sum()/2))
    forward, reverse = divergence(p, q), divergence(q, p)
    tv = float(min(1., np.sqrt(min(forward, reverse)/2)))
    analytic_kl = float(gamma**2*np.sum((r0-r1)**2)/(4*(1-gamma)**2)) if gamma < 1 else None
    return dict(kl_0_to_1=forward, kl_1_to_0=reverse, tv_upper=tv,
        analytic_kl_upper=analytic_kl, power_upper=float(min(1., alpha+tv)),
        simultaneous_membership_lower=float(max(0., 1-2*alpha-tv)),
        worst_error_probability_lower=float((1-tv)/2))
