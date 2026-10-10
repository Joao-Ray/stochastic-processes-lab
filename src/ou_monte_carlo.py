"""Stationary OU rank tests and conservative continuous-cell inversion.

Proofs assume exact arithmetic and independent exact Gaussian draws. Numerical
calculations use floating-point arithmetic and pseudorandom draws. Calibration
is unconditional over the stationary initial state, not conditional on X_0.
"""
from dataclasses import dataclass

import numpy as np
from scipy.signal import lfilter

from src.ou_uncertainty import ar_statistics, _batch_statistics


def _phi(phi, allow_one=False):
    if not np.isfinite(phi) or phi < 0 or phi > 1 or (phi == 1 and not allow_one):
        raise ValueError('phi must be in [0,1)')


def stationary_null_paths(phi, n, count, rng):
    """Canonical AR paths with stationary variance one and mean zero."""
    _phi(phi)
    if type(n) is not int or n < 3 or type(count) is not int or count < 1:
        raise ValueError('n >= 3 and count >= 1 must be integers')
    initial = rng.standard_normal(count)
    noise = np.sqrt((1-phi)*(1+phi))*rng.standard_normal((count, n))
    paths = np.empty((count, n+1))
    paths[:, 0] = initial
    paths[:, 1:] = lfilter([1.], [1., -phi], noise, axis=1)+initial[:, None]*phi**np.arange(1, n+1)
    return paths


def _lr(phi_hat, xx, sse, candidate, n):
    best = np.clip(phi_hat, 0., 1.)
    denominator = sse+xx*(best-phi_hat)**2
    numerator = sse+xx*(candidate-phi_hat)**2
    return np.maximum(n*np.log(numerator/denominator), 0.)


def profile_statistic(series, candidate):
    _phi(candidate)
    stats = ar_statistics(series)
    return float(_lr(stats.phi, stats.xx, stats.sse, candidate, stats.n))


def rank_pvalue(observed, simulated):
    """Plus-one upper-tail p-value; ties are retained conservatively."""
    simulated = np.asarray(simulated, dtype=float)
    if not np.isfinite(observed) or simulated.ndim != 1 or not simulated.size or not np.isfinite(simulated).all():
        raise ValueError('finite observed statistic and nonempty simulated vector required')
    exceedances = int(np.count_nonzero(simulated >= observed))
    return (1+exceedances)/(simulated.size+1), exceedances


def stationary_rank_test(series, candidate_phi, rng, simulations=399):
    """Stationary-null rank test without estimated nuisance parameters."""
    if type(simulations) is not int or simulations < 19:
        raise ValueError('simulations must be an integer >= 19')
    observed = profile_statistic(series, candidate_phi)
    n = len(series)-1
    paths = stationary_null_paths(candidate_phi, n, simulations, rng)
    slope, xx, sse = _batch_statistics(paths)
    null_statistics = _lr(slope, xx, sse, candidate_phi, n)
    pvalue, exceedances = rank_pvalue(observed, null_statistics)
    return dict(pvalue=float(pvalue), exceedances=exceedances,
                statistic=observed, simulations=simulations, candidate_phi=float(candidate_phi))


def stationary_kl(phi, reference_phi, n):
    """KL(P_phi || P_reference) for n stationary canonical AR transitions."""
    _phi(phi)
    _phi(reference_phi)
    if type(n) is not int or n < 1:
        raise ValueError('n must be a positive integer')
    q = (1-reference_phi)*(1+reference_phi)
    x = (reference_phi-phi)*(reference_phi+phi)/q
    if abs(x) < 1e-4:
        # Stable x-log(1+x), including an upper bound on the series remainder.
        variance_part = x*x/2-x**3/3+x**4/4-x**5/5+x**6/6
        variance_part += abs(x)**7/(7*(1-abs(x)))
    else:
        variance_part = x-np.log1p(x)
    return float(n/2*(variance_part+(phi-reference_phi)**2/q))


def cell_tv_bound(lower, upper, center, n):
    """Endpoint KL supremum/Pinsker bound with an outward numerical cushion."""
    _phi(lower)
    _phi(upper, allow_one=True)
    _phi(center)
    if type(n) is not int or n < 3 or not lower <= center <= upper or lower >= upper:
        raise ValueError('n >= 3 and ordered cell bounds containing center required')
    if upper == 1:
        return 1.
    kl = max(stationary_kl(lower, center, n), stationary_kl(upper, center, n))
    return float(min(1., np.sqrt(max(0., kl)/2)+1e-10))


@dataclass(frozen=True)
class PhiCell:
    lower: float
    upper: float
    center: float
    tv_bound: float


def phi_cells(n, resolution=.04, u_max=5.):
    """Data-independent grid in atanh(phi), including the final open-end tail."""
    if type(n) is not int or n < 3 or not np.isfinite([resolution, u_max]).all() or resolution <= 0 or u_max <= 0:
        raise ValueError('n >= 3, positive finite resolution and u_max required')
    if np.tanh(u_max) == 1:
        raise ValueError('u_max too large for floating-point resolution')
    count = int(np.ceil(u_max*np.sqrt(n)/resolution))
    u_edges = np.linspace(0, u_max, count+1)
    edges = np.tanh(u_edges)
    cells = []
    for index in range(count):
        center = 0. if index == 0 else float(np.tanh((u_edges[index]+u_edges[index+1])/2))
        lower, upper = float(edges[index]), float(edges[index+1])
        cells.append(PhiCell(lower, upper, center, cell_tv_bound(lower, upper, center, n)))
    cells.append(PhiCell(float(edges[-1]), 1., float(edges[-1]), 1.))
    return cells


def cell_at(phi, cells):
    _phi(phi)
    for index, cell in enumerate(cells):
        if cell.lower <= phi < cell.upper:
            return index, cell
    raise ValueError('phi outside the supplied grid')


def test_cell(series, cell, rng, simulations=399, alpha=.05):
    if type(simulations) is not int or simulations < 19:
        raise ValueError('simulations must be an integer >= 19')
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError('alpha must be in (0,1)')
    ar_statistics(series)
    required = cell_tv_bound(cell.lower, cell.upper, cell.center, len(series)-1)
    if not np.isfinite(cell.tv_bound) or cell.tv_bound < required:
        raise ValueError('cell must carry a sufficient TV bound for this n')
    effective_alpha = max(0., alpha-cell.tv_bound)
    if effective_alpha == 0:
        return dict(accepted=True, pvalue=None, effective_alpha=0.,
                    tv_bound=cell.tv_bound, simulated=False)
    result = stationary_rank_test(series, cell.center, rng, simulations)
    return dict(**result, accepted=bool(result['pvalue'] > effective_alpha),
                effective_alpha=float(effective_alpha), tv_bound=cell.tv_bound, simulated=True)


def invert_cells(series, dt, seed, simulations=399, alpha=.05, resolution=.04, u_max=5.):
    """Return the union of accepted cells, its enclosing hull, and all decisions.

    The final cell is accepted without testing. The hull is a conservative
    enlargement, not an exact two-sided interval.
    """
    ar_statistics(series)
    if not np.isfinite(dt) or dt <= 0 or type(seed) is not int or seed < 0:
        raise ValueError('positive finite dt and nonnegative integer seed required')
    cells = phi_cells(len(series)-1, resolution, u_max)
    decisions, accepted = [], []
    for index, cell in enumerate(cells):
        result = test_cell(series, cell, np.random.default_rng(np.random.SeedSequence([seed, index])), simulations, alpha)
        decisions.append(dict(index=index, lower=cell.lower, upper=cell.upper, center=cell.center, **result))
        if result['accepted']:
            if accepted and accepted[-1][1] == cell.lower:
                accepted[-1][1] = cell.upper
            else:
                accepted.append([cell.lower, cell.upper])
    components = []
    for lower, upper in reversed(accepted):
        components.append(dict(lower=float(max(0., -np.log(upper)/dt)),
            upper=float(-np.log(lower)/dt) if lower > 0 else None, unbounded_upper=lower == 0))
    hull = dict(lower=components[0]['lower'], upper=components[-1]['upper'],
                unbounded_upper=components[-1]['unbounded_upper'])
    return dict(components=components, hull=hull, decisions=decisions,
                resolution=resolution, u_max=u_max, simulations=simulations, seed=seed,
                alpha=alpha, cells=len(cells), always_accepted_tail=True)
