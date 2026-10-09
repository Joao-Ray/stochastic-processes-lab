"""Numerically verify the OU theory with recorded independent benchmarks."""
from pathlib import Path
import argparse
import json
import platform
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize
from scipy.stats import chi2, norm

from src.ou_inference import conditional_ou_nll, fit_ou, forecast_ou
from src.ou_process import ou_mean, ou_variance, simulate_ou_exact


def wilson(hits, n):
    p, z = hits/n, norm.ppf(.975)
    center = (p+z*z/(2*n))/(1+z*z/n)
    width = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return [float(center-width), float(center+width)]


def verification(config, paths=50_000, seed=20261009):
    if type(paths) is not int or paths < 2 or type(seed) is not int or seed < 0:
        raise ValueError('paths must be at least two and seed must be nonnegative')
    theta, mu, sigma = (float(config['truth'][k]) for k in ('theta','mu','sigma'))
    if not np.isfinite([theta,mu,sigma]).all() or theta <= 0 or sigma <= 0:
        raise ValueError('truth parameters must be finite with theta and sigma positive')
    checks = []

    def check(name, error, tolerance):
        checks.append(dict(name=name, error=float(error), tolerance=float(tolerance),
                           passed=bool(np.isfinite(error) and error <= tolerance)))

    integral_rows = []
    for speed in (.05, .7, 3.):
        for h in (.01, .2, 2.):
            integral, quadrature_error = quad(lambda u: sigma**2*np.exp(-2*speed*u), 0, h,
                                               epsabs=1e-13, epsrel=1e-13)
            formula = float(ou_variance(h, speed, sigma))
            error = abs(integral-formula)
            integral_rows.append(dict(theta=speed, h=h, quadrature=integral, formula=formula,
                absolute_error=error, relative_error=error/integral, quadrature_error_estimate=quadrature_error))
            check(f'ito_variance_theta_{speed}_h_{h}', error, 1e-11)

    horizons = np.array([.1,.5,1.,3.])
    dt, x0 = .1, mu-2.5
    indices = np.rint(horizons/dt).astype(int)
    simulations = simulate_ou_exact(theta,mu,sigma,x0,dt,int(indices.max()),paths,np.random.default_rng(seed))
    moment_rows = []
    for index, h in zip(indices,horizons):
        endpoint = simulations[:,index]
        mean, variance = float(ou_mean(h,theta,mu,x0)), float(ou_variance(h,theta,sigma))
        empirical_mean, empirical_variance = float(endpoint.mean()), float(endpoint.var(ddof=1))
        mean_se = np.sqrt(variance/paths)
        variance_se = variance*np.sqrt(2/(paths-1))
        _, lower, upper = forecast_ou(theta,mu,sigma,x0,np.array([h]))
        hits = int(((endpoint >= lower[0]) & (endpoint <= upper[0])).sum())
        oracle_mse = float(np.mean((endpoint-mean)**2))
        persistence_mse = float(np.mean((endpoint-x0)**2))
        displacement = mean-x0
        theory_persistence = variance+displacement**2
        oracle_se = variance*np.sqrt(2/paths)
        persistence_se = np.sqrt((2*variance**2+4*displacement**2*variance)/paths)
        moment_rows.append(dict(h=float(h), theory_mean=mean, empirical_mean=empirical_mean,
            mean_error=empirical_mean-mean, mean_mc_se=float(mean_se),
            mean_ci_95=[float(empirical_mean-norm.ppf(.975)*mean_se),
                        float(empirical_mean+norm.ppf(.975)*mean_se)],
            theory_variance=variance, empirical_variance=empirical_variance,
            variance_error=empirical_variance-variance, variance_mc_se=float(variance_se),
            variance_ci_95=[float((paths-1)*empirical_variance/chi2.ppf(.975,paths-1)),
                            float((paths-1)*empirical_variance/chi2.ppf(.025,paths-1))],
            oracle_coverage=hits/paths, coverage_ci_95=wilson(hits,paths), coverage_hits=hits,
            oracle_mse=oracle_mse, theory_oracle_mse=variance, oracle_mse_mc_se=float(oracle_se),
            persistence_mse=persistence_mse, theory_persistence_mse=theory_persistence,
            persistence_mse_mc_se=float(persistence_se)))
        for label, error, standard_error in [
            ('mean', empirical_mean-mean, mean_se),
            ('variance', empirical_variance-variance, variance_se),
            ('oracle_coverage', hits/paths-.95, np.sqrt(.95*.05/paths)),
            ('oracle_mse', oracle_mse-variance, oracle_se),
            ('persistence_mse', persistence_mse-theory_persistence, persistence_se)]:
            check(f'{label}_h_{h}', abs(error)/standard_error, 6.)
    del simulations

    likelihood_rows = []
    likelihood_seeds = list(range(seed+100,seed+108))
    starts = [[np.log(.4),.5,np.log(.5)], [np.log(1.3),2.5,np.log(1.)]]
    bounds = [(np.log(1e-4),np.log(20.)),(-10.,10.),(np.log(1e-4),np.log(10.))]
    for likelihood_seed in likelihood_seeds:
        series = simulate_ou_exact(theta,mu,sigma,mu,.2,1000,1,np.random.default_rng(likelihood_seed))[0]
        fitted = fit_ou(series,.2)

        def objective(parameters):
            return conditional_ou_nll(series,.2,np.exp(parameters[0]),parameters[1],np.exp(parameters[2]))

        optimized = [minimize(objective,start,method='L-BFGS-B',bounds=bounds,
            options=dict(ftol=1e-13,gtol=1e-8,maxiter=1000)) for start in starts]
        best = min(optimized,key=lambda result: result.fun)
        closed = np.array([fitted.theta,fitted.mu,fitted.sigma])
        numerical = np.array([np.exp(best.x[0]),best.x[1],np.exp(best.x[2])])
        gap = abs(best.fun-fitted.negative_log_likelihood)
        parameter_error = float(np.max(abs(numerical-closed)/np.maximum(abs(closed),1.)))
        residual = series[1:]-fitted.intercept-fitted.phi*series[:-1]
        normal_equations_error = float(max(abs(residual.mean()),abs(np.mean(series[:-1]*residual)),
                                           abs(np.mean(residual**2)-fitted.innovation_variance)))
        likelihood_rows.append(dict(seed=likelihood_seed,closed_parameters=closed.tolist(),
            numerical_parameters=numerical.tolist(),closed_nll=fitted.negative_log_likelihood,
            numerical_nll=float(best.fun),absolute_nll_gap=float(gap),scaled_parameter_error=parameter_error,
            normal_equations_error=normal_equations_error,
            optimizer_runs=[dict(success=bool(r.success),status=int(r.status),message=str(r.message),
                                 nll=float(r.fun),iterations=int(r.nit)) for r in optimized]))
        check(f'mle_objective_seed_{likelihood_seed}',gap,1e-6)
        check(f'mle_parameters_seed_{likelihood_seed}',parameter_error,1e-4)
        check(f'mle_normal_equations_seed_{likelihood_seed}',normal_equations_error,1e-10)

    phi = np.exp(-theta*dt)
    q = quad(lambda u: sigma**2*np.exp(-2*theta*u),0,dt)[0]
    mean, variance = x0, 0.
    composition_error = 0.
    for step in range(1,51):
        mean, variance = mu+phi*(mean-mu), phi**2*variance+q
        composition_error = max(composition_error, abs(mean-float(ou_mean(step*dt,theta,mu,x0))),
                                abs(variance-float(ou_variance(step*dt,theta,sigma))))
    check('transition_composition_50_steps',composition_error,1e-11)

    delta, shift_h = float(config['mean_shift']), 5.
    if not np.isfinite(delta):
        raise ValueError('mean_shift must be finite')
    base, shifted = np.full(paths,x0), np.full(paths,x0)
    rng = np.random.default_rng(seed+200)
    path_difference_error = 0.
    for step in range(1,51):
        noise = np.sqrt(q)*rng.standard_normal(paths)
        base = mu+phi*(base-mu)+noise
        shifted = mu+delta+phi*(shifted-mu-delta)+noise
        predicted_difference = delta*(1-np.exp(-theta*step*dt))
        path_difference_error = max(path_difference_error,float(np.max(abs(shifted-base-predicted_difference))))
    ode = solve_ivp(lambda t,d: theta*(delta-d),(0,shift_h),[0.],rtol=1e-10,atol=1e-12)
    theoretical_difference = delta*(1-np.exp(-theta*shift_h))
    ode_error = abs(ode.y[0,-1]-theoretical_difference)
    check('paired_mean_shift_50_steps',path_difference_error,1e-11)
    check('mean_shift_numerical_ode',ode_error,1e-8)
    shift_variance = float(ou_variance(shift_h,theta,sigma))
    _, old_lower, old_upper = forecast_ou(theta,mu,sigma,x0,np.array([shift_h]))
    sd, z = np.sqrt(shift_variance), norm.ppf(.975)
    analytic_coverage = float(norm.cdf(z-theoretical_difference/sd)-norm.cdf(-z-theoretical_difference/sd))
    hits = int(((shifted>=old_lower[0]) & (shifted<=old_upper[0])).sum())
    coverage_se = np.sqrt(analytic_coverage*(1-analytic_coverage)/paths)
    check('old_interval_shifted_coverage',abs(hits/paths-analytic_coverage)/coverage_se,6.)
    shifted_result = dict(delta=delta,h=shift_h,theory_difference=theoretical_difference,
        max_coupled_path_difference_error=path_difference_error,numerical_ode_error=float(ode_error),
        analytic_old_interval_coverage=analytic_coverage,empirical_old_interval_coverage=hits/paths,
        coverage_mc_se=float(coverage_se),coverage_ci_95=wilson(hits,paths),hits=hits)
    return dict(seed=seed,paths=paths,truth=dict(theta=theta,mu=mu,sigma=sigma),x0=x0,
        moment_horizons=horizons.tolist(),likelihood_seeds=likelihood_seeds,shift_seed=seed+200,
        likelihood_record_transitions=1000,likelihood_record_dt=.2,optimizer_starts=starts,optimizer_bounds=bounds,
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        variance_quadrature=integral_rows,moments_and_forecasts=moment_rows,likelihood=likelihood_rows,
        transition_composition_error=composition_error,shift=shifted_result,checks=checks,
        all_checks_passed=all(c['passed'] for c in checks))


def save(result, folder):
    folder.mkdir(parents=True,exist_ok=True)
    (folder/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    rows=result['moments_and_forecasts']
    fig,axes=plt.subplots(1,2,figsize=(11,4.2))
    for key,marker in [('oracle','o'),('persistence','s')]:
        axes[0].errorbar([r['h'] for r in rows],[r[key+'_mse'] for r in rows],
                        yerr=[1.96*r[key+'_mse_mc_se'] for r in rows],fmt=marker+'-',capsize=3,
                        label=f'{key}: estimate + 1.96 MC SE')
        axes[0].plot([r['h'] for r in rows],[r['theory_'+key+'_mse'] for r in rows],'--',label=f'{key}: theory')
    axes[0].set(xlabel='Forecast horizon',ylabel='Conditional MSE',title='Theory versus 50,000-path Monte Carlo' if result['paths']==50000 else f"Theory versus {result['paths']:,}-path Monte Carlo")
    for key,marker in [('mean','o'),('variance','s')]:
        axes[1].plot([r['h'] for r in rows],[r[key+'_error']/r[key+'_mc_se'] for r in rows],marker+'-',label=key)
    axes[1].axhline(0,color='black',lw=.8)
    axes[1].axhline(1.96,color='gray',ls='--')
    axes[1].axhline(-1.96,color='gray',ls='--',label='Pointwise normal reference +/-1.96')
    axes[1].set(xlabel='Horizon',ylabel='Error / theoretical MC SE',title='Sampling error is retained')
    for ax in axes:
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(folder/'forecast_mse_verification.png',dpi=140)
    plt.close(fig)
    lines=['# OU Numerical Verification: Executed Report','',
        'Derivations and assumptions: [THEORY.md](../../THEORY.md). Full results: [results.json](results.json).','',
        f"Seed {result['seed']}; {result['paths']:,} paths; eight independently seeded likelihood records.",'',
        '## Itô variance versus adaptive quadrature','',
        '| theta | h | Formula | Quadrature | Absolute error |','| --- | --- | --- | --- | --- |']
    for r in result['variance_quadrature']:
        lines.append(f"| {r['theta']} | {r['h']} | {r['formula']:.9f} | {r['quadrature']:.9f} | {r['absolute_error']:.2e} |")
    lines+=['','## Conditional moments and known-parameter prediction intervals','',
        '| h | Mean theory / empirical | Variance theory / empirical | Mean error / SE | Variance error / SE | Coverage [Wilson 95%] |',
        '| --- | --- | --- | --- | --- | --- |']
    for r in rows:
        ci=r['coverage_ci_95']
        lines.append(f"| {r['h']} | {r['theory_mean']:.6f} / {r['empirical_mean']:.6f} | "
            f"{r['theory_variance']:.6f} / {r['empirical_variance']:.6f} | {r['mean_error']/r['mean_mc_se']:.3f} | "
            f"{r['variance_error']/r['variance_mc_se']:.3f} | {r['oracle_coverage']:.4f} [{ci[0]:.4f}, {ci[1]:.4f}] |")
    lines+=['','Moments condition on the fixed initial value x0. Horizons share paths, so errors across rows are correlated. These are known-parameter intervals, not fitted plug-in intervals.','',
        '## Conditional MSE: optimal prediction versus persistence','',
        '| h | Oracle theory / empirical | Persistence theory / empirical |','| --- | --- | --- |']
    misses = [r for r in rows if not r['mean_ci_95'][0] <= r['theory_mean'] <= r['mean_ci_95'][1]]
    lines[-4:-4] = [f"Pointwise 95% mean intervals miss their theoretical target at {len(misses)} of {len(rows)} horizons: "
        + (', '.join(f"h={r['h']} ({r['mean_error']/r['mean_mc_se']:.2f} SE)" for r in misses) if misses else 'none')
        + '. These outcomes are retained; shared horizons are not independent tests.', '']
    for r in rows:
        lines.append(f"| {r['h']} | {r['theory_oracle_mse']:.6f} / {r['oracle_mse']:.6f} | {r['theory_persistence_mse']:.6f} / {r['persistence_mse']:.6f} |")
    lines+=['','![MSE verification and moment sampling errors](forecast_mse_verification.png)','',
        '## Conditional MLE: closed form versus numerical minimization','',
        'Each numerical reference uses two fixed starts in log(theta), mu, log(sigma), with wide bounds recorded in JSON. No start uses the closed-form estimate. Optimization statuses for both starts are retained; agreement is judged by objective and parameter errors, not the status flag alone.','',
        '| Seed | Absolute NLL gap | Max scaled parameter error | Normal equations error | Start convergence flags |',
        '| --- | --- | --- | --- | --- |']
    for r in result['likelihood']:
        flags=', '.join(str(o['success']) for o in r['optimizer_runs'])
        lines.append(f"| {r['seed']} | {r['absolute_nll_gap']:.2e} | {r['scaled_parameter_error']:.2e} | {r['normal_equations_error']:.2e} | {flags} |")
    s=result['shift']
    lines+=['','## Composition and paired mean shift','',
        f"- Maximum 50-step recursive-moment discrepancy: **{result['transition_composition_error']:.2e}**.",
        f"- Maximum coupled-path mean-shift identity error: **{s['max_coupled_path_difference_error']:.2e}**.",
        f"- Numerical ODE versus closed-form displacement error: **{s['numerical_ode_error']:.2e}**.",
        f"- Old known-parameter interval coverage after shift: analytic **{s['analytic_old_interval_coverage']:.4f}**, simulated **{s['empirical_old_interval_coverage']:.4f}**, MC SE **{s['coverage_mc_se']:.4f}**.",'',
        'The shifted-coverage formula is a known-parameter benchmark; it does not claim exact plug-in coverage.','',
        '## Verification gates and limitations','',
        f"Passed **{sum(c['passed'] for c in result['checks'])}/{len(result['checks'])}** recorded checks.",'',
        'Deterministic tolerances and six-standard-error distributional gates are fixed in verify.py. They detect implementation regressions. Six-SE gates are not 95% scientific significance tests. The reported 95% intervals and all seed outcomes are retained whether or not they contain the target.',
        '', 'The simulations are synthetic; agreement does not validate OU assumptions for observed data. Quadrature checks integration and formula implementation; likelihood optimization checks the interior estimator. Neither replaces the derivations.', '']
    (folder/'report.md').write_text('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('config.json'))
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('results')/'verification')
    parser.add_argument('--paths',type=int,default=50_000)
    parser.add_argument('--seed',type=int,default=20261009)
    args=parser.parse_args()
    result=verification(json.loads(args.config.read_text()),args.paths,args.seed)
    save(result,args.output)
    print(f"Verification: {sum(c['passed'] for c in result['checks'])}/{len(result['checks'])} checks passed.")
    print(f'Report: {(args.output/"report.md").resolve()}')
    if not result['all_checks_passed']:
        raise SystemExit('Numerical verification failed; inspect results.json for recorded discrepancies.')


if __name__=='__main__':
    main()
