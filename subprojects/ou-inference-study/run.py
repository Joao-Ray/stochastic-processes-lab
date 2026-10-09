"""Execute the fixed OU uncertainty protocol; save all outer attempts."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import norm

from src.ou_inference import OUFitError, fit_ou
from src.ou_process import simulate_ou_exact
from src.ou_uncertainty import ar_statistics, bootstrap_ou, profile_theta_interval, wald_theta_interval

METHODS = ('wald', 'profile', 'bootstrap_basic', 'bootstrap_profile')


def validate(config):
    for key, minimum in [('replications', 2), ('bootstrap_replications', 19), ('seed', 0)]:
        if type(config[key]) is not int or config[key] < minimum:
            raise ValueError(f'{key} must be an integer >= {minimum}')
    if not np.isfinite(config['level']) or not 0 < config['level'] < 1:
        raise ValueError('level must be in (0,1)')
    for key in ('theta', 'stationary_variance'):
        if not np.isfinite(config[key]) or config[key] <= 0:
            raise ValueError(f'{key} must be positive and finite')
    if not np.isfinite(config['mu']):
        raise ValueError('mu must be finite')
    lengths, steps = config['record_lengths'], config['sampling_steps']
    if not lengths or not steps or len(set(lengths)) != len(lengths) or len(set(steps)) != len(steps):
        raise ValueError('nonempty unique duration and sampling lists required')
    for value in lengths+steps:
        if not np.isfinite(value) or value <= 0:
            raise ValueError('durations and steps must be positive and finite')
    fine = min(steps)
    if any(not np.isclose(s/fine, round(s/fine), rtol=0, atol=1e-8) for s in steps):
        raise ValueError('all sampling steps must be nested on the finest grid')
    if any(r/a < 3 or not np.isclose(r/a, round(r/a), rtol=0, atol=1e-8) for r in lengths for a in steps):
        raise ValueError('durations must divide exactly with at least three transitions')
    if config['measurement_noise_sd']:
        if config['noise_record_length'] not in lengths or config['noise_sampling_step'] not in steps:
            raise ValueError('noise design must occur in clean designs')
        noise = config['measurement_noise_sd']
        if len(set(noise)) != len(noise) or any(not np.isfinite(x) or x <= 0 for x in noise):
            raise ValueError('noise SDs must be unique, positive and finite')
    return config


def coverage(hits):
    hits = np.asarray(hits, dtype=float)
    n = len(hits)
    if not n:
        return dict(n=0, coverage=None, mc_se=None, wilson_95=None)
    p, z = float(hits.mean()), float(norm.ppf(.975))
    center = (p+z*z/(2*n))/(1+z*z/n)
    radius = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return dict(n=n, coverage=p, mc_se=float(np.sqrt(p*(1-p)/n)),
                wilson_95=[float(center-radius), float(center+radius)])


def point_metrics(values, truth):
    errors = np.array([v-truth for v in values if v is not None])
    n = errors.size
    if not n:
        return dict(n=0, bias=None, bias_mc_se=None, rmse=None, rmse_mc_se=None)
    squares = errors**2
    rmse = float(np.sqrt(squares.mean()))
    return dict(n=int(n), bias=float(errors.mean()),
        bias_mc_se=float(errors.std(ddof=1)/np.sqrt(n)) if n > 1 else None,
        rmse=rmse, rmse_mc_se=float(squares.std(ddof=1)/(2*rmse*np.sqrt(n))) if n > 1 and rmse > 0 else None)


def difference(values):
    values = np.asarray(values, dtype=float)
    return dict(n=len(values), difference=float(values.mean()) if len(values) else None,
        mc_se=float(values.std(ddof=1)/np.sqrt(len(values))) if len(values) > 1 else None)


def record_interval(row, method, interval, theta):
    available = interval is not None
    row[method+'_available'] = available
    row[method+'_lower'] = interval.lower if available else None
    row[method+'_upper'] = interval.upper if available else None
    row[method+'_empty'] = interval.empty if available else None
    row[method+'_hit'] = interval.contains(theta) if available else False


def summarize(rows, scenario, theta):
    methods = {}
    for method in METHODS:
        available = [r for r in rows if r[method+'_available']]
        finite = [r for r in available if not r[method+'_empty'] and np.isfinite(r[method+'_upper'])]
        widths = [r[method+'_upper']-r[method+'_lower'] for r in finite]
        unbounded = sum(not np.isfinite(r[method+'_upper']) for r in available)
        methods[method] = dict(available=len(available), attempted=len(rows),
            unavailable=len(rows)-len(available), empty=sum(bool(r[method+'_empty']) for r in available),
            unbounded_upper=unbounded,
            conditional=coverage([r[method+'_hit'] for r in available]),
            success_and_coverage=coverage([r[method+'_hit'] for r in rows]),
            finite_nonempty_width_n=len(widths),
            finite_nonempty_mean_width=float(np.mean(widths)) if widths else None,
            finite_nonempty_width_mc_se=float(np.std(widths, ddof=1)/np.sqrt(len(widths))) if len(widths) > 1 else None,
            available_mean_width_is_infinite=bool(unbounded))
    common = [r for r in rows if r['theta'] is not None and r['corrected_theta'] is not None]
    paired = [r for r in rows if r['bootstrap_profile_available'] and r['profile_available']]
    correction_mse_difference = difference([(r['corrected_theta']-theta)**2-(r['theta']-theta)**2 for r in common])
    return dict(**scenario, attempted=len(rows), interior_fits=sum(r['theta'] is not None for r in rows),
        mle=point_metrics([r['theta'] for r in rows], theta),
        slope=point_metrics([r['phi'] for r in rows], np.exp(-scenario['a'])),
        corrected=point_metrics([r['corrected_theta'] for r in rows], theta),
        common_point_estimates=dict(mle=point_metrics([r['theta'] for r in common], theta),
                                    corrected=point_metrics([r['corrected_theta'] for r in common], theta)),
        paired_correction_minus_mle_mse=correction_mse_difference,
        bootstrap_attempted=sum(r['bootstrap_count'] for r in rows),
        bootstrap_inadmissible_slopes=sum(r['bootstrap_bad_slopes'] for r in rows), methods=methods,
        paired_bootstrap_minus_profile= difference([int(r['bootstrap_profile_hit'])-int(r['profile_hit']) for r in paired]),
        joint_bootstrap_minus_profile=difference([int(r['bootstrap_profile_hit'])-int(r['profile_hit']) for r in rows]))


def study(config):
    validate(config)
    theta, mu, v = (config[k] for k in ('theta', 'mu', 'stationary_variance'))
    sigma, fine = np.sqrt(2*theta*v), min(config['sampling_steps'])
    replications, seed = config['replications'], config['seed']
    records, scenarios, paired_grids = [], [], []
    scenario_index = 0
    for group, r in enumerate(config['record_lengths']):
        rng = np.random.default_rng(np.random.SeedSequence([seed, 0, group]))
        steps = int(round(r/fine))
        initial = mu+np.sqrt(v)*rng.standard_normal(replications)
        latent = simulate_ou_exact(theta, mu, sigma, mu, fine/theta, steps, replications, rng)
        latent += (initial-mu)[:, None]*np.exp(-fine*np.arange(steps+1))
        designs = [(a, 0.) for a in config['sampling_steps']]
        if r == config['noise_record_length']:
            designs += [(config['noise_sampling_step'], eta) for eta in config['measurement_noise_sd']]
        standardized_noise = None
        rows_by_grid = {}
        for a, eta in designs:
            dt = a/theta
            training = latent[:, ::int(round(a/fine))]
            if eta:
                if standardized_noise is None:
                    noise_rng = np.random.default_rng(np.random.SeedSequence([seed, 1, group]))
                    standardized_noise = noise_rng.standard_normal(training.shape)
                training = training+eta*standardized_noise
            scenario = dict(id=f'r={r:g},a={a:g},eta={eta:g}', r=r, a=a, noise_sd=eta,
                duration=r/theta, dt=dt, n_transitions=training.shape[1]-1,
                pseudo_theta=theta+np.log1p(eta**2/v)/dt)
            rows = []
            for index, series in enumerate(training):
                stats = ar_statistics(series)
                row = dict(scenario=scenario['id'], replication=index, phi=stats.phi,
                    theta=None, corrected_theta=None, failure='', correction_failure='',
                    bootstrap_count=0, bootstrap_bad_slopes=0, bootstrap_critical=None,
                    bootstrap_mean_phi=None, corrected_phi=None)
                profile = profile_theta_interval(series, dt, config['level'])
                intervals = dict(profile=profile, wald=None, bootstrap_basic=None, bootstrap_profile=None)
                try:
                    fitted = fit_ou(series, dt)
                    row['theta'] = fitted.theta
                    intervals['wald'] = wald_theta_interval(series, dt, config['level'])
                    bootstrap_rng = np.random.default_rng(np.random.SeedSequence([seed, 2, scenario_index, index]))
                    boot = bootstrap_ou(series, dt, bootstrap_rng, config['bootstrap_replications'], config['level'])
                    intervals.update(bootstrap_basic=boot['basic'], bootstrap_profile=boot['calibrated_profile'])
                    row.update(corrected_theta=boot['corrected_theta'], bootstrap_count=boot['replications'],
                        bootstrap_bad_slopes=boot['inadmissible_slopes'], bootstrap_critical=boot['critical_value'],
                        bootstrap_mean_phi=boot['bootstrap_mean_phi'], corrected_phi=boot['corrected_phi'])
                    if boot['corrected_theta'] is None:
                        row['correction_failure'] = 'corrected slope outside (0,1)'
                except OUFitError as exc:
                    row['failure'] = str(exc)
                for method in METHODS:
                    record_interval(row, method, intervals[method], theta)
                rows.append(row)
            scenarios.append(summarize(rows, scenario, theta))
            records.extend(rows)
            if eta == 0:
                rows_by_grid[a] = rows
            scenario_index += 1
            print(f"Completed {scenario['id']}: {len(rows)} attempts", flush=True)
        if len(rows_by_grid) > 1:
            reference = rows_by_grid[fine]
            for a, rows in rows_by_grid.items():
                if a == fine:
                    continue
                common = [(fine_row, coarse) for fine_row, coarse in zip(reference, rows)
                          if fine_row['theta'] is not None and coarse['theta'] is not None]
                paired_grids.append(dict(r=r, fine_a=fine, coarse_a=a,
                    squared_error_fine_minus_coarse=difference([
                        (f['theta']-theta)**2-(c['theta']-theta)**2 for f, c in common])))
    return dict(config=config, scenarios=scenarios, paired_sampling_comparisons=paired_grids,
        data_source='Synthetic stationary exact OU; optional independent Gaussian measurement error',
        seed_scheme='SeedSequence([seed, role, group_or_scenario, optional_replication]); role 0 latent, 1 noise, 2 bootstrap',
        bootstrap_initialization='Observed first value, fixed conditional on each training series',
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      matplotlib=matplotlib.__version__)), records


def plots(result, folder):
    cfg = result['config']
    clean = [s for s in result['scenarios'] if s['noise_sd'] == 0]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    colors = dict(zip(METHODS, plt.get_cmap('tab10').colors[:4]))
    for method in METHODS:
        reference = sorted((s for s in clean if s['a'] == min(cfg['sampling_steps'])), key=lambda s:s['r'])
        x = [s['r'] for s in reference]
        y = [s['methods'][method]['conditional']['coverage'] for s in reference]
        se = [s['methods'][method]['conditional']['mc_se'] for s in reference]
        axes[0].errorbar(x, y, yerr=np.array(se)*1.96, fmt='o-', capsize=3,
                         color=colors[method], label=method)
        reference = [s for s in result['scenarios'] if s['r'] == cfg['noise_record_length'] and s['a'] == cfg['noise_sampling_step']]
        x = [s['noise_sd'] for s in reference]
        y = [s['methods'][method]['conditional']['coverage'] for s in reference]
        se = [s['methods'][method]['conditional']['mc_se'] for s in reference]
        axes[1].errorbar(x, y, yerr=np.array(se)*1.96, fmt='o-', capsize=3,
                         color=colors[method], label=method)
    axes[0].set(xlabel='Dimensionless duration theta T', ylabel='Coverage among available sets',
                ylim=(0,1.02), title=f"Sampling a={min(cfg['sampling_steps'])}; 1.96 MC SE")
    axes[1].set(xlabel='Measurement-noise SD', ylabel='Coverage of latent theta',
                ylim=(0,1.02), title='Noise-free procedures under misspecification')
    for ax in axes:
        ax.axhline(cfg['level'], color='black', ls='--', lw=.8)
        ax.legend(fontsize=8)
    fig.tight_layout()
    for suffix in ('png', 'pdf'):
        fig.savefig(folder/f'coverage_and_noise.{suffix}', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for a in cfg['sampling_steps']:
        rows = sorted((s for s in clean if s['a'] == a), key=lambda s:s['r'])
        axes[0].plot([s['r'] for s in rows], [s['mle']['rmse']/cfg['theta'] for s in rows], 'o-', label=f'a={a}')
    rows = sorted((s for s in clean if s['a'] == min(cfg['sampling_steps'])), key=lambda s:s['r'])
    for key in ('mle', 'corrected'):
        axes[1].plot([s['r'] for s in rows], [s['common_point_estimates'][key]['rmse']/cfg['theta'] for s in rows], 'o-', label=key)
    axes[0].set(xlabel='theta T', ylabel='MLE RMSE / theta', title='Density versus duration: interior fits')
    axes[1].set(xlabel='theta T', ylabel='RMSE / theta', title='Bias correction: common admissible subset')
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    for suffix in ('png', 'pdf'):
        fig.savefig(folder/f'drift_estimation.{suffix}', dpi=180)
    plt.close(fig)


def report(result, folder):
    cfg = result['config']
    lines = ['# OU uncertainty: executed benchmark', '',
        '[Protocol](../../PROTOCOL.md) · [Derivations](../../THEORY.md) · [All attempts](replications.csv.gz) · [Full summary](summary.json)', '',
        f"Seed **{cfg['seed']}**; **{cfg['replications']:,}** outer paths per scenario; **B={cfg['bootstrap_replications']}**. All data are synthetic.", '',
        '## Point estimation', '',
        'Bias/RMSE are conditional on admissibility. The two common-subset RMSEs use exactly the same records. Their denominator can differ from the MLE total.', '',
        '| Scenario | Interior / attempted | MLE bias (MC SE) | MLE RMSE (MC SE) | Correction valid | Common MLE / corrected RMSE |',
        '| --- | --- | --- | --- | --- | --- |']
    def number(x): return 'undefined' if x is None else f'{x:.4f}'
    for s in result['scenarios']:
        m, c = s['mle'], s['common_point_estimates']
        lines.append(f"| {s['id']} | {s['interior_fits']} / {s['attempted']} | {number(m['bias'])} ({number(m['bias_mc_se'])}) | "
            f"{number(m['rmse'])} ({number(m['rmse_mc_se'])}) | {s['corrected']['n']} | {number(c['mle']['rmse'])} / {number(c['corrected']['rmse'])} (n={c['mle']['n']}) |")
    lines += ['', '### Paired squared-error change under slope correction', '',
        'Corrected minus MLE MSE on the common admissible subset; negative differences favor correction. Availability is still reported separately.', '',
        '| Scenario | Common n | Paired MSE change | Paired MC SE |',
        '| --- | --- | --- | --- |']
    for s in result['scenarios']:
        p = s['paired_correction_minus_mle_mse']
        lines.append(f"| {s['id']} | {p['n']} | {number(p['difference'])} | {number(p['mc_se'])} |")
    lines += ['', '## Confidence sets', '',
        'Conditional coverage uses available sets. Success-and-coverage uses all attempts, counting an unavailable procedure as unsuccessful. Empty sets remain available with coverage zero. Finite widths exclude empty/unbounded sets and have an explicit denominator.', '',
        '| Scenario | Method | Available | Conditional coverage [Wilson 95%] | MC SE | Success-and-coverage | Empty / unbounded | Finite mean width (n) |',
        '| --- | --- | --- | --- | --- | --- | --- | --- |']
    for s in result['scenarios']:
        for method in METHODS:
            m = s['methods'][method]
            c = m['conditional']
            ci = c['wilson_95']
            ci_text = f"{number(c['coverage'])} [{number(ci[0])}, {number(ci[1])}]" if ci else 'undefined'
            lines.append(f"| {s['id']} | {method} | {m['available']}/{m['attempted']} | {ci_text} | {number(c['mc_se'])} | "
                f"{number(m['success_and_coverage']['coverage'])} | {m['empty']} / {m['unbounded_upper']} | {number(m['finite_nonempty_mean_width'])} (n={m['finite_nonempty_width_n']}) |")
    lines += ['', '## Paired bootstrap-profile minus chi-square-profile coverage', '',
        'Same available records; differences and MC SE account for pairing. These descriptive comparisons are not multiplicity-adjusted superiority tests.', '',
        '| Scenario | Common n | Coverage difference | Paired MC SE | All-attempt joint-success difference |',
        '| --- | --- | --- | --- | --- |']
    for s in result['scenarios']:
        p = s['paired_bootstrap_minus_profile']
        lines.append(f"| {s['id']} | {p['n']} | {number(p['difference'])} | {number(p['mc_se'])} | {number(s['joint_bootstrap_minus_profile']['difference'])} |")
    lines += ['', '## Paired sampling-density comparisons', '',
        'Negative squared-error differences favor the fine grid. Fits must be admissible on both grids.', '',
        '| theta T | Fine a | Coarse a | Common n | MSE fine - coarse | Paired MC SE |',
        '| --- | --- | --- | --- | --- | --- |']
    for s in result['paired_sampling_comparisons']:
        p=s['squared_error_fine_minus_coarse']
        lines.append(f"| {s['r']} | {s['fine_a']} | {s['coarse_a']} | {p['n']} | {number(p['difference'])} | {number(p['mc_se'])} |")
    lines += ['', '## Observation-noise benchmark', '',
        '| Noise SD | True theta | Expanding-domain pseudo-theta | Recorded MLE mean |',
        '| --- | --- | --- | --- |']
    for s in result['scenarios']:
        if s['noise_sd']:
            lines.append(f"| {s['noise_sd']} | {cfg['theta']} | {s['pseudo_theta']:.4f} | {number(cfg['theta']+s['mle']['bias']) if s['mle']['bias'] is not None else 'undefined'} |")
    lines += ['', 'Pseudo-theta is a long-record limit of the misspecified regression, not a finite-sample target or replacement truth.', '',
        '## Figures', '', '![Coverage and observation noise](coverage_and_noise.png)', '',
        '![Drift estimation](drift_estimation.png)', '',
        '[Vector coverage figure](coverage_and_noise.pdf) · [Vector estimation figure](drift_estimation.pdf)', '',
        '## Limitations and reproducibility', '',
        '- This benchmarks established methods; it does not establish a new estimator or a publication-ready contribution.',
        '- Coverage is not guaranteed by the use of a chi-square quantile or bootstrap. Unavailable fits and corrected estimates are retained.',
        '- Interval endpoints at zero and infinity denote closures of the physical domain. CSV uses inf for infinite upper limits; JSON reports unbounded-set counts and flags without a finite surrogate width.',
        '- Noise stress cases deliberately fit the wrong observation model. No noise-aware state-space inference is implemented here.',
        '- The confirmation run uses a new seed and larger B, jointly; it cannot isolate either factor. Main/confirmation results are not pooled.',
        '- All bootstrap slopes, including inadmissible ones, are used for the calculations. The archive saves their counts and means per outer record; full inner paths are regenerated from the seed scheme.',
        '- Code/configuration/protocol SHA256 hashes and software versions are in summary.json. Runtime is machine-dependent.', '']
    (folder/'report.md').write_text('\n'.join(lines))


def save(result, records, folder):
    folder.mkdir(parents=True, exist_ok=True)
    csv_path = folder/'replications.csv.gz'
    with csv_path.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding='utf-8', newline='') as text:
                writer = csv.DictWriter(text, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)
    result['record_archive_sha256'] = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    (folder/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    plots(result, folder)
    report(result, folder)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('config.json'))
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('results')/'main')
    args = parser.parse_args()
    started = time.perf_counter()
    config = json.loads(args.config.read_text())
    result, records = study(config)
    source_paths = [Path(__file__), ROOT/'src/ou_uncertainty.py', ROOT/'src/ou_inference.py',
                    ROOT/'src/ou_process.py', Path(__file__).with_name('PROTOCOL.md'),
                    Path(__file__).with_name('THEORY.md')]
    result['source_sha256'] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths}
    result['config_sha256'] = hashlib.sha256(args.config.read_bytes()).hexdigest()
    result['elapsed_study_seconds'] = time.perf_counter()-started
    save(result, records, args.output)
    print(f"Saved {len(records):,} outer attempts to {args.output.resolve()}", flush=True)


if __name__ == '__main__':
    main()
