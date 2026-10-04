"""Reproduce the OU calibration subproject from synthetic exact OU paths."""
from pathlib import Path
import argparse
import csv
from dataclasses import asdict
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
from scipy.stats import norm

from src.ou_inference import OUFitError, fit_ou, forecast_ou
from src.ou_process import simulate_ou_exact


def coverage_summary(hits):
    hits = np.asarray(hits, dtype=bool)
    n = hits.size
    if n == 0:
        return dict(n=0, coverage=None, lower=None, upper=None)
    p, z = float(hits.mean()), norm.ppf(.975)
    center = (p+z*z/(2*n))/(1+z*z/n)
    width = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return dict(n=int(n), coverage=p, lower=float(center-width), upper=float(center+width))


def checked_config(config):
    truth = config['truth']
    if not all(np.isfinite(truth[k]) for k in ('theta', 'mu', 'sigma')) or truth['theta'] <= 0 or truth['sigma'] <= 0:
        raise ValueError('truth requires finite parameters, theta > 0 and sigma > 0')
    for key in ('replications', 'seed'):
        if type(config[key]) is not int or config[key] < (1 if key == 'replications' else 0):
            raise ValueError(f'{key} must be an integer in its valid range')
    base_dt = config['base_dt']
    if not np.isfinite(base_dt) or base_dt <= 0:
        raise ValueError('base_dt must be positive and finite')
    durations, steps, horizons = config['training_durations'], config['sampling_steps'], config['forecast_horizons']
    if not durations or not steps or not horizons:
        raise ValueError('duration, sampling-step and forecast lists must be nonempty')
    if len(set(durations)) != len(durations) or len(set(steps)) != len(steps):
        raise ValueError('duration and sampling-step lists must be unique')
    if config['sampling_duration'] not in durations:
        raise ValueError('sampling_duration must occur in training_durations')
    for value in durations + steps + horizons:
        if not np.isfinite(value) or value <= 0 or not np.isclose(value/base_dt, round(value/base_dt), atol=1e-8, rtol=0):
            raise ValueError('durations, steps and horizons must be positive grid multiples of base_dt')
    for step in steps:
        if step < base_dt or not np.isclose(config['sampling_duration']/step, round(config['sampling_duration']/step), atol=1e-8, rtol=0):
            raise ValueError('sampling steps must divide sampling_duration and be at least base_dt')
    if min(durations)/base_dt < 3 or min(config['sampling_duration']/s for s in steps) < 3:
        raise ValueError('every fit requires at least three transitions')
    if not np.isfinite(config['prediction_level']) or not 0 < config['prediction_level'] < 1:
        raise ValueError('prediction_level must be in (0, 1)')
    if not np.isfinite(config['mean_shift']):
        raise ValueError('mean_shift must be finite')
    return config


def study(config):
    theta, mu, sigma = (config['truth'][k] for k in ('theta', 'mu', 'sigma'))
    dt, replications = config['base_dt'], config['replications']
    horizons = np.asarray(config['forecast_horizons'], dtype=float)
    forecast_steps = np.rint(horizons/dt).astype(int)
    streams = np.random.SeedSequence(config['seed']).spawn(len(config['training_durations']))
    scenarios, records, streams_record, demo = [], [], [], None
    for duration, stream in zip(config['training_durations'], streams):
        train_steps = int(round(duration/dt))
        paths = simulate_ou_exact(theta, mu, sigma, mu, dt,
            train_steps+int(forecast_steps.max()), replications, np.random.default_rng(stream))
        streams_record.append(dict(duration=duration, entropy=int(stream.entropy), spawn_key=list(stream.spawn_key)))
        sampling_steps = sorted(set([dt]+config['sampling_steps'])) if duration == config['sampling_duration'] else [dt]
        for sampling_dt in sampling_steps:
            stride = int(round(sampling_dt/dt))
            training = paths[:, :train_steps+1:stride]
            future = paths[:, train_steps+forecast_steps]
            shifted_future = future + config['mean_shift']*(1-np.exp(-theta*horizons))
            name = f'T={duration:g}, dt={sampling_dt:g}'
            fits, reasons = [], {}
            hit_lists = {key: [[] for _ in horizons] for key in ('plugin', 'oracle', 'shift_plugin', 'shift_oracle')}
            square_errors = {key: [[] for _ in horizons] for key in ('plugin', 'persistence')}
            for index, series in enumerate(training):
                row = dict(scenario=name, replication=index, success=False, failure_reason='',
                           theta=None, mu=None, sigma=None, half_life=None)
                try:
                    fit = fit_ou(series, sampling_dt)
                except OUFitError as exc:
                    reason = str(exc).split(' slope ')[0]
                    reasons[reason] = reasons.get(reason, 0)+1
                    row['failure_reason'] = str(exc)
                    records.append(row)
                    continue
                fits.append(fit)
                row.update(success=True, theta=fit.theta, mu=fit.mu, sigma=fit.sigma, half_life=fit.half_life)
                records.append(row)
                plugin = forecast_ou(fit.theta, fit.mu, fit.sigma, series[-1], horizons, config['prediction_level'])
                oracle = forecast_ou(theta, mu, sigma, series[-1], horizons, config['prediction_level'])
                shift_oracle = forecast_ou(theta, mu+config['mean_shift'], sigma, series[-1], horizons, config['prediction_level'])
                for key, actual, prediction in [('plugin', future[index], plugin), ('oracle', future[index], oracle),
                        ('shift_plugin', shifted_future[index], plugin), ('shift_oracle', shifted_future[index], shift_oracle)]:
                    for j in range(horizons.size):
                        hit_lists[key][j].append(prediction[1][j] <= actual[j] <= prediction[2][j])
                for j in range(horizons.size):
                    square_errors['plugin'][j].append((plugin[0][j]-future[index,j])**2)
                    square_errors['persistence'][j].append((series[-1]-future[index,j])**2)
                if duration == config['sampling_duration'] and sampling_dt == dt and index == 0:
                    all_horizons = np.arange(int(forecast_steps.max())+1)*dt
                    demo_forecast = forecast_ou(fit.theta, fit.mu, fit.sigma, series[-1], all_horizons, config['prediction_level'])
                    demo = dict(fit=asdict(fit), times=(np.arange(paths.shape[1])*dt).tolist(),
                                values=paths[index].tolist(), train_steps=train_steps,
                                forecast_mean=demo_forecast[0].tolist(), forecast_lower=demo_forecast[1].tolist(),
                                forecast_upper=demo_forecast[2].tolist())
            truth = {**config['truth'], 'half_life': float(np.log(2)/theta)}
            parameters = {}
            for key in truth:
                values = np.array([f.half_life if key == 'half_life' else getattr(f, key) for f in fits])
                parameters[key] = (dict(mean=float(values.mean()), bias=float(values.mean()-truth[key]),
                    rmse=float(np.sqrt(np.mean((values-truth[key])**2))),
                    percentile_025=float(np.quantile(values,.025)), percentile_975=float(np.quantile(values,.975)))
                    if values.size else dict(mean=None, bias=None, rmse=None, percentile_025=None, percentile_975=None))
            forecasts = []
            for j, horizon in enumerate(horizons):
                forecasts.append(dict(horizon=float(horizon),
                    **{key: coverage_summary(hits[j]) for key, hits in hit_lists.items()},
                    **{key+'_rmse': float(np.sqrt(np.mean(errors[j]))) if errors[j] else None
                       for key, errors in square_errors.items()}))
            scenarios.append(dict(name=name, duration=duration, dt=sampling_dt,
                n_transitions=training.shape[1]-1, attempted=replications, successful=len(fits),
                failures=replications-len(fits), failure_reasons=reasons,
                parameters=parameters, forecasts=forecasts))
        del paths
    return dict(config=config, versions=dict(python=platform.python_version(), numpy=np.__version__,
        scipy=scipy.__version__, matplotlib=matplotlib.__version__), streams=streams_record,
        data_source='Synthetic exact OU transitions; deterministic initial value equals true mu',
        coverage_population='Successful interior OU fits only; oracle uses the same retained paths',
        shifted_future_design='Same future Gaussian innovations; mean shift begins at forecast origin',
        scenarios=scenarios), records, demo


def plots(summary, demo, folder):
    cfg = summary['config']
    base = [s for s in summary['scenarios'] if s['dt'] == cfg['base_dt']]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.1))
    for ax, key in zip(axes, ('theta', 'mu', 'sigma')):
        ax.plot([s['duration'] for s in base], [s['parameters'][key]['rmse'] for s in base], 'o-')
        ax.set(xscale='log', xlabel='Training duration T', ylabel='Empirical RMSE', title=f'Estimating {key}')
    fig.suptitle(f"{cfg['replications']} paths per design; dt={cfg['base_dt']} (successful fits)")
    fig.tight_layout()
    fig.savefig(folder/'parameter_recovery.png', dpi=140)
    plt.close(fig)
    sampling = [s for s in summary['scenarios'] if s['duration']==cfg['sampling_duration']]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.1))
    for ax, key in zip(axes, ('theta', 'sigma')):
        ax.plot([s['dt'] for s in sampling], [s['parameters'][key]['rmse'] for s in sampling], 'o-')
        ax.set(xlabel='Sampling interval dt', ylabel='Empirical RMSE', title=f'{key}: fixed duration T={cfg["sampling_duration"]}')
    fig.tight_layout()
    fig.savefig(folder/'sampling_frequency.png', dpi=140)
    plt.close(fig)
    reference = next(s for s in base if s['duration']==cfg['sampling_duration'])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for name, label, shift in [('plugin', 'Fitted OU: stable future', -.04),
                              ('oracle', 'True parameters: stable future', -.013),
                              ('shift_plugin', 'Fitted OU: shifted future', .013),
                              ('shift_oracle', 'True shifted parameters', .04)]:
        rows = [f[name] for f in reference['forecasts']]
        if not reference['successful']:
            continue
        y = np.array([r['coverage'] for r in rows])
        error = [y-np.array([r['lower'] for r in rows]), np.array([r['upper'] for r in rows])-y]
        axes[0].errorbar(np.array(cfg['forecast_horizons'])+shift, y, yerr=error, fmt='o-', capsize=3, label=label)
    axes[0].axhline(cfg['prediction_level'], color='black', ls='--', label='Nominal prediction level')
    axes[0].set(xlabel='Forecast horizon', ylabel='Coverage + Wilson 95% CI', ylim=(0, 1.02), title='Parameter uncertainty and mean shift')
    axes[0].legend(fontsize=7)
    for rows in base:
        if not rows['successful']:
            continue
        coverage_rows = [f['plugin'] for f in rows['forecasts']]
        y = np.array([r['coverage'] for r in coverage_rows])
        error = [y-np.array([r['lower'] for r in coverage_rows]), np.array([r['upper'] for r in coverage_rows])-y]
        axes[1].errorbar(cfg['forecast_horizons'], y, yerr=error, fmt='o-', capsize=3,
                         label=f'T={rows["duration"]}')
    axes[1].axhline(cfg['prediction_level'], color='black', ls='--')
    axes[1].set(xlabel='Forecast horizon', ylabel='Plug-in coverage', title='Training duration affects interval calibration')
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(folder/'forecast_coverage.png', dpi=140)
    plt.close(fig)
    if demo:
        fig, ax = plt.subplots(figsize=(10, 4.1))
        cut = demo['train_steps']
        times = np.array(demo['times'])
        ax.plot(times, demo['values'], lw=.9, label='Synthetic observations: replication 0')
        ax.plot(times[cut:], demo['forecast_mean'], label='Forecast mean fitted on training only')
        ax.fill_between(times[cut:], demo['forecast_lower'], demo['forecast_upper'], alpha=.25, label='Plug-in prediction interval')
        ax.axvline(times[cut], color='black', ls='--', label='Forecast origin')
        ax.set(xlabel='Time', ylabel='X', title='One illustrative held-out forecast; no seed selection')
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(folder/'example_forecast.png', dpi=140)
        plt.close(fig)


def write_report(summary, folder):
    cfg = summary['config']
    lines = ['# OU Calibration: Executed Results', '',
        'All observations are synthetic exact OU transitions. Results below are generated by `run.py`.', '',
        f"Seed: {cfg['seed']}; {cfg['replications']} independent paths per record-duration design.", '',
        '## Parameter recovery', '',
        'Errors below are conditional on successful interior fits. All attempts and failures remain in `replications.csv`.', '',
        '| T | dt | Successful / attempted | theta RMSE | mu RMSE | sigma RMSE |',
        '| --- | --- | --- | --- | --- | --- |']
    for s in summary['scenarios']:
        errors = [s['parameters'][key]['rmse'] for key in ('theta', 'mu', 'sigma')]
        cells = [f'{v:.5f}' if v is not None else 'undefined' for v in errors]
        lines.append(f"| {s['duration']} | {s['dt']} | {s['successful']} / {s['attempted']} | {' | '.join(cells)} |")
    lines += ['', '## Forecast coverage', '',
        'All methods use the same successful-fit paths. Oracle parameters are a simulation benchmark, not information available to a fitted model.', '',
        '| T | dt | Horizon | Fitted OU | True parameters | Fitted OU after mean shift | True shifted parameters |',
        '| --- | --- | --- | --- | --- | --- | --- |']
    for s in summary['scenarios']:
        for f in s['forecasts']:
            cells=[]
            for key in ('plugin', 'oracle', 'shift_plugin', 'shift_oracle'):
                row=f[key]
                cells.append(f"{row['coverage']:.3f} [{row['lower']:.3f}, {row['upper']:.3f}]" if row['n'] else 'undefined')
            lines.append(f"| {s['duration']} | {s['dt']} | {f['horizon']} | {' | '.join(cells)} |")
    lines += ['', 'Brackets are Wilson 95% intervals for the coverage experiment, not a prediction interval for one future value.', '',
        'The mean-shift future uses the same innovations as the stable future. Only the equilibrium mean changes at the forecast origin.', '',
        '## Interpretation and limits', '',
        '- Longer records provide additional time over which mean reversion can be observed. Dense sampling at fixed duration does not create a longer record.',
        '- Finite-sample conditional likelihood estimates can be biased. Parameter-distribution percentiles in summary.json are not parameter confidence intervals.',
        '- Plug-in prediction intervals exclude uncertainty in estimated parameters. Coverage must be measured rather than assumed.',
        '- A mean shift violates the fitted constant-mean model; the shifted oracle illustrates that lost coverage is model-dependent.',
        '- Replications estimate performance at the specified synthetic settings. These are not observations of SST or evidence that ocean anomalies follow an OU model.',
        '- Forecast horizons within one path are correlated. Coverage intervals are pointwise, with no multiplicity adjustment.', '',
        '## Figures', '',
        '![Parameter recovery](figures/parameter_recovery.png)', '',
        '![Sampling frequency](figures/sampling_frequency.png)', '',
        '![Forecast coverage](figures/forecast_coverage.png)', '',
        '']
    if (folder/'example_series.json').exists():
        lines += ['![Example forecast](figures/example_forecast.png)', '']
    lines += ['## Forecast RMSE versus persistence', '',
        'Persistence predicts the last training observation at every horizon. Both errors use the same successful-fit paths.', '',
        '| T | dt | Horizon | Fitted OU RMSE | Persistence RMSE |', '| --- | --- | --- | --- | --- |']
    for scenario in summary['scenarios']:
        for row in scenario['forecasts']:
            values = [row['plugin_rmse'], row['persistence_rmse']]
            cells = [f'{v:.5f}' if v is not None else 'undefined' for v in values]
            lines.append(f"| {scenario['duration']} | {scenario['dt']} | {row['horizon']} | {' | '.join(cells)} |")
    lines.append('')
    (folder/'report.md').write_text('\n'.join(lines), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('config.json'))
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('results'))
    args = parser.parse_args()
    config = checked_config(json.loads(args.config.read_text()))
    summary, records, demo = study(config)
    folder = args.output
    (folder/'figures').mkdir(parents=True, exist_ok=True)
    (folder/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    with (folder/'replications.csv').open('w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    if demo:
        (folder/'example_series.json').write_text(json.dumps(demo, indent=2, allow_nan=False)+'\n')
    else:
        # Do not retain an illustration from an earlier configuration.
        (folder/'example_series.json').unlink(missing_ok=True)
        (folder/'figures/example_forecast.png').unlink(missing_ok=True)
    plots(summary, demo, folder/'figures')
    write_report(summary, folder)
    for scenario in summary['scenarios']:
        print(f"{scenario['name']}: {scenario['successful']}/{scenario['attempted']} successful; "
              f"theta RMSE={scenario['parameters']['theta']['rmse']}")
    print(f'Report and results saved to {folder.resolve()}')


if __name__ == '__main__':
    main()
