"""Reproduce finite-sample stationary rank/cell calibration and stress tests."""
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

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy

from src.ou_monte_carlo import cell_at, invert_cells, phi_cells, stationary_rank_test, test_cell
from src.ou_process import simulate_ou_exact
from src.ou_uncertainty import ar_statistics, profile_theta_interval


def proportion(hits):
    values = np.asarray(hits, dtype=float)
    n = len(values)
    if n == 0:
        return dict(n=0, value=None, mc_se=None, wilson_95=None)
    p, z = float(values.mean()), 1.959963984540054
    center = (p+z*z/(2*n))/(1+z*z/n)
    radius = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return dict(n=n, value=p, mc_se=float(np.sqrt(p*(1-p)/n)),
                wilson_95=[float(center-radius), float(center+radius)])


def paired(a, b):
    delta = np.asarray(a, dtype=float)-np.asarray(b, dtype=float)
    return dict(n=len(delta), difference=float(delta.mean()),
                mc_se=float(delta.std(ddof=1)/np.sqrt(len(delta))))


def validate(config, base):
    for key, minimum in [('data_seed', 0), ('mc_seed', 0), ('replications', 2), ('simulations', 19)]:
        if type(config[key]) is not int or config[key] < minimum:
            raise ValueError(f'{key} must be an integer >= {minimum}')
    if not np.isfinite(config['alpha']) or not 0 < config['alpha'] < 1:
        raise ValueError('alpha must be in (0,1)')
    if type(config['examples']) is not bool or not config['scenarios']:
        raise ValueError('examples must be boolean and scenarios nonempty')
    seen = set()
    for s in config['scenarios']:
        identity = (s['r'], s['a'], s['noise_sd'], s['start'])
        if identity in seen or s['r'] not in base['record_lengths'] or s['a'] not in base['sampling_steps']:
            raise ValueError('unique scenarios on the parent design required')
        seen.add(identity)
        if not np.isfinite(s['noise_sd']) or s['noise_sd'] < 0 or s['start'] not in ('stationary', 'fixed5'):
            raise ValueError('nonnegative finite noise and known initialization required')
        phi_cells(int(round(s['r']/s['a'])), config['resolution'], config['u_max'])


def study(config):
    base_path = ROOT/'subprojects/ou-inference-study/config.json'
    base = json.loads(base_path.read_text())
    validate(config, base)
    theta, mu, v = (base[k] for k in ('theta', 'mu', 'stationary_variance'))
    count, data_seed, bank_seed = (config[k] for k in ('replications', 'data_seed', 'mc_seed'))
    fine, sigma = min(base['sampling_steps']), np.sqrt(2*theta*v)
    prior = {}
    baseline_digest = None
    if count == base['replications'] and data_seed == base['seed']:
        archive = ROOT/'subprojects/ou-inference-study/results/main/replications.csv.gz'
        baseline_digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        with gzip.open(archive, 'rt') as file:
            prior = {(r['scenario'], int(r['replication'])): r for r in csv.DictReader(file)}
    latent_cache = {}
    summaries, records, examples = [], [], []
    for design_index, design in enumerate(config['scenarios']):
        r, a, eta, start = (design[k] for k in ('r', 'a', 'noise_sd', 'start'))
        group = base['record_lengths'].index(r)
        if r not in latent_cache:
            rng = np.random.default_rng(np.random.SeedSequence([data_seed, 0, group]))
            initial = mu+np.sqrt(v)*rng.standard_normal(count)
            steps = int(round(r/fine))
            latent = simulate_ou_exact(theta, mu, sigma, mu, fine/theta, steps, count, rng)
            latent += (initial-mu)[:, None]*np.exp(-fine*np.arange(steps+1))
            latent_cache[r] = latent
        training = latent_cache[r][:, ::int(round(a/fine))].copy()
        if start == 'fixed5':
            training += (mu+5*np.sqrt(v)-training[:, :1])*np.exp(-a*np.arange(training.shape[1]))
        if eta:
            noise_rng = np.random.default_rng(np.random.SeedSequence([data_seed, 1, group]))
            training += eta*noise_rng.standard_normal(training.shape)
        n, dt = training.shape[1]-1, a/theta
        true_phi = float(np.exp(-a))
        grid = phi_cells(n, config['resolution'], config['u_max'])
        cell_index, cell = cell_at(true_phi, grid)
        name = f'r={r:g},a={a:g},eta={eta:g},start={start}'
        rows = []
        for index, series in enumerate(training):
            stats = ar_statistics(series)
            bank = lambda role: np.random.default_rng(np.random.SeedSequence([bank_seed, 20, design_index, index, role]))
            rank = stationary_rank_test(series, true_phi, bank(0), config['simulations'])
            buffered = test_cell(series, cell, bank(1), config['simulations'], config['alpha'])
            profile = profile_theta_interval(series, dt, 1-config['alpha'])
            row = dict(scenario=name, replication=index, phi_hat=stats.phi,
                rank_pvalue=rank['pvalue'], rank_hit=rank['pvalue'] > config['alpha'],
                rank_exceedances=rank['exceedances'], cell_pvalue=buffered['pvalue'],
                cell_hit=buffered['accepted'], profile_hit=profile.contains(theta),
                old_bootstrap_available=None, old_bootstrap_hit=None,
                power_half_reject=None, power_double_reject=None)
            key = (f'r={r:g},a={a:g},eta={eta:g}', index)
            if start == 'stationary' and key in prior:
                old = prior[key]
                if not np.isclose(float(old['phi']), stats.phi, rtol=0, atol=1e-10):
                    raise RuntimeError('paired baseline does not reproduce the same observed series')
                row.update(old_bootstrap_available=old['bootstrap_profile_available'] == 'True',
                           old_bootstrap_hit=old['bootstrap_profile_hit'] == 'True')
            if r == 10 and a == .1 and eta == 0 and start == 'stationary':
                for label, multiplier, role in [('half', .5, 2), ('double', 2., 3)]:
                    result = stationary_rank_test(series, float(np.exp(-a*multiplier)), bank(role), config['simulations'])
                    row['power_'+label+'_reject'] = result['pvalue'] <= config['alpha']
            rows.append(row)
        old_available = [row for row in rows if row['old_bootstrap_available'] is True]
        summary = dict(id=name, **design, n_transitions=n, dt=dt,
            theorem_assumptions_met=start == 'stationary' and eta == 0,
            rank_coverage=proportion([row['rank_hit'] for row in rows]),
            cell_coverage=proportion([row['cell_hit'] for row in rows]),
            profile_coverage=proportion([row['profile_hit'] for row in rows]),
            old_bootstrap_conditional=proportion([row['old_bootstrap_hit'] for row in old_available]),
            paired_rank_minus_profile=paired([row['rank_hit'] for row in rows], [row['profile_hit'] for row in rows]),
            paired_cell_minus_profile=paired([row['cell_hit'] for row in rows], [row['profile_hit'] for row in rows]),
            true_cell=dict(index=cell_index, lower=cell.lower, upper=cell.upper, center=cell.center,
                tv_bound=cell.tv_bound, effective_alpha=max(0., config['alpha']-cell.tv_bound)),
            power_half=proportion([row['power_half_reject'] for row in rows if row['power_half_reject'] is not None]),
            power_double=proportion([row['power_double_reject'] for row in rows if row['power_double_reject'] is not None]))
        summaries.append(summary)
        records.extend(rows)
        if config['examples'] and r in (2, 10) and eta == 0 and start == 'stationary':
            example = invert_cells(training[0], dt, bank_seed+1000+design_index, config['simulations'],
                config['alpha'], config['resolution'], config['u_max'])
            examples.append(dict(scenario=name, replication=0, values=training[0].tolist(),
                truth_theta=theta, dt=dt, **example))
        print(f'Completed {name}: {count} records', flush=True)
    return dict(config=config, truth=dict(theta=theta, mu=mu, sigma=float(sigma), stationary_variance=v),
        scenarios=summaries, parent_data_config=base,
        matched_baseline_archive_sha256=baseline_digest,
        bank_seed_scheme='SeedSequence([mc_seed,20,scenario_index,replication,role]); roles 0 point,1 cell,2 half,3 double',
        data_seed_scheme='Parent SeedSequence([data_seed,0,duration_group]); noise role 1',
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      matplotlib=matplotlib.__version__)), records, examples


def write_report(result, examples, folder):
    cfg = result['config']
    def fmt(value): return 'unavailable' if value is None else f'{value:.4f}'
    lines = ['# Executed finite-sample stationary calibration', '',
        '[Design](../../README.md) · [Proofs](../../../FINITE_SAMPLE_THEORY.md) · [Summary](summary.json) · [Every outer record](replications.csv.gz)', '',
        f"Data seed {cfg['data_seed']}; bank seed {cfg['mc_seed']}; {cfg['replications']} records per design; B={cfg['simulations']}.", '',
        'Coverage here is point-null nonrejection or membership in the true parameter cell. Full continuous inversions are computed only for the saved examples. It is not an empirical width study of full inversions.', '',
        '| Scenario | Assumptions met? | Profile | Stationary point-rank [Wilson 95%] | MC SE | TV-buffered cell [Wilson 95%] | MC SE | Prior bootstrap-profile (available n) |',
        '| --- | --- | --- | --- | --- | --- | --- | --- |']
    for s in result['scenarios']:
        rank, cell, old = s['rank_coverage'], s['cell_coverage'], s['old_bootstrap_conditional']
        def cell_text(m): return f"{fmt(m['value'])} [{fmt(m['wilson_95'][0])},{fmt(m['wilson_95'][1])}]"
        lines.append(f"| {s['id']} | {s['theorem_assumptions_met']} | {fmt(s['profile_coverage']['value'])} | {cell_text(rank)} | {fmt(rank['mc_se'])} | {cell_text(cell)} | {fmt(cell['mc_se'])} | {fmt(old['value'])} (n={old['n']}) |")
    lines += ['', 'The prior bootstrap values use audited matched records when available. They are conditional on an interior generating fit; the new methods retain every nondegenerate record. A changed data seed or reduced run has no matched archived bootstrap baseline.', '',
        '## Paired improvements over chi-square-profile', '',
        '| Scenario | Point-rank minus profile (MC SE) | Buffered-cell minus profile (MC SE) | TV margin | Effective alpha |',
        '| --- | --- | --- | --- | --- |']
    for s in result['scenarios']:
        p, c, t = s['paired_rank_minus_profile'], s['paired_cell_minus_profile'], s['true_cell']
        lines.append(f"| {s['id']} | {fmt(p['difference'])} ({fmt(p['mc_se'])}) | {fmt(c['difference'])} ({fmt(c['mc_se'])}) | {fmt(t['tv_bound'])} | {fmt(t['effective_alpha'])} |")
    lines += ['', '## Two illustrative power points in the central clean design', '']
    for s in result['scenarios']:
        if s['power_half']['n']:
            lines.append(f"- {s['id']}: reject false theta=0.5×truth **{fmt(s['power_half']['value'])}**, MC SE **{fmt(s['power_half']['mc_se'])}**; theta=2×truth **{fmt(s['power_double']['value'])}**, MC SE **{fmt(s['power_double']['mc_se'])}**.")
    lines += ['', 'Power is measured for point-null tests, not the more conservative full cell envelope.', '',
        '## Predetermined full inversions', '',
        '[Every grid decision and raw example series](examples.json). Replication zero is fixed before execution; examples are not selected for narrow sets.', '']
    for example in examples:
        h = example['hull']
        upper = 'infinity' if h['unbounded_upper'] else fmt(h['upper'])
        lines.append(f"- {example['scenario']}: **{example['cells']}** cells; **{len(example['components'])}** connected components; conservative hull **[{fmt(h['lower'])}, {upper}]**. The phi-to-one tail is retained without testing.")
    if not examples:
        lines += ['No full inversion examples were requested for this run; examples.json is an empty list.']
    lines += ['', '## Figures', '', '![Coverage and assumption stress](coverage.png)', '']
    if examples:
        lines += ['![Full accepted-cell sets](accepted_sets.png)', '',
            'The display is restricted to theta in [0.001,100]; the actual saved sets retain the full zero-limit tail and any infinite upper endpoint.', '']
    lines += ['## Interpretation and limitations', '',
        '- The proof gives a coverage lower bound under stationary noise-free Gaussian OU, over data and simulation randomness. It does not guarantee conditional coverage at a fixed X0.',
        '- Stress cases violate assumptions. Whether a particular stress run has high or low coverage does not extend the theorem.',
        '- Buffered sets trade power and precision for a bound on cell discretization. They can be disconnected; their enclosing hull always reaches zero because the final tail is included.',
        '- Independent banks per outer record support the reported MC SEs. Main data reuse is explicit; confirmation is separate and is not pooled.',
        '- Existing Monte Carlo rank-test theory and Pinsker bounds are credited. Novelty of the OU specialization is not established.',
        '- Proofs assume exact draws/arithmetic; this floating-point implementation is checked numerically but is not a universal rounding certificate.', '']
    (folder/'report.md').write_text('\n'.join(lines))


def figures(result, examples, folder):
    rows = result['scenarios']
    fig, ax = plt.subplots(figsize=(11, 4.7))
    for method, offset, label in [('profile_coverage', -.18, 'Chi-square profile'),
                                ('rank_coverage', 0., 'Stationary point-rank'),
                                ('cell_coverage', .18, 'TV-buffered true-cell')]:
        x = np.arange(len(rows))+offset
        ax.errorbar(x, [s[method]['value'] for s in rows],
            yerr=[1.96*s[method]['mc_se'] for s in rows], fmt='o', capsize=3, label=label)
    ax.axhline(1-result['config']['alpha'], color='black', ls='--', lw=.8)
    labels = [f"r={s['r']:g}, a={s['a']:g}\n"+('clean' if s['theorem_assumptions_met'] else ('noise' if s['noise_sd'] else 'fixed start')) for s in rows]
    ax.set(xticks=np.arange(len(rows)), xticklabels=labels, ylabel='Coverage + 1.96 MC SE',
           ylim=(0, 1.02), title='Finite-sample calibration: coverage and assumption stress')
    ax.legend(fontsize=8)
    fig.tight_layout()
    for suffix in ('png', 'pdf'): fig.savefig(folder/f'coverage.{suffix}', dpi=160)
    plt.close(fig)
    if not examples: return
    fig, axes = plt.subplots(len(examples), 1, figsize=(10, 2.5*len(examples)), squeeze=False)
    for ax, example in zip(axes[:, 0], examples):
        for component in example['components']:
            low = max(component['lower'], .001)
            high = min(component['upper'] if component['upper'] is not None else 100., 100.)
            if low <= high: ax.axvspan(low, high, color='tab:blue', alpha=.5)
        ax.axvline(example['truth_theta'], color='black', ls='--', label='True theta (illustration only)')
        ax.set(xscale='log', xlim=(.001, 100.), ylim=(0, 1), yticks=[], xlabel='Theta, displayed range [0.001,100]', title=example['scenario'])
        ax.legend(fontsize=8)
    fig.tight_layout()
    for suffix in ('png', 'pdf'): fig.savefig(folder/f'accepted_sets.{suffix}', dpi=160)
    plt.close(fig)


def save(result, records, examples, folder):
    folder.mkdir(parents=True, exist_ok=True)
    archive = folder/'replications.csv.gz'
    with archive.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding='utf-8', newline='') as text:
                writer = csv.DictWriter(text, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)
    result['record_archive_sha256'] = hashlib.sha256(archive.read_bytes()).hexdigest()
    (folder/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    (folder/'examples.json').write_text(json.dumps(examples, indent=2, allow_nan=False)+'\n')
    figures(result, examples, folder)
    write_report(result, examples, folder)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('config.json'))
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('results')/'main')
    args = parser.parse_args()
    started = time.perf_counter()
    result, records, examples = study(json.loads(args.config.read_text()))
    source_files = [Path(__file__), ROOT/'src/ou_monte_carlo.py', ROOT/'src/ou_uncertainty.py', ROOT/'src/ou_inference.py',
        ROOT/'src/ou_process.py', Path(__file__).parent/'README.md',
        Path(__file__).parent.parent/'FINITE_SAMPLE_THEORY.md', Path(__file__).parent.parent/'config.json']
    result['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    result['config_sha256'] = hashlib.sha256(args.config.read_bytes()).hexdigest()
    result['elapsed_study_seconds'] = time.perf_counter()-started
    save(result, records, examples, args.output)
    print(f'Saved {len(records)} records to {args.output.resolve()}', flush=True)


if __name__ == '__main__': main()
