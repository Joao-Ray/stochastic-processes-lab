"""Reproduce noisy-OU inference, fixed-budget replication and information limits."""
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
from scipy.optimize import differential_evolution
from scipy.stats import chi2
from src.ou_inference import OUFitError, fit_ou
from src.ou_uncertainty import profile_theta_interval
from src.ou_monte_carlo import stationary_null_paths
from src.ou_measurement_error import (fit_noisy_ou, profile_membership, profile_shape,
    information_matrix, information_limit, MIN_STEP, MAX_STEP)

METHODS = ('conditional', 'stationary_ignorant', 'ml', 'reml')


def proportion(values):
    n = len(values)
    if n == 0: return dict(n=0, value=None, mc_se=None, wilson_95=None)
    p, z = float(np.mean(values)), 1.959963984540054
    center = (p+z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return dict(n=n, value=p, mc_se=float(np.sqrt(p*(1-p)/n)), wilson_95=[float(center-half),float(center+half)])


def point_metrics(values, truth):
    x = np.array([v for v in values if v is not None], dtype=float)-truth
    if not x.size: return dict(available=0, bias=None, rmse=None, bias_mc_se=None, rmse_mc_se=None)
    mse = float(np.mean(x*x)); rmse = np.sqrt(mse)
    return dict(available=len(x), bias=float(x.mean()), rmse=float(rmse),
        bias_mc_se=float(x.std(ddof=1)/np.sqrt(len(x))) if len(x)>1 else None,
        rmse_mc_se=float((x*x).std(ddof=1)/(2*np.sqrt(len(x))*rmse)) if len(x)>1 and rmse>0 else None)


def validate(cfg):
    for key, minimum in [('seed',0), ('replications',2), ('budget',8), ('audit_records_per_scenario',0)]:
        if type(cfg[key]) is not int or cfg[key] < minimum: raise ValueError(f'invalid {key}')
    if cfg['audit_records_per_scenario'] > cfg['replications']: raise ValueError('too many audit records')
    if not np.isfinite([cfg['theta'],cfg['mu'],cfg['stationary_variance'],cfg['level']]).all() or min(cfg['theta'],cfg['stationary_variance'])<=0 or not 0<cfg['level']<1: raise ValueError('invalid truth/level')
    if not cfg['scenarios']: raise ValueError('nonempty scenarios required')
    seen=set()
    for s in cfg['scenarios']:
        key=(s['duration'],s['repeats'],s['noise_sd'],s['stress'])
        if key in seen: raise ValueError('duplicate scenario')
        seen.add(key)
        if type(s['repeats']) is not int or s['repeats']<1 or cfg['budget']%s['repeats'] or cfg['budget']//s['repeats']<4: raise ValueError('budget must support >=4 sites and exact repeat count')
        if not np.isfinite([s['duration'],s['noise_sd']]).all() or s['duration']<=0 or s['noise_sd']<0 or s['stress'] not in ('none','colored','fixed5'): raise ValueError('invalid scenario')
        if s['repeats']>1 and s['noise_sd']==0: raise ValueError('singular zero-noise replicated design')


def generate(cfg, scenario, index):
    n, k = cfg['replications'], scenario['repeats']
    m = cfg['budget']//k
    dt = scenario['duration']/(m-1)
    rng = np.random.default_rng(np.random.SeedSequence([cfg['seed'],300,index,0]))
    latent = cfg['mu']+np.sqrt(cfg['stationary_variance'])*stationary_null_paths(float(np.exp(-cfg['theta']*dt)),m-1,n,rng)
    if scenario['stress']=='fixed5':
        latent += (cfg['mu']+5*np.sqrt(cfg['stationary_variance'])-latent[:,:1])*np.exp(-cfg['theta']*dt*np.arange(m))
    noise_rng=np.random.default_rng(np.random.SeedSequence([cfg['seed'],300,index,1]))
    if scenario['stress']=='colored':
        noise = stationary_null_paths(.6,m-1,n*k,noise_rng).reshape(n,k,m).transpose(0,2,1)
    else: noise = noise_rng.standard_normal((n,m,k))
    return latent[:,:,None]+scenario['noise_sd']*noise,dt


def audit_fit(y, fit, seed):
    restricted=fit['restricted']
    upper=1. if y.shape[1]==1 else 1-1e-9
    result=differential_evolution(lambda z:profile_shape(y,float(np.exp(-z[0])),float(z[1]),restricted)['nll'],
        bounds=[(MIN_STEP,MAX_STEP),(0.,upper)],seed=seed,popsize=12,maxiter=160,tol=1e-9,polish=True)
    best=min(float(result.fun),profile_shape(y,fit['phi'],fit['gamma'],restricted)['nll'])
    return dict(seed=seed, objective_local=fit['nll'], objective_alternative=float(result.fun),
        gap=float(fit['nll']-best), discrepancy=bool(fit['nll']-best>1e-4),
        alternative_converged=bool(result.success), alternative_phi=float(np.exp(-result.x[0])),
        alternative_gamma=float(result.x[1]))


def study(cfg):
    validate(cfg)
    rows,summaries,audits=[],[],[]
    for design,s in enumerate(cfg['scenarios']):
        data,dt=generate(cfg,s,design)
        name=f"T={s['duration']:g},k={s['repeats']},eta={s['noise_sd']:g},stress={s['stress']}"
        block=[]
        for rep,y in enumerate(data):
            means=y.mean(axis=1)
            row=dict(scenario=name,replication=rep)
            try: naive=fit_ou(means,dt).theta
            except OUFitError: naive=None
            row.update(conditional_theta=naive,conditional_test_available=True,
                conditional_hit=profile_theta_interval(means,dt,cfg['level']).contains(cfg['theta']))
            try:
                ignored=fit_noisy_ou(means,dt,fixed_gamma=1.)
                lr=2*(profile_shape(means,float(np.exp(-cfg['theta']*dt)),1.)['nll']-ignored['nll'])
                ok=bool(ignored['optimizer_converged'] and lr>=-1e-6)
                row.update(stationary_ignorant_theta=ignored['theta'],stationary_ignorant_test_available=ok,
                    stationary_ignorant_hit=bool(ok and lr<=chi2.ppf(cfg['level'],1)))
            except (ValueError,FloatingPointError,np.linalg.LinAlgError):
                row.update(stationary_ignorant_theta=None,stationary_ignorant_test_available=False,stationary_ignorant_hit=False)
            for method,restricted in [('ml',False),('reml',True)]:
                try:
                    fit=fit_noisy_ou(y,dt,restricted)
                    member=profile_membership(y,dt,cfg['theta'],fit,cfg['level'])
                    row.update({method+'_theta':fit['theta'],method+'_test_available':member['available'],
                        method+'_phi':fit['phi'],method+'_nll':fit['nll'],method+'_mu':fit['mu'],method+'_scale':fit['scale'],
                        method+'_hit':member['hit'],method+'_lr':member['lr'],
                        method+'_noise_variance':fit['noise_variance'],method+'_latent_variance':fit['latent_variance'],
                        method+'_gamma':fit['gamma'],method+'_white_limit':fit['white_limit'],
                        method+'_search_cap':fit['search_cap'],method+'_noise_boundary':fit['noise_boundary'],
                        method+'_optimizer_converged':fit['optimizer_converged'],method+'_optimizer_gap':member['optimizer_gap'],
                        method+'_exception':None})
                    if rep<cfg['audit_records_per_scenario']:
                        seed=cfg['seed']+10000+design*20+rep*2+int(restricted)
                        audits.append(dict(scenario=name,replication=rep,method=method,**audit_fit(y,fit,seed)))
                except (ValueError,FloatingPointError,np.linalg.LinAlgError) as e:
                    row.update({method+'_theta':None,method+'_test_available':False,method+'_hit':False,
                        method+'_phi':None,method+'_nll':None,method+'_mu':None,method+'_scale':None,
                        method+'_lr':None,method+'_noise_variance':None,method+'_latent_variance':None,
                        method+'_gamma':None,method+'_white_limit':False,method+'_search_cap':False,
                        method+'_noise_boundary':False,method+'_optimizer_converged':False,
                        method+'_optimizer_gap':False,method+'_exception':type(e).__name__})
            block.append(row)
        metrics={}
        for method in METHODS:
            available=[r for r in block if r[method+'_test_available']]
            metrics[method]=dict(point=point_metrics([r[method+'_theta'] for r in block],cfg['theta']),
                test_available=len(available),conditional_coverage=proportion([r[method+'_hit'] for r in available]),
                success_and_coverage=proportion([r[method+'_hit'] for r in block]))
            if method in ('ml','reml'):
                metrics[method]['diagnostics']={flag:sum(bool(r[method+'_'+flag]) for r in block) for flag in ['white_limit','search_cap','noise_boundary','optimizer_gap','exception']}
                metrics[method]['diagnostics']['nonconvergence']=sum(not r[method+'_optimizer_converged'] for r in block)
                noise=[r[method+'_noise_variance'] for r in block if r[method+'_noise_variance'] is not None]
                metrics[method]['noise_variance_n']=len(noise)
                metrics[method]['noise_variance_mean']=float(np.mean(noise)) if noise else None
        paired={}
        for a,b in [('ml','conditional'),('reml','conditional'),('reml','ml')]:
            common=[r for r in block if r[a+'_theta'] is not None and r[b+'_theta'] is not None]
            diff=np.array([(r[a+'_theta']-cfg['theta'])**2-(r[b+'_theta']-cfg['theta'])**2 for r in common])
            paired[a+'_minus_'+b]=dict(n=len(common),mse_difference=float(diff.mean()) if len(common) else None,
                mc_se=float(diff.std(ddof=1)/np.sqrt(len(common))) if len(common)>1 else None)
        info=None
        if s['noise_sd']>0 and s['stress']=='none':
            info=information_matrix(cfg['theta'],cfg['stationary_variance'],s['noise_sd']**2,dt,y.shape[0],s['repeats'])
        summaries.append(dict(id=name,**s,dt=dt,sites=y.shape[0],attempted=len(block),
            model_assumptions_met=s['stress']=='none',methods=metrics,paired=paired,information=info))
        rows.extend(block)
        print(f'Completed {name}: {len(block)} records',flush=True)
    limits=[]
    for gamma in [1.,.8,.5,.2,.1,.05,.01,0.]:
        limits.append(dict(gamma=gamma,theta0=1.,theta1=2.,dt=.1,sites=101,
            **information_limit(float(np.exp(-.1)),float(np.exp(-.2)),gamma,101)))
    return dict(config=cfg,scenarios=summaries,information_limits=limits,
        data_seed_scheme='SeedSequence([seed,300,scenario_index,role]); role 0 latent,1 sensor',
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        optimizer_audit_discrepancies=sum(a['discrepancy'] for a in audits)),rows,audits


def write_report(result,folder):
    cfg=result['config']
    def fmt(x):return 'NA' if x is None else f'{x:.4f}'
    lines=['# Executed noisy-OU study','',
        '[Protocol](../../PROTOCOL.md) · [Derivations](../../THEORY.md) · [Summary](summary.json) · [Every record](replications.csv.gz) · [Optimizer audit](optimizer_audit.json)','',
        f"Seed {cfg['seed']}; {cfg['replications']} attempts/design; total scalar-reading budget {cfg['budget']}; nominal level {cfg['level']}.",'',
        'ML and REML estimate unknown sensor variance using the full joint Gaussian law. Profile coverage uses an asymptotic cutoff, not a finite-sample guarantee. Different designs have independent data; methods within each record are paired. All RMSE values below are conditional on an available point estimate.','',
        '| Design | Conditional OLS RMSE (n) | Stationary noise-ignored RMSE (n) | Noisy ML RMSE (n) | Noisy REML RMSE (n) |',
        '| --- | --- | --- | --- | --- |']
    for s in result['scenarios']:
        cells=[f"{fmt(s['methods'][m]['point']['rmse'])} ({s['methods'][m]['point']['available']})" for m in METHODS]
        lines.append('| '+s['id']+' | '+' | '.join(cells)+' |')
    lines+=['','Conditional bias, RMSE MC SE and paired squared-error differences are in summary.json. Methods with different availability populations cannot be ranked from conditional RMSE alone.','',
        '## Profile membership and availability','',
        '| Design | OLS profile, all | Noise-ignored profile, all | Noisy ML conditional / all (available n) | Noisy REML conditional / all (available n) |',
        '| --- | --- | --- | --- | --- |']
    for s in result['scenarios']:
        ms=s['methods']; vals=[fmt(ms[m]['success_and_coverage']['value']) for m in METHODS[:2]]
        for m in METHODS[2:]:
            p=ms[m]; vals.append(f"{fmt(p['conditional_coverage']['value'])} / {fmt(p['success_and_coverage']['value'])} ({p['test_available']})")
        lines.append('| '+s['id']+' | '+' | '.join(vals)+' |')
    lines+=['','All-attempt values count unavailable membership as failure. Conditional denominators, Wilson intervals and MC SEs are saved. White-noise limits can lack a finite point estimate while still admitting profile membership. Full-set endpoints and widths are not measured.','',
        '## Diagnostics','', '| Design | ML white / cap / nonconvergence / LR gap / exceptions | REML white / cap / nonconvergence / LR gap / exceptions |', '| --- | --- | --- |']
    for s in result['scenarios']:
        cells=[]
        for m in METHODS[2:]:
            d=s['methods'][m]['diagnostics']; cells.append(' / '.join(str(d[k]) for k in ['white_limit','search_cap','nonconvergence','optimizer_gap','exception']))
        lines.append('| '+s['id']+' | '+' | '.join(cells)+' |')
    lines += ['',f"Independent optimizer audit discrepancies >1e-4: **{result['optimizer_audit_discrepancies']}**. Audits are numerical comparisons, not global optimization certificates.",'',
        '## Equal-budget local information','', '| Design | Efficient theta information | Local unbiased SE bound |','| --- | --- | --- |']
    for s in result['scenarios']:
        if s['information']:
            lines.append(f"| {s['id']} | {fmt(s['information']['efficient_theta'])} | {fmt(s['information']['local_se_bound'])} |")
    lines += ['','This local bound applies to regular unbiased estimators at interior true parameters. It is not a universal bound on biased-estimator RMSE, and does not assert that replication wins for all designs.','',
        '## Weak-signal two-point illustration','', '| Signal fraction | KL(1 to 2) | TV upper bound | Size-0.05 power upper bound | Both-drift confidence membership lower bound |', '| --- | --- | --- | --- | --- |']
    for p in result['information_limits']:
        lines.append(f"| {fmt(p['gamma'])} | {fmt(p['kl_0_to_1'])} | {fmt(p['tv_upper'])} | {fmt(p['power_upper'])} | {fmt(p['simultaneous_membership_lower'])} |")
    lines += ['','Two alternatives theta=1 and 2, dt=0.1, 101 sites, same unit total variance and signal fraction. These population bounds illustrate a submodel and do not validate fitted intervals. The gamma=0 row is a white-noise limit.','',
        '![Coverage and conditional error](coverage_and_error.png)','',
        '## Scope','',
        '- Sensor noise is independent Gaussian and initial states stationary, except explicitly labeled diagnostics.',
        '- Correcting the likelihood does not by itself ensure short-window, boundary or weak-signal calibration.',
        '- Equal budget means scalar reading count; timing costs and real sensor constraints are excluded.',
        '- Established likelihood/REML/information tools are credited. Originality and journal suitability have not been established.','']
    (folder/'report.md').write_text('\n'.join(lines))


def save(result,rows,audits,folder):
    folder.mkdir(parents=True,exist_ok=True)
    archive=folder/'replications.csv.gz'
    with archive.open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as zipped:
            with io.TextIOWrapper(zipped,encoding='utf-8',newline='') as text:
                writer=csv.DictWriter(text,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    result['record_archive_sha256']=hashlib.sha256(archive.read_bytes()).hexdigest()
    (folder/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (folder/'optimizer_audit.json').write_text(json.dumps(audits,indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(2,1,figsize=(13,8),sharex=True)
    scenarios=result['scenarios']
    for method,offset in zip(METHODS,[-.24,-.08,.08,.24]):
        x=np.arange(len(scenarios))+offset
        axes[0].plot(x,[s['methods'][method]['point']['rmse'] for s in scenarios],'o',label=method)
        axes[1].errorbar(x,[s['methods'][method]['success_and_coverage']['value'] for s in scenarios],
            yerr=[1.96*s['methods'][method]['success_and_coverage']['mc_se'] for s in scenarios],fmt='o',capsize=2,label=method)
    axes[0].set(yscale='log',ylabel='Conditional drift RMSE',title='Unknown measurement error: availability-aware comparison')
    axes[0].legend(ncol=4,fontsize=8)
    axes[1].axhline(result['config']['level'],color='black',ls='--',lw=.8)
    axes[1].set(ylabel='All-attempt profile membership',ylim=(-.02,1.03))
    axes[1].set_xticks(np.arange(len(scenarios)),[f"T={s['duration']:g},k={s['repeats']}\neta={s['noise_sd']:g}\n{s['stress']}" for s in scenarios],rotation=45,ha='right',fontsize=8)
    fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(folder/f'coverage_and_error.{suffix}',dpi=160)
    plt.close(fig)
    write_report(result,folder)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('config.json'))
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('results')/'main')
    args=parser.parse_args();started=time.perf_counter()
    result,rows,audits=study(json.loads(args.config.read_text()))
    sources=[Path(__file__),ROOT/'src/ou_measurement_error.py',ROOT/'src/ou_inference.py',ROOT/'src/ou_uncertainty.py',ROOT/'src/ou_monte_carlo.py',Path(__file__).with_name('THEORY.md'),Path(__file__).with_name('PROTOCOL.md'),Path(__file__).with_name('AMENDMENTS.md')]
    result['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    result['config_sha256']=hashlib.sha256(args.config.read_bytes()).hexdigest()
    result['elapsed_study_seconds']=time.perf_counter()-started
    save(result,rows,audits,args.output)
    print(f'Saved {len(rows)} attempts to {args.output.resolve()}',flush=True)

if __name__=='__main__':main()
