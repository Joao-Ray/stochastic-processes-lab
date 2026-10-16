"""Post hoc audit of the worst finite REML error in each confirmation design.

Diagnostic only: this never changes point fits, calibration, seeds or records.
The displayed likelihood grid is not a certified confidence-set inversion.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,importlib.util,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2
from src.ou_measurement_error import profile_shape,profile_phi,MIN_STEP,MAX_STEP


def diagnose(folder,config_path):
    script=Path(__file__).with_name('run.py')
    spec=importlib.util.spec_from_file_location('measurement_study_runner',script)
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    cfg=json.loads(config_path.read_text())
    summary=json.loads((folder/'summary.json').read_text())
    if cfg!=summary['config']:raise ValueError('configuration does not match records')
    archive=folder/'replications.csv.gz';digest=hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest!=summary['record_archive_sha256']:raise ValueError('archive checksum mismatch')
    with gzip.open(archive,'rt') as file:rows=list(csv.DictReader(file))
    results=[]
    fig,axes=plt.subplots(len(cfg['scenarios']),1,figsize=(10,3.5*len(cfg['scenarios'])),squeeze=False)
    for design,(s,ax) in enumerate(zip(cfg['scenarios'],axes[:,0])):
        name=summary['scenarios'][design]['id']
        selected=[r for r in rows if r['scenario']==name and r['reml_theta']!='']
        if not selected:raise ValueError('no finite REML estimates to diagnose')
        worst=max(selected,key=lambda r:abs(float(r['reml_theta'])-cfg['theta']))
        rep=int(worst['replication'])
        data,dt=runner.generate(cfg,s,design);y=data[rep]
        fit={k:float(worst['reml_'+k]) for k in ['nll','phi','gamma']};fit['restricted']=True
        recalculated=profile_shape(y,fit['phi'],fit['gamma'],True)['nll']
        if not np.isclose(recalculated,fit['nll'],rtol=0,atol=1e-8):raise ValueError('record reconstruction failed')
        audit=runner.audit_fit(y,fit,cfg['seed']+950000+design)
        theta_grid=np.unique(np.r_[np.geomspace(MIN_STEP/dt,MAX_STEP/dt,151),float(worst['reml_theta']),cfg['theta']])
        lr=np.array([2*(profile_phi(y,float(np.exp(-theta*dt)),True)['nll']-fit['nll']) for theta in theta_grid])
        results.append(dict(scenario=name,replication=rep,theta_estimate=float(worst['reml_theta']),
            raw_observations=y.tolist(),dt=dt,truth_theta=cfg['theta'],audit=audit,
            theta_grid=theta_grid.tolist(),profile_lr=lr.tolist(),negative_profile_gap=bool(np.min(lr)<-1e-6)))
        ax.plot(theta_grid,lr,label='Evaluated REML profile grid')
        ax.axhline(chi2.ppf(cfg['level'],1),ls='--',color='grey',label='Asymptotic cutoff')
        ax.axvline(cfg['theta'],ls=':',color='black',label='True drift (annotation)')
        ax.axvline(float(worst['reml_theta']),ls='--',color='tab:red',label='Recorded estimate')
        ax.set(xscale='log',yscale='symlog',xlabel='Drift theta',ylabel='Profile LR',
            title=f'{name}; worst finite REML error, replication {rep}')
        ax.legend(fontsize=8)
    fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(folder/f'extreme_profiles.{suffix}',dpi=160)
    plt.close(fig)
    result=dict(posthoc=True,selection='Largest absolute finite REML error per supplied scenario; no trimming or refitting of scientific results',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),script,ROOT/'src/ou_measurement_error.py']},
        config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),archive_sha256=digest,examples=results)
    (folder/'extreme_profiles.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f'Audited {len(results)} post hoc extremes; DE discrepancies: {sum(e["audit"]["discrepancy"] for e in results)}')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=Path(__file__).with_name('results')/'confirmation')
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('confirmation.json'))
    args=parser.parse_args();diagnose(args.results,args.config)
