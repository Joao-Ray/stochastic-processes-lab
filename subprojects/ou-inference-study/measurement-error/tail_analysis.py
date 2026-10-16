"""Audit tails of reported finite drift estimates without trimming observations."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import numpy as np


def analyze(folder):
    summary=json.loads((folder/'summary.json').read_text())
    archive=folder/'replications.csv.gz'
    digest=hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest!=summary['record_archive_sha256']: raise ValueError('archive hash mismatch')
    with gzip.open(archive,'rt') as stream: records=list(csv.DictReader(stream))
    result=dict(archive_sha256=digest,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scenarios=[],definition='Finite available point estimates only; no records trimmed or replaced. Top-five share uses squared errors within that same availability population.')
    lines=['# Tail sensitivity of finite drift errors','',
        '[Complete report](report.md) · [Every record](replications.csv.gz) · [Full tail audit](tails.json)','',
        'No record is removed or winsorized. Each denominator is the finite point-estimate population of that method. Missing boundary/optimizer estimates remain missing and are counted in the main report. This does not estimate unconditional risk of a procedure that always returns a finite drift.','',
        '| Design | Method | Available n | Median abs error | 99th percentile | Max abs error | Top-five share of squared error |',
        '| --- | --- | --- | --- | --- | --- | --- |']
    for scenario in summary['scenarios']:
        block=[r for r in records if r['scenario']==scenario['id']]
        methods={}
        for method in ('conditional','stationary_ignorant','ml','reml'):
            selected=[r for r in block if r[method+'_theta']!='']
            errors=np.array([float(r[method+'_theta'])-summary['config']['theta'] for r in selected])
            if not len(errors):
                methods[method]=dict(available=0);continue
            squared=errors**2; largest=np.argsort(squared)[-5:][::-1]
            fraction=float(squared[largest].sum()/squared.sum()) if squared.sum()>0 else 0.
            values=dict(available=len(errors),median_absolute_error=float(np.median(np.abs(errors))),
                q99_absolute_error=float(np.quantile(np.abs(errors),.99)),max_absolute_error=float(np.max(np.abs(errors))),
                top_five_squared_error_share=fraction,
                largest_records=[dict(replication=int(selected[j]['replication']),theta=float(selected[j][method+'_theta']),absolute_error=float(abs(errors[j]))) for j in largest])
            methods[method]=values
            lines.append(f"| {scenario['id']} | {method} | {len(errors)} | {values['median_absolute_error']:.4f} | {values['q99_absolute_error']:.4f} | {values['max_absolute_error']:.4f} | {fraction:.4f} |")
        result['scenarios'].append(dict(id=scenario['id'],methods=methods))
    (folder/'tails.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (folder/'tail_report.md').write_text('\n'.join(lines)+'\n')
    print(f'Audited tails in {folder.resolve()}')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=Path(__file__).with_name('results')/'main')
    analyze(parser.parse_args().results)
