"""Audit supplied Juelich aggregates, not raw data or new fitted models."""
import argparse, json, math, hashlib, csv
from pathlib import Path
import numpy as np
import pandas as pd

def audit(dataset):
    root=Path(dataset); data=root/'data'; checks=[]
    def ck(name,cond):
        checks.append({'check':name,'pass':bool(cond)})
        if not cond:raise ValueError(name)
    def near(name,a,b):ck(name,math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-10))
    def read(p):return pd.read_csv(data/p)
    for x in csv.DictReader((root/'DATA_PROVENANCE.csv').open(encoding='utf-8-sig')):
        p=root/x['file'];ck('hash '+x['file'],p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==x['sha256'])
    j=read('04_ANALYSIS/JUELICH_2023_H30_FILE_EFFECTS.csv');s=read('04_ANALYSIS/JUELICH_2023_ALL_FILE_SUPPORT.csv')
    ck('2023 original support',len(j)==15 and j.n.sum()==6217 and s.time_rows.sum()==149889 and len(s)==20)
    for x in j.itertuples():
        near(f'{x.file_id} signed gain',x.rmse_H**2-x.rmse_A**2,x.mse_F+x.interaction)
        near(f'{x.file_id} residual',x.rmse_A**2,x.mse_R)
        near(f'{x.file_id} total',x.rmse_H**2,x.mse_total)
    gain=j.rmse_H**2-j.rmse_A**2;top=gain.nlargest(2);h=j.rmse_H.mean();a=j.rmse_A.mean();remaining=j.drop(top.index)
    summary=dict(files=len(j),origins=int(j.n.sum()),row_coverage_pct=float(100*j.n.sum()/s.time_rows.sum()),macro_H=float(h),macro_A=float(a),macro_reduction_pct=float(100*(1-a/h)),top_two_gain_share_pct=float(100*top.sum()/gain.sum()),remaining_median_relative_pct=float(np.median(100*(1-remaining.rmse_A/remaining.rmse_H))),improved=int((j.rmse_H>j.rmse_A).sum()),worsened=int((j.rmse_H<j.rmse_A).sum()))
    ck('direction counts',summary['improved']==11 and summary['worsened']==4)
    c=read('04_REPORTS/CANDIDATES.csv');f=read('04_REPORTS/FIT_MOMENTS.csv');s=read('04_REPORTS/SELECT_FILE_LOSS.csv');ck('archived candidate cardinality',len(c)==126 and len(f)==2772 and len(s)==7308)
    for x in c.itertuples(index=False):
        d=f[f['index']==x.index];v=s[s['index']==x.index];m=d[['mean_x2','mean_xy','mean_y2']].mean();beta=max(0,m.mean_xy/m.mean_x2)
        ck(f'{x.index} original fit/select support',len(d)==22 and len(v)==58)
        near(f'{x.index} beta',beta,x.beta);near(f'{x.index} fit MSE',beta**2*m.mean_x2-2*beta*m.mean_xy+m.mean_y2,x.fit_macro_mse);near(f'{x.index} select RMSE',v.rmse.mean(),x.select_macro_rmse)
    for family,g in c.groupby('family'):
        chosen=g[g.selected];ck(f'{family} archived winner',len(chosen)==1);near(f'{family} winner loss',chosen.select_macro_rmse.iloc[0],g.select_macro_rmse.min())
    for names,keys in [('FILE_METRICS MACRO_METRICS',['partition','model']),('STRATUM_FILE_METRICS STRATUM_MACRO_METRICS',['partition','model','stratum'])]:
        a,b=[read('04_REPORTS/'+n+'.csv') for n in names.split()]
        for _,x in b.iterrows():
            g=a
            for k in keys:g=g[g[k]==x[k]]
            tag=' / '.join(str(x[k]) for k in keys)
            ck(tag+' files and rows',len(g)==x.files and g.n.sum()==x.rows)
            for k in ['rmse','mse','mae','bias']:near(tag+' macro '+k,g[k].mean(),x['macro_'+k])
            near(tag+' pooled RMSE',math.sqrt(np.average(g.mse,weights=g.n)),x.pooled_rmse)
            for k in ['mae','bias']:near(tag+' pooled '+k,np.average(g[k],weights=g.n),x['pooled_'+k])
    return dict(status='PASS_JUELICH_AGGREGATE_ONLY',checks=len(checks),passed=sum(i['pass'] for i in checks),key_results=summary,scope='Archived summary identities and arithmetic only. No acquisition, signal access, new fit/selection/prediction, or reserved replay.',details=checks)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('dataset',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    if a.output.exists() or a.output.resolve().is_relative_to(a.dataset.resolve()):raise ValueError('use a new output directory outside the dataset')
    res=audit(a.dataset);a.output.mkdir(parents=True);(a.output/'AUDIT_RESULT.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='details'},indent=2))
