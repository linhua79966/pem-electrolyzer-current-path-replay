import math
from collections import defaultdict
import numpy as np
from timing_core import SPECS,PARTS,FAMILIES,fit_select,metrics,family,regular_grid_age
REFS=("M0_NOW_REFERENCE","M0_LAG_FIXED","M1_FIXED","M2_FIXED")
STRATA=("UP","DOWN","LOW_CHANGE")

def sufficient(g):
    return np.mean(g.X*g.X,axis=0),np.mean(g.X*g.y[:,None],axis=0),float(np.mean(g.y*g.y))

def evaluate(prepared,legacy_models):
    groups=[x['group'] for x in prepared]
    cand,winners,near=fit_select(groups)
    beta0=legacy_models['M0']['beta'][0]
    moments=[];losses=[];fm=[];sm=[];diff=[]
    for it in prepared:
        g=it['group'];X=g.X;y=g.y;meta={'file_id':g.file_id,'partition':g.partition}
        if g.partition=='FIT':
            xx,xy,yy=sufficient(g)
            moments.extend(dict(**meta,index=i,mean_x2=float(xx[i]),mean_xy=float(xy[i]),mean_y2=yy,n=len(y)) for i in range(len(SPECS)))
        if g.partition=='SELECT':
            for c in cand:
                if c['eligible']:
                    i=c['index'];losses.append(dict(**meta,index=i,rmse=metrics(y,X[:,i]*c['beta'])['rmse'],n=len(y)))
        pred={f:X[:,winners[f]['index']]*winners[f]['beta'] for f in FAMILIES}
        pred.update(dict(zip((REFS[0],REFS[2],REFS[3]),it['legacy_predictions'].T)))
        pred[REFS[1]]=X[:,SPECS.index((0.,0.))]*beta0
        base=metrics(y,pred['M1_FIXED']);base_lag=metrics(y,pred['M0_LAG_FIXED'])
        for m,p in pred.items():
            r=metrics(y,p);fm.append(dict(**meta,model=m,**r))
            diff.append(dict(**meta,model=m,delta_rmse_vs_M1_fixed=base['rmse']-r['rmse'],delta_rmse_vs_M0_lag_fixed=base_lag['rmse']-r['rmse']))
            for s in STRATA:
                mask=it['strata']==s
                if mask.any():sm.append(dict(**meta,model=m,stratum=s,**metrics(y[mask],p[mask])))
    def aggregate(rs,keys):
        out=[];by=defaultdict(list)
        for r in rs:by[tuple(r[k] for k in keys)].append(r)
        for k,rr in sorted(by.items()):
            n=sum(x['n'] for x in rr);v=dict(zip(keys,k));out.append(dict(**v,files=len(rr),rows=n,macro_rmse=float(np.mean([x['rmse'] for x in rr])),macro_mse=float(np.mean([x['mse'] for x in rr])),macro_mae=float(np.mean([x['mae'] for x in rr])),macro_bias=float(np.mean([x['bias'] for x in rr])),pooled_rmse=math.sqrt(math.fsum(x['n']*x['mse'] for x in rr)/n),pooled_mae=math.fsum(x['n']*x['mae'] for x in rr)/n,pooled_bias=math.fsum(x['n']*x['bias'] for x in rr)/n))
        return out
    for c in cand:
        c['near_family_min']=c['index'] in near[c['family']]
        c['selected']=c['index']==winners[c['family']]['index']
        c['tau_edge']=c['tau_s'] in (.25,60) if c['tau_s']>0 else False
        c['delay_edge']=c['delay_s'] in (1,30) if c['delay_s']>0 else False
        c['regular_2s_effective_age_s']=regular_grid_age(2,c['tau_s'],c['delay_s'])
    return {'CANDIDATES':cand,'WINNERS':winners,'NEAR_OPTIMAL':near,'FIT_MOMENTS':moments,'SELECT_FILE_LOSS':losses,'FILE_METRICS':fm,'STRATUM_FILE_METRICS':sm,'MACRO_METRICS':aggregate(fm,('partition','model')),'STRATUM_MACRO_METRICS':aggregate(sm,('partition','model','stratum')),'PAIRED_FILE_DIFFERENCES':diff}
