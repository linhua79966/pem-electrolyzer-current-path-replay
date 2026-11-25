"""Development-only numerical timing grid. Inputs are caller-supplied arrays. No file access."""
from dataclasses import dataclass
import math
import numpy as np
TAUS=(0,.25,.5,1,2,3,5,8,10,15,20,30,45,60)
DELAYS=(0,1,2,4,6,10,15,20,30)
SPECS=tuple((float(t),float(d)) for t in TAUS for d in DELAYS)
PARTS=('FIT','SELECT','DEV_DIAGNOSTIC')
FAMILIES=('LAG','DELAY','FOP','FOPDT')

def family(tau,delay):
    return ('LAG' if delay==0 else 'DELAY') if tau==0 else ('FOP' if delay==0 else 'FOPDT')

def check_segment(t,u):
    t=np.asarray(t,dtype=float);u=np.asarray(u,dtype=float)
    if t.ndim!=1 or u.shape!=t.shape or len(t)==0 or not np.all(np.isfinite(t)) or not np.all(np.isfinite(u)) or np.any(u<0) or np.any(np.diff(t)<=0):
        raise ValueError('invalid segment: require increasing finite time and finite nonnegative current')
    return t-t[0],u

def base_state(t,u,tau):
    if not math.isfinite(tau) or tau<=0:raise ValueError('positive tau required')
    z=np.empty_like(u);z[0]=u[0]
    a=np.exp(-np.diff(t)/tau);om=-np.expm1(-np.diff(t)/tau)
    for k in range(1,len(t)): z[k]=a[k-1]*z[k-1]+om[k-1]*u[k-1]
    return z

def evaluate_basis(t,u,target_idx,tau,delay,*,base=None):
    """Exact ZOH at recorded knots, constant prehistory, left-limit static.
    t must already be normalized; no samples after target enter the feature.
    A positive tau uses the continuous state; tau=0 uses u((t-delay)^-).
    """
    if not math.isfinite(tau) or not math.isfinite(delay) or tau<0 or delay<0: raise ValueError('invalid parameters')
    ix=np.asarray(target_idx,dtype=int)
    if ix.ndim!=1 or np.any(ix<0) or np.any(ix>=len(t)):raise ValueError('target index')
    s=t[ix]-delay
    if tau==0:
        j=np.maximum(0,np.searchsorted(t,s,side='left')-1)
        return u[j]
    z=base_state(t,u,tau) if base is None else base
    j=np.maximum(0,np.searchsorted(t,s,side='right')-1)
    age=np.maximum(0,s-t[j]);a=np.exp(-age/tau)
    return a*z[j]-np.expm1(-age/tau)*u[j]

def features(t,u,target_idx,specs=SPECS):
    t,u=check_segment(t,u)
    cache={};out=[]
    for tau,d in specs:
        if tau>0 and tau not in cache:cache[tau]=base_state(t,u,tau)
        out.append(evaluate_basis(t,u,target_idx,tau,d,base=cache.get(tau)))
    return np.column_stack(out)

def consumed_change_10s(t,u,idx):
    t,u=check_segment(t,u);idx=np.asarray(idx,dtype=int)
    return evaluate_basis(t,u,idx,0,0)-evaluate_basis(t,u,idx,0,10)

def change_strata(d):
    d=np.asarray(d,dtype=float)
    return np.where(d>=30,'UP',np.where(d<=-30,'DOWN','LOW_CHANGE'))

def metrics(y,p):
    y=np.asarray(y,dtype=float);p=np.asarray(p,dtype=float)
    if y.shape!=p.shape or y.ndim!=1 or len(y)==0 or not np.all(np.isfinite(y)) or not np.all(np.isfinite(p)):raise ValueError('metric input')
    e=y-p;mse=float(np.mean(e*e))
    return dict(n=len(e),mse=mse,rmse=math.sqrt(mse),mae=float(np.mean(np.abs(e))),bias=float(np.mean(e)))
@dataclass(frozen=True)
class Group:
    file_id:int
    partition:str
    X:np.ndarray
    y:np.ndarray

def validate_groups(groups,ncols):
    seen=set();cnt={p:0 for p in PARTS}
    for g in groups:
        if g.file_id in seen or g.partition not in PARTS:raise ValueError('duplicate file or forbidden partition')
        seen.add(g.file_id);cnt[g.partition]+=1
        if g.X.shape!=(len(g.y),ncols) or g.y.ndim!=1 or len(g.y)<1 or not np.all(np.isfinite(g.y)) or not np.all(np.isfinite(g.X)) or np.any(g.X<0):raise ValueError('invalid design group')
    if any(n<2 for n in cnt.values()):raise ValueError('at least two files per partition')

def fit_select(groups,specs=SPECS):
    """2019 FIT gain only, 2020 SELECT only, 2021 target excluded from fitting/selection.
    No global winner chosen across families; all four family winners are reported.
    """
    validate_groups(groups,len(specs))
    fit=[g for g in groups if g.partition=='FIT'];sel=[g for g in groups if g.partition=='SELECT']
    denom=np.mean([np.mean(g.X*g.X,axis=0) for g in fit],axis=0)
    numer=np.mean([np.mean(g.X*g.y[:,None],axis=0) for g in fit],axis=0)
    beta=np.divide(numer,denom,out=np.zeros_like(numer),where=denom>0);beta=np.maximum(beta,0)
    fitmse=np.mean([np.mean((g.y[:,None]-g.X*beta)**2,axis=0) for g in fit],axis=0)
    select=np.mean([np.sqrt(np.mean((g.y[:,None]-g.X*beta)**2,axis=0)) for g in sel],axis=0)
    rows=[]
    for i,(tau,d) in enumerate(specs):
        ok=denom[i]>0 and np.isfinite(beta[i]) and np.isfinite(select[i])
        rows.append(dict(index=i,family=family(tau,d),tau_s=tau,delay_s=d,beta=float(beta[i]) if ok else None,fit_macro_mse=float(fitmse[i]) if ok else None,select_macro_rmse=float(select[i]) if ok else None,eligible=bool(ok)))
    winners={};near={}
    for f in FAMILIES:
        valid=[r for r in rows if r['family']==f and r['eligible']]
        if not valid: raise ValueError('no eligible candidate: '+f)
        low=min(r['select_macro_rmse'] for r in valid)
        tied=[r for r in valid if r['select_macro_rmse']<=low+1e-12]
        winners[f]=min(tied,key=lambda r:(r['delay_s'],r['tau_s']))
        near[f]=[r['index'] for r in valid if r['select_macro_rmse']<=low*1.01+1e-12]
    return rows,winners,near

def evaluate_winners(groups,winners):
    out=[]
    for g in groups:
        for f,r in winners.items():
            out.append(dict(file_id=g.file_id,partition=g.partition,model=f,**metrics(g.y,g.X[:,r['index']]*r['beta'])))
    return out

def summarize(rows):
    out=[]
    for p in PARTS:
        for f in FAMILIES:
            r=[x for x in rows if x['partition']==p and x['model']==f];n=sum(x['n'] for x in r)
            out.append(dict(partition=p,model=f,files=len(r),rows=n,macro_rmse=float(np.mean([x['rmse'] for x in r])),macro_mae=float(np.mean([x['mae'] for x in r])),macro_bias=float(np.mean([x['bias'] for x in r])),pooled_rmse=math.sqrt(sum(x['n']*x['mse'] for x in r)/n)))
    return out

def regular_grid_age(dt,tau,delay):
    """Steady ramp effective age for left ZOH; time constants in seconds."""
    if dt<=0 or tau<0 or delay<0:raise ValueError('age domain')
    m=math.floor(delay/dt);phi=delay-m*dt
    if tau==0:return (m+1)*dt
    return (m+1)*dt+dt*math.exp(-(dt-phi)/tau)/(-math.expm1(-dt/tau))
