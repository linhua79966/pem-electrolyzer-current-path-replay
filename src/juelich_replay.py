import numpy as np,math

def simulate(s,models):
    # s: row,counter,relative seconds,absolute seconds,total current,target
    t=np.array([r[2] for r in s],float);u=np.array([r[4] for r in s],float);q=np.array([np.nan if r[5] is None else r[5] for r in s]);rn=np.array([r[0] for r in s],int)
    if not len(t) or not np.isfinite(t).all() or not np.isfinite(u).all() or np.any(u<0) or np.any(np.diff(t)<=0):raise ValueError('invalid segment')
    z={}
    for tau in sorted({tau for f in models.values() for tau in f['tau_s']}):
        a=np.exp(-np.diff(t)/tau);x=np.empty(len(t));x[0]=u[0]
        for i in range(1,len(t)):x[i]=a[i-1]*x[i-1]+(1-a[i-1])*u[i-1]
        z[tau]=x
    pp={k:(v['beta'][0]*u if not v['tau_s'] else sum(b*z[tau] for b,tau in zip(v['beta'],v['tau_s']))) for k,v in models.items()}
    ii=np.flatnonzero((t-t[0]>=300)&np.isfinite(q));ix=np.flatnonzero(t-t[0]>=300);ok=np.ones(len(ix),bool);js=[]
    for h in (10,30):
        v=t[ix]+h;r=np.minimum(np.searchsorted(t,v),len(t)-1);l=np.maximum(np.searchsorted(t,v)-1,0);j=np.where(abs(t[l]-v)<=abs(t[r]-v),l,r);ok&=(j>ix)&(abs(t[j]-v)<=.5)&np.isfinite(q[j]);js.append(j)
    pairs=[];beta=models['M1']['beta'][0];tau=models['M1']['tau_s'][0]
    for h,j in zip((10,30),js):
        i=ix[ok];k=j[ok];elapsed=t[k]-t[i];hp=beta*(u[i]+(z[tau][i]-u[i])*np.exp(-elapsed/tau));ap=pp['M1'][k];pairs.append((h,i,k,elapsed,q[k],hp,ap))
    return rn,t,u,q,ii,pp,pairs

def metric(y,p):
    e=y-p;ms=float(np.mean(e*e));v=float(np.var(y));return dict(n=len(y),rmse=math.sqrt(ms),mae=float(np.mean(abs(e))),bias=float(np.mean(e)),r2=None if v==0 else 1-ms/v,negative_predictions=int((p<0).sum()))
