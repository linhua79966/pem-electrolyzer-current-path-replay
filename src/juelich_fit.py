import numpy as np,math
TAUS=(1,2,3,5,8,10,15,20,30,45,60)

def design(segments):
    zs=[];ys=[];refs=[];span=[]
    for sn,s in enumerate(segments,1):
        t=np.array([x[2] for x in s]);u=np.array([x[4] for x in s]);y=np.array([x[5] if x[5] is not None else np.nan for x in s]);idx=np.flatnonzero((t-t[0]>=300)&np.isfinite(y))
        if not len(idx):continue
        X=np.empty((len(idx),len(TAUS)+1));X[:,0]=u[idx];delta=np.diff(t)
        for j,tau in enumerate(TAUS,1):
            a=np.exp(-delta/tau);state=np.empty(len(u));state[0]=u[0]
            for k in range(1,len(u)):state[k]=a[k-1]*state[k-1]+(1-a[k-1])*u[k-1]
            X[:,j]=state[idx]
        zs.append(X);ys.append(y[idx]);refs.extend((sn,s[k][0],float(t[k]),float(u[k])) for k in idx);span.extend((float(u[idx].min()),float(u[idx].max())))
    if not zs:return None
    X=np.concatenate(zs);y=np.concatenate(ys)
    if len(y)<30:return None
    return X,y,refs,min(span),max(span)

def gram_fit(groups):
    # Equal total experiment weights, not pooled sample count.
    return sum((X.T@X)/len(y) for X,y,*_ in groups)/len(groups),sum((X.T@y)/len(y) for X,y,*_ in groups)/len(groups)

def coeff(G,g,idx):
    A=G[np.ix_(idx,idx)];b=g[idx];diag=np.diag(A)
    if np.any(diag<=0):return None,0.,'ZERO_DESIGN'
    if len(idx)==1:return np.maximum(0,b/diag),1.,'ELIGIBLE'
    D=np.sqrt(diag);eig=np.linalg.eigvalsh(A/np.outer(D,D));ratio=max(0.,float(eig[0]/eig[-1]))
    if ratio<1e-4:return None,ratio,'COLLINEAR'
    # exact 2D NNLS: compare unconstrained interior, two boundary solutions and origin
    cand=[np.zeros(2),np.array([max(0,b[0]/diag[0]),0]),np.array([0,max(0,b[1]/diag[1])])];u=np.linalg.solve(A,b)
    if np.all(u>=0):cand.append(u)
    q=min(cand,key=lambda x:float(x@A@x-2*b@x))
    if q.sum()<=0 or q.min()/q.sum()<.01:return None,ratio,'DEGENERATE_GAIN'
    return q,ratio,'ELIGIBLE'

def metric(y,p):
    e=y-p;mse=float(np.mean(e*e));var=float(np.var(y));return dict(n=len(y),rmse=math.sqrt(mse),mae=float(np.mean(abs(e))),bias=float(np.mean(e)),r2=1-mse/var if var>0 else None,negative_predictions=int((p<0).sum()))

def candidates(groups,selection):
    G,g=gram_fit(groups);rows=[]
    specs=[('M0',(0,))]+[('M1',(i,)) for i in range(1,12)]+[('M2',(i,j)) for i in range(1,12) for j in range(i+1,12) if TAUS[j-1]/TAUS[i-1]>=2]
    for family,idx in specs:
        b,ratio,status=coeff(G,g,list(idx));r=dict(family=family,indices=list(idx),tau_s=[TAUS[i-1] for i in idx if i],beta=None if b is None else b.tolist(),gram_ratio=ratio,status=status,select_rmse=None)
        if b is not None:r['select_rmse']=float(np.mean([metric(y,X[:,idx]@b)['rmse'] for X,y,*_ in selection]))
        rows.append(r)
    winners={}
    for family in ('M0','M1','M2'):
        ok=[r for r in rows if r['family']==family and r['status']=='ELIGIBLE']
        if ok:winners[family]=min(ok,key=lambda r:(r['select_rmse'],r['tau_s']))
    if not {'M0','M1'}<=winners.keys():raise ValueError('no eligible M0/M1')
    return rows,winners
