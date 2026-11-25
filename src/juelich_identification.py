"""In-memory equal-file Gram and nonnegative coefficient routines."""
import numpy as np
import math

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
