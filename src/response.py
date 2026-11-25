"""Explicit, parameterized left-endpoint ZOH and matched path diagnostics.
Release adapter of the response mathematics in Phase2D/E and Phase3A; no input loader.
"""
import numpy as np

def zoh(t, u, tau):
    t,u=np.asarray(t,float),np.asarray(u,float)
    if t.ndim!=1 or u.shape!=t.shape or len(t)==0 or not np.isfinite(t).all() or not np.isfinite(u).all() or np.any(np.diff(t)<=0) or np.any(u<0) or not np.isfinite(tau) or tau<=0:
        raise ValueError('invalid single continuous segment or time constant')
    z=np.empty(len(u));z[0]=u[0]
    for k in range(1,len(u)):
        a=np.exp(-(t[k]-t[k-1])/tau)
        z[k]=a*z[k-1]+(1-a)*u[k-1]
    return z

def matched_paths(t,actual,command,origin,horizon,tau,beta):
    """Integer indices within one caller-validated segment; h is steps, not seconds.
    H holds command after consuming actual at origin; C consumes command thereafter;
    A consumes actual; U holds actual. No acquisition or experiment selection.
    """
    t,a,c=map(lambda x:np.asarray(x,float),(t,actual,command))
    if a.shape!=c.shape or not np.isfinite(c).all() or np.any(c<0) or not 0<=origin<origin+horizon<len(a) or not np.isfinite(beta) or beta<0:
        raise ValueError('invalid paths or indices')
    initial=beta*zoh(t,a,tau)[origin]; states=dict.fromkeys(('H','C','A','U'),initial)
    for j in range(origin,origin+horizon):
        alpha=np.exp(-(t[j+1]-t[j])/tau)
        inputs={'H':a[origin] if j==origin else c[origin],'C':a[origin] if j==origin else c[j],'A':a[j],'U':a[origin]}
        for k in states:states[k]=alpha*states[k]+beta*(1-alpha)*inputs[k]
    return states

def signed_terms(y,h,c,a):
    """Equal-weight matched-origin squared-error identity, native squared units."""
    y,h,c,a=map(lambda x:np.asarray(x,float),(y,h,c,a))
    if not y.size or any(x.shape!=y.shape or not np.isfinite(x).all() for x in (y,h,c,a)):
        raise ValueError('finite common support required')
    r,t,f=y-a,a-c,c-h
    m=lambda x:float(np.mean(x))
    return dict(mse_R=m(r*r),mse_T=m(t*t),mse_F=m(f*f),cross_R_T=2*m(r*t),cross_R_F=2*m(r*f),cross_T_F=2*m(t*f),mse_total=m((y-h)**2))
