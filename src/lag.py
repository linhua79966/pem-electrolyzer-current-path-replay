"""Archived common-support lag calculations, in-memory only. Positive lag: command leads measured current."""
import numpy as np
K=np.arange(-10,11)

def corr(x,y):
    if len(x)<100: return np.nan
    x=x-x.mean(); y=y-y.mean(); a=np.linalg.norm(x); b=np.linalg.norm(y)
    if a==0 or b==0: return np.nan
    return float(np.clip(np.dot(x,y)/(a*b),-1,1))

def pairs(segments, kind):
    out=[]
    for c,a in segments:
        if kind=='increment': c,a=np.diff(c),np.diff(a)
        if len(c)>20: out.append((c,a))
    x=np.concatenate([c[10:len(c)-10] for c,a in out]) if out else np.empty(0)
    ys=[np.concatenate([a[10+k:len(a)-10+k] for c,a in out]) if out else np.empty(0) for k in K]
    return x,ys

def peak(rs):
    if not np.all(np.isfinite(rs)): return None,[]
    ties=K[np.abs(rs-np.max(rs))<=1e-12].tolist()
    return min(ties,key=lambda k:(abs(k),k)),ties

def selftests():
    rng=np.random.default_rng(171); c=rng.normal(size=2048); a=np.r_[np.zeros(3),c[:-3]]
    tests=[]
    for label,seg,expected in [('positive_3',[(c,a)],3),('negative_3',[(a,c)],-3),('zero',[(c,c)],0)]:
        for kind in ['level','increment']:
            x,ys=pairs(seg,kind); r=np.array([corr(x,y) for y in ys]); p,_=peak(r)
            assert p==expected; tests.append([label+'_'+kind,'PASS'])
    x,ys=pairs([(c[:71],a[:71]),(c[71:170],a[71:170])],'level'); assert len(x)==130 and all(len(y)==130 for y in ys); tests.append(['no_cross_segment_common_support','PASS'])
    x,ys=pairs([(c[:100],a[:100])],'level'); assert all(np.isnan(corr(x,y)) for y in ys); tests.append(['minimum_n','PASS'])
    assert np.isnan(corr(np.ones(101),np.arange(101))); tests.append(['zero_variance','PASS'])
    r=np.zeros(21); r[8]=1;r[12]=1; assert peak(r)==(-2,[-2,2]); r[10]=1; assert peak(r)==(0,[-2,0,2]);tests.append(['tie_rule','PASS'])
    r=np.zeros(21);r[0]=1; assert peak(r)[0]==-10; tests.append(['edge_peak','PASS'])
    r=np.zeros(21);r[10]=np.nan; assert peak(r)==(None,[]);tests.append(['undefined_peak','PASS'])
    return tests
