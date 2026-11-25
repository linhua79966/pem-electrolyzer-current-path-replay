"""In-memory NLR methods; caller supplies data and block indices."""
import math
MIN_EIG_RATIO=1e-4
MIN_FRACTION=.01

def state(u,tau,blocks):
    """Exact ZOH update on [k-1,k), using left endpoint. Each block reinitializes."""
    a=math.exp(-1/tau);b=-math.expm1(-1/tau);z=[0.]*len(u)
    for _,l,r in blocks:
        z[l]=u[l]
        for k in range(l+1,r):z[k]=a*z[k-1]+b*u[k-1]
    return z

def dot(x,y,ix):return math.fsum(x[i]*y[i] for i in ix)

def fit1(x,y,ix):
    s=dot(x,x,ix)
    if s<=0:raise ValueError('zero input energy')
    return max(0.,dot(x,y,ix)/s)

def loss(xs,bs,y,ix):
    return math.fsum((y[i]-sum(b*x[i] for x,b in zip(xs,bs)))**2 for i in ix)

def fit2(x,z,y,ix):
    aa=dot(x,x,ix);bb=dot(z,z,ix);ab=dot(x,z,ix)
    if aa<=0 or bb<=0:return None,'ZERO_ENERGY',None
    rho=max(-1.,min(1.,ab/math.sqrt(aa*bb)));er=(1-abs(rho))/(1+abs(rho))
    if er<MIN_EIG_RATIO:return None,'COLLINEAR_DESIGN',er
    ay=dot(x,y,ix);by=dot(z,y,ix);det=aa*bb-ab*ab
    a=(ay*bb-by*ab)/det;b=(by*aa-ay*ab)/det
    options=[(0.,0.),(max(0.,ay/aa),0.),(0.,max(0.,by/bb))]
    if a>=0 and b>=0:options.append((a,b))
    best=min(options,key=lambda v:(loss((x,z),v,y,ix),v))
    if sum(best)<=0 or min(best)/sum(best)<MIN_FRACTION:return None,'BOUNDARY_OR_SMALL_WEIGHT',er
    return best,'ELIGIBLE',er

def metrics(y,p,ix):
    if not ix:raise ValueError('empty metric support')
    e=[y[i]-p[i] for i in ix];m=math.fsum(y[i] for i in ix)/len(ix)
    s=math.fsum(x*x for x in e);v=math.fsum((y[i]-m)**2 for i in ix)
    return {'n':len(ix),'mae_kg_h':math.fsum(abs(x) for x in e)/len(ix),'rmse_kg_h':math.sqrt(s/len(ix)),'bias_observed_minus_predicted_kg_h':math.fsum(e)/len(ix),'r2':1-s/v if v>0 else None,'negative_predictions':sum(p[i]<0 for i in ix)}

def pred(xs,bs):return [sum(b*x for b,x in zip(bs,v)) for v in zip(*xs)]
