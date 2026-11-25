import math,itertools
N=18000; BURN=300; TAUS=(1.,2.,3.,5.,8.,10.,15.,20.,30.,45.,60.)
BLOCKS=(("FIT",0,9000),("SELECT",9000,13500),("DEV_DIAGNOSTIC",13500,18000))
MIN_EIG_RATIO=1e-4; MIN_FRACTION=.01

def support(u,c,q):
    rows=[];indices={}
    for name,l,r in BLOCKS:
        ix=[i for i in range(l+BURN,r) if q[i] is not None]
        indices[name]=ix
        row={'block':name,'first_row_1based':l+1,'last_row_1based':r,'burn_in_rows':BURN,'eligible':len(ix),'potential':r-l-BURN,'nonfinite_targets':sum(q[i] is None for i in range(l,r)),'negative_targets':sum(q[i] is not None and q[i]<0 for i in range(l,r)),'actual_current_min_A':min(u[l:r]),'actual_current_max_A':max(u[l:r]),'command_min_A':min(c[l:r]),'command_max_A':max(c[l:r]),'command_change_edges':sum(abs(c[i]-c[i-1])>1e-9 for i in range(l+1,r)),'actual_change_edges':sum(abs(u[i]-u[i-1])>1e-9 for i in range(l+1,r))}
        rows.append(row)
        if len(ix)<.9*(r-l-BURN):raise ValueError('insufficient common target support: '+name)
        if name in ('FIT','SELECT') and max(u[i] for i in ix)-min(u[i] for i in ix)<30:raise ValueError('insufficient actual-current range (30 A) in '+name)
    return rows,indices

def state(u,tau,blocks=BLOCKS):
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

def identify(u,q,ix):
    # Fit only on FIT. Candidate selection only on SELECT. No refit before DEV_DIAGNOSTIC.
    tr=ix['FIT'];va=ix['SELECT'];zs={t:state(u,t) for t in TAUS};cand=[];models={};preds={}
    beta=fit1(u,q,tr);models['M0']={'tau_s':[],'beta_kg_h_per_A':[beta]};preds['M0']=pred([u],[beta]);
    viable=[]
    for t in TAUS:
        b=fit1(zs[t],q,tr);v=math.sqrt(loss([zs[t]],[b],q,va)/len(va));row={'model':'M1','tau_s':[t],'beta_kg_h_per_A':[b],'selection_rmse_kg_h':v,'status':'ELIGIBLE'};cand.append(row);viable.append((v,t,b))
    _,t,b=min(viable);models['M1']={'tau_s':[t],'beta_kg_h_per_A':[b]};preds['M1']=pred([zs[t]],[b]);viable=[]
    for t1,t2 in itertools.combinations(TAUS,2):
        if t2/t1<2:continue
        b,status,er=fit2(zs[t1],zs[t2],q,tr);row={'model':'M2','tau_s':[t1,t2],'train_normalized_gram_eigen_ratio':er,'status':status}
        if b is not None:
            v=math.sqrt(loss([zs[t1],zs[t2]],b,q,va)/len(va));row.update(beta_kg_h_per_A=list(b),selection_rmse_kg_h=v);viable.append((v,t1,t2,b))
        cand.append(row)
    if viable:
        _,t1,t2,b=min(viable);models['M2']={'tau_s':[t1,t2],'beta_kg_h_per_A':list(b)};preds['M2']=pred([zs[t1],zs[t2]],b)
    else:models['M2']={'status':'NO_ELIGIBLE_TWO_TIMESCALE_CANDIDATE'}
    # Model choices are now immutable. Calculate last-block diagnostics once.
    results={m:{k:metrics(q,p,v) for k,v in ix.items()} for m,p in preds.items()}
    return models,cand,preds,results
