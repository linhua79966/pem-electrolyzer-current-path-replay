import math
PRED=("H","C","A","U","M0_H","M0_ORACLE_ENDPOINT")
CROSS=((0,1),(0,2),(1,2)); ENAMES=("R","T","F","total"); EPS=1e-9
def need(ok,msg):
    if not ok: raise ValueError(msg)

def valid_input(r):return r[2] is not None and r[3] is not None and min(r[2],r[3])>=0

def split_segments(g):
    segs=[];seg=[];last=None
    for r in g:
        if not valid_input(r):
            if seg:segs.append(seg);seg=[]
            last=None;continue
        if last is not None and (r[1]-last).total_seconds()!=1:
            segs.append(seg);seg=[]
        seg.append(r);last=r[1]
    if seg:segs.append(seg)
    return segs

class Scores:
    def __init__(self):self.n=0;self.mean=0.;self.m2=0.;self.ss=[0.]*6;self.abs=[0.]*6;self.bias=[0.]*6;self.es=[0.]*4;self.em=[0.]*4;self.cross=[0.]*3
    def add(self,q,p,e):
        self.n+=1;delta=q-self.mean;self.mean+=delta/self.n;self.m2+=delta*(q-self.mean)
        for j,x in enumerate(p):
            d=q-x;self.ss[j]+=d*d;self.abs[j]+=abs(d);self.bias[j]+=d
        for j,x in enumerate(e):self.es[j]+=x*x;self.em[j]+=x
        for j,(k,l) in enumerate(CROSS):self.cross[j]+=2*e[k]*e[l]
    def report(self):
        r={'n':self.n}
        if not self.n:return r
        n=self.n
        for j,k in enumerate(PRED):
            r.update({'rmse_'+k:math.sqrt(self.ss[j]/n),'mae_'+k:self.abs[j]/n,'bias_'+k:self.bias[j]/n,'r2_'+k:(1-self.ss[j]/self.m2 if self.m2>0 else None)})
        for j,k in enumerate(ENAMES):r.update({'mse_e_'+k:self.es[j]/n,'mean_e_'+k:self.em[j]/n})
        for j,(k,l) in enumerate(CROSS):r['cross_'+ENAMES[k]+'_'+ENAMES[l]]=self.cross[j]/n
        r['mse_identity_gap']=(self.es[3]-sum(self.es[:3])-sum(self.cross))/n
        need(abs(r['mse_identity_gap'])<=EPS,'aggregate MSE identity mismatch');return r
