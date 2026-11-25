"""Local, read-only acquisition and reproducibility adapter. No automatic download.
Use separately obtained legitimate source copies; this is not a production controller.
Results are local outputs, not a redistribution license. No historical execution gates.
"""
import argparse,collections,csv,datetime as dt,hashlib,io,json,math,re,sys,zipfile
from pathlib import Path
import numpy as np
import nlr_fit as nf
import nlr_core as nc
from nlr_clock import parse_clock
import lag
import juelich_reader as jr
import juelich_fit as jf
import juelich_replay as jp
import timing_core as tc
import timing_reports as tr
ROOT=Path(__file__).resolve().parents[1]
HORIZONS=(1,5,10,30); PARTS=('FIT','SELECT','DEV_DIAGNOSTIC')
SUBSETS=('ALL','CHANGE','NO_CHANGE','SMALL_CHANGE_LT30A','LARGE_CHANGE_GE30A')
ACT='H2E_f_PSU_A_Current'; CMD='H2E_f_NLR_CurrentCmd (A)'; FLOW='IVAL_f_FM011_Flow'

def require(ok,msg):
    if not ok: raise ValueError(msg)
def sha_file(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def put_json(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def put_csv(p,rows):
    rows=list(rows)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def num(s):
    try:v=float(s);return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None
def check_path(base,name,size,digest):
    require(isinstance(name,str) and name not in ('','.','..') and '/' not in name and '\\' not in name and ':' not in name,'invalid filename')
    root=Path(base).resolve();p=root/name
    require(p.is_file() and p.resolve().is_relative_to(root),'missing or external source: '+name)
    require(p.stat().st_size==size and sha_file(p)==digest,'source identity differs: '+name)
    return p

def preflight(spec,source,stage):
    if stage=='nlr-fit': items=[spec['characterization']]
    elif stage in ('nlr-replay','nlr-lag'):items=spec['nlr_combined']+spec['nlr_archives']
    elif stage in ('juelich-fit','juelich-timing'):items=spec['juelich_development']
    elif stage=='juelich-replay':items=spec['juelich_evaluation']
    else:raise ValueError('unknown stage')
    names=[x['filename'] for x in items];require(len(set(names))==len(names),'duplicate input filename')
    return {x['filename']:check_path(source,x['filename'],x.get('bytes',x.get('archive_bytes')),x.get('sha256',x.get('archive_sha256'))) for x in items}

def csv_records(data,header,cols):
    r=csv.reader(io.StringIO(data.decode('utf-8-sig'),newline=''),strict=True)
    h=next(r);require(h==header,'source header differs');require(all(h.count(k)==1 for k in cols),'mapped field missing/ambiguous')
    ix=[h.index(k) for k in cols];out=[]
    for row in r:
        require(len(row)==len(h),'ragged CSV record');out.append([row[i] for i in ix])
    require(bool(out),'empty source CSV');return out

def zip_member(p,member,size,digest,metadata=None):
    with zipfile.ZipFile(p) as z:
        names=z.namelist();require(len(names)==len(set(names)),'duplicate archive member')
        def read(n,s,h):
            info=z.getinfo(n);require(info.file_size==s,'member size differs');b=z.read(n);require(hashlib.sha256(b).hexdigest()==h,'member identity differs');return b
        b=read(member,size,digest)
        if metadata:
            n,h=metadata;require(hashlib.sha256(z.read(n)).hexdigest()==h,'metadata identity differs')
        return b

def fit_nlr(spec,paths):
    s=spec['characterization'];p=paths[s['filename']];mem={x['member']:x for x in s['members']};a=mem['characterization_200.csv']
    with zipfile.ZipFile(p) as z:
        require(len(z.namelist())==len(mem) and set(z.namelist())==set(mem),'characterization members differ')
        for x in mem.values():
            require(z.getinfo(x['member']).file_size==x['bytes'] and hashlib.sha256(z.read(x['member'])).hexdigest()==x['sha256'],'characterization member identity')
        data=z.read(a['member'])
    rows=csv_records(data,s['header'],['',ACT,CMD,FLOW]);t=[];u=[];c=[];q=[]
    for ts,aa,cc,yy in rows:
        x=dt.datetime.fromisoformat(ts.strip());require(x.tzinfo is None and x.microsecond==0,'characterization clock');t.append(x);u.append(num(aa));c.append(num(cc));q.append(num(yy))
    require(len(t)==nf.N==s['rows'] and all((b-a).total_seconds()==1 for a,b in zip(t,t[1:])),'characterization time/count differs')
    require(all(v is not None and v>=0 for v in u+c),'invalid characterization current')
    coverage,ix=nf.support(u,c,q);models,candidates,_,metrics=nf.identify(u,q,ix)
    return {'models':models,'candidates':candidates,'metrics':metrics,'support':coverage,'source_sha256':s['sha256'],'fit_clock':[str(t[0]),str(t[-1])]}

def strict_second(x):
    require(re.fullmatch(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}',x) is not None,'combined clock syntax');return dt.datetime.fromisoformat(x)
def combined_runs(p,s):
    cols=[s[k] for k in ('time_column','experiment_column','actual_column','command_column','target_column')];rows=csv_records(p.read_bytes(),s['header'],cols)
    groups={};current=None;seen=set()
    for i,(literal,ex,u,c,q) in enumerate(rows,1):
        require(bool(ex),'empty experiment label');t=strict_second(literal)
        if ex!=current:
            require(ex not in seen,'non-contiguous experiment');seen.add(ex);groups[ex]=[];current=ex
        g=groups[ex];require(not g or t>g[-1][1],'duplicate or reversed time');g.append((i,t,num(u),num(c),num(q),u,c,q))
    expected={g['experiment']:g for g in s['groups']};require(set(groups)==set(expected),'combined experiment identities differ')
    for ex,g in groups.items():
        a=expected[ex];gaps=sum((y[1]-x[1]).total_seconds()!=1 for x,y in zip(g,g[1:]));
        if 'rows' in a:require(len(g)==a['rows'] and g[0][1]==strict_second(a['first']) and g[-1][1]==strict_second(a['last']) and gaps==a['gaps'],'combined support differs: '+ex)
    return [(s['source'],ex,g) for ex,g in groups.items()]
def archive_run(p,s):
    b=zip_member(p,s['csv_member'],s['csv_bytes'],s['csv_sha256'],(s['metadata_member'],s['metadata_sha256']))
    rows=csv_records(b,s['header'],['',ACT,CMD,FLOW]);times,_=parse_clock([x[0] for x in rows]);delta=[(y-x).total_seconds() for x,y in zip(times,times[1:])];require(all(d>0 for d in delta) and 1 in delta,'archive clock invalid')
    return s['source'],s['experiment'],[(i,t,*map(num,row[1:]),*row[1:]) for i,(t,row) in enumerate(zip(times,rows),1)]
def nlr_groups(spec,paths,fit_clock=None):
    groups=[]
    for s in spec['nlr_combined']:groups.extend(combined_runs(paths[s['filename']],s))
    for s in spec['nlr_archives']:groups.append(archive_run(paths[s['filename']],s))
    require(len(groups)==23 and len({x[1] for x in groups})==23,'expected 23 unique trajectories')
    seen=set();period=[];fit0,fit1=map(dt.datetime.fromisoformat,fit_clock) if fit_clock else (None,None)
    for src,ex,g in groups:
        # Historical 2D and 2E each check within their own source collection; do not assert cross-collection physical independence.
        cluster='archive' if src in ('SIMULATED_WIND','SIMULATED_WAVE') else 'combined'
        lo,hi=g[0][1],g[-1][1];require(fit0 is None or hi<fit0 or lo>fit1,'fit time overlap');require(all(cluster!=v[0] or hi<v[1] or lo>v[2] for v in period),'recording time overlap within collection');period.append((cluster,lo,hi));digest=hashlib.sha256(''.join(repr(r[2:5])+'\n' for r in g).encode()).hexdigest();require((cluster,digest) not in seen,'duplicate exact sequence');seen.add((cluster,digest))
    return groups

def lag_reports(groups):
    summary=[];curves=[]
    for src,ex,g in groups:
        segs=[(np.array([r[3] for r in s]),np.array([r[2] for r in s])) for s in nc.split_segments(g)]
        for kind in ('level','increment'):
            x,ys=lag.pairs(segs,kind);rr=np.array([lag.corr(x,y) for y in ys]);p,ties=lag.peak(rr)
            def finite(v):return float(v) if np.isfinite(v) else None
            summary.append(dict(source=src,experiment=ex,kind=kind,n_common=len(x),defined_lags=int(np.isfinite(rr).sum()),r_zero=finite(rr[10]),peak_lag_s=p,r_peak=None if p is None else float(rr[p+10]),peak_separation=None if p is None else float(np.sort(rr)[-1]-np.sort(rr)[-2]),tied_lags=';'.join(map(str,ties)),boundary_peak=p in (-10,10)))
            for k,r in zip(lag.K,rr):curves.append(dict(source=src,experiment=ex,kind=kind,lag_s=int(k),n_common=len(x),r=finite(r)))
    return summary,curves

def validate_models(models,kind):
    for f in ('M0','M1'):
        require(f in models,'missing model '+f)
    if kind=='nlr':
        b0=models['M0']['beta_kg_h_per_A'];m=models['M1'];b=m['beta_kg_h_per_A'];t=m['tau_s'];require(not models['M0']['tau_s'] and len(b0)==len(b)==len(t)==1 and all(math.isfinite(x) and x>=0 for x in b0+b) and math.isfinite(t[0]) and t[0]>0,'invalid NLR models');return b0[0],b[0],t[0]
    for f,m in models.items():
        b=m['beta'];t=m['tau_s'];require(f in ('M0','M1','M2') and len(b)=={'M0':1,'M1':1,'M2':2}[f] and len(t)==(0 if f=='M0' else len(b)) and all(math.isfinite(x) and x>=0 for x in b) and all(math.isfinite(x) and x>0 for x in t),'invalid Juelich model')
    return models

def nlr_replay(groups,models,burn=300):
    b0,beta,tau=validate_models(models,'nlr');aa=math.exp(-1/tau);bb=beta*(-math.expm1(-1/tau));step=lambda y,u:aa*y+bb*u
    scores=collections.defaultdict(nc.Scores);coverage=[]
    for src,ex,g in groups:
        levels=(('EXPERIMENT',src,ex),('SOURCE_POOLED',src,''));n=0;segs=nc.split_segments(g)
        for level in levels:
            for h in HORIZONS:
                for subset in SUBSETS:scores[level+(h,subset)]
        for s in segs:
            z=[beta*s[0][2]]
            for r in s[:-1]:z.append(step(z[-1],r[2]))
            for i in range(burn,len(s)-30):
                if any(s[i+h][4] is None for h in HORIZONS):continue
                n+=1;u0,c0=s[i][2:4];a=c=hold=u=z[i];exc=0.
                for j in range(30):
                    x,y=s[i+j][2:4]
                    if j:exc=max(exc,abs(y-c0))
                    a=step(a,x);c=step(c,u0 if j==0 else y);hold=step(hold,u0 if j==0 else c0);u=step(u,u0);h=j+1
                    if h not in HORIZONS:continue
                    q=s[i+h][4];e=(q-a,a-c,c-hold,q-hold);p=(hold,c,a,u,b0*u0,b0*s[i+h][2]);require(abs(a-z[i+h])<=1e-9 and abs(e[3]-sum(e[:3]))<=1e-9,'replay identity')
                    if h==1:require(a==c==hold==u,'H1 information boundary')
                    subs=('ALL','CHANGE','LARGE_CHANGE_GE30A' if exc>=30 else 'SMALL_CHANGE_LT30A') if exc>0 else ('ALL','NO_CHANGE')
                    for lev in levels:
                        for sub in subs:scores[lev+(h,sub)].add(q,p,e)
        coverage.append(dict(source=src,experiment=ex,rows=len(g),segments=len(segs),anchors=n))
    rows=[dict(level=lev,source=src,experiment=ex,h_s=h,subset=sub,**v.report()) for (lev,src,ex,h,sub),v in scores.items()]
    macro=[]
    for src in dict.fromkeys(s for s,_,_ in groups):
        for h in HORIZONS:
            for sub in SUBSETS:
                r=[x for x in rows if x['level']=='EXPERIMENT' and x['source']==src and x['h_s']==h and x['subset']==sub and x['n']>0]
                keys=set(k for x in r for k in x if k.startswith(('rmse_','mae_','bias_','r2_','mse_','mean_','cross_')));m=dict(source=src,h_s=h,subset=sub,eligible_experiments=len(r),total_experiments=sum(x[0]==src for x in groups),n=sum(x['n'] for x in r));m.update({k:math.fsum(x[k] for x in r)/len(r) if all(x[k] is not None for x in r) else None for k in sorted(keys)});macro.append(m)
    return rows,macro,coverage

class DiscardTrace:
    def writerow(self,row): pass

def load_juelich(items,paths):
    for m in items:
        yield m,*jr.extract(paths[m['filename']],m,DiscardTrace())
def prepare_juelich(items,paths):
    groups={p:[] for p in PARTS};coverage=[];seen=set()
    for m,segs,counts in load_juelich(items,paths):
        require(m['partition'] in PARTS,'non-development partition');d=jf.design(segs);n=0 if d is None else len(d[1]);require(n==m['eligible_rows'],'development support differs '+str(m['file_id']));coverage.append(dict(file_id=m['file_id'],partition=m['partition'],eligible_rows=n,segments=len(segs),**counts))
        if d is not None:
            X,y,*_=d;sig=hashlib.sha256(np.column_stack([X[:,0],y]).astype('<f8').tobytes()).hexdigest();require(sig not in seen,'duplicate development sequence');seen.add(sig);groups[m['partition']].append((m,d))
    for p,g in groups.items():require(len(g)>=2 and (p=='DEV_DIAGNOSTIC' or max(d[4] for _,d in g)-min(d[3] for _,d in g)>=30),'insufficient group count/excitation '+p)
    return groups,coverage

def fit_juelich(items,paths):
    groups,coverage=prepare_juelich(items,paths);candidates,models=jf.candidates([d for _,d in groups['FIT']],[d for _,d in groups['SELECT']]);file_metrics=[]
    for part,gg in groups.items():
        for m,(X,y,*_) in gg:
            for k,v in models.items():file_metrics.append(dict(file_id=m['file_id'],partition=part,model=k,**jf.metric(y,X[:,v['indices']]@np.asarray(v['beta']))))
    return models,candidates,file_metrics,coverage

def juelich_fixed_replay(items,paths,models):
    validate_models(models,'juelich');fm=[];rm=[];coverage=[];seen=set()
    for m,segs,counts in load_juelich(items,paths):
        records=[jp.simulate(s,models) for s in segs if s[-1][2]-s[0][2]>=300];eligible=sum(len(v[4]) for v in records);n=eligible if eligible>=30 else 0;nr=0;fid=m['file_id'];camp=m['campaign'];meta=dict(file_id=fid,campaign=camp,abort_name=m['abort_name'])
        if n:
            y=np.concatenate([v[3][v[4]] for v in records]);u=np.concatenate([v[2][v[4]] for v in records]);sig=hashlib.sha256(np.column_stack([u,y]).astype('<f8').tobytes()).hexdigest();require(sig not in seen,'duplicate eligible evaluation sequence');seen.add(sig)
            for k in models:fm.append(dict(**meta,model=k,**jp.metric(y,np.concatenate([v[5][k][v[4]] for v in records]))))
            for h in (10,30):
                rr=[(y,hp,ap) for v in records for hh,i,j,el,y,hp,ap in v[-1] if hh==h and len(i)]
                if rr:
                    y,hp,ap=(np.concatenate([r[j] for r in rr]) for j in range(3));er=y-ap;ef=ap-hp;et=y-hp;require(np.max(np.abs(et-er-ef))<=1e-9,'Juelich replay identity');nr=len(y) if h==30 else nr;rm.append(dict(**meta,h=h,n=len(y),rmse_H=float(np.sqrt(np.mean(et**2))),rmse_A=float(np.sqrt(np.mean(er**2))),mse_R=float(np.mean(er**2)),mse_F=float(np.mean(ef**2)),interaction=float(2*np.mean(er*ef)),mse_total=float(np.mean(et**2))))
        coverage.append(dict(**meta,segments=len(segs),eligible_rows=n,paired_origins=nr,**counts))
    return fm,rm,coverage

def timing_prepared(items,paths,models):
    validate_models(models,'juelich');require('M2' in models,'legacy M2 reference required');prepared=[];seen=set()
    for m,segs,_ in load_juelich(items,paths):
        require(m['partition'] in PARTS,'timing stage development only');X=[];Y=[];refs=[];changes=[];legacy=[];actual=[]
        for sn,s in enumerate(segs,1):
            t=np.asarray([r[2] for r in s]);u=np.asarray([r[4] for r in s]);q=np.asarray([np.nan if r[5] is None else r[5] for r in s]);ix=np.flatnonzero((t-t[0]>=300)&np.isfinite(q))
            if not len(ix):continue
            a=tc.features(t,u,ix);X.append(a);Y.append(q[ix]);actual.append(u[ix]);changes.append(tc.consumed_change_10s(t,u,ix));refs.extend((sn,s[i][0]) for i in ix)
            m0=models['M0']['beta'][0]*u[ix];m1=sum(b*a[:,tc.SPECS.index((float(tau),0.))] for b,tau in zip(models['M1']['beta'],models['M1']['tau_s']));m2=sum(b*a[:,tc.SPECS.index((float(tau),0.))] for b,tau in zip(models['M2']['beta'],models['M2']['tau_s']));legacy.append(np.column_stack([m0,m1,m2]))
        n=len(refs) if len(refs)>=30 else 0;require(n==m['eligible_rows'],'timing support differs '+str(m['file_id']))
        if n:
            y=np.concatenate(Y);sig=hashlib.sha256(np.column_stack([np.concatenate(actual),y]).astype('<f8').tobytes()).hexdigest();require(sig not in seen,'duplicate timing sequence');seen.add(sig);prepared.append(dict(group=tc.Group(m['file_id'],m['partition'],np.concatenate(X),y),keys=refs,strata=tc.change_strata(np.concatenate(changes)),legacy_predictions=np.concatenate(legacy)))
    return prepared

def macro_juelich(rows,keys,metrics):
    by=collections.defaultdict(list)
    for r in rows:by[tuple(r[k] for k in keys)].append(r)
    return [dict(zip(keys,k),files=len(v),rows=sum(x['n'] for x in v),**{m:float(np.mean([r[m] for r in v])) if all(r[m] is not None for r in v) else None for m in metrics}) for k,v in by.items()]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=('nlr-fit','nlr-replay','nlr-lag','juelich-fit','juelich-replay','juelich-timing'));p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--models',type=Path);p.add_argument('--preflight-only',action='store_true');p.add_argument('--local-reproduction',action='store_true',help='explicitly request calculations on independent lawful copies; not a publication approval');a=p.parse_args();spec=read_json(ROOT/'acquisition/SOURCE_SPEC.json');source=a.source.resolve();out=a.out.resolve();require(source.is_dir() and not out.exists() and not out.is_relative_to(source) and not source.is_relative_to(out),'separate source and NEW output directory required');require(a.preflight_only or a.local_reproduction,'explicit local-reproduction flag required');paths=preflight(spec,source,a.stage)
    if a.preflight_only:print(json.dumps(dict(status='INPUT_IDENTITY_ONLY',stage=a.stage,files=len(paths))));return
    models=None;fit_clock=None
    if a.stage in ('nlr-replay','juelich-replay','juelich-timing'):
        require(a.models is not None,'explicit model JSON required');bundle=read_json(a.models);models=bundle.get('models',bundle);validate_models(models,'nlr' if a.stage=='nlr-replay' else 'juelich')
        if a.stage=='nlr-replay':
            require(bundle.get('source_sha256')==spec['characterization']['sha256'] and isinstance(bundle.get('fit_clock'),list) and len(bundle['fit_clock'])==2,'use MODELS.json from the independent NLR characterization stage');fit_clock=bundle['fit_clock']
    out.mkdir(parents=True,exist_ok=False)
    # No output traces, raw data copies, per-origin prediction dumps or downloads. Reports local only.
    try:
        if a.stage=='nlr-fit':
            d=fit_nlr(spec,paths);put_json(out/'MODELS.json',{k:d[k] for k in ('models','source_sha256','fit_clock')});put_json(out/'IDENTIFICATION.json',d)
        elif a.stage in ('nlr-replay','nlr-lag'):
            g=nlr_groups(spec,paths,fit_clock)
            if a.stage=='nlr-lag':s,c=lag_reports(g);put_csv(out/'TRAJECTORY_LAG_SUMMARY.csv',s);put_csv(out/'LAG_CURVES_21.csv',c)
            else:r,m,c=nlr_replay(g,models);put_csv(out/'DECOMPOSITION_METRICS.csv',r);put_csv(out/'DOMAIN_MACRO.csv',m);put_csv(out/'COVERAGE.csv',c)
        elif a.stage=='juelich-fit':
            m,c,r,v=fit_juelich(spec['juelich_development'],paths);put_json(out/'MODELS.json',m);put_json(out/'CANDIDATES.json',c);put_csv(out/'FILE_METRICS.csv',r);put_csv(out/'COVERAGE.csv',v);put_csv(out/'MACRO_METRICS.csv',macro_juelich(r,('partition','model'),('rmse','mae','bias','r2')))
        elif a.stage=='juelich-replay':
            f,r,c=juelich_fixed_replay(spec['juelich_evaluation'],paths,models);put_csv(out/'CONDITIONAL_METRICS.csv',f);put_csv(out/'REPLAY_METRICS.csv',r);put_csv(out/'COVERAGE.csv',c);put_csv(out/'CONDITIONAL_MACRO.csv',macro_juelich(f,('campaign','model'),('rmse','mae','bias','r2')));put_csv(out/'REPLAY_MACRO.csv',macro_juelich(r,('campaign','h'),('rmse_H','rmse_A','mse_R','mse_F','interaction','mse_total')))
        else:
            prepared=timing_prepared(spec['juelich_development'],paths,models);reports=tr.evaluate(prepared,models)
            for name,value in reports.items():
                if isinstance(value,list):put_csv(out/(name+'.csv'),value)
                else:put_json(out/(name+'.json'),value)
        put_json(out/'STATUS.json',dict(status='LOCAL_CALCULATION_COMPLETED_UNCOMPARED',stage=a.stage,input_count=len(paths),public_release=False,reference_equivalence_certified=False))
    except Exception as e:
        put_json(out/'STATUS.json',dict(status='FAILED',stage=a.stage,error=str(e)));raise
    print(out/'STATUS.json')
if __name__=='__main__':main()
