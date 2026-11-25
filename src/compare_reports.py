"""Compare locally regenerated report rows against lawfully held reference CSVs.
The reference is not downloaded or distributed by this tool. Numeric tolerance is explicit.
"""
import argparse,csv,json,math
from pathlib import Path

def load(path,keys,columns):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f);required=set(keys+columns)
        if not required<=set(r.fieldnames or []):raise ValueError('missing columns in '+str(path))
        out={}
        for row in r:
            key=tuple(row[k] for k in keys)
            if key in out:raise ValueError('duplicate key: '+str(key))
            out[key]=row
    return out

def compare(actual,reference,keys,columns,atol=1e-9,rtol=1e-9):
    if not keys or not columns or atol<0 or rtol<0:raise ValueError('invalid comparison settings')
    a=load(actual,keys,columns);b=load(reference,keys,columns);miss=sorted(set(b)-set(a));extra=sorted(set(a)-set(b));errors=[];n=0
    for k in sorted(set(a)&set(b)):
        for c in columns:
            x,y=a[k][c],b[k][c];n+=1
            if x==y:continue
            try:
                z,w=float(x),float(y);ok=math.isfinite(z) and math.isfinite(w) and (z==w if c in ('n','rows','files','peak_lag_s','n_common','eligible_rows','paired_origins','index','defined_lags') else math.isclose(z,w,abs_tol=atol,rel_tol=rtol))
            except (ValueError,TypeError):ok=False
            if not ok:errors.append(dict(key=k,column=c,actual=x,reference=y))
    return dict(status='PASS_SELECTED_COLUMNS' if n and not(miss or extra or errors) else 'FAIL',compared_cells=n,missing_keys=miss,extra_keys=extra,mismatches=errors,keys=keys,columns=columns,atol=atol,rtol=rtol)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('actual');p.add_argument('reference');p.add_argument('--keys',required=True);p.add_argument('--columns',required=True);p.add_argument('--atol',type=float,default=1e-9);p.add_argument('--rtol',type=float,default=1e-9);p.add_argument('--output',type=Path,required=True);a=p.parse_args();result=compare(a.actual,a.reference,a.keys.split(','),a.columns.split(','),a.atol,a.rtol);f=a.output.open('x',encoding='utf-8');json.dump(result,f,indent=2);f.close();print(result['status']);raise SystemExit(0 if result['status']=='PASS_SELECTED_COLUMNS' else 2)
