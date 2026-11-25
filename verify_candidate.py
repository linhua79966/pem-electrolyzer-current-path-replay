"""Safe aggregate/synthetic acceptance; does not call any raw-data workflow."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
import numpy,pandas,matplotlib
ROOT=Path(__file__).resolve().parent

def check_manifest(root):
    root=root.resolve();seen=set();lines=(root/'SHA256SUMS.txt').read_text(encoding='utf-8-sig').splitlines()
    for line in lines:
        if not line.strip():continue
        h,n=line.split(maxsplit=1);n=n.strip();p=root/n
        if n in seen or Path(n).is_absolute() or not p.resolve().is_relative_to(root) or not p.is_file():raise ValueError('invalid manifest entry '+n)
        if hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise ValueError('SHA256 mismatch '+n)
        seen.add(n)
    if not seen:raise ValueError('empty manifest')
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and not any(c in ('.git','.venv','__pycache__') for c in p.relative_to(root).parts) and p.name!='SHA256SUMS.txt'}
    if actual!=seen:raise ValueError('manifest coverage differs: '+str(sorted(actual^seen)))
    return len(seen)

def run(dataset,out):
    dataset=dataset.resolve();out=out.resolve()
    if out.exists() or out.is_relative_to(dataset) or out.is_relative_to(ROOT) or dataset.is_relative_to(out) or ROOT.is_relative_to(out):raise ValueError('new output directory outside software and dataset required')
    ns=check_manifest(ROOT);nd=check_manifest(dataset);out.mkdir(parents=True)
    status=dict(status='FAILED',software_manifest_entries=ns,dataset_manifest_entries=nd,python=sys.version.split()[0],numpy=numpy.__version__,pandas=pandas.__version__,matplotlib=matplotlib.__version__,original_signal_access=False,study_fit_or_replay=False)
    def call(args,name):
        r=subprocess.run([sys.executable,'-B',*map(str,args)],cwd=ROOT,capture_output=True,text=True)
        (out/(name+'.txt')).write_text(r.stdout+'\n'+r.stderr,encoding='utf-8')
        if r.returncode:raise RuntimeError(name+' failed: see log')
        return r.stdout+'\n'+r.stderr
    try:
        log=call(['-m','unittest','discover','-s','tests','-v'],'SYNTHETIC_TESTS');status['synthetic_tests_passed']=True
        call([ROOT/'src/audit_juelich.py',dataset,out/'aggregate'],'AGGREGATE_TEST');a=json.loads((out/'aggregate/AUDIT_RESULT.json').read_text());
        if a['checks']!=1363 or a['passed']!=1363:raise ValueError('unexpected aggregate check count')
        status.update(aggregate_checks=a['checks'],aggregate_passed=a['passed'])
        call([ROOT/'src/plot_timing_grid_v015.py',dataset/'data/04_REPORTS',out/'plots'],'PLOT_TIMING');call([ROOT/'src/plot_evidence_v015.py',dataset/'data/TIMING_DERIVED',out/'plots'],'PLOT_COMPOSITION')
        names=['F_S3_timing_identifiability.png','F_S4_timing_composition.png'];status['plot_sha256']={n:hashlib.sha256((out/'plots'/n).read_bytes()).hexdigest() for n in names};status['status']='PASS_SYNTHETIC_AND_AGGREGATE_ONLY';status['full_raw_reproduction']=False
    except Exception as e:status['error']=str(e);raise
    finally:(out/'ACCEPTANCE.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(status,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('dataset',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.dataset,a.output)
