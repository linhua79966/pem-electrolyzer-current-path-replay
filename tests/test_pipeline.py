"""Synthetic file-format and core workflow tests; no real data, no network."""
import csv,datetime as dt,hashlib,io,json,math,shutil,subprocess,sys,tempfile,unittest,zipfile
from pathlib import Path
import xml.sax.saxutils as sx
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import reproduce as r
import timing_core as tc
from juelich_identification import metric
from response import zoh,matched_paths
from nlr_fit import state,BLOCKS
from nlr_clock import parse_clock
from lag import selftests
from compare_reports import compare
from plot_local_reports import make_figures

def sha(b):return hashlib.sha256(b).hexdigest()
def ident(p):return dict(filename=p.name,bytes=p.stat().st_size,sha256=r.sha_file(p))
def csv_bytes(header,rows):
    f=io.StringIO(newline='');w=csv.writer(f);w.writerow(header);w.writerows(rows);return f.getvalue().encode()
def signals(n,offset=0):
    t=np.arange(n);a=150+30*np.sin((t+offset)/7)+20*np.cos((t+offset)/17)+offset*.1;c=150+30*np.sin((t+offset+2)/7)+20*np.cos((t+offset+2)/17)+offset*.1;q=.007*zoh(t,a,8);return a,c,q

def xlsx(p,offset=0,n=210,bad=False):
    ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main';rel='http://schemas.openxmlformats.org/officeDocument/2006/relationships';pkg='http://schemas.openxmlformats.org/package/2006/relationships';rows=[]
    headers={'A':'Messnummer','B':'Datum','C':'Zeit','D':'t','E':'t','H':'I_Stack-1','K':'I_Stack-2','AK':'FI-24'}
    def st(col,rr,value):return f'<c r="{col}{rr}" t="inlineStr"><is><t>{sx.escape(value)}</t></is></c>'
    rows.append('<row r="51">'+''.join(st(k,51,v) for k,v in headers.items())+'</row>');rows.append('<row r="52">'+''.join(st(k,52,v) for k,v in {'D':'s','E':'h','H':'[A]','K':'[A]','AK':'[m³/h]'}.items())+'</row>')
    t=np.arange(n)*2.;u=150+40*np.sin((t+offset)/8)+offset*.15;q=.01*zoh(t,u,3)
    for i in range(n):
        rn=i+53;vv={'A':i+1,'B':45000.,'C':(3600+t[i])/86400.,'D':t[i],'E':t[i]/3600,'H':u[i]*.6,'K':u[i]*.4,'AK':q[i]}
        if bad and i==2:vv['H']=-1
        rows.append(f'<row r="{rn}">'+''.join(f'<c r="{k}{rn}"><v>{v:.17g}</v></c>' for k,v in vv.items())+'</row>')
    with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('xl/workbook.xml',f'<workbook xmlns="{ns}" xmlns:r="{rel}"><sheets><sheet name="Data" sheetId="1" r:id="rId1"/></sheets></workbook>');z.writestr('xl/_rels/workbook.xml.rels',f'<Relationships xmlns="{pkg}"><Relationship Id="rId1" Type="{rel}/worksheet" Target="worksheets/sheet1.xml"/></Relationships>');z.writestr('xl/worksheets/sheet1.xml',f'<worksheet xmlns="{ns}"><sheetData>'+''.join(rows)+'</sheetData></worksheet>')
    return dict(**ident(p),file_id=offset+1,sheet='Data',partition='FIT',fields={'I_Stack-1':'H','I_Stack-2':'K','FI-24':'AK'},time_rows=n,mode_s=2.,eligible_rows=n-150,abort_name=False,campaign='SYNTHETIC')

class Pipeline(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def test_identity_and_paths(self):
        f=self.p/'a.csv';f.write_text('a');d=ident(f);self.assertEqual(r.check_path(self.p,f.name,d['bytes'],d['sha256']),f);f.write_text('b');self.assertRaises(ValueError,r.check_path,self.p,f.name,d['bytes'],d['sha256']);self.assertRaises(ValueError,r.check_path,self.p,'../a.csv',1,d['sha256'])
    def test_csv_header_and_clock(self):
        self.assertRaises(ValueError,r.csv_records,b'a,a\n1,2\n',['a','a'],['a']);self.assertRaises(ValueError,r.csv_records,b'a,b\n1\n',['a','b'],['a']);self.assertRaises(ValueError,parse_clock,['2/3/2025 1:02:03']);self.assertRaises(ValueError,r.strict_second,'2025-01-01 00:00');self.assertEqual(parse_clock(['2/13/2025 1:02:03','2/13/2025 1:02:04'])[1],'SLASH_MDY_DISAMBIGUATED')
    def test_nlr_characterization_fit(self):
        n=18000;t=np.arange(n);a=180+60*np.sin(t/17)+20*np.cos(t/9);c=a.copy();q=np.array(state(a,8,BLOCKS))*.007;ts=dt.datetime(2025,7,11,10,44,28);head=['',r.ACT,r.CMD,r.FLOW];b=csv_bytes(head,[(str(ts+dt.timedelta(seconds=i)),a[i],c[i],q[i]) for i in range(n)]);p=self.p/'char.zip';meta=b'synthetic metadata';pic=b'synthetic image placeholder';members={'characterization_200.csv':b,'electrolyzer_metadata.xlsx':meta,'characterization_200-timeseries.png':pic}
        with zipfile.ZipFile(p,'w') as z:
            for k,v in members.items():z.writestr(k,v)
        s=dict(**ident(p),header=head,rows=n,members=[dict(member=k,bytes=len(v),sha256=sha(v)) for k,v in members.items()]);res=r.fit_nlr({'characterization':s},{p.name:p});self.assertEqual(res['models']['M1']['tau_s'],[8.]);self.assertAlmostEqual(res['models']['M1']['beta_kg_h_per_A'][0],.007,places=12);self.assertEqual(res['support'][0]['eligible'],8700)
    def test_combined_and_archive_then_replay_lag(self):
        n=420;a,c,q=signals(n);ts=dt.datetime(2025,8,13);times=[str(ts+dt.timedelta(seconds=i)) for i in range(n)];head=['','Experiment',r.ACT,r.CMD,r.FLOW];p=self.p/'combined.csv';p.write_bytes(csv_bytes(head,zip(times,['ex1']*n,a,c,q)));s=dict(**ident(p),source='HISTORICAL_WIND',header=head,time_column='',experiment_column='Experiment',actual_column=r.ACT,command_column=r.CMD,target_column=r.FLOW,groups=[dict(experiment='ex1',rows=n,first=times[0],last=times[-1],gaps=0)]);g=r.combined_runs(p,s)
        ts+=dt.timedelta(days=1);head2=['',r.ACT,r.CMD,r.FLOW];b=csv_bytes(head2,zip([str(ts+dt.timedelta(seconds=i)) for i in range(n)],a,c,q));p2=self.p/'run.zip';meta=b'metadata'
        with zipfile.ZipFile(p2,'w') as z:z.writestr('run.csv',b);z.writestr('meta.xlsx',meta)
        s2=dict(source='SIMULATED_WAVE',experiment='ex2',header=head2,csv_member='run.csv',csv_bytes=len(b),csv_sha256=sha(b),metadata_member='meta.xlsx',metadata_sha256=sha(meta));g.append(r.archive_run(p2,s2));models={'M0':{'tau_s':[],'beta_kg_h_per_A':[.007]},'M1':{'tau_s':[8.],'beta_kg_h_per_A':[.007]}}
        rows,macro,cov=r.nlr_replay(g,models);self.assertEqual([v['anchors'] for v in cov],[90,90]);h1=[x for x in rows if x['h_s']==1 and x['subset']=='ALL'];self.assertTrue(all(x['rmse_A']==x['rmse_C']==x['rmse_U']==x['rmse_H'] for x in h1));self.assertTrue(all(abs(x['mse_identity_gap'])<1e-9 for x in rows if x['n']));s3,curves=r.lag_reports(g);self.assertEqual(len(curves),84);self.assertEqual(len(s3),4);self.assertTrue(all(x['peak_lag_s']==2 for x in s3));self.assertEqual(len(selftests()),12)
        out=self.p/'replay';out.mkdir();r.put_csv(out/'DOMAIN_MACRO.csv',macro);make_figures('nlr-replay',out,self.p/'plots');self.assertEqual(len(list((self.p/'plots').glob('*.png'))),2)
        lagout=self.p/'lag';lagout.mkdir();r.put_csv(lagout/'LAG_CURVES_21.csv',curves);make_figures('nlr-lag',lagout,self.p/'lagplots');self.assertEqual(len(list((self.p/'lagplots').glob('*.png'))),2)
    def test_invalid_row_and_common_support(self):
        a,c,q=signals(420);ts=dt.datetime(2025,8,13);g=[(i,ts+dt.timedelta(seconds=i),float(a[i]),float(c[i]),float(q[i])) for i in range(420)];g[100]=(*g[100][:2],None,g[100][3],g[100][4]);self.assertEqual([len(s) for s in r.nc.split_segments(g)],[100,319]);m={'M0':{'tau_s':[],'beta_kg_h_per_A':[.007]},'M1':{'tau_s':[8.],'beta_kg_h_per_A':[.007]}};self.assertEqual(r.nlr_replay([('s','e',g)],m)[2][0]['anchors'],0)
    def test_juelich_headers_segment_fit_and_timing(self):
        spec=[];paths={}
        for i,part in enumerate(['FIT','FIT','SELECT','SELECT','DEV_DIAGNOSTIC','DEV_DIAGNOSTIC']):
            p=self.p/f'{i}.xlsx';s=xlsx(p,offset=i*10);s['partition']=part;spec.append(s);paths[p.name]=p
        m,c,rows,cov=r.fit_juelich(spec,paths);self.assertTrue({'M0','M1'}<=set(m));self.assertEqual(m['M1']['tau_s'],[3]);self.assertTrue(all(x['eligible_rows']==60 for x in cov));self.assertEqual(len(c),len({(x['family'],tuple(x['tau_s'])) for x in c}));self.assertAlmostEqual(metric(np.array([1.,2.]),np.array([1.,2.]))['rmse'],0.)
        # Two-timescale reference provided explicitly for diagnostic even if synthetic fit rejects M2.
        m['M2']={'tau_s':[1,5],'beta':[.003,.007]};v=r.timing_prepared(spec,paths,m);out=r.tr.evaluate(v,m);self.assertEqual(len(out['CANDIDATES']),126);self.assertEqual(len(out['FIT_MOMENTS']),252);self.assertEqual(len(out['WINNERS']),4);self.assertEqual(len(out['FILE_METRICS']),48)
        # Verify diagnostics do not affect fitting or selection.
        groups=[x['group'] for x in v];_,w,_=tc.fit_select(groups);groups2=[tc.Group(x.file_id,x.partition,x.X,x.y+100 if x.partition=='DEV_DIAGNOSTIC' else x.y) for x in groups];_,w2,_=tc.fit_select(groups2);self.assertEqual(w,w2)
    def test_juelich_fixed_replay(self):
        p=self.p/'x.xlsx';s=xlsx(p);s['partition']='RESERVED_2023';m={'M0':{'tau_s':[],'beta':[.01]},'M1':{'tau_s':[3],'beta':[.01]},'M2':{'tau_s':[1,5],'beta':[.003,.007]}};f,z,c=r.juelich_fixed_replay([s],{p.name:p},m);self.assertEqual(c[0]['eligible_rows'],60);self.assertEqual(c[0]['paired_origins'],45);self.assertEqual([x['n'] for x in z],[45,45]);self.assertTrue(all(abs(x['mse_total']-x['mse_R']-x['mse_F']-x['interaction'])<1e-12 for x in z))
        p2=self.p/'bad.xlsx';s2=xlsx(p2,bad=True);segs,counts=r.jr.extract(p2,s2,r.DiscardTrace());self.assertEqual(counts['invalid_input'],1);self.assertEqual(len(segs),2)
    def test_report_comparison(self):
        a=self.p/'a.csv';b=self.p/'b.csv';a.write_text('id,n,rmse\n1,3,0.123456789\n');b.write_text('id,n,rmse\n1,3,0.1234567889\n');self.assertEqual(compare(a,b,['id'],['n','rmse'])['status'],'PASS_SELECTED_COLUMNS');b.write_text('id,n,rmse\n1,4,0.1234567889\n');self.assertEqual(compare(a,b,['id'],['n','rmse'])['status'],'FAIL');b.write_text('id,n,rmse\n2,3,0.1234567889\n');self.assertEqual(compare(a,b,['id'],['n','rmse'])['status'],'FAIL')
    def test_cli_preflight_and_no_overwrite(self):
        # Only in a disposable synthetic software copy; real source spec untouched.
        copy=self.p/'software';shutil.copytree(Path(r.ROOT),copy,ignore=shutil.ignore_patterns('__pycache__'));source=self.p/'source';source.mkdir();p=source/'fake.xlsx';s=xlsx(p);s['partition']='RESERVED_2023';spec={'juelich_evaluation':[s]};(copy/'acquisition/SOURCE_SPEC.json').write_text(json.dumps(spec));model=self.p/'model.json';model.write_text(json.dumps({'M0':{'tau_s':[],'beta':[.01]},'M1':{'tau_s':[3],'beta':[.01]}}));out=self.p/'out';args=[sys.executable,'-B',str(copy/'src/reproduce.py'),'juelich-replay','--source',str(source),'--out',str(out),'--models',str(model)]
        self.assertNotEqual(subprocess.run(args,capture_output=True).returncode,0);self.assertFalse(out.exists());self.assertEqual(subprocess.run(args+['--preflight-only'],capture_output=True).returncode,0);self.assertFalse(out.exists());done=subprocess.run(args+['--local-reproduction'],capture_output=True);self.assertEqual(done.returncode,0,done.stderr.decode());self.assertEqual(json.loads((out/'STATUS.json').read_text())['status'],'LOCAL_CALCULATION_COMPLETED_UNCOMPARED');self.assertNotEqual(subprocess.run(args+['--local-reproduction'],capture_output=True).returncode,0)
if __name__=='__main__':unittest.main(verbosity=2)
