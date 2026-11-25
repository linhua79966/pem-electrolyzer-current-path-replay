import collections,csv,math,zipfile,re
import xml.etree.ElementTree as ET
import xlsx_headers as H
import xlsx_cells as C
NS=H.NS

def clock_ok(a,b,c,d,e):
    return all(x is not None for x in (a,b,c,d,e)) and abs(b-round(b))<1e-8 and 0<=c<1 and abs(d-3600*e)<=.1

def edge_ok(old,new,mode):
    # row, counter, relative seconds, absolute serial seconds, current sum, target
    return new[0]==old[0]+1 and new[1]-old[1]==1 and 0<new[2]-old[2]<=1.5*mode and abs((new[3]-old[3])-(new[2]-old[2]))<=.1

def extract(p,m,trace):
    sh=[s for s in H.xlsx_headers(p,[50,51,52]) if s['sheet']==m['sheet']]
    if len(sh)!=1:raise ValueError('sheet header mismatch')
    vals={x['cell']:str(x['value']).strip() for x in sh[0]['cells']}
    if [vals.get(k+'51') for k in 'ABCDE']!=['Messnummer','Datum','Zeit','t','t'] or vals.get('D52')!='s' or vals.get('E52')!='h':raise ValueError('clock header mismatch')
    cols=m['fields'];
    for tag,col in cols.items():
        if vals.get(col+'51')!=tag or vals.get(col+'52')!=('[m³/h]' if tag=='FI-24' else '[A]'):raise ValueError('signal header mismatch '+tag)
    if len(set(cols.values()))!=3 or any(col in 'ABCDE' for col in cols.values()):raise ValueError('invalid mapping')
    wanted=set('ABCDE')|set(cols.values());segments=[];cur=[];counts=collections.Counter();last=None;mode=float(m['mode_s'])
    def close():
        nonlocal cur
        if cur:segments.append(cur);cur=[]
    with zipfile.ZipFile(p) as z,z.open(C.sheet_path(z,m['sheet'])) as f:
        for _,row in ET.iterparse(f,events=('end',)):
            if row.tag!=NS+'row':continue
            rn=int(row.get('r','0'))
            if rn>=53:
                cells={re.sub(r'\d+$','',c.get('r','')):c for c in row.findall(NS+'c') if re.sub(r'\d+$','',c.get('r','')) in wanted}
                if any(C.present(cells.get(k)) for k in 'ABCDE'):
                    counts['time_rows']+=1;a,b,c,d,e=[C.num(cells.get(k)) for k in 'ABCDE'];i1,i2,q=[C.num(cells.get(cols[k])) for k in ('I_Stack-1','I_Stack-2','FI-24')]
                    ck=clock_ok(a,b,c,d,e);iv=i1 is not None and i2 is not None and i1>=0 and i2>=0;counts['invalid_clock']+=not ck;counts['invalid_input']+=not iv;counts['missing_target']+=q is None;counts['negative_target']+=q is not None and q<0
                    if not ck or not iv:close();last=None;sn=None;pos=None
                    else:
                        v=(rn,a,d,(b+c)*86400,i1+i2,q)
                        if last is None or not edge_ok(last,v,mode):close()
                        cur.append(v);last=v;sn=len(segments)+1;pos=len(cur)-1
                    trace.writerow([m['file_id'],m['partition'],rn,a,b,c,d,e,i1,i2,q,int(ck),int(iv),sn,pos])
            row.clear()
    close()
    if counts['time_rows']!=m['time_rows']:raise ValueError('time row count differs')
    return segments,dict(counts)
