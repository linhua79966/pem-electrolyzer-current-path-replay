import zipfile,posixpath
import xml.etree.ElementTree as ET
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

def header_text(el):
    return ''.join(x.text or '' for x in el.iter(NS+'t'))

def xlsx_headers(p, wanted):
    """Read only XML header rows; shared-string table may be scanned to resolve header strings."""
    with zipfile.ZipFile(p) as z:
        names=z.namelist()
        if len(names)!=len(set(names)): raise ValueError('duplicate ZIP member')
        wb=ET.fromstring(z.read('xl/workbook.xml'))
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        relmap={e.get('Id'):e for e in rels}
        sheets=[]; needed=set()
        for s in wb.iter(NS+'sheet'):
            rel=relmap.get(s.get(REL+'id'))
            if rel is None: raise ValueError('unresolved worksheet relationship')
            if not rel.get('Type','').endswith('/worksheet'): continue
            if rel.get('TargetMode')=='External': raise ValueError('external worksheet')
            target=rel.get('Target','')
            path=posixpath.normpath(target.lstrip('/') if target.startswith('/') else posixpath.join('xl',target))
            if not path.startswith('xl/') or path not in names: raise ValueError('unsafe or missing worksheet')
            data=[]; seen=set()
            with z.open(path) as f:
                for event,e in ET.iterparse(f,events=('end',)):
                    if e.tag!=NS+'row': continue
                    rn=int(e.get('r','0'))
                    if rn>max(wanted): break
                    if rn in wanted:
                        seen.add(rn)
                        for c in e.findall(NS+'c'):
                            ref=c.get('r',''); typ=c.get('t','n'); v=c.find(NS+'v'); formula=c.find(NS+'f')
                            val=v.text if v is not None else None
                            if formula is not None: typ='formula_literal'; val='='+str(formula.text or '')
                            elif typ=='inlineStr': val=header_text(c)
                            elif typ=='s' and val is not None: val=int(val); needed.add(val)
                            if val is not None: data.append({'row':rn,'cell':ref,'type':typ,'value':val})
                    e.clear()
            sheets.append({'sheet':s.get('name'),'state':s.get('state','visible'),'present_rows':sorted(seen),'cells':data})
        vals={}
        if needed:
            if min(needed)<0 or 'xl/sharedStrings.xml' not in names: raise ValueError('invalid shared strings')
            with z.open('xl/sharedStrings.xml') as f:
                i=0
                for _,el in ET.iterparse(f,events=('end',)):
                    if el.tag==NS+'si':
                        if i in needed: vals[i]=header_text(el)
                        i+=1; el.clear()
                        if len(vals)==len(needed): break
            if set(vals)!=needed: raise ValueError('unresolved shared string')
        for s in sheets:
            for c in s['cells']:
                if c['type']=='s': c['value']=vals[c['value']]
        return sheets
