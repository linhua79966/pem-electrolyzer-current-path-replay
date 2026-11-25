import math,posixpath
import xml.etree.ElementTree as ET
from xlsx_headers import NS,REL

def sheet_path(z,name):
    wb=ET.fromstring(z.read('xl/workbook.xml')); rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    s=[x for x in wb.iter(NS+'sheet') if x.get('name')==name]
    if len(s)!=1:raise ValueError('Frozen worksheet missing or ambiguous')
    rel=[x for x in rels if x.get('Id')==s[0].get(REL+'id')]
    if len(rel)!=1 or rel[0].get('TargetMode')=='External' or not rel[0].get('Type','').endswith('/worksheet'):raise ValueError('Invalid worksheet link')
    t=rel[0].get('Target',''); p=posixpath.normpath(t.lstrip('/') if t.startswith('/') else posixpath.join('xl',t))
    if not p.startswith('xl/') or p not in z.namelist():raise ValueError('Invalid worksheet path')
    return p

def num(c):
    if c is None or c.get('t','n')!='n':return None
    v=c.find(NS+'v')
    try:x=float(v.text) if v is not None else math.nan;return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None

def present(c):return c is not None and (c.find(NS+'v') is not None or c.find(NS+'is') is not None or c.find(NS+'f') is not None)
