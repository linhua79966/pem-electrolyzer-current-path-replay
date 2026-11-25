import datetime as dt,re
TIME=re.compile(r"\d{4}[-/]\d{2}[-/]\d{2}[ T]\d{2}:\d{2}:\d{2}(?:\.\d{1,6})?\Z")
SLASH=re.compile(r"(\d{1,2})/(\d{1,2})/(\d{4}) (\d{1,2}):(\d{2}):(\d{2})(?:\.\d{1,6})?\Z")
def need(ok,msg):
    if not ok: raise ValueError(msg)

def parse_clock(literals):
    # Year-first timestamps carry an unambiguous calendar date. Slash-date formats are accepted only if at least one row disambiguates date order and all rows agree. No guessed seconds, no reordering.
    if all(TIME.fullmatch(x) for x in literals):return [dt.datetime.fromisoformat(x.replace('/','-')) for x in literals],'YEAR_FIRST_EXPLICIT_SECONDS'
    ss=[SLASH.fullmatch(x) for x in literals]
    if all(ss):
        mdy=any(int(m[2])>12 and int(m[1])<=12 for m in ss);dmy=any(int(m[1])>12 and int(m[2])<=12 for m in ss)
        need(mdy!=dmy,'ambiguous or mixed month/day order')
        dates=[]
        for x,m in zip(literals,ss):
            mo,day=(int(m[1]),int(m[2])) if mdy else (int(m[2]),int(m[1]));suffix=x.split(' ',1)[1];dates.append(dt.datetime.fromisoformat(f'{int(m[3]):04}-{mo:02}-{day:02} {int(m[4]):02}:{m[5]}:{m[6]}'+('.'+suffix.split('.',1)[1] if '.' in suffix else '')))
        return dates,'SLASH_MDY_DISAMBIGUATED' if mdy else 'SLASH_DMY_DISAMBIGUATED'
    raise ValueError('time syntax lacks explicit seconds or unsupported/mixed date format')
