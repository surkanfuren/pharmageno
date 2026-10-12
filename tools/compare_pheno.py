import json, re, sys, urllib.request, urllib.parse
sys.path.insert(0, '.')
from xlsx_read import read
from collections import Counter
rows = read('../samples/real/getrm_consolidated.xlsx')['PGx HLA Genotypes']
hdr = [h.strip() for h in rows[0]]
print('DPYD-like headers:', [h for h in hdr if 'DPYD' in h.upper()])
calls = json.load(open('data/calls_1000g.json'))
cache = {}
def pheno(g, d):
    d = re.sub(r'\s+', ' ', d).strip()
    if (g, d) in cache: return cache[g, d]
    res = None
    parts = d.split('/')
    for x in ([d, '/'.join(parts[::-1])] if len(parts) == 2 else [d]):
        q = urllib.parse.urlencode({'genesymbol': f'eq.{g}', 'diplotype': f'eq.{x}', 'select': 'generesult'})
        r = json.load(urllib.request.urlopen(f'https://api.cpicpgx.org/v1/diplotype?{q}'))
        if r: res = r[0]['generesult']; break
    cache[g, d] = res; return res
genes = [g for g in ['CYP2C19', 'CYP2C9', 'SLCO1B1', 'TPMT', 'NUDT15'] if g in hdr]
tot = Counter()
for g in genes:
    ci = hdr.index(g); s = Counter(); bad = []
    for r in rows[1:]:
        sid = r[1].strip() if len(r) > 1 else ''
        if sid not in calls or ci >= len(r) or not r[ci].strip(): continue
        truth_d = r[ci].strip(); ours_d = calls[sid][g]
        tp, op = pheno(g, truth_d), pheno(g, ours_d)
        if tp is None: s['truth not mappable to current CPIC names'] += 1; continue
        s['compared'] += 1
        if tp == op: s['same phenotype'] += 1
        else:
            s['different phenotype'] += 1
            falsenormal = op and op.startswith('Normal') and not tp.startswith('Normal')
            if falsenormal: s['  of which: we said NORMAL, truth abnormal'] += 1
            bad.append((sid, truth_d, tp, ours_d, op))
    print(f'\n{g}:', dict(s))
    for b in Counter((t, tp, o, op) for _, t, tp, o, op in bad).most_common(12): print('   ', b)
    tot += s
print('\nTOTAL', dict(tot))
