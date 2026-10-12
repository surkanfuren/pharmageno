import json, re, sys
sys.path.insert(0, '.')
from xlsx_read import read
rows = read('../samples/real/getrm_consolidated.xlsx')['PGx HLA Genotypes']
hdr = [h.strip() for h in rows[0]]
col = {g: hdr.index(g) for g in ['CYP2C19', 'CYP2C9', 'SLCO1B1', 'TPMT', 'NUDT15', 'DPYD'] if g in hdr}
calls = json.load(open('data/calls_1000g.json'))
def norm(g, d):
    d = re.sub(r'\s+', '', d or '')
    if not d: return None
    return d
def canon(d):  # order-insensitive comparison
    return tuple(sorted(d.split('/'))) if d and '/' in d else (d,)
res = {g: {'match': 0, 'mismatch': [], 'n': 0} for g in col}
nsamp = 0
for r in rows[1:]:
    sid = (r[1] if len(r) > 1 else '').strip()
    if sid not in calls: continue
    nsamp += 1
    for g, ci in col.items():
        truth = norm(g, r[ci] if ci < len(r) else '')
        if not truth: continue
        ours = calls[sid][g].replace(' ', '')
        res[g]['n'] += 1
        # DPYD: GeT-RM notation may differ; report it as written
        if canon(ours) == canon(truth): res[g]['match'] += 1
        else: res[g]['mismatch'].append((sid, truth, ours))
print('GeT-RM samples also in 1000G:', nsamp)
for g, x in res.items():
    print(f"\n{g}: {x['match']} / {x['n']} match")
    from collections import Counter
    c = Counter((t, o) for _, t, o in x['mismatch'])
    for (t, o), n in c.most_common(25): print(f"   GeT-RM {t!r:40} ours {o!r:45} x{n}")
