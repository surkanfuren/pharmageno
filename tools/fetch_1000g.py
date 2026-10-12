"""1000 Genomes (phase 3) genotypes for our rsIDs via Ensembl REST -> data/geno_1000g.json
Usage (from tools/): python3 fetch_1000g.py   (15 rsIDs, ~1-2 min)"""
import json, time, urllib.request
m = json.load(open('../allele_map.json'))
rsids = sorted({r for g in m['genes'].values() for r in g['rsids']})
geno = {}
for rs in rsids:
    for attempt in range(4):
        try:
            req = urllib.request.Request(f'https://rest.ensembl.org/variation/human/{rs}?genotypes=1', headers={'Accept': 'application/json'})
            d = json.load(urllib.request.urlopen(req, timeout=60)); break
        except Exception:
            time.sleep(2 * (attempt + 1))
    for g in d.get('genotypes', []):
        if g['sample'].startswith('1000GENOMES:phase_3:'):
            geno.setdefault(g['sample'].split(':')[-1], {})[rs] = g['genotype'].replace('|', '').replace('/', '')
    print(rs, 'ok')
json.dump(geno, open('data/geno_1000g.json', 'w'))
print('samples', len(geno))
