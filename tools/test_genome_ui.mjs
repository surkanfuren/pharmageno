// Execute the single-file UI with a minimal DOM; no browser dependency.
import fs from 'node:fs';
import vm from 'node:vm';
import assert from 'node:assert/strict';
const html = fs.readFileSync(new URL('../index.html', import.meta.url), 'utf8');
const source = html.match(/<script>([\s\S]*?)<\/script>/)[1];
const elements = new Map();
const element = id => {
  if (!elements.has(id)) elements.set(id, { innerHTML: '', textContent: '', disabled: false, dataset: {}, style: {}, insertAdjacentHTML() {}, scrollIntoView() {}, focus() {}, classList: { add() {}, remove() {}, toggle() {} }, addEventListener() {} });
  return elements.get(id);
};
const doc = JSON.parse(fs.readFileSync(new URL('../samples/real/giab/HG001.research-passport.json', import.meta.url), 'utf8'));
const map = JSON.parse(fs.readFileSync(new URL('../allele_map.json', import.meta.url), 'utf8'));
let requestCount = 0, cpicHit = false;
const context = vm.createContext({ document: { getElementById: element, querySelectorAll: () => [] }, performance, setTimeout,
  URL, URLSearchParams, Blob, fetch: async url => {
    if (url === '/api/cpic/status') return {status:404};
    if (url === 'allele_map.json') return { ok: true, json: async () => map };
    if (url.includes('/diplotype?')) {
      const q = new URL(url, 'http://x').searchParams, gene = q.get('genesymbol').slice(3), dip = q.get('diplotype').slice(3);
      return { ok: true, json: async () => cpicHit ? [{ diplotype: dip, generesult: 'Normal Metabolizer', lookupkey: { [gene]: 'Normal Metabolizer' } }] : [] };
    }
    if (url.includes('/pair_view?')) { requestCount++; return { ok: true, json: async () => [{genesymbol:'CYP2C19'}] }; }
    throw Error('Unexpected API request: ' + url);
  }
});
vm.runInContext(source, context);
context.data = doc;
await vm.runInContext('showGenome(data)', context);
assert.equal(element('passportName').textContent, 'HG001');
assert.match(element('identityBadge').textContent, /not clinically validated/);
assert.equal(vm.runInContext('Object.values(passport).every(p => p.status === "unknown")', context), true);
const verdict = await vm.runInContext('checkDrug("clopidogrel", "clopidogrel")', context);
assert.equal(verdict.color, 'yellow');
assert.deepEqual(Array.from(verdict.missing), ['CYP2C19']);
// Even a tampered research JSON cannot promote genes to clinical status.
context.data = structuredClone(doc);
context.data.genes.CYP2C19.status = 'ok';
context.data.genes.CYP2C19.lookup = 'Normal Metabolizer';
await vm.runInContext('showGenome(data)', context);
assert.equal(vm.runInContext('passport.CYP2C19.status', context), 'unknown');
assert.equal(requestCount, 1);
// With CPIC matches, genes are called from the genome evidence (all key positions read), others stay unknown.
cpicHit = true;
context.data = structuredClone(doc);
await vm.runInContext('showGenome(data)', context);
assert.equal(vm.runInContext('passport.CYP2C19.status + "/" + passport.CYP2C19.source + "/" + passport.CYP2C19.diplotype', context), 'ok/genome/*1/*2');
assert.equal(vm.runInContext('passport.TPMT.status', context), 'unknown');      // rs1142345 not covered
assert.equal(vm.runInContext('passport.CYP2D6.status', context), 'unknown');    // not callable from SNP evidence
// Removing the evidence makes the gene unknown again, whatever the file's own status field says.
context.data = structuredClone(doc);
context.data.evidence = context.data.evidence.filter(site => !(site.rsids || []).includes('rs4244285'));
context.data.genes.CYP2C19.status = 'ok';
await vm.runInContext('showGenome(data)', context);
assert.equal(vm.runInContext('passport.CYP2C19.status', context), 'unknown');
console.log('PASS: full UI script parses; HG001 renders; without CPIC matches genes stay unknown and verdict is yellow; with matches genes are called only from complete evidence; status fields in a file cannot promote a gene.');
console.log('HG001 limited-panel candidates:', Object.entries(map.genes).map(([gene, def]) => {
 context.def = def; context.gene = gene;
 return [gene, vm.runInContext('callGene(gene, def, Object.fromEntries(data.evidence.filter(s=>s.genotype).flatMap(s=>(s.rsids||[]).map(rs=>[rs,s.genotype]))))',context)];
}));
