// Runs the actual single-file UI script with deterministic resource API responses.
import fs from 'node:fs';
import vm from 'node:vm';
import assert from 'node:assert/strict';
const source = fs.readFileSync(new URL('../index.html', import.meta.url), 'utf8').match(/<script>([\s\S]*?)<\/script>/)[1];
const elements = new Map();
const element = id => {
  if (!elements.has(id)) elements.set(id, { innerHTML:'', textContent:'', disabled:false, value:'Plavix',
    handlers:{}, dataset:{}, style:{}, insertAdjacentHTML(){}, scrollIntoView(){}, focus(){}, classList:{add(){},remove(){},toggle(){}}, addEventListener(event,callback){ this.handlers[event]=callback; } });
  return elements.get(id);
};
const A = 'a'.repeat(64), B = 'b'.repeat(64);
const manifest = id => ({snapshot_id:id, retrieved_at:'2026-10-12T08:00:00Z', source:'https://api.cpicpgx.org/v1',
  table_hashes:{diplotype:id}, release:{tag:'vTEST',verified:true}, consistency:'two matching reads',validation:{checks:['schema']}});
let state = {freshness:'current_as_of_check',checked_at:'2026-10-12T09:00:00Z', active:manifest(A),candidate:null,latest_release:{tag:'vTEST'}};
let mappingChanged = false;
const calls = [];
const context = vm.createContext({document:{getElementById:element,querySelectorAll:()=>[]}, performance,URL,URLSearchParams,Blob,setTimeout,
  fetch:async (path,options={})=>{
    calls.push(path);
    const url = new URL(path,'http://localhost');
    const table = url.pathname.split('/').at(-1);
    let data;
    if (table==='status') data=state;
    else if (table==='passport-impact') {
      const passport = JSON.parse(options.body);
      data = passport.cpic_provenance?.snapshot_id===state.active.snapshot_id ? {status:'same_snapshot'} : {
        status:'changed_snapshot',genome_reanalysis_required:false,gene_result_review_required:mappingChanged,
        drug_recheck_required:true,passport_affected_genes:['CYP2C19'],new_genes:['NEWGENE']};
    } else if (table==='manifest') data=manifest(url.searchParams.get('snapshot'));
    else if (table==='diplotype') {
      const gene=url.searchParams.get('genesymbol').slice(3);
      data=[{diplotype:url.searchParams.get('diplotype').slice(3),generesult:'Normal Metabolizer',lookupkey:{[gene]:'Normal Metabolizer'}}];
    } else if (table==='pair_view') data=[{genesymbol:'CYP2C19',cpiclevel:'A'}];
    else if (table==='recommendation_view') data=[{drugrecommendation:'Use standard dose.',phenotypes:{CYP2C19:'Normal Metabolizer'}}];
    else if (table==='activate') {
      assert.equal(JSON.parse(options.body).snapshot_id,B);
      state={...state,active:manifest(B),candidate:null,freshness:'current_as_of_check'};
      data=state;
    } else throw Error('Unexpected request: '+path);
    return {ok:true,status:200,json:async()=>structuredClone(data)};
  }
});
vm.runInContext(source,context);
// Build a verified passport through the real CPIC lookup path (what the genome import does per gene)
await vm.runInContext(`withResources(async () => {
  passport = {};
  for (const [gene, dip] of [["CYP2C19", "*2/*2"], ["SLCO1B1", "*1/*1"]]) passport[gene] = { ...(await lookupDiplotype(gene, dip)), source: "genome", status: "ok" };
  passportProvenance = currentProvenance();
  renderPassport();
})`,context);
assert.equal(vm.runInContext('passportProvenance.snapshot_id',context),A);
assert.match(element('passportVersion').textContent,/aaaaaaaaaaaa/);
const normal = await vm.runInContext('withResources(async () => {await reviewPassportVersion(); return checkDrug("Plavix", "clopidogrel");})',context);
assert.equal(normal.color,'green');
assert(calls.filter(path=>/\/(diplotype|pair_view|recommendation_view)\?/.test(path)).every(path=>new URL(path,'http://localhost').searchParams.get('snapshot')===A));
// A failed live check must remain explicit even with an available saved snapshot.
state={...state,freshness:'unverified',error:'offline'};
await vm.runInContext('prepareResources()',context);
assert.match(element('resourceStatus').textContent,/Could not verify/);
assert.doesNotMatch(element('resourceStatus').textContent,/matched the live source/);
// Reviewing and activating a mapping update invalidates affected old clinical calls.
state={...state,freshness:'update_available',error:null,candidate:manifest(B),impact:{changed_tables:['diplotype'],affected_genes:['CYP2C19'],affected_drugs:['clopidogrel'],gene_result_review_required:true,genome_reanalysis_required:false,drug_recheck_required:true}};
mappingChanged=true;
await vm.runInContext('refreshResources()',context);
assert.equal(element('resourceActivate').disabled,false);
await element('resourceActivate').handlers.click();
assert.equal(vm.runInContext('passport.CYP2C19.status',context),'unknown');
assert.equal(vm.runInContext('passportProvenance.snapshot_id',context),A); // history is never rewritten
assert.equal(vm.runInContext('cpicContext.snapshot_id',context),B);
assert.match(element('passportImpact').textContent,/New guideline genes/);
const uncertain=await vm.runInContext('checkDrug("Plavix","clopidogrel")',context);
assert.equal(uncertain.color,'yellow');
assert.deepEqual(Array.from(uncertain.missing),['CYP2C19']);
console.log('PASS: UI pins requests and provenance, reports offline freshness, reviews updates, preserves history, and blocks affected old calls from green verdicts.');
