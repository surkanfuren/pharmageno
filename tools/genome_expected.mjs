// Runs every test genome in genomes/manifest.json through the real upload endpoint and the real
// callGene() from index.html, then records file sizes, SHA-256 and the expected passport in the manifest.
// usage: node tools/genome_expected.mjs [server]   (default http://127.0.0.1:8020; genomes/data must be downloaded)
import fs from "fs";
import crypto from "crypto";

const server = process.argv[2] || "http://127.0.0.1:8020";
const root = new URL("..", import.meta.url);
const html = fs.readFileSync(new URL("index.html", root), "utf8");
const grab = name => { const i = html.indexOf(`function ${name}(`); let d = 0, j = html.indexOf("{", i);
  for (let k = j; k < html.length; k++) { if (html[k] === "{") d++; else if (html[k] === "}" && --d === 0) return html.slice(i, k + 1); } };
const callGene = eval("(" + grab("callGene") + ")");
const map = JSON.parse(fs.readFileSync(new URL("allele_map.json", root), "utf8"));
const manifestPath = new URL("genomes/manifest.json", root);
const manifest = JSON.parse(fs.readFileSync(manifestPath, "utf8"));

async function phenotype(gene, diplotype) {
  const ask = async d => (await (await fetch(`https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.${gene}&diplotype=eq.${encodeURIComponent(d)}&select=generesult`)).json())[0]?.generesult;
  const [a, b] = diplotype.split("/");
  return await ask(diplotype) || await ask(`${b}/${a}`) || "no CPIC match";
}

for (const person of manifest.genomes) {
  const buffers = person.files.map(file => {
    const data = fs.readFileSync(new URL(`genomes/data/${person.id}/${file.name}`, root));
    file.bytes = data.length;
    file.sha256 = crypto.createHash("sha256").update(data).digest("hex");
    return data;
  });
  const t0 = Date.now();
  const response = await fetch(`${server}/api/genome`, { method: "POST", headers: { "X-File-Sizes": buffers.map(b => b.length).join(",") }, body: Buffer.concat(buffers) });
  const doc = await response.json();
  if (!response.ok) throw new Error(`${person.id}: ${doc.error}`);
  const geno = {};
  for (const site of doc.evidence) if (["observed", "reference_in_benchmark"].includes(site.status) && /^[ACGT]{2}$/.test(site.genotype || "")) for (const rs of site.rsids || []) geno[rs] = site.genotype;
  person.expected = {};
  for (const [gene, def] of Object.entries(map.genes)) {
    const call = callGene(gene, def, geno);
    person.expected[gene] = call.diplotype ? `${call.diplotype} (${await phenotype(gene, call.diplotype)})` : "unknown";
  }
  person.variants = doc.records_scanned;
  console.log(person.id.padEnd(8), `${doc.records_scanned} variants, ${((Date.now() - t0) / 1000).toFixed(1)} s`, JSON.stringify(person.expected));
}
fs.writeFileSync(manifestPath, JSON.stringify(manifest, null, 2) + "\n");
console.log("manifest updated");
