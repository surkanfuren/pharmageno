// Runs the REAL callGene() from index.html (not a copy) on 2,504 1000 Genomes people
import fs from "fs";
const html = fs.readFileSync(process.argv[2], "utf8");
const grab = name => { const i = html.indexOf(`function ${name}(`); let d = 0, j = html.indexOf("{", i);
  for (let k = j; k < html.length; k++) { if (html[k] === "{") d++; else if (html[k] === "}" && --d === 0) return html.slice(i, k + 1); } };
const callGene = eval("(" + grab("callGene") + ")");
const map = JSON.parse(fs.readFileSync(process.argv[3], "utf8"));
const geno = JSON.parse(fs.readFileSync(process.argv[4], "utf8"));
const out = {}, stats = {};
for (const [sample, g] of Object.entries(geno)) {
  // Same genotype map the genome path builds from VCF evidence (lower-case rsID -> two letters)
  const parsed = Object.fromEntries(Object.entries(g).map(([rs, gt]) => [rs.toLowerCase(), String(gt).toUpperCase()]));
  out[sample] = {};
  for (const [gene, def] of Object.entries(map.genes)) {
    const c = callGene(gene, def, parsed);
    const v = c.unknown ? "UNKNOWN: " + c.unknown : c.diplotype;
    out[sample][gene] = v;
    const k = c.unknown ? "UNKNOWN" : c.diplotype;
    stats[gene] ??= {}; stats[gene][k] = (stats[gene][k] || 0) + 1;
  }
}
fs.writeFileSync("data/calls_1000g.json", JSON.stringify(out));
for (const [gene, s] of Object.entries(stats)) console.log(gene, JSON.stringify(Object.entries(s).sort((a, b) => b[1] - a[1]).slice(0, 8)));
const unk = Object.values(out).flatMap(o => Object.values(o)).filter(v => v.startsWith("UNKNOWN"));
const reasons = {}; unk.forEach(u => { const r = u.replace(/rs\d+/g, "rsX").replace(/\(.*?\)/g, ""); reasons[r] = (reasons[r] || 0) + 1; });
console.log("unknown reasons:", reasons);
for (const s of ["NA12878", "NA18519", "NA19143", "HG00436", "NA18484", "NA10831"]) if (out[s]) console.log(s, JSON.stringify(out[s]));
