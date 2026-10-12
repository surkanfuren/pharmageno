"""
Builds the rsID -> allele map used by the genome reader, from CPIC.

SOURCE: CPIC API (https://api.cpicpgx.org/v1) tables:
  allele_definition  -> allele_location_value -> sequence_location (dbsnpid, chromosomelocation GRCh38)
rsIDs are never typed by hand; this file only chooses which alleles are supported (the main clinically tested ones).

An allele's "core" SNPs: positions where the allele has a single base (A/C/G/T) that differs from the reference allele (*1 / Reference).
IUPAC ambiguity codes (R, Y...) mean "either base" and do not tell alleles apart.

Run:  python3 build_allele_map.py   ->  allele_map.json
"""
import json
import re
import urllib.parse
import urllib.request
from datetime import date

API = "https://api.cpicpgx.org/v1"

# Supported genes: reference allele name + main alleles defined by one or a few SNPs
GENES = {
    "CYP2C19": ("*1", ["*2", "*3", "*17"]),
    "CYP2C9": ("*1", ["*2", "*3"]),
    "SLCO1B1": ("*1", ["*5"]),
    "NUDT15": ("*1", ["*3"]),
    "TPMT": ("*1", ["*3A", "*2", "*3B", "*3C"]),  # multi-SNP *3A first: it "consumes" the same SNPs before *3B and *3C
    "DPYD": ("Reference", ["c.1905+1G>A (*2A)", "c.1679T>G (*13)", "c.2846A>T", "c.1129-5923C>G, c.1236G>A (HapB3)"]),
}


def get(path, params):
    url = f"{API}/{path}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=60) as r:
        return json.load(r)


def alleles(gene, names):
    rows = get("allele_definition", {
        "genesymbol": f"eq.{gene}",
        "name": "in.(" + ",".join(f'"{n}"' for n in names) + ")",
        "select": "name,allele_location_value(variantallele,sequence_location(dbsnpid,chromosomelocation))",
    })
    out = {}
    for a in rows:
        out[a["name"]] = {
            v["sequence_location"]["dbsnpid"]: (v["variantallele"], v["sequence_location"]["chromosomelocation"])
            for v in a["allele_location_value"] if v["sequence_location"].get("dbsnpid")
        }
    return out


def ref_base(chromosomelocation):
    # "g.94781859G>A" -> "G" (GRCh38 reference base)
    m = re.search(r"\d([ACGT])>[ACGT]$", chromosomelocation)
    return m.group(1) if m else None


def build():
    result = {"_source": f"CPIC API {API} allele_definition/allele_location_value/sequence_location",
              "_generated": date.today().isoformat(), "genes": {}}
    for gene, (ref_name, targets) in GENES.items():
        defs = alleles(gene, [ref_name] + targets)
        ref = defs.get(ref_name, {})
        gene_out = {"reference": ref_name, "alleles": {}, "rsids": {}}
        for name in targets:
            core = {}
            for rsid, (base, loc) in defs.get(name, {}).items():
                if base not in "ACGT" or len(base) != 1:
                    continue  # IUPAC ambiguity or indel: not distinctive / not readable as a single base
                ref_here = ref.get(rsid, (ref_base(loc), loc))[0]
                if base != ref_here:
                    core[rsid] = base
                    gene_out["rsids"][rsid] = ref_here
            gene_out["alleles"][name] = core
        result["genes"][gene] = gene_out
    return result


if __name__ == "__main__":
    m = build()
    with open("allele_map.json", "w") as f:
        json.dump(m, f, indent=1)
    for g, d in m["genes"].items():
        print(g, {a: c for a, c in d["alleles"].items()}, "| referans:", d["rsids"])
