"""Download the 10 public test genomes listed in manifest.json into genomes/data/.

Usage:  python3 genomes/download.py            (all 10 people, ~1.5 GB)
        python3 genomes/download.py HG003      (one person)
        python3 genomes/download.py --zip      (also pack each person into genomes/zips/<ID>.zip)

A person's VCF and confident-region BED belong together: the VCF lists where the person differs
from the reference genome, the BED says where the sequencing was read reliably. The zip holds both,
so it can be dropped on the site as a single file.

Standard library only. Every file is checked against the SHA-256 in manifest.json.
NA12879 and NA12881 are not published as single-person files: they are split out of
the public Platinum Pedigree family VCF (only that person's variant rows are kept and a
"##reference=GRCh38" line is added), so the result is reproducible from the public source.
"""
import gzip
import hashlib
import json
import os
import shutil
import sys
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
ZIPS = os.path.join(HERE, "zips")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url, path):
    if os.path.exists(path):
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"  downloading {url}")
    with urllib.request.urlopen(url, timeout=120) as response, open(path + ".part", "wb") as output:
        shutil.copyfileobj(response, output, 1024 * 1024)
    os.replace(path + ".part", path)


def split_sample(family_vcf, sample, path):
    """Write one person's variant rows from a multi-sample VCF as a single-sample VCF."""
    with gzip.open(family_vcf, "rt") as source, gzip.open(path + ".part", "wt", compresslevel=6) as output:
        column = None
        for line in source:
            if line.startswith("##"):
                output.write(line)
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                column = fields.index(sample)
                output.write("##reference=GRCh38\n")
                output.write(f"##source=Gene Passport genomes/download.py: sample {sample} split from {os.path.basename(family_vcf)}\n")
                output.write("\t".join(fields[:9] + [sample]) + "\n")
                continue
            fields = line.rstrip("\n").split("\t")
            genotype = fields[column].split(":")[0]
            if any(allele not in ("0", ".") for allele in genotype.replace("|", "/").split("/")):
                output.write("\t".join(fields[:9] + [fields[column]]) + "\n")
    os.replace(path + ".part", path)


def pack(person):
    """One zip per person: the already-compressed VCF is stored as is, the BED is deflated."""
    os.makedirs(ZIPS, exist_ok=True)
    path = os.path.join(ZIPS, f"{person['id']}.zip")
    with zipfile.ZipFile(path + ".part", "w") as archive:
        for item in person["files"]:
            source = os.path.join(DATA, person["id"], item["name"])
            compress = zipfile.ZIP_STORED if item["name"].endswith(".gz") else zipfile.ZIP_DEFLATED
            archive.write(source, f"{person['id']}/{item['name']}", compress_type=compress)
    os.replace(path + ".part", path)
    print(f"  zip {os.path.relpath(path, os.path.dirname(HERE))} ({os.path.getsize(path) // 1048576} MB)")


def main(wanted, make_zip=False):
    manifest = json.load(open(os.path.join(HERE, "manifest.json"), encoding="utf-8"))
    for person in manifest["genomes"]:
        if wanted and person["id"] not in wanted:
            continue
        print(person["id"], "-", person["label"])
        for item in person["files"]:
            path = os.path.join(DATA, person["id"], item["name"])
            if "url" in item:
                fetch(item["url"], path)
            else:
                source = item["split_from"]
                family = os.path.join(DATA, "_sources", os.path.basename(source))
                fetch(source, family)
                if not os.path.exists(path):
                    print(f"  splitting {item['sample']} from {os.path.basename(family)}")
                    os.makedirs(os.path.dirname(path), exist_ok=True)
                    split_sample(family, item["sample"], path)
            if item.get("sha256") and sha256(path) != item["sha256"]:
                sys.exit(f"  CHECKSUM MISMATCH: {path}")
            print(f"  ok  {item['name']}")
        if make_zip:
            pack(person)


if __name__ == "__main__":
    args = sys.argv[1:]
    main({a for a in args if not a.startswith("--")}, make_zip="--zip" in args)
