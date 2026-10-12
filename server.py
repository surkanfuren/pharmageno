"""
Gene Passport: small server (Python standard library only)

- Serves the static files in this folder (index.html).
- POST /api/llm  {"prompt": "..."}  ->  {"text": "...", "backend": "...", "ms": 1234}
  Backend order:
    1) OPENAI_API_KEY set: OpenAI Responses API (model: OPENAI_MODEL, else DEFAULT_MODEL)
    2) otherwise a local `codex exec`, if installed and signed in
    3) neither: a clear error that the page shows
  The key is never printed or logged.

Run:  python3 server.py   ->  http://localhost:8000
"""
import bisect
import gzip
import hashlib
import re
import secrets
import sqlite3
import urllib.parse
from contextlib import closing
from datetime import datetime, timezone
import json
import os
import shutil
import subprocess
import tempfile
import zipfile
import time
import threading
from cpic_resources import ResourceStore, CALLER_VERSION
import urllib.error
import urllib.request
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

PORT = int(os.environ.get("PORT", "8000"))
HOST = os.environ.get("HOST", "127.0.0.1")  # 0.0.0.0 = reachable from the network
DEFAULT_MODEL = "gpt-6-luna"  # fast, low-cost tier; override with OPENAI_MODEL
CODEX_TIMEOUT = 120
HERE = os.path.dirname(os.path.abspath(__file__))
RESOURCES = ResourceStore(os.environ.get("CPIC_DATA_DIR", os.path.join(HERE, ".cpic-data")))
SHARE_DB = os.environ.get("SHARE_DB", os.path.join(HERE, ".data", "passports.db"))
SHARE_FIELDS = ("diplotype", "result", "lookup", "source", "status", "note", "reason")


def share_db():
    os.makedirs(os.path.dirname(SHARE_DB), exist_ok=True)
    db = sqlite3.connect(SHARE_DB, timeout=10)
    db.execute("CREATE TABLE IF NOT EXISTS passports (id TEXT PRIMARY KEY, created_at TEXT NOT NULL, data TEXT NOT NULL)")
    return db


def clean_passport(value):
    # Keep only what the UI needs. Shared results are re-verified against CPIC when the link is opened.
    genes = value.get("genes")
    if not isinstance(genes, dict) or not 0 < len(genes) <= 64:
        raise ValueError("A passport needs 1-64 genes")
    out = {}
    for gene, entry in genes.items():
        if not re.fullmatch(r"[A-Za-z0-9-]{2,20}", str(gene)) or not isinstance(entry, dict):
            raise ValueError("Invalid gene entry")
        item = {k: str(entry[k])[:500] for k in SHARE_FIELDS if entry.get(k) not in (None, "")}
        if item.get("status") not in ("ok", "unknown") or item.get("source") not in ("lab", "consumer", "genome"):
            raise ValueError("Invalid gene status or source")
        out[gene] = item
    identity = value.get("identity") if isinstance(value.get("identity"), dict) else {}
    provenance = value.get("cpic_provenance")
    data = {"schema": "gene-passport-share-v1", "genes": out,
            "identity": {k: str(identity.get(k, ""))[:120] for k in ("name", "badge")},
            "cpic_provenance": provenance if isinstance(provenance, dict) else None}
    blob = json.dumps(data, ensure_ascii=False)
    if len(blob) > 256 * 1024:
        raise ValueError("Passport is too large to share")
    return blob


def via_openai(prompt):
    body = json.dumps({"model": os.environ.get("OPENAI_MODEL", DEFAULT_MODEL), "input": prompt}).encode()
    req = urllib.request.Request("https://api.openai.com/v1/responses", data=body, headers={
        "Authorization": "Bearer " + os.environ["OPENAI_API_KEY"], "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        d = json.load(r)
    if d.get("output_text"):
        return d["output_text"]
    # without the output_text convenience field, join the message contents
    return "".join(c.get("text", "") for o in d.get("output", []) if o.get("type") == "message"
                   for c in o.get("content", []))


def via_codex(prompt):
    # Run read-only in an empty folder, leaving nothing behind; read only the final message from a file.
    # Fast model (~7 s); retry once on transient errors such as "model at capacity".
    for attempt in range(2):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "out.txt")
            try:
                subprocess.run(["codex", "exec", "-s", "read-only", "--skip-git-repo-check", "--ephemeral",
                                "--ignore-user-config", "-m", os.environ.get("CODEX_MODEL", DEFAULT_MODEL),
                                "-C", tmp, "-o", out, "-"],
                               input=prompt, text=True, capture_output=True, timeout=CODEX_TIMEOUT, check=True)
                with open(out, encoding="utf-8") as f:
                    return f.read()
            except subprocess.CalledProcessError:
                if attempt == 1:
                    raise


CACHE = {}  # answer a repeated question instantly


def ask_llm(prompt):
    if prompt in CACHE:
        return CACHE[prompt], CACHE["_backend"] + " (cache)"
    if os.environ.get("OPENAI_API_KEY"):
        text, backend = via_openai(prompt), "openai"
    elif shutil.which("codex"):
        text, backend = via_codex(prompt), "codex"
    else:
        raise RuntimeError("No AI backend: OPENAI_API_KEY is not set and codex was not found.")
    CACHE[prompt], CACHE["_backend"] = text, backend
    return text, backend



def genome_catalog(snapshot=None):
    if snapshot is None:
        snapshot, _ = RESOURCES.ensure()
    tables = snapshot["tables"]
    pairs = [{"genesymbol": r["genesymbol"], "drugname": r["drugname"]}
             for r in tables["pair_view"] if r["cpiclevel"] in ("A", "B", "A/B")]
    genes = {r["genesymbol"] for r in pairs}
    chromosomes = {r["symbol"]: r["chr"] for r in tables["gene"]}
    locations = tables["sequence_location"]
    targets, unsupported = {}, {g: 0 for g in genes}
    for row in locations:
        gene = row["genesymbol"]
        if gene not in genes:
            continue
        # A SNV may list several alternates at one position ("g.18130687T>C; g.18130687T>G", TPMT*3C)
        parts = [re.fullmatch(r"g\.(\d+)([ACGT])>([ACGT])", part.strip()) for part in (row.get("chromosomelocation") or "").split(";")]
        match = parts[0] if all(parts) and len({(m[1], m[2]) for m in parts}) == 1 else None
        chromosome = chromosomes.get(gene)
        if not match or not chromosome:
            unsupported[gene] += 1
            continue
        key = (chromosome.removeprefix("chr"), int(match[1]))
        item = targets.setdefault(key, {"chrom": key[0], "position": key[1], "ref": match[2], "genes": [], "rsids": [], "versions": []})
        if item["ref"] != match[2]:
            raise ValueError("Conflicting CPIC reference bases")
        for field, value in [("genes", gene), ("rsids", row.get("dbsnpid")), ("versions", row["version"])]:
            if value is not None and value not in item[field]:
                item[field].append(value)
    return pairs, targets, unsupported


def open_genome(path):
    with open(path, "rb") as stream:
        compressed = stream.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if compressed else open(path, encoding="utf-8-sig")


def sniff_genome_file(path):
    """Classify an uploaded file by its content (not its name): vcf, bed, zip or other."""
    with open(path, "rb") as stream:
        if stream.read(2) == b"PK":
            return "zip"
    try:
        with open_genome(path) as stream:
            for line in stream:
                if not line.strip() or line.startswith(("track", "browser")):
                    continue
                if line.startswith("##fileformat=VCF"):
                    return "vcf"
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                return "bed" if len(fields) >= 3 and fields[1].isdigit() and fields[2].isdigit() else "other"
    except (UnicodeDecodeError, OSError, EOFError, gzip.BadGzipFile):
        return "other"
    return "other"


def split_genome_upload(paths, folder):
    """One upload may hold a VCF, an optional confident-region BED, and/or a zip with both. Returns (vcf, bed, ignored)."""
    found, ignored = {"vcf": [], "bed": []}, 0
    for path in paths:
        kind = sniff_genome_file(path)
        if kind == "zip":
            with zipfile.ZipFile(path) as archive:
                members = [m for m in archive.infolist() if not m.is_dir() and not m.filename.startswith("__MACOSX")]
                if sum(m.file_size for m in members) > 4 * 1024 ** 3:
                    raise ValueError("The zip file is too large once unpacked")
                for index, member in enumerate(members):
                    target = os.path.join(folder, f"zip{index}")
                    with archive.open(member) as source, open(target, "wb") as output:
                        shutil.copyfileobj(source, output, 1024 * 1024)
                    kind = sniff_genome_file(target)
                    if kind in found:
                        found[kind].append(target)
                    else:
                        ignored += 1
        elif kind in found:
            found[kind].append(path)
        else:
            ignored += 1
    if len(found["vcf"]) != 1:
        raise ValueError("Add exactly one genome file (.vcf or .vcf.gz)" if not found["vcf"] else "Add one person's genome at a time (found several VCF files)")
    if len(found["bed"]) > 1:
        raise ValueError("Add at most one confident-region (.bed) file")
    return found["vcf"][0], (found["bed"] or [None])[0], ignored


def analyze_genome(vcf_path, bed_path=None, catalog=None, trusted_bed=False):
    provenance = None
    if catalog is None:
        snapshot, freshness = RESOURCES.ensure()
        catalog = genome_catalog(snapshot)
        provenance = RESOURCES.provenance(snapshot, freshness)
        with open(os.path.join(HERE, "allele_map.json"), "rb") as panel_file:
            provenance["limited_panel_sha256"] = hashlib.sha256(panel_file.read()).hexdigest()
    pairs, targets, unsupported = catalog
    evidence = {key: dict(item, status="missing", genotype=None) for key, item in targets.items()}
    by_chrom = {}
    for chromosome, position in targets:
        by_chrom.setdefault(chromosome, []).append(position)
    for positions in by_chrom.values():
        positions.sort()
    digest = hashlib.sha256()
    with open(vcf_path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    sample, records, build_seen, chromosomes_seen = None, 0, False, set()
    with open_genome(vcf_path) as stream:
        for line in stream:
            if line.startswith("##"):
                if "GRCh38" in line or "hg38" in line:
                    build_seen = True
                if "GRCh37" in line or "hg19" in line:
                    raise ValueError("GRCh37 is not supported. Supply a GRCh38 VCF; coordinates are never guessed.")
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                if len(fields) != 10:
                    raise ValueError("Supply a single-sample VCF (exactly one sample column).")
                sample = fields[9]
                continue
            if line.startswith("#") or not line.strip():
                continue
            if sample is None:
                raise ValueError("Missing VCF sample header")
            f = line.rstrip().split("\t")
            if len(f) != 10:
                raise ValueError("Malformed VCF record")
            records += 1
            chrom, start, ref = f[0].removeprefix("chr"), int(f[1]), f[3]
            chromosomes_seen.add(chrom)
            positions = by_chrom.get(chrom, [])
            # A deletion or complex record overlapping a target blocks reference inference.
            end = start + len(ref) - 1
            end_match = re.search(r"(?:^|;)END=(\d+)(?:;|$)", f[7])
            if end_match:
                end = max(end, int(end_match[1]))
            for index in range(bisect.bisect_left(positions, start), bisect.bisect_right(positions, end)):
                key = (chrom, positions[index])
                item = evidence[key]
                if item["status"] != "missing":
                    item.update(status="conflicting_records", genotype=None)
                    continue
                item.update(status="unsupported_record", genotype=None, vcf_ref=ref, vcf_alt=f[4], filter=f[6])
                if start != key[1] or len(ref) != 1 or ref != item["ref"] or f[6] != "PASS":
                    continue
                values = dict(zip(f[8].split(":"), f[9].split(":")))
                gt = values.get("GT", "")
                indices = re.split(r"[/|]", gt)
                alleles = [ref] + f[4].split(",")
                item["gt"] = gt
                item["phase_set"] = values.get("PS")
                if len(indices) != 2 or any(not i.isdigit() or int(i) >= len(alleles) for i in indices):
                    item["status"] = "no_call"
                    continue
                bases = [alleles[int(i)] for i in indices]
                if any(b not in "ACGT" or len(b) != 1 for b in bases):
                    continue
                item.update(status="observed", genotype="".join(bases), phased="|" in gt)
    if sample is None or not build_seen:
        raise ValueError("VCF must declare GRCh38 in its reference or contig headers.")
    bed_digest = None
    if bed_path and trusted_bed:
        with open(bed_path, "rb") as bed_stream:
            bed_digest = hashlib.sha256(bed_stream.read()).hexdigest()
        with open_genome(bed_path) as stream:
            for line in stream:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.split()
                chrom, start, end = f[0].removeprefix("chr"), int(f[1]), int(f[2])
                if start < 0 or end <= start:
                    raise ValueError("Invalid BED interval")
                positions = by_chrom.get(chrom, [])
                for index in range(bisect.bisect_right(positions, start), bisect.bisect_right(positions, end)):
                    item = evidence[(chrom, positions[index])]
                    if item["status"] == "missing":
                        item.update(status="reference_in_benchmark", genotype=item["ref"] * 2)
    genes = {}
    for gene in sorted({r["genesymbol"] for r in pairs}):
        sites = [item for item in evidence.values() if gene in item["genes"]]
        known = sum(item["genotype"] is not None for item in sites)
        genes[gene] = {"status": "unknown", "source": "genome", "diplotype": "", "note":
            f"{known}/{len(sites)} supported positions read; {unsupported[gene]} complex definitions require additional analysis. A validated gene call is needed before this result can guide drug checks.",
            "positions_read": known, "positions_total": len(sites), "unsupported_definitions": unsupported[gene],
            "chromosome": sites[0]["chrom"] if sites else None}
    return {"schema": "gene-passport-genome-v1", "sample": sample, "assembly": "GRCh38",
            "chromosomes": sorted(chromosomes_seen, key=lambda c: (not c.isdigit(), int(c) if c.isdigit() else 0, c)),
            "created_at": datetime.now(timezone.utc).isoformat(), "source_sha256": digest.hexdigest(),
            "bed_sha256": bed_digest, "records_scanned": records, "genes": genes,
            "evidence": list(evidence.values()), "cpic_pairs": pairs,
            "catalog_sha256": hashlib.sha256(json.dumps([list(targets.values()), unsupported, pairs], sort_keys=True).encode()).hexdigest(),
            "cpic_source": "https://api.cpicpgx.org/v1", "cpic_provenance": provenance, "caller_version": CALLER_VERSION, "clinical_ready": False,
            "limitations": "Research evidence passport, not a clinical gene call. SNP evidence does not resolve rare alleles, phase, copy number, HLA typing, or all complex variants. Keep the original genome and coverage files for reanalysis when CPIC definitions change."}


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=HERE, **kw)

    def json_response(self, result, code=200):
        data = json.dumps(result, ensure_ascii=False).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def request_json(self):
        size = int(self.headers.get("Content-Length", "0"))
        if not 0 < size <= 2 * 1024 * 1024:
            raise ValueError("Invalid resource request size")
        value = json.loads(self.rfile.read(size))
        if not isinstance(value, dict):
            raise ValueError("Expected a JSON object")
        return value

    def do_GET(self):
        parsed = urllib.parse.urlsplit(self.path)
        if parsed.path.startswith("/api/cpic/"):
            try:
                params = dict(urllib.parse.parse_qsl(parsed.query))
                if parsed.path == "/api/cpic/status":
                    return self.json_response(RESOURCES.status())
                snapshot_id = params.pop("snapshot", None)
                if parsed.path == "/api/cpic/manifest":
                    return self.json_response(RESOURCES.load(snapshot_id)["manifest"])
                table = parsed.path.removeprefix("/api/cpic/")
                return self.json_response(RESOURCES.query(snapshot_id, table, params))
            except Exception as error:
                return self.json_response({"error": str(error)}, 400)
        if parsed.path.startswith("/api/passports/"):
            share_id = parsed.path.removeprefix("/api/passports/")
            if not re.fullmatch(r"[A-Za-z0-9_-]{6,32}", share_id):
                return self.json_response({"error": "Shared passport not found"}, 404)
            with closing(share_db()) as db:
                row = db.execute("SELECT created_at, data FROM passports WHERE id = ?", (share_id,)).fetchone()
            if not row:
                return self.json_response({"error": "Shared passport not found"}, 404)
            return self.json_response({**json.loads(row[1]), "id": share_id, "shared_at": row[0]})
        # Never serve hidden folders (.cpic-data, .data with the share database, .env)
        if parsed.path.startswith("/genomes/data/") or any(part.startswith(".") for part in urllib.parse.unquote(parsed.path).split("/")):
            return self.send_error(404)
        return super().do_GET()

    def do_POST(self):
        if self.path == "/api/passports":
            try:
                blob = clean_passport(self.request_json())
                share_id, created = secrets.token_urlsafe(9), datetime.now(timezone.utc).isoformat()
                with closing(share_db()) as db, db:
                    db.execute("INSERT INTO passports VALUES (?, ?, ?)", (share_id, created, blob))
                return self.json_response({"id": share_id, "shared_at": created})
            except Exception as error:
                return self.json_response({"error": str(error)}, 400)
        if self.path.startswith("/api/cpic/"):
            try:
                if self.path == "/api/cpic/check":
                    result = RESOURCES.check(force=True)
                elif self.path == "/api/cpic/activate":
                    result = RESOURCES.activate(self.request_json().get("snapshot_id"))
                elif self.path == "/api/cpic/passport-impact":
                    result = RESOURCES.passport_impact(self.request_json())
                else:
                    return self.send_error(404)
                return self.json_response(result)
            except Exception as error:
                return self.json_response({"error": str(error)}, 400)
        if self.path == "/api/genome":
            try:
                size = int(self.headers.get("Content-Length", "0"))
                if not 0 < size <= 1024 ** 3:
                    raise ValueError("Supply a VCF up to 1 GB")
                # Body = the dropped files back to back; X-File-Sizes says where each one ends
                sizes = [int(x) for x in (self.headers.get("X-File-Sizes") or str(size)).split(",")]
                if len(sizes) > 20 or any(x <= 0 for x in sizes) or sum(sizes) != size:
                    raise ValueError("Invalid upload")
                with tempfile.TemporaryDirectory() as folder:
                    paths = []
                    for index, length in enumerate(sizes):
                        paths.append(os.path.join(folder, f"part{index}"))
                        with open(paths[-1], "wb") as output:
                            while length:
                                chunk = self.rfile.read(min(1024 * 1024, length))
                                if not chunk:
                                    raise ValueError("Incomplete upload")
                                output.write(chunk)
                                length -= len(chunk)
                    vcf_path, bed_path, ignored = split_genome_upload(paths, folder)
                    result = analyze_genome(vcf_path, bed_path, trusted_bed=bool(bed_path))
                    result["upload"] = {"files": len(sizes), "confident_regions": bool(bed_path), "ignored_files": ignored}
                data, code = json.dumps(result).encode(), 200
            except Exception as error:
                data, code = json.dumps({"error": str(error)}).encode(), 400
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)
            return
        if self.path != "/api/llm":
            return self.send_error(404)
        t0 = time.time()
        try:
            prompt = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))))["prompt"]
            text, backend = ask_llm(prompt)
            code, res = 200, {"text": text, "backend": backend}
        except subprocess.TimeoutExpired:
            code, res = 504, {"error": "The AI did not answer in time."}
        except subprocess.CalledProcessError as e:
            code, res = 502, {"error": "codex error: " + (e.stderr or "")[-300:]}
        except urllib.error.HTTPError as e:
            code, res = 502, {"error": f"OpenAI API error: HTTP {e.code}"}
        except Exception as e:  # noqa: BLE001 (return a readable error to the page)
            code, res = 500, {"error": str(e)}
        res["ms"] = int((time.time() - t0) * 1000)
        data = json.dumps(res, ensure_ascii=False).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, fmt, *args):
        # Path and status only; bodies and headers (keys) are never logged
        print(f"{self.command} {self.path} {args[1] if len(args) > 1 else ''}")


def resource_schedule():
    # Daemon ends with this server; no OS-level scheduled task is installed.
    while True:
        try:
            RESOURCES.check()
        except Exception as error:
            print("CPIC resource check unavailable:", str(error))
        time.sleep(60)


if __name__ == "__main__":
    threading.Thread(target=resource_schedule, daemon=True).start()
    print(f"Gene Passport: http://{HOST}:{PORT}  (LLM: "
          f"{'OpenAI' if os.environ.get('OPENAI_API_KEY') else 'codex' if shutil.which('codex') else 'NONE'})")
    ThreadingHTTPServer((HOST, PORT), Handler).serve_forever()
