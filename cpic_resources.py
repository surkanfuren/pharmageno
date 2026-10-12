"""Versioned public CPIC API snapshots and update impact analysis (stdlib only).

A local snapshot is reproducible. Two matching API reads detect observable drift,
but do not provide the transactional guarantee of a CPIC database export.
No patient genomes or report text are stored here.
"""
from contextlib import contextmanager
import fcntl
import concurrent.futures
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import threading
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone

API = 'https://api.cpicpgx.org/v1'
RELEASE_URL = 'https://api.github.com/repos/cpicpgx/cpic-data/releases/latest'
CHECK_INTERVAL = 24 * 60 * 60
CALLER_VERSION = 'genome-evidence-1.1'
TABLES = {
    'pair_view': ('pairid', 'drugname', 'genesymbol', 'cpiclevel'),
    'gene': ('symbol', 'chr', 'lookupmethod'),
    'sequence_location': ('id', 'genesymbol', 'chromosomelocation', 'version'),
    'allele_definition': ('id', 'genesymbol', 'name', 'version'),
    'allele_location_value': ('alleledefinitionid', 'locationid', 'variantallele'),
    'allele': ('id', 'genesymbol', 'name', 'definitionid', 'clinicalfunctionalstatus'),
    'gene_result': ('id', 'genesymbol', 'result', 'activityscore'),
    'diplotype': ('genesymbol', 'diplotype', 'generesult', 'lookupkey'),
    'recommendation_view': ('recommendationid', 'drugname', 'lookupkey', 'drugrecommendation', 'population'),
    'guideline': ('id', 'name', 'genes', 'version'),
    'publication': ('id', 'guidelineid', 'title', 'version'),
    'change_log_view': ('date', 'type', 'entityname', 'note', 'deployedrelease'),
}
KEYS = {table: fields[:1] for table, fields in TABLES.items()}
KEYS.update(allele_location_value=('alleledefinitionid', 'locationid'), diplotype=('genesymbol', 'diplotype'))
ORDER = {table: ','.join(keys) for table, keys in KEYS.items()}
ORDER['change_log_view'] = 'date,type,entityname,note'
GENOTYPE_TABLES = {'gene', 'sequence_location', 'allele_definition', 'allele_location_value', 'allele'}
MAPPING_TABLES = {'diplotype', 'gene_result'}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def normalize_tables(tables):
    return {table: sorted(rows, key=canonical) for table, rows in tables.items()}


def fetch_json(url, count=False):
    headers = {'User-Agent': 'GenePassport/1.1', 'Accept-Encoding': 'gzip'}
    if count:
        headers['Prefer'] = 'count=exact'
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=60) as response:
        body = response.read()
        if response.headers.get('Content-Encoding') == 'gzip':
            body = gzip.decompress(body)
        return json.loads(body), dict(response.headers)


def fetch_table(table):
    rows, offset, expected = [], 0, None
    while offset < 500000:
        params = {'order': ORDER[table], 'limit': 5000, 'offset': offset}
        if table == 'diplotype':
            params['select'] = ','.join(TABLES[table])
        page, headers = fetch_json(API + '/' + table + '?' + urllib.parse.urlencode(params), count=True)
        if not isinstance(page, list):
            raise ValueError('CPIC returned a non-list for ' + table)
        content_range = next((v for k, v in headers.items() if k.lower() == 'content-range'), '')
        match = re.search(r'/(\d+)$', content_range)
        if not match:
            raise ValueError('CPIC did not supply an exact row count for ' + table)
        total = int(match[1])
        if expected is not None and expected != total:
            raise ValueError('CPIC changed during pagination: ' + table)
        expected = total
        rows.extend(page)
        offset += len(page)
        if offset == total:
            return rows
        if not page or offset > total:
            raise ValueError('Incomplete CPIC download: ' + table)
    raise ValueError('CPIC table exceeds supported size: ' + table)


def fetch_tables():
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = {table: pool.submit(fetch_table, table) for table in TABLES}
        # Observe every result, including failures, before leaving the pool.
        results, errors = {}, []
        for table, future in futures.items():
            try:
                results[table] = future.result()
            except Exception as error:
                errors.append(f'{table}: {error}')
        if errors:
            raise ValueError('; '.join(errors))
    return normalize_tables(results)


def latest_release():
    try:
        data, _ = fetch_json(RELEASE_URL)
        return {'tag': data['tag_name'], 'published_at': data['published_at'], 'url': data['html_url'], 'verified': True}
    except Exception as error:
        return {'tag': None, 'verified': False, 'error': str(error)}


def validate_tables(tables):
    checks = []
    for table, required in TABLES.items():
        rows = tables.get(table)
        if not isinstance(rows, list) or not rows:
            raise ValueError('Empty or absent CPIC table: ' + table)
        if any(not isinstance(row, dict) or any(key not in row for key in required) for row in rows):
            raise ValueError('Unsupported CPIC schema: ' + table)
        if table != 'change_log_view':
            ids = [tuple(str(row[key]) for key in KEYS[table]) for row in rows]
            if len(set(ids)) != len(ids):
                raise ValueError('Duplicate CPIC row identifiers: ' + table)
    gene_names = {r['symbol'] for r in tables['gene']}
    location_ids = {r['id'] for r in tables['sequence_location']}
    definition_ids = {r['id'] for r in tables['allele_definition']}
    if any(r['genesymbol'] not in gene_names for r in tables['pair_view'] + tables['diplotype'] + tables['allele_definition']):
        raise ValueError('CPIC references a missing gene')
    if any(r['alleledefinitionid'] not in definition_ids or r['locationid'] not in location_ids for r in tables['allele_location_value']):
        raise ValueError('CPIC allele definitions have broken links')
    if any(not isinstance(r['lookupkey'], dict) for r in tables['diplotype'] + tables['recommendation_view']):
        raise ValueError('Unsupported CPIC recommendation lookup keys')
    if not any(r['cpiclevel'] in ('A', 'B') for r in tables['pair_view']):
        raise ValueError('No actionable CPIC gene-drug pairs')
    checks.extend(['required schemas', 'unique identifiers', 'gene and location links', 'lookup-key shapes', 'actionable pairs'])
    # Safety properties rather than freezing historical phenotype values.
    for row in tables['diplotype']:
        if not row['lookupkey'] or row['generesult'] is None:
            raise ValueError('Incomplete CPIC diplotype mapping')
    checks.append('complete diplotype mappings')
    return {'passed': True, 'checks': checks, 'clinical_validation': False}


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix='.write-')
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(canonical(data))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def table_changes(old, new):
    return [table for table in TABLES if digest(old.get(table, [])) != digest(new.get(table, []))]


def impact(old_tables, new_tables):
    changed = table_changes(old_tables, new_tables)
    genes, drugs = set(), set()
    definitions = {row['id']: row['genesymbol'] for tables in (old_tables, new_tables) for row in tables.get('allele_definition', [])}
    guidelines = {row['id']: row.get('genes', []) for tables in (old_tables, new_tables) for row in tables.get('guideline', [])}
    def collect(row):
        if row.get('genesymbol'):
            genes.add(row['genesymbol'])
        if row.get('symbol'):
            genes.add(row['symbol'])
        genes.update(row.get('genes') or [])
        genes.update((row.get('lookupkey') or {}).keys())
        if row.get('alleledefinitionid') in definitions:
            genes.add(definitions[row['alleledefinitionid']])
        if row.get('guidelineid') in guidelines:
            genes.update(guidelines[row['guidelineid']])
        if row.get('drugname'):
            drugs.add(row['drugname'])
    row_counts = {}
    for table in changed:
        before = {canonical(row): row for row in old_tables.get(table, [])}
        after = {canonical(row): row for row in new_tables.get(table, [])}
        delta = before.keys() ^ after.keys()
        row_counts[table] = len(delta)
        for key in delta:
            collect(before.get(key) or after[key])
    for tables in (old_tables, new_tables):
        for pair in tables.get('pair_view', []):
            if pair.get('genesymbol') in genes:
                drugs.add(pair['drugname'])
    return {'changed_tables': changed, 'changed_rows': row_counts, 'affected_genes': sorted(genes), 'affected_drugs': sorted(drugs),
            'genome_reanalysis_required': bool(set(changed) & GENOTYPE_TABLES),
            'gene_result_review_required': bool(set(changed) & MAPPING_TABLES),
            'drug_recheck_required': bool(set(changed) - {'change_log_view'}),
            'scope': 'conservative; metadata-only changes can also trigger review'}


class ResourceStore:
    def __init__(self, root):
        self.root = Path(root)
        self.lock = threading.RLock()
        self.check_lock = threading.Lock()
        self.cache = {}

    @contextmanager
    def update_lock(self):
        # Coordinate the CLI and server as well as concurrent HTTP requests.
        with self.check_lock:
            self.root.mkdir(parents=True, exist_ok=True)
            with open(self.root / "update.lock", "a") as stream:
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
                try:
                    yield
                finally:
                    fcntl.flock(stream.fileno(), fcntl.LOCK_UN)

    def state(self):
        with self.lock:
            path = self.root / 'state.json'
            return json.loads(path.read_text()) if path.exists() else {'active': None, 'candidate': None, 'last_attempt': None, 'last_success': None, 'error': None}

    def write_state(self, state):
        with self.lock:
            atomic_json(self.root / 'state.json', state)

    def load(self, snapshot_id):
        if not re.fullmatch(r'[a-f0-9]{64}', snapshot_id or ''):
            raise ValueError('Invalid CPIC snapshot identifier')
        with self.lock:
            if snapshot_id not in self.cache:
                folder = self.root / 'snapshots' / snapshot_id
                manifest = json.loads((folder / 'manifest.json').read_text())
                with gzip.open(folder / 'tables.json.gz', 'rt') as stream:
                    tables = json.load(stream)
                if digest(tables) != snapshot_id or manifest['snapshot_id'] != snapshot_id:
                    raise ValueError('CPIC snapshot integrity check failed')
                self.cache[snapshot_id] = {'manifest': manifest, 'tables': tables}
            return self.cache[snapshot_id]

    def save_snapshot(self, tables, release, fetched_from, checks):
        snapshot_id = digest(tables)
        folder = self.root / 'snapshots' / snapshot_id
        manifest = {'schema': 'cpic-snapshot-v1', 'snapshot_id': snapshot_id, 'source': API,
                    'retrieved_at': utc_now(), 'retrieval_started_at': fetched_from, 'release': release,
                    'api_interface': 'v1', 'consistency': 'two matching complete API reads; not a database transaction',
                    'table_hashes': {t: digest(rows) for t, rows in tables.items()},
                    'row_counts': {t: len(rows) for t, rows in tables.items()}, 'validation': checks,
                    'source_registry': [
                        {'source': 'CPIC', 'role': 'prescribing and translation', 'used_for_verdict': True},
                        {'source': 'PharmVar', 'role': 'allele nomenclature; not connected', 'used_for_verdict': False},
                        {'source': 'DPWG', 'role': 'additional guidance; not connected', 'used_for_verdict': False},
                        {'source': 'ClinPGx', 'role': 'curated evidence; not independent CPIC confirmation', 'used_for_verdict': False},
                        {'source': 'regulatory labels', 'role': 'jurisdiction-specific guidance; not connected', 'used_for_verdict': False},
                    ]}
        with self.lock:
            if not folder.exists():
                folder.parent.mkdir(parents=True, exist_ok=True)
                temporary = Path(tempfile.mkdtemp(dir=folder.parent, prefix='.snapshot-'))
                try:
                    with gzip.open(temporary / 'tables.json.gz', 'wb') as stream:
                        stream.write(canonical(tables))
                    atomic_json(temporary / 'manifest.json', manifest)
                    os.replace(temporary, folder)
                finally:
                    if temporary.exists():
                        import shutil
                        shutil.rmtree(temporary)
            # Never rewrite a historical snapshot manifest with a new fetch date.
        return snapshot_id

    def check(self, force=False, fetcher=fetch_tables, release_fetcher=latest_release):
        with self.update_lock():
            state = self.state()
            if not force and state.get('last_attempt') and time.time() - state['last_attempt'] < CHECK_INTERVAL:
                return self.status()
            state.update(last_attempt=time.time(), error=None)
            self.write_state(state)
            try:
                started = utc_now()
                release_before = release_fetcher()
                first = normalize_tables(fetcher())
                checks = validate_tables(first)
                second = normalize_tables(fetcher())
                if digest(first) != digest(second):
                    raise ValueError('CPIC changed while reading. Retry before using these data.')
                release_after = release_fetcher()
                if release_before.get('verified') and release_after.get('verified') and release_before['tag'] != release_after['tag']:
                    raise ValueError('CPIC published a release during this check. Retry.')
                snapshot_id = self.save_snapshot(first, release_after, started, checks)
                state.update(last_success=time.time(), latest_snapshot=snapshot_id, latest_release=release_after, error=None)
                if not state['active']:
                    state.update(active=snapshot_id, candidate=None)
                elif state['active'] != snapshot_id:
                    state['candidate'] = snapshot_id
                else:
                    state['candidate'] = None
                self.write_state(state)
            except Exception as error:
                state['error'] = str(error)
                self.write_state(state)
            return self.status()

    def status(self):
        state = self.state()
        active = self.load(state['active'])['manifest'] if state['active'] else None
        candidate = self.load(state['candidate'])['manifest'] if state.get('candidate') else None
        fresh = bool(state.get('last_success') and time.time() - state['last_success'] < CHECK_INTERVAL)
        if state.get('error'):
            freshness = 'unverified'
        elif not fresh:
            freshness = 'stale' if active else 'uninitialized'
        elif candidate:
            freshness = 'update_available'
        else:
            freshness = 'current_as_of_check'
        result = {'freshness': freshness, 'checked_at': datetime.fromtimestamp(state['last_success'], timezone.utc).isoformat() if state.get('last_success') else None,
                  'last_attempt_at': datetime.fromtimestamp(state['last_attempt'], timezone.utc).isoformat() if state.get('last_attempt') else None,
                  'error': state.get('error'), 'active': active, 'candidate': candidate,
                  'latest_release': state.get('latest_release'), 'check_interval_seconds': CHECK_INTERVAL,
                  'checking': self.check_lock.locked()}
        if candidate and active:
            result['impact'] = impact(self.load(active['snapshot_id'])['tables'], self.load(candidate['snapshot_id'])['tables'])
        return result

    def ensure(self):
        status = self.check()
        if not status['active']:
            raise ValueError('No verified CPIC snapshot is available. Check resources and retry.')
        return self.load(status['active']['snapshot_id']), status

    def activate(self, snapshot_id):
        with self.update_lock():
            state = self.state()
            if state.get('candidate') != snapshot_id:
                raise ValueError('This is no longer the staged CPIC update. Refresh the review.')
            if state.get('error') or not state.get('last_success') or time.time() - state['last_success'] >= CHECK_INTERVAL:
                raise ValueError('Check CPIC successfully again before activating this update.')
            snapshot = self.load(snapshot_id)
            validate_tables(snapshot['tables'])
            report = impact(self.load(state['active'])['tables'], snapshot['tables'])
            history = state.setdefault('activation_history', [])
            history.append({'from': state['active'], 'to': snapshot_id, 'activated_at': utc_now(), 'impact': report})
            state.update(active=snapshot_id, candidate=None)
            self.write_state(state)
        return self.status()

    def provenance(self, snapshot, status):
        return {'snapshot_id': snapshot['manifest']['snapshot_id'], 'source': API,
                'retrieved_at': snapshot['manifest']['retrieved_at'], 'release': snapshot['manifest']['release'],
                'checked_at': status['checked_at'], 'freshness_at_analysis': status['freshness'],
                'caller_version': CALLER_VERSION, 'table_hashes': snapshot['manifest']['table_hashes'],
                'consistency': snapshot['manifest']['consistency']}

    def passport_impact(self, passport):
        state = self.state()
        if not state['active']:
            return {'status': 'unverified', 'message': 'No active CPIC snapshot'}
        previous = (passport.get('cpic_provenance') or {}).get('snapshot_id')
        current = state['active']
        if previous == current:
            return {'status': 'same_snapshot', 'active_snapshot': current, 'message': 'Passport uses the active CPIC snapshot.'}
        if not re.fullmatch(r'[a-f0-9]{64}', previous or '') or not (self.root / 'snapshots' / previous / 'manifest.json').exists():
            return {'status': 'untracked', 'active_snapshot': current, 'genome_reanalysis_required': True,
                    'message': 'The original CPIC snapshot is unavailable. Reanalysis is needed to establish version history.'}
        report = impact(self.load(previous)['tables'], self.load(current)['tables'])
        passport_genes = set((passport.get('genes') or {}).keys())
        relevant = sorted(passport_genes & set(report['affected_genes']))
        report.update(status='changed_snapshot', previous_snapshot=previous, active_snapshot=current,
                      passport_affected_genes=relevant,
                      new_genes=sorted({r['genesymbol'] for r in self.load(current)['tables']['pair_view'] if r['cpiclevel'] in ('A', 'B')} - passport_genes))
        # This endpoint only reports impact, never rewrites or upgrades patient calls.
        return report

    def query(self, snapshot_id, table, params):
        if table not in {'pair_view', 'diplotype', 'gene_result', 'recommendation_view'}:
            raise ValueError('Unsupported CPIC snapshot query')
        rows = self.load(snapshot_id)['tables'][table]
        for field, expression in params.items():
            if field == 'select':
                continue
            if expression.startswith('eq.'):
                value = expression[3:]
                rows = [row for row in rows if str(row.get(field)) == value]
            elif expression.startswith('in.(') and expression.endswith(')'):
                values = expression[4:-1].split(',')
                rows = [row for row in rows if str(row.get(field)) in values]
            elif expression.startswith('cs.') and field == 'lookupkey':
                key = json.loads(expression[3:])
                if not isinstance(key, dict):
                    raise ValueError('Invalid lookup key')
                rows = [row for row in rows if all((row.get(field) or {}).get(k) == v for k, v in key.items())]
            else:
                raise ValueError('Unsupported CPIC query filter')
        if 'select' in params:
            fields = params['select'].split(',')
            if any(not re.fullmatch(r'[a-zA-Z_]+', field) for field in fields):
                raise ValueError('Unsupported CPIC selection')
            rows = [{field: row.get(field) for field in fields} for row in rows]
        return rows
