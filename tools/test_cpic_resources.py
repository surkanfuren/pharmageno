"""Deterministic update/failure/impact tests. Never touches the real CPIC store."""
import copy
import gzip
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cpic_resources as cpic


def fixture():
    return {
        'pair_view': [{'pairid': 1, 'drugname': 'clopidogrel', 'genesymbol': 'CYP2C19', 'cpiclevel': 'A'}],
        'gene': [{'symbol': 'CYP2C19', 'chr': 'chr10', 'lookupmethod': 'PHENOTYPE'}],
        'sequence_location': [{'id': 1, 'genesymbol': 'CYP2C19', 'chromosomelocation': 'g.10G>A', 'version': 1}],
        'allele_definition': [{'id': 2, 'genesymbol': 'CYP2C19', 'name': '*2', 'version': 1}],
        'allele_location_value': [{'alleledefinitionid': 2, 'locationid': 1, 'variantallele': 'A'}],
        'allele': [{'id': 3, 'genesymbol': 'CYP2C19', 'name': '*2', 'definitionid': 2, 'clinicalfunctionalstatus': 'No function'}],
        'gene_result': [{'id': 4, 'genesymbol': 'CYP2C19', 'result': 'Poor Metabolizer', 'activityscore': 'n/a'}],
        'diplotype': [{'genesymbol': 'CYP2C19', 'diplotype': '*2/*2', 'generesult': 'Poor Metabolizer', 'lookupkey': {'CYP2C19': 'Poor Metabolizer'}}],
        'recommendation_view': [{'recommendationid': 5, 'drugname': 'clopidogrel', 'lookupkey': {'CYP2C19': 'Poor Metabolizer'}, 'drugrecommendation': 'Avoid clopidogrel.', 'population': 'general'}],
        'guideline': [{'id': 6, 'name': 'CYP2C19 and clopidogrel', 'genes': ['CYP2C19'], 'version': 1}],
        'publication': [{'id': 7, 'guidelineid': 6, 'title': 'Guideline', 'version': 1}],
        'change_log_view': [{'date': '2026-10-12', 'type': 'GUIDELINE', 'entityname': 'CYP2C19', 'note': 'Test', 'deployedrelease': 'vTEST'}],
    }


def release():
    return {'tag': 'vTEST', 'verified': True, 'published_at': '2026-10-12T00:00:00Z', 'url': 'https://example.test/release'}


class ResourceSafety(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.store = cpic.ResourceStore(self.directory.name)
        self.data = fixture()

    def check(self, data=None, **kwargs):
        return self.store.check(force=True, fetcher=lambda: copy.deepcopy(data or self.data), release_fetcher=kwargs.pop('release_fetcher', release), **kwargs)

    def test_initial_verified_baseline_and_pinned_queries(self):
        result = self.check()
        self.assertEqual(result['freshness'], 'current_as_of_check')
        snapshot = result['active']['snapshot_id']
        self.assertFalse(result['active']['validation']['clinical_validation'])
        rows = self.store.query(snapshot, 'diplotype', {'genesymbol': 'eq.CYP2C19', 'diplotype': 'eq.*2/*2', 'select': 'generesult,lookupkey'})
        self.assertEqual(rows[0]['generesult'], 'Poor Metabolizer')
        self.assertEqual(self.store.query(snapshot, 'pair_view', {'cpiclevel': 'in.(A,B)'}), self.data['pair_view'])
        self.assertEqual(len(self.store.query(snapshot, 'recommendation_view', {'lookupkey': 'cs.{"CYP2C19":"Poor Metabolizer"}'})), 1)
        self.assertEqual(self.store.query(snapshot, 'recommendation_view', {'lookupkey': 'cs.{"CYP2C19":"Normal Metabolizer"}'}), [])

    def test_changed_data_staged_even_when_release_is_identical(self):
        original = self.check()['active']['snapshot_id']
        changed = copy.deepcopy(self.data)
        changed['recommendation_view'][0]['drugrecommendation'] = 'New recommendation'
        result = self.check(changed)
        self.assertEqual(result['freshness'], 'update_available')
        self.assertEqual(result['active']['snapshot_id'], original)
        self.assertTrue(result['impact']['drug_recheck_required'])
        self.assertFalse(result['impact']['genome_reanalysis_required'])
        self.assertEqual(result['impact']['affected_drugs'], ['clopidogrel'])
        candidate = result['candidate']['snapshot_id']
        self.store.activate(candidate)
        self.assertEqual(self.store.state()['active'], candidate)
        # Historical lookup remains available after activation.
        self.assertEqual(self.store.query(original, 'recommendation_view', {})[0]['drugrecommendation'], 'Avoid clopidogrel.')

    def test_allele_changes_require_reanalysis(self):
        self.check()
        changed = copy.deepcopy(self.data)
        changed['allele_location_value'][0]['variantallele'] = 'T'
        result = self.check(changed)
        self.assertTrue(result['impact']['genome_reanalysis_required'])
        self.assertEqual(result['impact']['affected_genes'], ['CYP2C19'])

    def test_mapping_changes_are_distinguished(self):
        self.check()
        changed = copy.deepcopy(self.data)
        changed['diplotype'][0]['generesult'] = 'Indeterminate'
        result = self.check(changed)
        self.assertTrue(result['impact']['gene_result_review_required'])
        self.assertFalse(result['impact']['genome_reanalysis_required'])

    def test_failure_preserves_active_and_never_claims_current(self):
        original = self.check()['active']['snapshot_id']
        def fail():
            raise OSError('offline')
        result = self.store.check(force=True, fetcher=fail, release_fetcher=release)
        self.assertEqual(result['freshness'], 'unverified')
        self.assertEqual(result['active']['snapshot_id'], original)
        self.assertIn('offline', result['error'])

    def test_failure_without_baseline(self):
        result = self.check(data={'invalid': []})
        self.assertIsNone(result['active'])
        self.assertEqual(result['freshness'], 'unverified')

    def test_partial_or_broken_links_cannot_be_activated(self):
        self.check()
        for table in ('diplotype', 'recommendation_view', 'allele_location_value'):
            changed = copy.deepcopy(self.data)
            changed[table] = []
            result = self.check(changed)
            self.assertEqual(result['freshness'], 'unverified')
            self.assertIsNone(result['candidate'])
        changed = copy.deepcopy(self.data)
        changed['allele_location_value'][0]['locationid'] = 999
        self.assertEqual(self.check(changed)['freshness'], 'unverified')

    def test_two_different_reads_rejected(self):
        changed = copy.deepcopy(self.data)
        changed['guideline'][0]['version'] = 2
        reads = iter([self.data, changed])
        result = self.store.check(force=True, fetcher=lambda: next(reads), release_fetcher=release)
        self.assertIsNone(result['active'])
        self.assertIn('changed while reading', result['error'])

    def test_release_change_during_download_rejected(self):
        versions = iter([release(), dict(release(), tag='vNEW')])
        result = self.check(release_fetcher=lambda: next(versions))
        self.assertIn('published a release', result['error'])
        self.assertIsNone(result['active'])

    def test_release_failure_separate_from_live_data_check(self):
        result = self.check(release_fetcher=lambda: {'verified': False, 'tag': None, 'error': 'GitHub unavailable'})
        self.assertEqual(result['freshness'], 'current_as_of_check')
        self.assertFalse(result['latest_release']['verified'])

    def test_row_order_does_not_create_an_update(self):
        self.data['publication'].append(dict(self.data['publication'][0], id=8))
        original = self.check()['active']['snapshot_id']
        reordered = copy.deepcopy(self.data)
        reordered['publication'].reverse()
        result = self.check(reordered)
        self.assertEqual(result['active']['snapshot_id'], original)
        self.assertIsNone(result['candidate'])

    def test_stale_candidate_cannot_be_activated(self):
        self.check()
        changed = copy.deepcopy(self.data)
        changed['guideline'][0]['version'] = 2
        candidate = self.check(changed)['candidate']['snapshot_id']
        state = self.store.state()
        state['last_success'] = time.time() - cpic.CHECK_INTERVAL - 1
        self.store.write_state(state)
        self.assertEqual(self.store.status()['freshness'], 'stale')
        with self.assertRaises(ValueError):
            self.store.activate(candidate)
        with self.assertRaises(ValueError):
            self.store.activate('0' * 64)

    def test_passport_impact_and_new_missing_gene(self):
        old = self.check()['active']['snapshot_id']
        changed = copy.deepcopy(self.data)
        changed['gene'].append({'symbol': 'NEWGENE', 'chr': 'chr2', 'lookupmethod': 'PHENOTYPE'})
        changed['pair_view'].append({'pairid': 20, 'drugname': 'newdrug', 'genesymbol': 'NEWGENE', 'cpiclevel': 'A'})
        new = self.check(changed)['candidate']['snapshot_id']
        self.store.activate(new)
        passport = {'cpic_provenance': {'snapshot_id': old}, 'genes': {'CYP2C19': {'status': 'unknown'}}}
        result = self.store.passport_impact(passport)
        self.assertEqual(result['new_genes'], ['NEWGENE'])
        self.assertEqual(result['status'], 'changed_snapshot')
        self.assertEqual(passport['genes']['CYP2C19']['status'], 'unknown')
        self.assertEqual(self.store.passport_impact({'genes': {}})['status'], 'untracked')

    def test_restarted_store_detects_corruption(self):
        snapshot = self.check()['active']['snapshot_id']
        path = Path(self.directory.name) / 'snapshots' / snapshot / 'tables.json.gz'
        with gzip.open(path, 'wt') as stream:
            json.dump({'tampered': []}, stream)
        restarted = cpic.ResourceStore(self.directory.name)
        with self.assertRaisesRegex(ValueError, 'integrity'):
            restarted.load(snapshot)

    def test_capped_pagination_and_count_change(self):
        # The server may cap pages below our requested 5,000 rows.
        pages = iter([([{'symbol': 'A'}], {'Content-Range': '0-0/2'}), ([{'symbol': 'B'}], {'Content-Range': '1-1/2'})])
        with patch.object(cpic, 'fetch_json', side_effect=lambda *a, **kw: next(pages)):
            self.assertEqual(len(cpic.fetch_table('gene')), 2)
        pages = iter([([{'symbol': 'A'}], {'Content-Range': '0-0/2'}), ([{'symbol': 'B'}], {'Content-Range': '1-1/3'})])
        with patch.object(cpic, 'fetch_json', side_effect=lambda *a, **kw: next(pages)):
            with self.assertRaisesRegex(ValueError, 'pagination'):
                cpic.fetch_table('gene')

    def test_daily_check_reuses_snapshot_without_network(self):
        self.check()
        def unexpected():
            self.fail('Unnecessary daily download')
        result = self.store.check(fetcher=unexpected, release_fetcher=unexpected)
        self.assertEqual(result['freshness'], 'current_as_of_check')


if __name__ == '__main__':
    unittest.main()
