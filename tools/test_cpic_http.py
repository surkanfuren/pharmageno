"""HTTP update lifecycle against isolated synthetic CPIC resources."""
import copy
import json
import pathlib
import sys
import tempfile
import threading
import unittest
import urllib.error
import urllib.parse
import urllib.request
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import server
from test_cpic_resources import fixture, release
from cpic_resources import ResourceStore


class UpdateHTTP(unittest.TestCase):
    def test_update_lifecycle_and_historical_query(self):
        with tempfile.TemporaryDirectory() as root:
            store = ResourceStore(root)
            old_store = server.RESOURCES
            server.RESOURCES = store
            httpd = server.ThreadingHTTPServer(('127.0.0.1', 0), server.Handler)
            thread = threading.Thread(target=httpd.serve_forever, daemon=True)
            thread.start()
            base = 'http://127.0.0.1:' + str(httpd.server_port)
            def request(path, body=None):
                data = None if body is None else json.dumps(body).encode()
                req = urllib.request.Request(base + path, data=data, headers={'Content-Type': 'application/json'})
                with urllib.request.urlopen(req, timeout=10) as response:
                    self.assertEqual(response.headers['Cache-Control'], 'no-store')
                    return json.load(response)
            try:
                data = fixture()
                baseline = store.check(force=True, fetcher=lambda: copy.deepcopy(data), release_fetcher=release)
                sid = baseline['active']['snapshot_id']
                self.assertEqual(request('/api/cpic/status')['active']['snapshot_id'], sid)
                query = '/api/cpic/diplotype?' + urllib.parse.urlencode({'snapshot': sid, 'genesymbol': 'eq.CYP2C19', 'diplotype': 'eq.*2/*2'})
                self.assertEqual(request(query)[0]['generesult'], 'Poor Metabolizer')
                changed = copy.deepcopy(data)
                changed['diplotype'][0]['generesult'] = 'Indeterminate'
                staged = store.check(force=True, fetcher=lambda: copy.deepcopy(changed), release_fetcher=release)
                candidate = staged['candidate']['snapshot_id']
                self.assertEqual(request('/api/cpic/status')['freshness'], 'update_available')
                activated = request('/api/cpic/activate', {'snapshot_id': candidate})
                self.assertEqual(activated['active']['snapshot_id'], candidate)
                report = request('/api/cpic/passport-impact', {'cpic_provenance': {'snapshot_id': sid}, 'genes': {'CYP2C19': {}}})
                self.assertTrue(report['gene_result_review_required'])
                self.assertEqual(report['passport_affected_genes'], ['CYP2C19'])
                self.assertEqual(request(query)[0]['generesult'], 'Poor Metabolizer')
                for path, body in [('/api/cpic/activate', {'snapshot_id': sid}), ('/api/cpic/manifest?snapshot=invalid', None)]:
                    with self.assertRaises(urllib.error.HTTPError) as failed:
                        request(path, body)
                    self.assertEqual(failed.exception.code, 400)
                    failed.exception.close()
            finally:
                httpd.shutdown()
                httpd.server_close()
                server.RESOURCES = old_store


if __name__ == '__main__':
    unittest.main()
