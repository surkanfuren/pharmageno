"""Safety checks for the stdlib full-genome evidence importer (no network)."""
import importlib.util
import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

spec = importlib.util.spec_from_file_location('passport_server', pathlib.Path(__file__).resolve().parents[1] / 'server.py')
server = importlib.util.module_from_spec(spec)
spec.loader.exec_module(server)


class GenomeSafety(unittest.TestCase):
    def analyze(self, records='', bed=None, header='GRCh38', samples='SAMPLE'):
        sites = {('1', p): {'chrom': '1', 'position': p, 'ref': 'A', 'genes': ['TEST'], 'rsids': [f'rs{p}'], 'versions': [1]} for p in (10, 11, 12)}
        catalog = ([{'genesymbol': 'TEST', 'drugname': 'example'}], sites, {'TEST': 0})
        with tempfile.TemporaryDirectory() as root:
            vcf = pathlib.Path(root) / 'input.vcf'
            vcf.write_text(f'##fileformat=VCFv4.2\n##reference={header}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n' + records)
            path = pathlib.Path(root) / 'coverage.bed'
            if bed is not None:
                path.write_text(bed)
            return server.analyze_genome(str(vcf), str(path) if bed is not None else None, catalog, trusted_bed=True)

    def test_missing_is_unknown(self):
        result = self.analyze()
        self.assertTrue(all(i['genotype'] is None for i in result['evidence']))
        self.assertEqual(result['genes']['TEST']['status'], 'unknown')

    def test_bed_zero_based_half_open(self):
        result = self.analyze(bed='chr1\t10\t11\n')
        self.assertEqual([i['genotype'] for i in result['evidence']], [None, 'AA', None])

    def test_multiallelic_and_phased(self):
        result = self.analyze('chr1\t10\t.\tA\tG,C\t100\tPASS\t.\tGT:PS\t1|2:10\n')
        item = result['evidence'][0]
        self.assertEqual(item['genotype'], 'GC')
        self.assertTrue(item['phased'])
        self.assertEqual(item['phase_set'], '10')
        self.assertEqual(result['genes']['TEST']['status'], 'unknown')

    def test_no_call_and_failed_filter_block_reference(self):
        for gt, filt in [('./.', 'PASS'), ('0/1', 'LowQual'), ('0', 'PASS'), ('5/5', 'PASS')]:
            result = self.analyze(f'chr1\t10\t.\tA\tG\t10\t{filt}\t.\tGT\t{gt}\n', bed='1\t0\t20\n')
            self.assertIsNone(result['evidence'][0]['genotype'])

    def test_deletion_blocks_spanned_sites(self):
        result = self.analyze('1\t9\t.\tAAAA\tA\t100\tPASS\t.\tGT\t0/1\n', bed='1\t0\t20\n')
        self.assertTrue(all(i['genotype'] is None for i in result['evidence']))

    def test_duplicate_blocks_reference(self):
        row = '1\t10\t.\tA\tG\t100\tPASS\t.\tGT\t0/1\n'
        result = self.analyze(row + row, bed='1\t0\t20\n')
        self.assertEqual(result['evidence'][0]['status'], 'conflicting_records')
        self.assertIsNone(result['evidence'][0]['genotype'])

    def test_wrong_reference(self):
        result = self.analyze('1\t10\t.\tT\tG\t100\tPASS\t.\tGT\t0/1\n', bed='1\t0\t20\n')
        self.assertIsNone(result['evidence'][0]['genotype'])

    def test_unsupported_build_and_multiple_samples(self):
        for kwargs in ({'header': 'GRCh37'}, {'header': 'unknown'}, {'samples': 'A\tB'}):
            with self.assertRaises(ValueError):
                self.analyze(**kwargs)


if __name__ == '__main__':
    unittest.main()
