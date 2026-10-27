"""Safety checks use synthetic fixtures, never real credentials or scientific data."""
from pathlib import Path
import importlib.util
import json
import struct
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('github_preflight', ROOT/'scripts/project/github_preflight.py')
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)


class UploadSafetyTests(unittest.TestCase):
    def setUp(self):
        base = ROOT/'docs/audit/2026-10-27/github-sync/local/preflight-tests'
        base.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=base)
        self.root = Path(self.temp.name).resolve()
        assert self.root.is_relative_to(base.resolve())
        (self.root/'docs').mkdir()
        self.original_root = tool.ROOT
        tool.ROOT = self.root

    def tearDown(self):
        tool.ROOT = self.original_root
        assert self.root.is_relative_to((ROOT/'docs/audit/2026-10-27/github-sync/local/preflight-tests').resolve())
        self.temp.cleanup()

    def check(self, name, content):
        path = self.root/'docs'/name
        path.write_bytes(content)
        return tool.inspect_file('docs/'+name)

    def test_refuses_escape_without_reading_outside_project(self):
        row, flags = tool.inspect_file('../outside.md')
        self.assertIsNone(row)
        self.assertTrue(any(f['severity'] == 'block' for f in flags))

    def test_synthetic_credential_is_blocked_and_value_not_reported(self):
        field, value = 'password', 'a'*32
        _, flags = self.check('fixture.txt', f'{field} = {json.dumps(value)}\n'.encode())
        self.assertTrue(any(f['kind'] == 'credential_literal' for f in flags))
        self.assertNotIn(value, json.dumps(flags))

    def test_direct_identity_blocked_but_deidentified_fields_require_source_review(self):
        _, flags = self.check('direct.tsv', b'patient_name\tage\nSynthetic\t42\n')
        self.assertTrue(any(f['kind'] == 'direct_identifier_header' for f in flags))
        row, flags = self.check('summary.tsv', b'patient_id\tmedian_UMI\nPublicCase\t1000\n')
        self.assertEqual(row['classification'], 'case_or_clinical_fields_review_source')
        self.assertFalse(flags)

    def test_unknown_binary_format_is_blocked_even_when_small(self):
        _, flags = self.check('fixture.h5ad', b'tiny fixture')
        self.assertTrue(any(f['kind'] == 'unapproved_extension' for f in flags))

    def test_png_metadata_is_flagged_for_manual_review(self):
        data = b'Software\0Synthetic renderer'
        png = b'\x89PNG\r\n\x1a\n'+struct.pack('>I', len(data))+b'tEXt'+data+b'\0'*4
        _, flags = self.check('figure.png', png)
        self.assertTrue(any(f['kind'] == 'png_metadata_tEXt' for f in flags))

    def test_non_utf8_text_is_blocked(self):
        _, flags = self.check('fixture.md', b'\xff\xfe\x80')
        self.assertTrue(any(f['kind'] == 'non_utf8_text' for f in flags))


if __name__ == '__main__':
    unittest.main()
