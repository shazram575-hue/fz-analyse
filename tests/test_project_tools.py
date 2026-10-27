import unittest,tempfile,csv,hashlib,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TOOL=ROOT/'scripts/project/restore_legacy_layout.py'
class RestoreSafetyTests(unittest.TestCase):
 def setUp(self):
  base=ROOT/'docs/audit/2026-10-26/after/tool_tests';base.mkdir(parents=True,exist_ok=True)
  self.temp=tempfile.TemporaryDirectory(dir=base);self.root=Path(self.temp.name).resolve()
  assert self.root.is_relative_to(base.resolve()) # Verify cleanup target is confined to new audit fixtures.
  (self.root/'docs').mkdir();(self.root/'analysis').mkdir();self.code=self.root/'analysis/example.py';self.code.write_bytes(b'print("test")\n')
  self.write_map('outputs/stage/scripts/example.py')
 def tearDown(self):
  assert self.root.is_relative_to((ROOT/'docs/audit/2026-10-26/after/tool_tests').resolve())
  self.temp.cleanup()
 def write_map(self,target,sha=None):
  with (self.root/'docs/PATH_MAP.tsv').open('w',encoding='utf-8',newline='') as f:
   w=csv.DictWriter(f,fieldnames=['source','repository_path','sha256'],delimiter='\t');w.writeheader();w.writerow(dict(source=target,repository_path='analysis/example.py',sha256=sha or hashlib.sha256(self.code.read_bytes()).hexdigest()))
 def run_tool(self,*flags):return subprocess.run([sys.executable,str(TOOL),'--root',str(self.root),*flags],capture_output=True,text=True)
 def test_dry_run_never_writes(self):
  r=self.run_tool('--dry-run');self.assertEqual(r.returncode,0,r.stderr);self.assertFalse((self.root/'outputs').exists())
 def test_apply_identical_bytes_and_idempotent(self):
  r=self.run_tool('--apply');self.assertEqual(r.returncode,0,r.stderr);self.assertEqual((self.root/'outputs/stage/scripts/example.py').read_bytes(),self.code.read_bytes());self.assertEqual(self.run_tool('--apply').returncode,0)
 def test_refuses_different_existing_file(self):
  t=self.root/'outputs/stage/scripts/example.py';t.parent.mkdir(parents=True);t.write_bytes(b'KEEP\n');r=self.run_tool('--apply');self.assertNotEqual(r.returncode,0);self.assertEqual(t.read_bytes(),b'KEEP\n')
 def test_rejects_path_escape(self):
  self.write_map('../escape.py');r=self.run_tool('--apply');self.assertNotEqual(r.returncode,0);self.assertFalse((self.root.parent/'escape.py').exists())
 def test_bad_source_hash_rejected(self):
  self.write_map('outputs/stage/scripts/example.py','0'*64);self.assertNotEqual(self.run_tool('--apply').returncode,0);self.assertFalse((self.root/'outputs').exists())
if __name__=='__main__':unittest.main()
