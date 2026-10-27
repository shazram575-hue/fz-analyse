"""Check archived delivery hashes without changing them; classify mutable-document drift."""
from pathlib import Path
import csv,hashlib,json
root=Path(__file__).resolve().parents[2];out=root/'docs/audit/2026-10-26';rows=[];cache={}
manifests=[
 '第一阶段_2026-10-23/audit/Final_artifact_SHA256_manifest.tsv',
 '单细胞基础补全_2026-10-23/audit/Deliverable_SHA256_manifest.tsv',
 '身份深核与独立交集_2026-10-24/audit/Deliverable_SHA256_manifest.tsv',
 'reverse_docking47_2026-10-24/audit/Deliverable_SHA256_manifest.tsv',
 'MFSD10_表达分布摸底_2026-10-26/audit/Deliverable_SHA256.tsv',
 'MFSD10_程序发现与独立验证_2026-10-26/audit/Deliverable_SHA256.tsv',
 'MFSD10_稀疏与聚合探索_2026-10-26/audit/Deliverable_SHA256.tsv']
def sha(p):
 key=str(p)
 if key in cache:return cache[key]
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 cache[key]=h.hexdigest();return cache[key]
for rel in manifests:
 m=root/'outputs'/rel;base=m.parent.parent
 with m.open(encoding='utf-8-sig') as f:entries=list(csv.DictReader(f,delimiter='\t'))
 for r in entries:
  s=r.get('relative_path',r.get('path',''));value=Path(s)
  p=value if value.is_absolute() else ((root/value) if (root/value).exists() else (base/value))
  if not p.exists():status='missing';actual=''
  else:
   actual=sha(p);status='match' if actual.lower()==r['sha256'].lower() else 'hash_mismatch'
  kind='scientific_or_stage_file'
  if p.parent.resolve()==root:kind='mutable_root_document'
  if 'snapshot' in s.lower() or '快照' in s:kind='historical_snapshot'
  rows.append(dict(manifest='outputs/'+rel,path=s,resolved_path=str(p),status=status,category=kind,expected_sha256=r['sha256'],actual_sha256=actual))
 print('MANIFEST_CHECKED',rel,len(entries),flush=True)
with (out/'STAGE_MANIFEST_CHECKS.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
summary=dict(records=len(rows),matched=sum(r['status']=='match' for r in rows),missing=sum(r['status']=='missing' for r in rows),mismatches=sum(r['status']=='hash_mismatch' for r in rows),nonroot_mismatches=sum(r['status']!='match' and r['category']=='scientific_or_stage_file' for r in rows),root_document_drift=sum(r['status']!='match' and r['category']=='mutable_root_document' for r in rows),manifests_modified=0)
(out/'STAGE_MANIFEST_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
