"""Finalize only after numerical audit and explicit figure review have completed."""
import hashlib
import json
import pathlib
import re
import pandas as pd
ROOT=pathlib.Path(__file__).resolve().parents[1]
def main():
    summary=json.loads((ROOT/'audit'/'Independent_validation_summary.json').read_text())
    assert summary['n_checks']==summary['n_pass'] and not summary['failed']
    assert (ROOT/'audit'/'Figure_visual_QA.md').exists()
    links=[]
    for document in [ROOT/'README_当前状态.md',ROOT/'身份核对结论.md',ROOT/'独立预测与疾病交集_结果报告.md']:
        for target in re.findall(r'\]\(([^)]+)\)',document.read_text(encoding='utf-8')):
            if target.startswith(('https:','http:','#')):continue
            path=document.parent/target
            links.append({'document':document.name,'target':target,'exists':path.exists()})
    pd.DataFrame(links).to_csv(ROOT/'audit'/'Document_links.tsv',sep='\t',index=False)
    assert all(row['exists'] for row in links),[r for r in links if not r['exists']]
    files=[]
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file() or any(part in {'cache','__pycache__'} for part in path.relative_to(ROOT).parts):continue
        if path.name in {'Deliverable_SHA256_manifest.tsv','INDEPENDENT_INTERSECTION_STAGE_COMPLETED.txt'}:continue
        files.append({'relative_path':path.relative_to(ROOT).as_posix(),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    pd.DataFrame(files).to_csv(ROOT/'audit'/'Deliverable_SHA256_manifest.tsv',sep='\t',index=False)
    (ROOT/'audit'/'INDEPENDENT_INTERSECTION_STAGE_COMPLETED.txt').write_text('Completed independent reviewed human proteome prediction (20431 x 8 = 163448), 1/5/10-percent drug sets and twelve disease1881 intersections; all numerical checks, source structure comparison, figure review and document links completed. No single main threshold selected.\n',encoding='utf-8')
    print(json.dumps({'validation_passed':summary['n_pass'],'document_links_checked':len(links),'manifest_files':len(files),'status':'completed'}))
if __name__=='__main__':main()
