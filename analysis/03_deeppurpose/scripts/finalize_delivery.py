"""Check local delivery links and record final artifact checksums after visual QA."""
import hashlib
import json
import pathlib
import re
import pandas as pd

ROOT=pathlib.Path(__file__).resolve().parents[1]

def main():
    audit=json.loads((ROOT/"audit"/"Final_validation_summary.json").read_text(encoding="utf-8"))
    assert audit["n_checks"]==audit["n_pass"]==27
    assert (ROOT/"audit"/"Figure_visual_QA.md").exists()
    links=[]
    for path in [ROOT/"DeepPurpose_结果报告.md",ROOT/"README_当前状态.md",ROOT/"nomenclature_review"/"README_命名核对与采用决定.md"]:
        for target in re.findall(r"\]\(([^)]+)\)",path.read_text(encoding="utf-8")):
            if target.startswith(("https://","http://")):continue
            target_path=(path.parent/target).resolve()
            links.append({"document":str(path.relative_to(ROOT)),"target":target,"exists":target_path.exists()})
    pd.DataFrame(links).to_csv(ROOT/"audit"/"Final_document_links.tsv",sep="\t",index=False)
    assert all(item["exists"] for item in links)
    paths=[ROOT/"DeepPurpose_结果报告.md",ROOT/"README_当前状态.md",
           ROOT/"predictions"/"Frozen_prediction_design.json",ROOT/"predictions"/"All_model_run_summary.json",
           ROOT/"structure_verification_CID209254"/"Adopted_computational_input.json"]
    for directory in ["results","figures","audit","scripts"]:
        paths.extend(p for p in (ROOT/directory).glob("*") if p.is_file() and p.name not in ["Final_delivery_SHA256.tsv","PREDICTION_STAGE_COMPLETED.txt"])
    paths=sorted(set(paths))
    manifest=[{"relative_path":str(path.relative_to(ROOT)),"bytes":path.stat().st_size,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()} for path in paths]
    pd.DataFrame(manifest).to_csv(ROOT/"audit"/"Final_delivery_SHA256.tsv",sep="\t",index=False)
    marker=("Completed 2026-10-24: user-approved fixed CID209254 structure, main1881 plus76 diagnostic-only genes, "
            "1993 protein pairs, 10 pretrained models, 19930 raw scores. 27 consistency checks passed; figures reviewed; "
            "frozen stage1 files preserved. Experimental accuracy, direct binding and final target selection are not established.\n")
    (ROOT/"audit"/"PREDICTION_STAGE_COMPLETED.txt").write_text(marker,encoding="utf-8")
    print(json.dumps({"delivery_files":len(manifest),"local_document_links":len(links),"all_links_exist":True,"completed_marker":str(ROOT/"audit"/"PREDICTION_STAGE_COMPLETED.txt")}),flush=True)

if __name__=="__main__":main()
