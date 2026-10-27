"""Independent reconstruction from per-type significant tables and preservation audit."""
import hashlib
import json
import pathlib
import pandas as pd

ROOT=pathlib.Path(__file__).resolve().parents[1]
WORKSPACE=ROOT.parents[1]
STAGE1=WORKSPACE/"outputs"/"第一阶段_2026-10-23"
FOUNDATION=WORKSPACE/"outputs"/"单细胞基础补全_2026-10-23"

def main():
    out=ROOT/"audit"
    out.mkdir(parents=True,exist_ok=True)
    bulk=pd.read_csv(STAGE1/"results"/"Bulk_exploratory_DE_protein_coding.tsv",sep="\t")
    biological=["Endothelial","Epithelial","Fibroblast","Immune","Mast","Perivascular_SMC"]
    diagnostic=["Mixed_epithelial_immune","Unresolved"]
    definitions={"celllevel_6biological":(biological,"celllevel_DE_Bonf005_absFC05"),
                 "celllevel_all8_including_diagnostic":(biological+diagnostic,"celllevel_DE_Bonf005_absFC05"),
                 "pseudobulk_4biological_estimable":(biological[:4],"pseudobulk_DE_BH005_absFC05"),
                 "pseudobulk_all6_estimable_including_diagnostic":(biological[:4]+diagnostic,"pseudobulk_DE_BH005_absFC05")}
    checks=[]
    for name,(groups,suffix) in definitions.items():
        ids=set()
        for group in groups:
            table=pd.read_csv(FOUNDATION/"results"/"within_celltype"/(group+"_"+suffix+".tsv"),sep="\t")
            ids.update(table.loc[table.gene_type.eq("protein_coding"),"canonical_id"])
        ids &= set(bulk.ensembl_stable)
        published=pd.read_csv(ROOT/"candidate_scope_review"/(name+"_genes.tsv"),sep="\t")
        checks.append({"check":name+"_independent_reconstruction","pass":ids==set(published.canonical_id),"n":len(ids)})
    frozen=pd.read_csv(FOUNDATION/"audit"/"Frozen_stage1_before.tsv",sep="\t")
    preserved=[]
    for row in frozen.itertuples(index=False):
        path=STAGE1/pathlib.Path(row.relative_path)
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        preserved.append({"relative_path":row.relative_path,"unchanged":actual.lower()==row.sha256.lower(),"current_sha256":actual})
    pd.DataFrame(preserved).to_csv(out/"Frozen_stage1_preservation.tsv",sep="\t",index=False)
    checks.append({"check":"all_7_frozen_stage1_files_unchanged","pass":all(r["unchanged"] for r in preserved),"n":len(preserved)})
    runtimes=json.loads((ROOT/"public_validation"/"Runtime_checks.json").read_text(encoding="utf-8"))
    checks.append({"check":"all_10_models_technical_runtime_only","pass":len(runtimes)==10 and all(r["status"]=="technical_checks_pass" for r in runtimes),"n":len(runtimes)})
    public=json.loads((ROOT/"public_human_proteome"/"Download_audit.json").read_text(encoding="utf-8"))
    checks.append({"check":"candidate_list_not_transmitted_in_public_download","pass":public["candidate_ids_transmitted"] is False,"n":public["n_entries"]})
    pd.DataFrame(checks).to_csv(out/"Preparation_validation.tsv",sep="\t",index=False)
    print(json.dumps(checks),flush=True)
    if not all(item["pass"] for item in checks):
        raise SystemExit(1)

if __name__=="__main__":
    main()
