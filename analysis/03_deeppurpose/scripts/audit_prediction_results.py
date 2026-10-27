"""Independent score/consensus/provenance validation, distinct from model accuracy."""
import hashlib
import json
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from rdkit import Chem

ROOT=pathlib.Path(__file__).resolve().parents[1]
WORKSPACE=ROOT.parents[1]

def sha256(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    out=ROOT/"audit";out.mkdir(parents=True,exist_ok=True)
    results=ROOT/"results";prediction=ROOT/"predictions"
    checks=[]
    def check(name,value,detail=""):
        checks.append({"check":name,"pass":bool(value),"detail":str(detail)})
    design=json.loads((prediction/"Frozen_prediction_design.json").read_text(encoding="utf-8"))
    main_genes=pd.read_csv(ROOT/"main_input"/"Disease_matching_candidates_1881.tsv",sep="\t")
    extra_genes=pd.read_csv(ROOT/"main_input"/"Additional_diagnostic_only_76genes.tsv",sep="\t")
    check("confirmed_main1881_unique",len(main_genes)==1881 and main_genes.canonical_id.is_unique)
    check("diagnostic_extra76_disjoint",len(extra_genes)==76 and not set(extra_genes.canonical_id)&set(main_genes.canonical_id))
    catalog=pd.read_csv(prediction/"Protein_sequence_catalog.tsv",sep="\t",keep_default_na=False)
    raw=pd.read_csv(results/"All_protein_model_predictions.tsv",sep="\t")
    check("1993_protein_pairs",len(catalog)==1993 and catalog.prediction_row_id.is_unique)
    check("all1881_main_genes_have_scores",set(raw.loc[raw.input_scope.eq("main"),"canonical_id"])==set(main_genes.canonical_id))
    check("all76_extra_genes_have_scores",set(raw.loc[raw.input_scope.eq("diagnostic_extra"),"canonical_id"])==set(extra_genes.canonical_id))
    check("19930_complete_finite_raw_scores",len(raw)==19930 and np.isfinite(raw.predicted_log_affinity).all() and not raw.duplicated(["prediction_row_id","model"]).any())
    grouped=raw.groupby("model").prediction_row_id.nunique()
    check("10_models_cover_same_pairs",len(grouped)==10 and grouped.eq(1993).all())
    check("8_primary_2_MPNN_sensitivity",raw.loc[raw.role.eq("primary_encoding_ensemble"),"model"].nunique()==8 and raw.loc[raw.role.eq("sensitivity_MPNN"),"model"].nunique()==2)
    source_keys=dict(zip(catalog.prediction_row_id,catalog.sequence_sha256))
    check("raw_pair_sequence_identity",all(source_keys[pair]==key for pair,key in zip(raw.prediction_row_id,raw.sequence_sha256)))
    check("actual_structure_CID209254",Chem.MolToInchiKey(Chem.MolFromSmiles(design["actual_model_input_smiles"]))=="SLGWGPQWJRVPAD-PZOAWPBASA-N")
    check("drug_SMILES_avoids_CNN_truncation",len(design["actual_model_input_smiles"])<=100, len(design["actual_model_input_smiles"]))
    centers=Chem.FindMolChiralCenters(Chem.MolFromSmiles(design["actual_model_input_smiles"]),includeUnassigned=True,useLegacyImplementation=False)
    check("all9_drug_centers_assigned",len(centers)==9 and all(cip!="?" for _,cip in centers))
    runs=json.loads((prediction/"All_model_run_summary.json").read_text(encoding="utf-8"))
    check("all_weight_and_score_hashes_match",all(sha256(ROOT/"models"/r["model"]/"model.pt")==r["model_pt_sha256"] and sha256(prediction/(r["model"]+"_protein_scores.tsv"))==r["score_file_sha256"] for r in runs))
    check("actual_input_singleton_checks",all(r["actual_input_singleton_check_max_abs_difference"]<=1e-4 for r in runs))
    check("MPNN_single_pair_inference",all(r["inference_batch_size"]==1 for r in runs if r["role"]=="sensitivity_MPNN"))
    for endpoint in ["Kd","IC50"]:
        selected=raw[raw.input_scope.eq("main") & raw.endpoint.eq(endpoint) & raw.role.eq("primary_encoding_ensemble")]
        matrix=selected.groupby(["canonical_id","model"]).predicted_log_affinity.median().unstack().sort_index()
        percentile=rankdata(matrix.to_numpy(),method="average",axis=0)/len(matrix)
        consensus=np.median(percentile,axis=1)
        ranks=rankdata(-consensus,method="min").astype(int)
        table=pd.read_csv(results/("Main_"+endpoint+"_consensus_ranked.tsv"),sep="\t").set_index("canonical_id").reindex(matrix.index)
        check(endpoint+"_4model1881_independent_consensus",matrix.shape==(1881,4) and np.allclose(consensus,table.consensus_percentile,rtol=0,atol=1e-12))
        check(endpoint+"_independent_rank_matches",np.array_equal(ranks,table.consensus_rank.to_numpy()))
        extra=pd.read_csv(results/("Diagnostic_extended1957_"+endpoint+"_ranked.tsv"),sep="\t")
        check(endpoint+"_diagnostic_scope1957_complete",len(extra)==1957 and extra.canonical_id.nunique()==1957)
    frozen=pd.read_csv(WORKSPACE/"outputs"/"单细胞基础补全_2026-10-23"/"audit"/"Frozen_stage1_before.tsv",sep="\t")
    preserved=[]
    for row in frozen.itertuples(index=False):
        path=WORKSPACE/"outputs"/"第一阶段_2026-10-23"/row.relative_path
        actual=sha256(path)
        preserved.append({"relative_path":row.relative_path,"pass":actual.lower()==row.sha256.lower(),"sha256":actual})
    pd.DataFrame(preserved).to_csv(out/"Final_stage1_preservation.tsv",sep="\t",index=False)
    check("all7_frozen_stage1_preserved",all(item["pass"] for item in preserved))
    scope=json.loads((ROOT/"main_input"/"Confirmed_scope_and_provenance.json").read_text(encoding="utf-8"))
    check("confirmed_disease_input_files_unchanged",all(sha256(ROOT/"main_input"/name)==expected for name,expected in scope["files_sha256"].items()))
    public=json.loads((ROOT/"public_human_proteome"/"Download_audit.json").read_text(encoding="utf-8"))
    check("candidate_list_not_transmitted",public["candidate_ids_transmitted"] is False)
    for stem in ["Model_rank_agreement","Kd_vs_IC50_consensus","Top20_encoding_agreement"]:
        check(stem+"_PNG_SVG_present",all((ROOT/"figures"/(stem+"."+ext)).stat().st_size>1000 for ext in ["png","svg"]))
    pd.DataFrame(checks).to_csv(out/"Final_prediction_validation.tsv",sep="\t",index=False)
    summary={"n_checks":len(checks),"n_pass":sum(item["pass"] for item in checks),"failed":[item for item in checks if not item["pass"]],
             "validation_scope":"data/score/model-configuration consistency and frozen-result preservation, not predictive accuracy or experimental validation"}
    (out/"Final_validation_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary),flush=True)
    if not all(item["pass"] for item in checks):raise SystemExit(1)

if __name__=="__main__":main()
