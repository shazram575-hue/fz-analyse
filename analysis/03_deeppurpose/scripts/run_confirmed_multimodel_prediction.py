"""Local inference only, with frozen structure, cached encodings and explicit model roles."""
import contextlib
import hashlib
import json
import pathlib
import time
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from DeepPurpose import utils
from validate_model_runtime import load_model

ROOT=pathlib.Path(__file__).resolve().parents[1]

def sha256(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    cache=ROOT/"inference_cache"; cache.mkdir(parents=True,exist_ok=True)
    out=ROOT/"predictions";out.mkdir(parents=True,exist_ok=True)
    logs=ROOT/"logs";logs.mkdir(parents=True,exist_ok=True)
    verify=json.loads((ROOT/"structure_verification_CID209254"/"Verification_report.json").read_text(encoding="utf-8"))
    original=verify["cleaned_user_smiles"]
    molecule=Chem.MolFromSmiles(original)
    smiles=Chem.MolToSmiles(molecule,isomericSmiles=True,canonical=True)
    assert Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))=="SLGWGPQWJRVPAD-PZOAWPBASA-N"
    assert len(smiles)<=utils.MAX_SEQ_DRUG
    main_table=pd.read_csv(ROOT/"main_input"/"Resolved_main_protein_sequences.tsv",sep="\t",keep_default_na=False)
    extra_table=pd.read_csv(ROOT/"main_input"/"Resolved_diagnostic_extra_protein_sequences.tsv",sep="\t",keep_default_na=False)
    catalog=pd.concat([main_table,extra_table],ignore_index=True).sort_values(["input_scope","canonical_id","protein_input_id","sequence_sha256"]).reset_index(drop=True)
    assert main_table.canonical_id.nunique()==1881 and extra_table.canonical_id.nunique()==76
    assert not catalog.duplicated(["canonical_id","sequence_sha256"]).any()
    assert catalog.protein_sequence.str.len().eq(catalog.sequence_length).all()
    assert all(hashlib.sha256(seq.encode()).hexdigest()==key for seq,key in zip(catalog.protein_sequence,catalog.sequence_sha256))
    catalog["prediction_row_id"]=[f"pair_{index:05d}" for index in range(len(catalog))]
    catalog.to_csv(out/"Protein_sequence_catalog.tsv",sep="\t",index=False,na_rep="NA")
    sequences=list(catalog.protein_sequence)
    sequence_order_hash=hashlib.sha256("\n".join(catalog.sequence_sha256).encode()).hexdigest()
    panel=json.loads((ROOT/"metadata"/"Model_panel_archive_audit.json").read_text(encoding="utf-8"))
    design={"scope":"1881_main_genes_plus76_diagnostic_sensitivity_only","structure_CID":209254,
            "structure_InChIKey":verify["InChIKey"],"original_PubChem_smiles":original,
            "actual_model_input_smiles":smiles,"original_smiles_length":len(original),"actual_smiles_length":len(smiles),
            "smiles_standardization":"RDKit canonical isomeric SMILES, same InChIKey; avoids truncating 101-character source string at CNN limit100",
            "user_authorization":"2026-10-24: 接受，按 CID 209254 计算并保留命名依据与限制",
            "n_protein_pairs":len(catalog),"sequence_order_sha256":sequence_order_hash,
            "primary_encodings":["CNN_CNN","Morgan_CNN","Morgan_AAC","Daylight_AAC"],
            "MPNN_role":"sensitivity_only; original batch dependency; single-pair inference; accuracy not restored by repeatability check",
            "endpoints":["Kd","IC50"],"scores":"raw predicted log-scale per trained endpoint; not calibrated probabilities or measured affinities",
            "gene_aggregation":"median predicted score across all retained unique protein sequences; per-protein scores and ranges retained",
            "ensemble_aggregation":"separate for each endpoint: median of four per-model gene percentile ranks, no raw endpoint mixing",
            "missing_amino_acids":"native author encoders without substitution; U/X explicitly flagged",
            "sequence_truncation":"native CNN first1000 residues, flagged; AAC uses full sequence",
            "training":False,"fine_tuning":False,"candidate_list_external_transfer":False,
            "topN_selection":None,"binding_threshold":None}
    (out/"Frozen_prediction_design.json").write_text(json.dumps(design,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"stage":"input_frozen","n_pairs":len(catalog),"smiles_length":len(smiles)}),flush=True)
    aac_path=cache/"AAC_features_float32.npy"
    aac_meta=cache/"AAC_features_metadata.json"
    if aac_path.exists() and aac_meta.exists():
        metadata=json.loads(aac_meta.read_text(encoding="utf-8"))
        assert metadata["sequence_order_sha256"]==sequence_order_hash
        aac=np.load(aac_path,allow_pickle=False)
    else:
        aac=np.empty((len(sequences),8420),dtype=np.float32)
        started=time.monotonic()
        for index,sequence in enumerate(sequences):
            feature=utils.target2aac(sequence).astype(np.float32)
            if feature.shape!=(8420,) or not np.isfinite(feature).all() or not np.any(feature):
                raise ValueError(f"Invalid/zero AAC features for {catalog.prediction_row_id.iloc[index]}")
            aac[index]=feature
            if index%100==0 or index==len(sequences)-1:
                print(json.dumps({"stage":"AAC_encoding","done":index+1,"total":len(sequences),"seconds":round(time.monotonic()-started,1)}),flush=True)
        np.save(aac_path,aac,allow_pickle=False)
        aac_meta.write_text(json.dumps({"sequence_order_sha256":sequence_order_hash,"dtype":"float32","shape":list(aac.shape),"encoder":"DeepPurpose0.1.5 target2aac native full sequence","sha256":sha256(aac_path)},indent=2),encoding="utf-8")
    assert aac.shape==(len(sequences),8420) and np.isfinite(aac).all()
    cnn=[utils.trans_protein(sequence) for sequence in sequences]
    model_results=[]
    with (logs/"inference_library_stdout.log").open("a",encoding="utf-8") as lib_log:
        for item in panel:
            name=item["model"]
            directory=ROOT/"models"/name
            assert sha256(directory/"model.pt")==item["model_pt_sha256"]
            model_path=out/(name+"_protein_scores.tsv")
            manifest_path=out/(name+"_run.json")
            if model_path.exists() and manifest_path.exists():
                old=json.loads(manifest_path.read_text(encoding="utf-8"))
                if old.get("sequence_order_sha256")==sequence_order_hash and old.get("actual_model_input_smiles")==smiles and old.get("model_pt_sha256")==item["model_pt_sha256"] and old.get("status")=="complete":
                    model_results.append(old)
                    print(json.dumps({"stage":"reuse_complete_model","model":name}),flush=True)
                    continue
            started=time.monotonic()
            net=load_model(directory)
            is_mpnn=net.config["drug_encoding"]=="MPNN"
            net.config["batch_size"]=1 if is_mpnn else 32
            with contextlib.redirect_stdout(lib_log):
                drug=utils.smiles2morgan(smiles) if net.drug_encoding=="Morgan" else utils.smiles2daylight(smiles) if net.drug_encoding=="Daylight" else utils.trans_drug(smiles) if net.drug_encoding=="CNN" else utils.smiles2mpnnfeature(smiles)
            if net.drug_encoding in ["Morgan","Daylight"]:
                assert np.isfinite(drug).all() and np.any(drug)
            target_features=list(aac) if net.target_encoding=="AAC" else cnn
            frame=pd.DataFrame({"SMILES":[smiles]*len(catalog),"Target Sequence":sequences,"Label":[-1]*len(catalog),
                                "drug_encoding":[drug]*len(catalog),"target_encoding":target_features})
            predictions=[]
            for start in range(0,len(frame),128):
                batch=frame.iloc[start:start+128].reset_index(drop=True)
                with contextlib.redirect_stdout(lib_log):values=net.predict(batch)
                predictions.extend(float(value) for value in values)
                if start%512==0 or start+128>=len(frame):
                    print(json.dumps({"stage":"predicting","model":name,"done":min(start+128,len(frame)),"total":len(frame)}),flush=True)
            values=np.asarray(predictions)
            assert len(values)==len(catalog) and np.isfinite(values).all()
            probe=[0,len(frame)//2,len(frame)-1]
            singleton=[]
            for index in probe:
                with contextlib.redirect_stdout(lib_log):singleton.extend(net.predict(frame.iloc[[index]].reset_index(drop=True)))
            max_difference=float(np.max(np.abs(values[probe]-np.asarray(singleton))))
            assert max_difference<=1e-4
            scores=catalog[["prediction_row_id","input_scope","canonical_id","gene_symbol","protein_input_id","sequence_sha256"]].copy()
            scores["model"]=name;scores["endpoint"]=item["endpoint"]
            scores["role"]="sensitivity_MPNN" if is_mpnn else "primary_encoding_ensemble"
            scores["predicted_p_value"]=values
            scores.to_csv(model_path,sep="\t",index=False)
            result={"model":name,"endpoint":item["endpoint"],"role":"sensitivity_MPNN" if is_mpnn else "primary_encoding_ensemble",
                    "status":"complete","n_scores":len(values),"all_finite":True,"actual_model_input_smiles":smiles,
                    "sequence_order_sha256":sequence_order_hash,"model_pt_sha256":item["model_pt_sha256"],
                    "inference_batch_size":net.config["batch_size"],"actual_input_singleton_check_max_abs_difference":max_difference,
                    "score_file_sha256":sha256(model_path),"elapsed_seconds":round(time.monotonic()-started,3),
                    "predictive_accuracy_not_established":True}
            manifest_path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
            model_results.append(result)
            print(json.dumps({"stage":"model_complete",**{k:result[k] for k in ["model","n_scores","elapsed_seconds","actual_input_singleton_check_max_abs_difference"]}}),flush=True)
            del net,frame
    assert len(model_results)==10
    (out/"All_model_run_summary.json").write_text(json.dumps(model_results,indent=2),encoding="utf-8")
    (out/"RAW_PREDICTIONS_COMPLETE.txt").write_text("Ten pretrained models completed local inference on fixed CID209254 structure and all confirmed main/sensitivity protein sequences. Results remain computational predictions; analysis and audit are separate.\n",encoding="utf-8")

if __name__=="__main__":main()
