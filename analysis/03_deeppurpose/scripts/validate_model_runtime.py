"""Technical QA on public author examples; not a predictive-accuracy benchmark."""
import csv
import hashlib
import io
import json
import pathlib
import pickle
import time
import traceback
import numpy as np
import torch
from rdkit import Chem
from DeepPurpose import DTI
from DeepPurpose.utils import data_process

ROOT = pathlib.Path(__file__).resolve().parents[1]

class PlainConfigUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError(f"Configuration class not allowed: {module}.{name}")

def load_model(directory):
    config = PlainConfigUnpickler(io.BytesIO((directory / "config.pkl").read_bytes())).load()
    config.update(num_workers=0, batch_size=1 if config["drug_encoding"] == "MPNN" else 4)
    config["result_folder"] = str(ROOT / "public_validation" / "model_internal_results" / directory.name)
    pathlib.Path(config["result_folder"]).mkdir(parents=True, exist_ok=True)
    net = DTI.model_initialize(**config)
    state = torch.load(directory / "model.pt", map_location="cpu", weights_only=True)
    if next(iter(state)).startswith("module."):
        state = {key[7:]: value for key, value in state.items()}
    net.model.load_state_dict(state, strict=True)
    net.binary = config["binary"]
    return net

def main():
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    validation = ROOT / "public_validation"
    validation.mkdir(parents=True, exist_ok=True)
    lines = (validation / "dti_author_toy_data.txt").read_text(encoding="utf-8").splitlines()[:2]
    original_records = [line.split() for line in lines]
    assert len(original_records) == 2 and all(len(row) == 3 for row in original_records)
    records = original_records + original_records
    drugs = [row[0] for row in records]
    proteins = [row[1] for row in records]
    assert all(Chem.MolFromSmiles(s) is not None for s in drugs)
    assert all(set(s) <= set("ACDEFGHIKLMNPQRSTVWY") for s in proteins)
    qa_input = [{"row": i, "smiles": d, "smiles_length": len(d), "protein_length": len(p),
                 "protein_sha256": hashlib.sha256(p.encode()).hexdigest(), "CNN_drug_truncation": len(d)>100,
                 "CNN_protein_truncation": len(p)>1000} for i,(d,p) in enumerate(zip(drugs,proteins))]
    (validation / "QA_input_audit.json").write_text(json.dumps(qa_input, indent=2), encoding="utf-8")
    panel = json.loads((ROOT / "metadata" / "Model_panel_archive_audit.json").read_text(encoding="utf-8"))
    checks = []
    scores = []
    for record in panel:
        model_name = record["model"]
        start = time.monotonic()
        check = {"model": model_name, "endpoint": record["endpoint"], "test_role": "technical_runtime_only_not_accuracy_validation"}
        try:
            net = load_model(ROOT / "models" / model_name)
            drug_encoding = net.config["drug_encoding"]
            target_encoding = net.config["target_encoding"]
            frame = data_process(drugs, proteins, [-1]*4, drug_encoding, target_encoding, split_method="no_split")
            pred = np.asarray(net.predict(frame), dtype=float)
            repeated = np.asarray(net.predict(frame), dtype=float)
            singleton = np.asarray(net.predict(frame.iloc[:1].copy()), dtype=float)
            assert pred.shape == (4,) and np.isfinite(pred).all()
            repeat_difference = float(np.max(np.abs(pred-repeated)))
            singleton_difference = float(abs(singleton[0]-pred[0]))
            duplicate_difference = float(np.max(np.abs(pred[:2]-pred[2:])))
            assert repeat_difference <= 1e-5
            assert singleton_difference <= 1e-4
            assert duplicate_difference <= 1e-4
            check.update(status="technical_checks_pass", strict_state_dict_loaded=True, n_scores=4,
                         repeated_prediction_max_abs_difference=repeat_difference,
                         batch_vs_singleton_max_abs_difference=singleton_difference,
                         duplicate_pair_max_abs_difference=duplicate_difference,
                         drug_encoding=drug_encoding, target_encoding=target_encoding,
                         inference_batch_size=net.config["batch_size"],
                         MPNN_accuracy_still_unvalidated=(drug_encoding=="MPNN"),
                         elapsed_seconds=round(time.monotonic()-start,3))
            scores.extend({"model":model_name,"public_example_row":i,"raw_score":float(value),
                           "score_scale":"predicted_p"+record["endpoint"],
                           "role":"technical_QA_not_Spiropachysine_A_prediction"} for i,value in enumerate(pred))
        except Exception as exc:
            check.update(status="failed", error=str(exc), traceback=traceback.format_exc())
        checks.append(check)
        (validation / "Runtime_checks.json").write_text(json.dumps(checks,indent=2), encoding="utf-8")
        print(json.dumps(check,ensure_ascii=True),flush=True)
    if scores:
        with (validation / "Public_example_scores.tsv").open("w",encoding="utf-8",newline="") as handle:
            writer=csv.DictWriter(handle,fieldnames=list(scores[0]),delimiter="\t")
            writer.writeheader();writer.writerows(scores)
    if not all(item["status"]=="technical_checks_pass" for item in checks):
        raise SystemExit(1)
    (validation / "TECHNICAL_RUNTIME_CHECKS_PASSED.txt").write_text("10 pretrained models loaded strictly; finite, repeatable and batch-consistent outputs on two public author pairs repeated to four rows. Not predictive accuracy validation; no study compound predictions.\n",encoding="utf-8")

if __name__ == "__main__":
    main()
