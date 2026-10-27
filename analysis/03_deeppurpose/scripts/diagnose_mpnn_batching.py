"""Record batch-dependent outputs of unmodified author MPNN implementation."""
import json
import pathlib
import numpy as np
import torch
from DeepPurpose.utils import data_process
from validate_model_runtime import load_model

ROOT = pathlib.Path(__file__).resolve().parents[1]

def main():
    torch.set_num_threads(2)
    pairs = [line.split() for line in (ROOT / "public_validation" / "dti_author_toy_data.txt").read_text(encoding="utf-8").splitlines()]
    pairs = pairs + pairs
    report = []
    for name in ["MPNN_CNN_BindingDB", "MPNN_CNN_BindingDB_IC50"]:
        net = load_model(ROOT / "models" / name)
        frame = data_process([p[0] for p in pairs], [p[1] for p in pairs], [-1]*len(pairs), "MPNN", "CNN", split_method="no_split")
        predictions = {}
        for batch_size in [1,2,4]:
            net.config["batch_size"] = batch_size
            predictions[str(batch_size)] = [float(x) for x in net.predict(frame)]
        difference = float(np.max(np.abs(np.array(predictions["1"])-np.array(predictions["4"]))))
        report.append({"model":name,"batch_predictions":predictions,"batch1_vs_batch4_max_abs_difference":difference,
                       "batch1_duplicate_difference":float(np.max(np.abs(np.array(predictions["1"][:2])-np.array(predictions["1"][2:])))),
                       "implementation":"unmodified_DeepPurpose_0.1.5","decision":"use_batch_size1_for_MPNN_to_avoid_cross_pair_batch_dependency"})
    (ROOT / "public_validation" / "MPNN_batch_dependency_diagnostic.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report),flush=True)

if __name__ == "__main__":
    main()
