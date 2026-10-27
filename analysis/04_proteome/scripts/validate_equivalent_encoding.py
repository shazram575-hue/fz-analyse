"""Validate accelerated encoding against all previous native features and new public examples."""
import hashlib
import json
import pathlib
import time
import numpy as np
import pandas as pd
from DeepPurpose import utils
from fast_equivalent_encoding import fast_aac

ROOT=pathlib.Path(__file__).resolve().parents[1]
PRIOR=ROOT.parent/"DeepPurpose_预测_2026-10-24"

def main():
    out=ROOT/"audit";out.mkdir(parents=True,exist_ok=True)
    previous=pd.read_csv(PRIOR/"predictions"/"Protein_sequence_catalog.tsv",sep="\t")
    old_path=PRIOR/"inference_cache"/"AAC_features_float32.npy"
    metadata=json.loads((PRIOR/"inference_cache"/"AAC_features_metadata.json").read_text(encoding="utf-8"))
    assert hashlib.sha256(old_path.read_bytes()).hexdigest()==metadata["sha256"]
    native=np.load(old_path,mmap_mode="r",allow_pickle=False)
    started=time.monotonic()
    for index,sequence in enumerate(previous.protein_sequence):
        result=fast_aac(sequence)
        if not np.array_equal(result,native[index]):raise ValueError(f"Native cached AAC mismatch at row {index}")
        if index%500==0:print(json.dumps({"stage":"native_cache_equivalence","done":index+1,"total":len(previous)}),flush=True)
    background=pd.read_csv(ROOT/"background"/"All_reviewed_human_protein_background.tsv",sep="\t")
    ordered=background.sort_values("sequence_length")
    indices=np.linspace(0,len(ordered)-1,48,dtype=int)
    examples=list(ordered.iloc[indices].protein_sequence)
    examples += ["AAAAAA","ABABABABAB","ACACACACA","UUUXAXAXX","A","", "ACDEFGHIKLMNPQRSTVWY"*3]
    for index,sequence in enumerate(examples):
        expected=utils.target2aac(sequence).astype(np.float32)
        if not np.array_equal(fast_aac(sequence),expected):raise ValueError(f"Fresh native AAC mismatch at example {index}")
    categories=list(utils.enc_protein.categories_[0]);char_map={char:index for index,char in enumerate(categories)}
    unknown=char_map["?"]
    for sequence in examples:
        normalized=utils.trans_protein(sequence)
        indices=np.array([char_map.get(c,unknown) for c in normalized],dtype=np.uint8)
        optimized=np.eye(26,dtype=np.float32)[indices].T
        expected=utils.protein_2_embed(normalized)
        if not np.array_equal(optimized,expected):raise ValueError("CNN one-hot category order mismatch")
    report={"all1993_previous_native_AAC_vectors_exactly_equal":True,"previous_AAC_elements_compared":int(native.size),
            "fresh_native_AAC_examples":len(examples),"fresh_native_AAC_exactly_equal":True,
            "CNN_onehot_category_order_exactly_equal":True,"nonoverlap_semantics":"leftmost greedy per literal kmer, matching native str.count/re.findall",
            "feature_order_and_rounding":"same native AALetter/product order and Python round; float32 is native torch input precision",
            "changes_to_model_weights":False,"elapsed_seconds":round(time.monotonic()-started,3)}
    (out/"Equivalent_encoding_validation.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report),flush=True)

if __name__=="__main__":main()
