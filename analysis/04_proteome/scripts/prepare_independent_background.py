"""Freeze all reviewed human proteins independently of any disease-selected gene list."""
import gzip
import hashlib
import json
import pathlib
from collections import defaultdict
import pandas as pd

ROOT=pathlib.Path(__file__).resolve().parents[1]
WORKSPACE=ROOT.parents[1]
PRIOR=ROOT.parent/"DeepPurpose_预测_2026-10-24"

def main():
    out=ROOT/"background";out.mkdir(parents=True,exist_ok=True)
    path=PRIOR/"public_human_proteome"/"Reviewed_human_UniProt.json.gz"
    public=json.loads((PRIOR/"public_human_proteome"/"Download_audit.json").read_text(encoding="utf-8"))
    assert hashlib.sha256(path.read_bytes()).hexdigest()==public["sha256"]
    records=json.load(gzip.open(path,"rt",encoding="utf-8"))["results"]
    annotation=pd.read_csv(WORKSPACE/"outputs"/"第一阶段_2026-10-23"/"metadata"/"genome_reference"/"GENCODE_v36_gene_annotation.tsv",sep="\t")
    annotation=annotation[annotation.gene_type.eq("protein_coding")]
    reference_ids=set(annotation.ensembl_stable)
    symbols=annotation.groupby("gene_name").ensembl_stable.apply(lambda x:sorted(set(x))).to_dict()
    rows=[];mapping=[]
    for protein in records:
        assert protein["organism"]["taxonId"]==9606 and "Swiss-Prot" in protein["entryType"]
        accession=protein["primaryAccession"];sequence=protein["sequence"]["value"]
        names=[gene.get("geneName",{}).get("value","") for gene in protein.get("genes",[])]
        ids=set()
        for item in protein.get("uniProtKBCrossReferences",[]):
            if item.get("database")=="Ensembl":
                ids.update(prop["value"].split(".")[0] for prop in item.get("properties",[]) if prop.get("key")=="GeneId")
        exact=ids & reference_ids
        for gene in sorted(exact):mapping.append({"uniprot_accession":accession,"ensembl_stable":gene,"mapping_method":"exact_Ensembl_GeneId"})
        rescued=set()
        for name in names:
            candidates=symbols.get(name,[])
            if len(candidates)==1 and candidates[0] not in exact:
                rescued.add(candidates[0]);mapping.append({"uniprot_accession":accession,"ensembl_stable":candidates[0],"mapping_method":"unique_primary_gene_symbol_in_full_GENCODEv36_reference"})
        rows.append({"uniprot_accession":accession,"primary_gene_names":";".join(names),"protein_sequence":sequence,
                     "sequence_length":len(sequence),"sequence_sha256":hashlib.sha256(sequence.encode()).hexdigest(),
                     "GENCODEv36_gene_ids":";".join(sorted(exact|rescued)),"direct_current_Ensembl_ids":";".join(sorted(ids)),
                     "CNN_truncates_after1000":len(sequence)>1000,"nonstandard_amino_acids":";".join(sorted(set(sequence)-set("ACDEFGHIKLMNPQRSTVWY")))})
    table=pd.DataFrame(rows).sort_values("uniprot_accession").reset_index(drop=True)
    assert len(table)==20431 and table.uniprot_accession.is_unique
    table.to_csv(out/"All_reviewed_human_protein_background.tsv",sep="\t",index=False)
    pd.DataFrame(mapping).drop_duplicates().sort_values(["uniprot_accession","ensembl_stable","mapping_method"]).to_csv(out/"Protein_to_full_reference_gene_mapping.tsv",sep="\t",index=False)
    protocol={"date":"2026-10-24","drug_side_background_selection":"all reviewed human UniProt canonical entries in archived release2026_03",
              "n_proteins":len(table),"n_unique_sequences":table.sequence_sha256.nunique(),
              "disease_gene_list_used_for_background_selection":False,"disease_gene_list_used_for_gene_mapping":False,
              "mapping_reference":"entire GENCODEv36 protein_coding annotation, not the DE/candidate set",
              "structure_CID":209254,"structure_InChIKey":"SLGWGPQWJRVPAD-PZOAWPBASA-N",
              "identity_confirmation":"published A Figure1 reconstructed from nine wedge/dash bonds; all9 centers and full InChIKey match CID209254",
              "primary_models":"4 non-MPNN encodings each for Kd and IC50; MPNN excluded because pre-existing batch defect",
              "selection_unit":"reviewed protein accession; each endpoint uses median of4 per-model percentiles in entire background",
              "candidate_thresholds_percent":[1,5,10],"boundary_ties":"include tied proteins at the cutoff using minimum rank <= ceil(N*percent)",
              "intersection_order":"select drug proteins without disease information; map to reference gene sets; compare endpoint gene union/intersection; intersect disease1881 last",
              "single_primary_threshold_selected":False,"user_choice":"完整报告1%/5%/10%比较，暂不指定唯一主阈值",
              "background_sha256":hashlib.sha256((out/"All_reviewed_human_protein_background.tsv").read_bytes()).hexdigest(),
              "mapping_sha256":hashlib.sha256((out/"Protein_to_full_reference_gene_mapping.tsv").read_bytes()).hexdigest()}
    (out/"Frozen_independent_drug_side_protocol.json").write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(protocol,ensure_ascii=True),flush=True)

if __name__=="__main__":main()
