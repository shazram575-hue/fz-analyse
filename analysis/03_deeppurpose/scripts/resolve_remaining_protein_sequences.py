"""Local sequence reconciliation using validated public GENCODE and UniProt files."""
import csv
import gzip
import hashlib
import json
import pathlib
from collections import defaultdict
import pandas as pd

ROOT=pathlib.Path(__file__).resolve().parents[1]

def fasta_records(path):
    with gzip.open(path,"rt",encoding="utf-8") as handle:
        header=None; parts=[]
        for line in handle:
            line=line.strip()
            if line.startswith(">"):
                if header is not None:yield header,"".join(parts)
                header=line[1:];parts=[]
            else:parts.append(line)
        if header is not None:yield header,"".join(parts)

def main():
    out=ROOT/"protein_mapping_preparation"
    genes=pd.read_csv(ROOT/"candidate_scope_review"/"Union_for_preparation_only.tsv",sep="\t")
    direct=pd.read_csv(out/"Direct_mapping_with_symbol_checks.tsv",sep="\t",keep_default_na=False)
    missing=set(genes.canonical_id)-set(direct.canonical_id)
    reference=ROOT/"public_human_proteome"/"gencode.v36.pc_translations.fa.gz"
    assert hashlib.md5(reference.read_bytes()).hexdigest()=="20e8a4eb0f54be8e7e85b293021da08a"
    peptides=defaultdict(lambda:defaultdict(list))
    for header,sequence in fasta_records(reference):
        fields=header.split("|")
        gene=fields[2].split(".")[0]
        if gene in missing:
            sequence=sequence.removesuffix("*")
            peptides[gene][sequence].append(fields[0])
    entries=json.load(gzip.open(ROOT/"public_human_proteome"/"Reviewed_human_UniProt.json.gz","rt",encoding="utf-8"))["results"]
    by_sequence=defaultdict(list)
    for protein in entries:
        by_sequence[protein["sequence"]["value"]].append(protein)
    rows=[]; sequence_audit=[]
    name_map=dict(zip(genes.canonical_id,genes.gene_symbol))
    for gene in sorted(missing):
        exact=defaultdict(lambda:{"protein":None,"translations":set()})
        for sequence,translations in peptides[gene].items():
            for protein in by_sequence.get(sequence,[]):
                key=protein["primaryAccession"]
                exact[key]["protein"]=protein;exact[key]["translations"].update(translations)
        if exact:
            for accession,hit in exact.items():
                protein=hit["protein"];sequence=protein["sequence"]["value"]
                symbols=[g.get("geneName",{}).get("value","") for g in protein.get("genes",[])]
                rows.append({"canonical_id":gene,"uniprot_accession":accession,"uniprot_gene_names":";".join(symbols),
                             "organism_id":9606,"entry_type":protein.get("entryType",""),"sequence_length":len(sequence),
                             "protein_sequence":sequence,"sequence_sha256":hashlib.sha256(sequence.encode()).hexdigest(),
                             "CNN_truncates_after1000":len(sequence)>1000,
                             "nonstandard_amino_acids":";".join(sorted(set(sequence)-set("ACDEFGHIKLMNPQRSTVWY"))),
                             "mapping_method":"GENCODEv36_Ensembl_translation_exact_sequence_to_reviewed_human_UniProt",
                             "scope_status":"sequence_preparation_not_drug_prediction","gene_symbol":name_map[gene],
                             "primary_symbol_exact_match":name_map[gene] in symbols,
                             "mapping_tier":"exact_GENCODEv36_protein_sequence","gencode_translation_ids":";".join(sorted(hit["translations"]))})
            status="rescued_exact_reviewed_protein_sequence"
        elif peptides[gene]:
            for index,(sequence,translations) in enumerate(sorted(peptides[gene].items()),1):
                rows.append({"canonical_id":gene,"uniprot_accession":"","uniprot_gene_names":"","organism_id":9606,
                             "entry_type":"GENCODEv36_translation_no_exact_reviewed_UniProt_match","sequence_length":len(sequence),
                             "protein_sequence":sequence,"sequence_sha256":hashlib.sha256(sequence.encode()).hexdigest(),
                             "CNN_truncates_after1000":len(sequence)>1000,
                             "nonstandard_amino_acids":";".join(sorted(set(sequence)-set("ACDEFGHIKLMNPQRSTVWY"))),
                             "mapping_method":"GENCODEv36_translation_fallback_all_unique_protein_sequences_preserved",
                             "scope_status":"sequence_preparation_not_drug_prediction","gene_symbol":name_map[gene],
                             "primary_symbol_exact_match":False,"mapping_tier":"GENCODEv36_translation_only",
                             "gencode_translation_ids":";".join(sorted(translations))})
            status="GENCODE_translation_only_no_reviewed_match"
        else:status="no_GENCODEv36_pc_translation_found"
        sequence_audit.append({"canonical_id":gene,"gene_symbol":name_map[gene],"status":status,
                               "n_unique_GENCODE_sequences":len(peptides[gene]),"n_exact_reviewed_matches":len(exact)})
    rescued=pd.DataFrame(rows)
    rescued.to_csv(out/"Sequence_reconciled_additional_mappings.tsv",sep="\t",index=False,na_rep="NA")
    pd.DataFrame(sequence_audit).to_csv(out/"Remaining_genes_sequence_reconciliation.tsv",sep="\t",index=False)
    direct["gencode_translation_ids"]=""
    all_rows=pd.concat([direct,rescued],ignore_index=True)
    all_rows["protein_input_id"]=[accession if accession else "GENCODE36_"+hash_[:16] for accession,hash_ in zip(all_rows.uniprot_accession,all_rows.sequence_sha256)]
    all_rows=all_rows.drop_duplicates(["canonical_id","sequence_sha256"])
    all_rows.to_csv(out/"Resolved_all_candidate_protein_sequences.tsv",sep="\t",index=False,na_rep="NA")
    scope={"main":pd.read_csv(ROOT/"main_input"/"Disease_matching_candidates_1881.tsv",sep="\t"),
           "diagnostic_extra":pd.read_csv(ROOT/"main_input"/"Additional_diagnostic_only_76genes.tsv",sep="\t")}
    summaries=[]
    for name,pool in scope.items():
        selected=all_rows[all_rows.canonical_id.isin(pool.canonical_id)].copy()
        selected["input_scope"]=name
        selected.to_csv(ROOT/"main_input"/("Resolved_"+name+"_protein_sequences.tsv"),sep="\t",index=False,na_rep="NA")
        unavailable=pool.loc[~pool.canonical_id.isin(selected.canonical_id),["canonical_id","gene_symbol"]]
        unavailable.to_csv(ROOT/"main_input"/(name+"_genes_without_available_sequence.tsv"),sep="\t",index=False)
        summaries.append({"scope":name,"input_genes":len(pool),"genes_with_sequences":selected.canonical_id.nunique(),
                          "gene_protein_rows":len(selected),"unique_protein_sequences":selected.sequence_sha256.nunique(),
                          "remaining_no_sequence":len(unavailable),"n_reviewed_sequence_rescued_rows":int(selected.mapping_tier.eq("exact_GENCODEv36_protein_sequence").sum()),
                          "n_GENCODE_translation_only_rows":int(selected.mapping_tier.eq("GENCODEv36_translation_only").sum()),
                          "CNN_truncates_rows":int(selected.CNN_truncates_after1000.astype(str).str.lower().eq("true").sum())})
    (out/"Resolved_scope_summary.json").write_text(json.dumps(summaries,indent=2),encoding="utf-8")
    print(json.dumps(summaries),flush=True)

if __name__=="__main__":main()
