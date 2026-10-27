"""Test feature sensitivity without predicting study-compound protein scores."""
import json
import pathlib
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors,Draw
from DeepPurpose.utils import smiles2morgan,smiles2daylight,trans_drug

ROOT=pathlib.Path(__file__).resolve().parents[1]

def main():
    source=ROOT/"structure_sources"
    a=json.loads((source/"PubChem_209254_properties.json").read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
    b=json.loads((source/"PubChem_500127_properties.json").read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
    stereo=Chem.MolFromSmiles(a["SMILES"])
    plain=Chem.Mol(stereo)
    Chem.RemoveStereochemistry(plain)
    connectivity=Chem.MolToSmiles(plain,isomericSmiles=False)
    assert connectivity==Chem.MolToSmiles(Chem.MolFromSmiles(b["ConnectivitySMILES"]),isomericSmiles=False)
    inverted=Chem.Mol(stereo)
    for atom in inverted.GetAtoms():
        if atom.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:atom.InvertChirality()
    variants={"PubChem209254_stereo_not_yet_confirmed_as_study_sample":Chem.MolToSmiles(stereo,isomericSmiles=True),
              "connectivity_only_from_same_graph":connectivity,
              "all_assigned_centers_inverted_technical_probe":Chem.MolToSmiles(inverted,isomericSmiles=True)}
    morgan={name:smiles2morgan(smiles) for name,smiles in variants.items()}
    daylight={name:smiles2daylight(smiles) for name,smiles in variants.items()}
    reference=list(variants)[0]
    audit={"variants":variants,"formula":rdMolDescriptors.CalcMolFormula(plain),"formal_charge":Chem.GetFormalCharge(plain),
           "n_atoms":plain.GetNumAtoms(),"unspecified_chiral_centers":Chem.FindMolChiralCenters(plain,includeUnassigned=True),
           "morgan_features_identical_across_variants":all(np.array_equal(x,morgan[reference]) for x in morgan.values()),
           "daylight_features_identical_across_variants":all(np.array_equal(x,daylight[reference]) for x in daylight.values()),
           "CNN_character_features_identical_across_variants":all(trans_drug(s)==trans_drug(variants[reference]) for s in variants.values()),
           "smiles_lengths":{name:len(s) for name,s in variants.items()},
           "study_protein_predictions_performed":False,
           "identity_status":"connectivity_compared_to_Fang2022_Fig1A; absolute_configuration_not_assigned_from_Fig1A"}
    (source/"Stereochemistry_feature_sensitivity.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    Draw.MolToFile(plain,str(source/"Connectivity_only_for_visual_verification.png"),size=(1000,650))
    print(json.dumps(audit),flush=True)

if __name__=="__main__":main()
