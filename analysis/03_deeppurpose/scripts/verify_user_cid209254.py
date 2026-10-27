"""Verify the user-supplied structure against PubChem; do not run target prediction."""
import hashlib
import io
import json
import pathlib
import urllib.request
from datetime import datetime, timezone
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/"structure_verification_CID209254"
RAW=r"C[C@@H]\([C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CC[C@@]5(C4)C6=CC=CC=C6C(=O)N5C)C)C)N(C)C"
URLS={"properties":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/209254/property/CanonicalSMILES,IsomericSMILES,InChI,InChIKey,MolecularFormula,IUPACName/JSON",
      "synonyms":"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/209254/synonyms/JSON"}

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    retrieval=[]
    bodies={}
    for name,url in URLS.items():
        record={"name":name,"url":url}
        try:
            cache_path=OUT/("PubChem_"+name+".json")
            if cache_path.exists():
                body=cache_path.read_bytes()
                bodies[name]=json.loads(body)
                record.update(status="reused_successful_API_response_from_this_verification",bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
                retrieval.append(record)
                continue
            req=urllib.request.Request(url,headers={"User-Agent":"Workflow3-structure-verification/1.0"})
            with urllib.request.urlopen(req,timeout=25) as response:body=response.read(2_000_000)
            bodies[name]=json.loads(body)
            (OUT/("PubChem_"+name+".json")).write_bytes(body)
            record.update(status="retrieved",bytes=len(body),sha256=hashlib.sha256(body).hexdigest())
        except Exception as exc:
            record.update(status="not_retrieved",error=str(exc))
        retrieval.append(record)
    if "properties" in bodies:
        reference=bodies["properties"]["PropertyTable"]["Properties"][0]
        reference_status="official_API_response_retrieved_during_this_verification"
    else:
        cached=ROOT/"structure_sources"/"PubChem_209254_properties.json"
        reference=json.loads(cached.read_text(encoding="utf-8"))["PropertyTable"]["Properties"][0]
        reference_status="previously_downloaded_official_API_response_current_refresh_failed"
    raw_molecule=Chem.MolFromSmiles(RAW)
    cleaned=RAW.replace(r"\(","(")
    molecule=Chem.MolFromSmiles(cleaned)
    if molecule is None:raise ValueError("Corrected user SMILES did not parse")
    Chem.AssignStereochemistry(molecule,cleanIt=True,force=True)
    canonical=Chem.MolToSmiles(molecule,isomericSmiles=True,canonical=True)
    inchi=Chem.MolToInchi(molecule)
    key=Chem.MolToInchiKey(molecule)
    centers=Chem.FindMolChiralCenters(molecule,includeUnassigned=True,useLegacyImplementation=False)
    assert key==reference["InChIKey"]
    assert inchi==reference["InChI"]
    assert rdMolDescriptors.CalcMolFormula(molecule)==reference["MolecularFormula"]
    assert len(centers)==9 and not any(cip=="?" for _,cip in centers)
    synonyms=bodies.get("synonyms",{}).get("InformationList",{}).get("Information",[{}])[0].get("Synonym",[])
    AllChem.Compute2DCoords(molecule)
    Chem.WedgeMolBonds(molecule,molecule.GetConformer())
    mol_path=OUT/"PubChem_CID209254_verified_2D.mol"
    sdf_path=OUT/"PubChem_CID209254_verified_2D.sdf"
    mol_path.write_text(Chem.MolToMolBlock(molecule),encoding="utf-8")
    sdf_buffer=io.StringIO()
    with Chem.SDWriter(sdf_buffer) as writer:
        molecule.SetProp("_Name","Spiropachysine | PubChem CID 209254 | study identity correspondence pending")
        molecule.SetProp("PubChem_CID","209254")
        molecule.SetProp("InChIKey",key)
        molecule.SetProp("Source","https://pubchem.ncbi.nlm.nih.gov/compound/209254")
        writer.write(molecule)
        writer.flush()
        sdf_path.write_text(sdf_buffer.getvalue(),encoding="utf-8")
    roundtrip=Chem.MolFromMolBlock(mol_path.read_text(encoding="utf-8"),removeHs=False)
    assert Chem.MolToInchiKey(roundtrip)==key
    for format_name in ["svg","png"]:
        drawer=rdMolDraw2D.MolDraw2DSVG(1100,720) if format_name=="svg" else rdMolDraw2D.MolDraw2DCairo(1100,720)
        drawer.drawOptions().addStereoAnnotation=True
        rdMolDraw2D.PrepareAndDrawMolecule(drawer,molecule,legend="PubChem CID 209254 | stereochemistry verified against the database")
        drawer.FinishDrawing()
        image=drawer.GetDrawingText()
        if format_name=="svg":(OUT/"CID209254_stereochemistry.svg").write_text(image,encoding="utf-8")
        else:(OUT/"CID209254_stereochemistry.png").write_bytes(image)
    (OUT/"verified_input_smiles.smi").write_text(cleaned+"\tPubChem_CID209254\n",encoding="utf-8")
    (OUT/"rdkit_canonical_isomeric_smiles.smi").write_text(canonical+"\tPubChem_CID209254\n",encoding="utf-8")
    report={"checked_at_utc":datetime.now(timezone.utc).isoformat(),"raw_user_smiles":RAW,"raw_parses":raw_molecule is not None,
            "formatting_correction":"removed only the backslash immediately preceding the branch-opening parenthesis; retained all stereochemical @/@@ markers",
            "cleaned_user_smiles":cleaned,"matches_official_SMILES_text":cleaned==reference.get("SMILES",reference.get("IsomericSMILES")),
            "matches_official_InChI":True,"matches_official_InChIKey":True,"InChIKey":key,"reference_status":reference_status,
            "molecular_formula":rdMolDescriptors.CalcMolFormula(molecule),"molecular_weight_g_mol_rdkit":Descriptors.MolWt(molecule),
            "formal_charge":Chem.GetFormalCharge(molecule),"n_disconnected_components":len(Chem.GetMolFrags(molecule)),
            "n_assigned_chiral_centers":len(centers),"chiral_centers_RDKit_zero_based_atom_index_CIP":centers,
            "n_unassigned_potential_chiral_centers":sum(cip=="?" for _,cip in centers),"mol_file_stereochemical_roundtrip_pass":True,
            "first_PubChem_synonyms":synonyms[:12],"Spiropachysine_A_exact_synonym_present":any(s.casefold()=="spiropachysine a" for s in synonyms),
            "conclusion":"User structure matches CID209254 after formatting correction, with all nine centers assigned. This does not independently establish identity with Fang2022 experimental material named Spiropachysine A.",
            "study_drug_identity_approved":False,"drug_target_predictions_run":False,"retrieval":retrieval}
    (OUT/"Verification_report.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(report,ensure_ascii=True),flush=True)

if __name__=="__main__":main()
