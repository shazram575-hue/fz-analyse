"""Read Figure1 wedge/dash assignments independently of the CID stereochemical tags.

Coordinates are manually transcribed from the public primary author's Figure1.
Only the constitution is reused from the already checked molecular graph.
All CID stereochemical tags are cleared before assigning figure-derived stereochemistry.
"""
import json
import pathlib
from rdkit import Chem
from rdkit.Geometry import Point3D
from rdkit.Chem.Draw import rdMolDraw2D

ROOT=pathlib.Path(__file__).resolve().parents[1]
PRIOR=ROOT.parent/"DeepPurpose_预测_2026-10-24"

def main():
    output=ROOT/"structure_comparison";output.mkdir(parents=True,exist_ok=True)
    prior=json.loads((PRIOR/"structure_verification_CID209254"/"Verification_report.json").read_text(encoding="utf-8"))
    reference=Chem.MolFromSmiles(prior["cleaned_user_smiles"])
    base=Chem.Mol(reference)
    Chem.RemoveStereochemistry(base)
    center_hydrogens=[2,5,9,10,13]
    base=Chem.AddHs(base,onlyOnAtoms=center_hydrogens)
    # Pixel coordinates within the enlarged source Figure1 detail (y downward).
    coords={0:(1034,113),1:(1055,218),2:(978,287),3:(1030,368),4:(979,458),5:(879,420),
            6:(879,316),7:(786,264),8:(694,316),9:(694,420),10:(786,474),11:(786,580),
            12:(694,635),13:(600,580),14:(600,474),15:(505,420),16:(413,474),17:(413,580),
            18:(506,636),19:(310,529),20:(208,430),21:(109,449),22:(81,555),23:(165,622),
            24:(270,611),25:(307,704),26:(277,793),27:(405,680),28:(492,727),29:(600,370),
            30:(879,212),31:(1117,239),32:(1190,185),33:(1160,333)}
    h_coords={2:(978,186),5:(879,524),9:(694,527),10:(786,370),13:(600,682)}
    explicit_h={}
    for center in center_hydrogens:
        hydrogen=next(atom.GetIdx() for atom in base.GetAtomWithIdx(center).GetNeighbors() if atom.GetAtomicNum()==1)
        explicit_h[center]=hydrogen;coords[hydrogen]=h_coords[center]
    # Nine source wedges/dashes: C20-Me dash; C13-Me wedge; C10-Me wedge;
    # C3-N dash; C17-H dash; C14-H dash; C9-H dash; C8-H wedge; C5-H dash.
    directions={(1,0):Chem.BondDir.BEGINDASH,(6,30):Chem.BondDir.BEGINWEDGE,
                (14,29):Chem.BondDir.BEGINWEDGE,(17,27):Chem.BondDir.BEGINDASH,
                (2,explicit_h[2]):Chem.BondDir.BEGINDASH,(5,explicit_h[5]):Chem.BondDir.BEGINDASH,
                (9,explicit_h[9]):Chem.BondDir.BEGINDASH,(10,explicit_h[10]):Chem.BondDir.BEGINWEDGE,
                (13,explicit_h[13]):Chem.BondDir.BEGINDASH}
    rw=Chem.RWMol()
    for atom in base.GetAtoms():rw.AddAtom(Chem.Atom(atom))
    for bond in base.GetBonds():
        begin,end=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if (end,begin) in directions:begin,end=end,begin
        rw.AddBond(begin,end,bond.GetBondType())
        if (begin,end) in directions:rw.GetBondBetweenAtoms(begin,end).SetBondDir(directions[(begin,end)])
    derived=rw.GetMol();Chem.SanitizeMol(derived)
    conformer=Chem.Conformer(derived.GetNumAtoms());conformer.Set3D(False)
    for index,(x,y) in coords.items():conformer.SetAtomPosition(index,Point3D(x/100.,-y/100.,0))
    derived.AddConformer(conformer,assignId=True)
    Chem.AssignChiralTypesFromBondDirs(derived,replaceExistingTags=True)
    Chem.AssignStereochemistry(derived,cleanIt=True,force=True)
    smiles=Chem.MolToSmiles(Chem.RemoveHs(derived),isomericSmiles=True)
    key=Chem.MolToInchiKey(derived)
    original_centers=dict(Chem.FindMolChiralCenters(reference,includeUnassigned=True,useLegacyImplementation=False))
    derived_centers=dict(Chem.FindMolChiralCenters(derived,includeUnassigned=True,useLegacyImplementation=False))
    comparison=[{"atom_index":index,"CID_CIP":cip,"published_A_figure_CIP":derived_centers.get(index),"same":derived_centers.get(index)==cip} for index,cip in original_centers.items()]
    report={"source":"Halaczkiewicz et al primary author manuscript, Figure1 Spiropachysine A; published DOI10.1002/chem.202404223",
            "method":"manual transcription of 2D coordinates and all nine wedge/dash bonds; clear all CID stereo tags before assignment",
            "source_coordinates":coords,"source_wedge_bonds":[{"begin_atom":a,"end_atom":b,"direction":str(direction)} for (a,b),direction in directions.items()],
            "source_connectivity":"visually matches Fang2022 Fig1A and the CID constitutional graph",
            "CID_InChIKey":prior["InChIKey"],"figure_derived_InChIKey":key,
            "figure_derived_isomeric_SMILES":smiles,"stereocenter_comparison":comparison,
            "all9_centers_match":len(comparison)==9 and all(row["same"] for row in comparison),
            "complete_InChIKey_match":key==prior["InChIKey"]}
    (output/"Published_A_Figure1_reconstruction.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    (output/"Published_A_Figure1_reconstructed.mol").write_text(Chem.MolToMolBlock(derived),encoding="utf-8")
    drawer=rdMolDraw2D.MolDraw2DCairo(1300,850);drawer.drawOptions().addStereoAnnotation=True
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,derived,legend="Reconstructed from published Figure 1: Spiropachysine A")
    drawer.FinishDrawing();(output/"Published_A_Figure1_reconstructed.png").write_bytes(drawer.GetDrawingText())
    print(json.dumps({k:v for k,v in report.items() if k not in ["source_coordinates","source_wedge_bonds"]},ensure_ascii=True),flush=True)

if __name__=="__main__":main()
