"""Independent proteome inference; no disease file is read by this program."""
import concurrent.futures
import contextlib
import hashlib
import io
import json
import pathlib
import pickle
import time
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from DeepPurpose import DTI, utils
from fast_equivalent_encoding import fast_aac

ROOT=pathlib.Path(__file__).resolve().parents[1]
PRIOR=ROOT.parent/'DeepPurpose_预测_2026-10-24'
SMILES='C[C@@H]([C@H]1CC[C@H]2[C@@H]3CC[C@H]4C[C@]5(CC[C@]4(C)[C@H]3CC[C@]12C)c1ccccc1C(=O)N5C)N(C)C'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
class PlainConfigUnpickler(pickle.Unpickler):
    def find_class(self,module,name):raise pickle.UnpicklingError('Non-plain model configuration')

def load_model(name):
    directory=PRIOR/'models'/name
    config=PlainConfigUnpickler(io.BytesIO((directory/'config.pkl').read_bytes())).load()
    config.update(num_workers=0,batch_size=16,result_folder=str(ROOT/'cache'/'internal'/name))
    pathlib.Path(config['result_folder']).mkdir(parents=True,exist_ok=True)
    net=DTI.model_initialize(**config)
    state=torch.load(directory/'model.pt',map_location='cpu',weights_only=True)
    if next(iter(state)).startswith('module.'):state={k[7:]:v for k,v in state.items()}
    net.model.load_state_dict(state,strict=True);net.binary=config['binary'];net.model.eval()
    return net

def prepare_cache(table):
    cache=ROOT/'cache';cache.mkdir(parents=True,exist_ok=True)
    order_hash=hashlib.sha256('\n'.join(table.sequence_sha256).encode()).hexdigest()
    meta=cache/'Features.json'
    if meta.exists():
        old=json.loads(meta.read_text());assert old['sequence_order_sha256']==order_hash
        for file,key in [('AAC.npy','AAC_sha256'),('CNN_indices.npy','CNN_sha256')]:assert sha(cache/file)==old[key]
        return order_hash
    assert json.loads((ROOT/'audit'/'Equivalent_encoding_validation.json').read_text())['all1993_previous_native_AAC_vectors_exactly_equal']
    aac=np.lib.format.open_memmap(cache/'AAC.npy',mode='w+',dtype=np.float32,shape=(len(table),8420))
    cats=list(utils.enc_protein.categories_[0]);lookup={char:i for i,char in enumerate(cats)}
    cnn=np.lib.format.open_memmap(cache/'CNN_indices.npy',mode='w+',dtype=np.uint8,shape=(len(table),1000))
    for index,seq in enumerate(table.protein_sequence):
        aac[index]=fast_aac(seq)
        chars=utils.trans_protein(seq);cnn[index]=[lookup[c] for c in chars]
        if index%4000==0:print(json.dumps({'stage':'features','done':index+1,'total':len(table)}),flush=True)
    aac.flush();cnn.flush();del aac,cnn
    meta.write_text(json.dumps({'sequence_order_sha256':order_hash,'AAC_sha256':sha(cache/'AAC.npy'),'CNN_sha256':sha(cache/'CNN_indices.npy'),'CNN_categories':cats},indent=2),encoding='utf-8')
    return order_hash

def run_model(item,order_hash):
    torch.set_num_threads(2)
    if torch.get_num_interop_threads()!=1:torch.set_num_interop_threads(1)
    name=item['model'];out=ROOT/'predictions';out.mkdir(parents=True,exist_ok=True)
    manifest=out/(name+'_run.json');scorepath=out/(name+'_scores.tsv')
    if manifest.exists():
        old=json.loads(manifest.read_text())
        assert old['sequence_order_sha256']==order_hash and old['model_sha256']==item['model_pt_sha256'] and old['score_sha256']==sha(scorepath)
        return old
    start=time.monotonic();assert sha(PRIOR/'models'/name/'model.pt')==item['model_pt_sha256']
    table=pd.read_csv(ROOT/'background'/'All_reviewed_human_protein_background.tsv',sep='\t',keep_default_na=False)
    net=load_model(name);model=net.model
    drug=utils.smiles2morgan(SMILES) if net.drug_encoding=='Morgan' else utils.smiles2daylight(SMILES) if net.drug_encoding=='Daylight' else utils.trans_drug(SMILES)
    drug_tensor=torch.as_tensor(utils.drug_2_embed(drug)[None],dtype=torch.float32) if net.drug_encoding=='CNN' else torch.as_tensor(drug[None],dtype=torch.float32)
    features=np.load(ROOT/'cache'/('AAC.npy' if net.target_encoding=='AAC' else 'CNN_indices.npy'),mmap_mode='r')
    def target_tensor(indices):
        array=np.array(features[indices],copy=True)
        if net.target_encoding=='AAC':return torch.from_numpy(array)
        return torch.nn.functional.one_hot(torch.from_numpy(array.astype(np.int64)),26).permute(0,2,1).float()
    with torch.no_grad():drug_embedding=model.model_drug(drug_tensor)
    def forward(indices):
        with torch.no_grad():
            protein=model.model_protein(target_tensor(indices))
            value=torch.cat((drug_embedding.expand(len(indices),-1),protein),1)
            for index,layer in enumerate(model.predictor):
                value=layer(value)
                if index<len(model.predictor)-1:value=torch.nn.functional.relu(model.dropout(value))
            return value.reshape(-1).numpy().astype(float)
    probe=np.unique(np.rint(np.linspace(0,len(table)-1,16)).astype(int))
    frame=pd.DataFrame({'SMILES':[SMILES]*len(probe),'Target Sequence':table.protein_sequence.iloc[probe].tolist(),'Label':[-1]*len(probe),'drug_encoding':[drug]*len(probe),'target_encoding':[np.array(features[i],copy=True) if net.target_encoding=='AAC' else utils.trans_protein(table.protein_sequence.iloc[i]) for i in probe]})
    with contextlib.redirect_stdout(io.StringIO()):native=np.asarray(net.predict(frame),float)
    net.model.eval()
    optimized=forward(probe);diff=float(np.max(np.abs(native-optimized)));assert diff<1e-4,(name,diff)
    values=np.empty(len(table),float)
    for offset in range(0,len(table),64):
        indices=np.arange(offset,min(offset+64,len(table)))
        values[indices]=forward(indices)
        if offset%4096==0:print(json.dumps({'stage':'predicting','model':name,'done':int(indices[-1])+1,'total':len(table)}),flush=True)
    assert np.isfinite(values).all()
    scores=table[['uniprot_accession','sequence_sha256']].copy();scores['predicted_p_value']=values
    scores.to_csv(scorepath,sep='\t',index=False)
    report={'model':name,'endpoint':item['endpoint'],'n_scores':len(table),'status':'complete','sequence_order_sha256':order_hash,'model_sha256':item['model_pt_sha256'],'score_sha256':sha(scorepath),'native_vs_cached_drug_forward_max_abs_difference':diff,'native_validation_n':len(probe),'elapsed_seconds':round(time.monotonic()-start,3),'actual_smiles':SMILES,'batch_size':64,'disease_list_read':False}
    manifest.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps({'stage':'model_complete',**report}),flush=True)
    return report

def main():
    assert Chem.MolToInchiKey(Chem.MolFromSmiles(SMILES))=='SLGWGPQWJRVPAD-PZOAWPBASA-N'
    protocol=json.loads((ROOT/'background'/'Frozen_independent_drug_side_protocol.json').read_text(encoding='utf-8'))
    assert sha(ROOT/'background'/'All_reviewed_human_protein_background.tsv')==protocol['background_sha256']
    table=pd.read_csv(ROOT/'background'/'All_reviewed_human_protein_background.tsv',sep='\t',keep_default_na=False)
    order=prepare_cache(table)
    panel=json.loads((PRIOR/'metadata'/'Model_panel_archive_audit.json').read_text())
    panel=[row for row in panel if row['config']['drug_encoding']!='MPNN'];assert len(panel)==8
    results=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        futures=[executor.submit(run_model,item,order) for item in panel]
        for future in concurrent.futures.as_completed(futures):results.append(future.result())
    out=ROOT/'predictions';(out/'All_model_run_summary.json').write_text(json.dumps(sorted(results,key=lambda r:r['model']),indent=2),encoding='utf-8')
    (out/'INDEPENDENT_RAW_PREDICTIONS_COMPLETE.txt').write_text('20431 independent reviewed human proteins x 8 models = 163448 scores. Disease lists not read during background selection or inference.\n',encoding='utf-8')
if __name__=='__main__':main()
