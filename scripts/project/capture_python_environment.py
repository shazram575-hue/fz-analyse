from pathlib import Path
import argparse,json,subprocess,sys,platform,importlib.util,importlib.metadata as im
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--name',required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
mods=['numpy','pandas','scipy','sklearn','anndata','scanpy','cnmf','statsmodels','torch','rdkit','DeepPurpose','h5py','matplotlib','seaborn','requests','pytest']
snap=dict(name=a.name,python=sys.version,platform=platform.platform(),executable=sys.executable,modules={m:bool(importlib.util.find_spec(m)) for m in mods})
(a.output/(a.name+'.runtime.json')).write_text(json.dumps(snap,indent=2),encoding='utf-8')
for suffix,cmd in [('pip_list.json',['list','--local','--format=json']),('freeze.txt',['freeze','--local'])]:
 r=subprocess.run([sys.executable,'-m','pip',*cmd],capture_output=True,text=True,encoding='utf-8');(a.output/(a.name+'.'+suffix)).write_text(r.stdout,encoding='utf-8');assert r.returncode==0,r.stderr
print(json.dumps(snap,indent=2))
