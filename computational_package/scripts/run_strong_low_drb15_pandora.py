"""PANDORA expected-versus-adjacent-register screen for DR15 strong-low pairs."""
from __future__ import annotations
import argparse, csv, json, shutil, sys
from pathlib import Path
from statistics import median
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'processed/high_priority_handoff_hunt_2026-09-07'
REC=ROOT/'processed/charge_reversal_computational_recovery_2026-09-07'
OUT=BASE/'pandora_strong_low_drb15_2026-09-07'
sys.path.insert(0,str(ROOT/'src'))
from run_charge_reversal_pandora import run_manifest

def read(p):
 with p.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def write(p,r,fields=None):
 r=list(r);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(r[0]),extrasaction='ignore');w.writeheader();w.writerows(r)
def main(offset=0, limit=None, finalize=False):
 OUT.mkdir(parents=True,exist_ok=True); (OUT/'pandora/templates/pmhc_only').mkdir(parents=True,exist_ok=True)
 for x in ['1BX2','6CQQ']: shutil.copy2(REC/f'pandora/templates/pmhc_only/{x}_MNP_pandora.pdb',OUT/f'pandora/templates/pmhc_only/{x}_MNP_pandora.pdb')
 shutil.copy2(REC/'pandora/results/CALIBRATION_GATE.json',OUT/'pandora/inherited_drb15_calibration_gate.json')
 arms={x['arm_id']:x for x in read(BASE/'prepared_inputs/peptide_arm_registry.csv')}
 runs=[]
 for cid in ['HP43','HP49']:
  for side in ['ebv','self']:
   a=arms[f'{cid}__{side}']; start=int(a['declared_core_start_1_based'])
   for label,st in [('minus_1',start-1),('expected',start),('plus_1',start+1)]:
    if st<1 or st+8>len(a['sequence']): continue
    core=a['sequence'][st-1:st+8]
    for template in ['1BX2','6CQQ']:
     runs.append({'run_id':f'{cid.lower()}_{side}_{label}_{template.lower()}','candidate_id':cid,'side':side,'peptide':a['sequence'],'declared_core':a['declared_core'],'register_label':label,'core_start_1_based':st,'core':core,'anchors_1_based':f'{st},{st+3},{st+5},{st+8}','forced_template_pdb':template,'requested_models':20})
 manifest=OUT/'pandora/candidate_runs_24_requesting_480_models.csv';write(manifest,runs)
 run_manifest(OUT,manifest,jobs_per_run=4,resume=True,offset=offset,limit=limit)
 if not finalize: return
 summary=[]; results=OUT/'pandora/results'
 for r in runs:
  scores=[]
  for line in (results/r['run_id']/'molpdf_DOPE.tsv').read_text().splitlines():
   n,m,d=line.split('\t');scores.append((n,float(m),float(d)))
  summary.append({**r,'produced_models':len(scores),'top_five_molpdf_median':round(median(x[1] for x in sorted(scores,key=lambda x:x[1])[:5]),5),'exact_sequence_qc':'deferred_to_model_parser'})
 write(results/'pandora_run_summary_24.csv',summary)
 armsum=[]
 for cid in ['HP43','HP49']:
  for side in ['ebv','self']:
   rr=[x for x in summary if x['candidate_id']==cid and x['side']==side]; pref={}
   for t in ['1BX2','6CQQ']:
    q=[x for x in rr if x['forced_template_pdb']==t];pref[t]=min(q,key=lambda x:x['top_five_molpdf_median'])['register_label']
   armsum.append({'candidate_id':cid,'side':side,'preferred_1BX2':pref['1BX2'],'preferred_6CQQ':pref['6CQQ'],'pandora_register_status':'supports_declared_register' if pref=={'1BX2':'expected','6CQQ':'expected'} else 'register_or_template_confounded','note':'Molpdf ranks PANDORA models only; it is not affinity.'})
 write(OUT/'pandora_arm_register_summary.csv',armsum)
 pair=[]
 for cid in ['HP43','HP49']:
  q=[x for x in armsum if x['candidate_id']==cid]; pair.append({'candidate_id':cid,'pandora_pair_status':'supports_declared_register_both_arms' if all(x['pandora_register_status']=='supports_declared_register' for x in q) else 'register_or_template_confounded','promotion':'weak_medium_provisional only if PANDORA supports both arms; otherwise remain strong_low'})
 write(OUT/'pandora_pair_decision.csv',pair)
 (OUT/'README.md').write_text('# DR15 strong-low PANDORA register review\n\nUses inherited passed DR15 calibration, two templates, expected and adjacent registers, and 20 models/run. PANDORA is structural register support, not an affinity or presentation assay.\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--offset',type=int,default=0);p.add_argument('--limit',type=int);p.add_argument('--finalize',action='store_true');a=p.parse_args();main(a.offset,a.limit,a.finalize)
