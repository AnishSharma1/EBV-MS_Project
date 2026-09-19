"""Build deliberately non-high provisional tiers from the 49-pair screen.

These labels are routing aids for a future binding/register handoff.  They do
not convert prediction into evidence of presentation or T-cell biology.
"""
from __future__ import annotations
import csv, json, hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/'processed/high_priority_handoff_hunt_2026-09-07'
OUT=IN/'secondary_provisional_tiers_2026-09-07'

def rows(p):
    with p.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def write(p,rs,fields=None):
    rs=list(rs);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rs[0]),extrasaction='ignore');w.writeheader();w.writerows(rs)
def sha(p):
    h=hashlib.sha256();h.update(p.read_bytes());return h.hexdigest()
def strict(x): return x['binding_pass'] and x['register_pass'] and x['register_matches_declared']
def near(x): return x['register_pass'] and x['register_matches_declared'] and sum(float(v)<=20 for v in x['predictor_percentiles'].split(';'))>=2

def main():
    if OUT.exists():
        for p in OUT.iterdir():
            if p.is_file():p.unlink()
    OUT.mkdir(exist_ok=True)
    matrix=rows(IN/'candidate_evidence_matrix.csv'); arms={r['arm_id']:r for r in rows(IN/'prepared_inputs/peptide_arm_registry.csv')}
    evidence=[{'candidate_id':'HP47','arm_id':'HP47__ebv','allele':'HLA-DRB1*15:01','sequence':'TGGVYHFVKKHVHES','iedb_assay_id':'1774674','reference':'Hansen et al., Tissue Antigens (2007)','method':'purified MHC/competitive/radioactivity; qualitative binding','outcome':'Positive','access_route':'IEDB website MHC Ligand Assays tab, verified 2026-09-07','decision_use':'exact-sequence/exact-allele direct-binding evidence for the BALF5 arm only'}]
    write(OUT/'iedb_ui_verified_direct_binding.csv',evidence)
    pandora_path=IN/'pandora_strong_low_drb15_2026-09-07/pandora_pair_decision.csv'
    pandora={r['candidate_id']:r for r in rows(pandora_path)} if pandora_path.exists() else {}
    out=[]
    for r in matrix:
        e=json.loads(r['ebv_predictor_summary']);s=json.loads(r['self_predictor_summary'])
        tier='not_tiered'; reason='fails_more_than_one_three-predictor binding/register condition'
        if r['candidate_id']=='HP47':
            tier='strong_medium_provisional'; reason='both arms pass all three predictor/register conditions; exact DRB1*15:01 direct binding is verified for BALF5, but inherited structural work remains register-confounded'
        elif strict(e) and strict(s):
            tier='weak_medium_provisional'; reason='both arms pass all three predictor/register conditions; no exact-allele direct binding/ligand record verified'
        elif pandora.get(r['candidate_id'],{}).get('pandora_pair_status')=='supports_declared_register_both_arms':
            tier='weak_medium_provisional'; reason='one predictor narrowly failed on one arm, but both arms retained the declared register across two calibrated DR15 PANDORA templates and expected/adjacent-register decoys'
        elif (strict(e) and near(s)) or (strict(s) and near(e)):
            tier='strong_low_provisional'; reason='one arm passes all three predictor/register conditions; the other has registered-core agreement and two of three predictor ranks <=20'
        out.append({'candidate_id':r['candidate_id'],'tier':tier,'rationale':reason,'allele':r['allele'],'pair_id':r['pair_id'],'ebv_protein':arms[r['ebv_arm_id']]['protein'],'ebv_sequence':arms[r['ebv_arm_id']]['sequence'],'ebv_core':arms[r['ebv_arm_id']]['declared_core'],'self_protein':arms[r['self_arm_id']]['protein'],'self_sequence':arms[r['self_arm_id']]['sequence'],'self_core':arms[r['self_arm_id']]['declared_core'],'ebv_percentiles':e['predictor_percentiles'],'self_percentiles':s['predictor_percentiles'],'claim_boundary':'Provisional routing tier only; not evidence of natural presentation, affinity measurement, T-cell recognition, cross-reactivity, molecular mimicry, or MS mechanism.'})
    out.sort(key=lambda r:({'strong_medium_provisional':0,'weak_medium_provisional':1,'strong_low_provisional':2,'not_tiered':3}[r['tier']],r['candidate_id']))
    write(OUT/'all_49_provisional_tiers.csv',out)
    write(OUT/'handoff_shortlist.csv',[r for r in out if r['tier']!='not_tiered'])
    summary=['# Secondary provisional handoff tiers','', '* Strong medium: HP47 only. It has both-arm three-predictor/register support plus current IEDB direct binding for its BALF5/DRB1*15:01 arm; it remains below high because the prior cross-method structural recovery was register-confounded.', '* Weak medium: HP36 has both-arm predictor/register support but lacks exact DRB1*13:03 public binding/ligand evidence. HP43 additionally passed the calibrated two-template DR15 PANDORA expected-versus-adjacent-register screen, but has a narrow single-predictor binding discrepancy.', '* Strong low: remaining pairs retain a failed predictor gate or PANDORA register/template confounding.', '', 'These are not biological claims and do not replace the frozen high-priority gate.']
    (OUT/'README.md').write_text('\n'.join(summary)+'\n',encoding='utf-8')
    write(OUT/'checksums.csv',[{'relative_path':str(p.relative_to(OUT)),'sha256':sha(p)} for p in sorted(OUT.iterdir()) if p.is_file()],['relative_path','sha256'])
if __name__=='__main__':main()
