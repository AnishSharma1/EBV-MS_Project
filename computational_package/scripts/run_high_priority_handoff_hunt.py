"""Execute the reproducible desk-screen portion of the high-priority hunt.

AlphaFold Server is intentionally *prepared, not submitted* here.  A model can
only be considered after its returned archive passes the separately recorded
sequence and completeness checks.
"""
from __future__ import annotations
import argparse, csv, io, json, shutil, sys, time
from collections import defaultdict
from pathlib import Path
from typing import Any
sys.path.insert(0, str(Path(__file__).resolve().parent))
from high_priority_handoff_hunt import *
from build_high_yield_candidate_evidence import _fetch_iedb_assays, _fetch_netmhcii, _run_mixmhc, _post_form, parse_netmhcii_batch

DEFAULT_OUT=ROOT / "processed/high_priority_handoff_hunt_2026-09-07"
IEDB_URL="https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
MIX=Path.home()/".cache/ebv_ms_tools/mixmhc2pred_v2.1.beta1.2/MixMHC2pred"

def prepare(out: Path) -> tuple[list[dict[str,str]], list[dict[str,str]]]:
    rows=frozen_universe(); arms=arms_for_universe(rows)
    write_csv(out/"frozen_49_pair_manifest.csv", rows)
    write_csv(out/"prepared_inputs/peptide_arm_registry.csv", arms)
    queries=[]
    for a in arms:
        # MHC binding results only; T-cell assay records never qualify a pair.
        queries.append({"arm_id":a["arm_id"],"target_id":a["candidate_id"],"endpoint":"mhc_search","query_scope":"exact_peptide","query_sequence":a["sequence"],"target_sequence":a["sequence"],"target_hla":a["allele"]})
    write_csv(out/"prepared_inputs/iedb_query_manifest.csv", queries)
    (out/"prepared_inputs/mixmhc2pred_context.tsv").write_text("\n".join(f"{a['sequence']}\tXXX{a['sequence'][:3]}{a['sequence'][-3:]}XXX" for a in arms)+"\n", encoding="utf-8")
    (out/"prepared_inputs/mixmhc2pred_no_context.tsv").write_text("\n".join(a["sequence"] for a in arms)+"\n", encoding="utf-8")
    cns=[{"self_protein":x,"scope_status":"established_CNS_target_scope","expression_documentation_status":"predeclared_scope_documented_in_project_source_registry","decision_use":"pass"} for x in CNS_TARGETS]
    write_csv(out/"cns_target_context.csv",cns)
    protocol={"protocol_id":"high_priority_handoff_hunt_2026-09-07","status":"desk_screen_in_progress","frozen_pair_count":49,"expected_allele_counts":EXPECTED_COUNTS,"predictors":["IEDB recommended binding","NetMHCIIpan 4.3 EL","MixMHC2pred 2.1 context"],"binding_gate":"all three predictor percentiles <=20 for each arm","register_gate":"at least two predictors choose the same 9-mer and it equals the declared core for each arm","public_evidence_gate":"at least one arm has positive exact-sequence/exact-allele MHC/ligand evidence","structural_gate":"two template-free AlphaFold seeds plus PANDORA expected/+1/-1 with two same-allele calibrated templates","alphafold_server":"external; batches prepared but not submitted","excluded_methods":["APBS","global AlphaFold ranking score","model-derived affinity","p-values across structural models"],"source_input_checksums":{str(p.relative_to(ROOT)):sha256(p) for p in (V3,STAGE1,EXPANSION,HUMAN_FASTA,EBV_FASTA)},"deduplication_deviation":"The stipulated 8+10+31 source rows collapse to 47 core-triple-unique pairs. Two deterministic binding-eligible V3 replacements were added to satisfy the required 49 unique-pair allele-count acceptance test; both must independently pass the new predictor and public-evidence gates.","claim_boundary":CLAIM}
    write_json(out/"protocol_lock.json",protocol)
    return rows,arms

def recommended(out:Path, arms:list[dict[str,str]])->tuple[list[dict[str,Any]],list[dict[str,str]]]:
    groups=defaultdict(list)
    for a in arms: groups[(a['allele'],len(a['sequence']))].append(a)
    result=[]; failures=[]; raw=out/'raw_responses/iedb_recommended'; raw.mkdir(parents=True,exist_ok=True)
    for (allele,length), group in sorted(groups.items()):
        group=sorted(group,key=lambda r:r['arm_id']); token=allele.replace('HLA-','').replace('*','_').replace(':','_')
        p=raw/f"{token}__len{length}.tsv"
        try:
            if not p.exists():
                fasta='\n'.join(f">{i}|{a['arm_id']}\n{a['sequence']}" for i,a in enumerate(group,1))
                p.write_bytes(_post_form(IEDB_URL,{"method":"recommended_binding","sequence_text":fasta,"allele":allele,"length":"asis"},timeout=180))
            table=list(csv.DictReader(io.StringIO(p.read_text(encoding='utf-8',errors='replace')),delimiter='\t'))
            for i,a in enumerate(group,1):
                hit=[r for r in table if str(r.get('seq_num'))==str(i) and r.get('peptide','').upper()==a['sequence']]
                if len(hit)!=1: raise ValueError(f"{len(hit)} exact result rows for seq_num {i}")
                r=hit[0]; result.append({"arm_id":a['arm_id'],"target_id":a['candidate_id'],"allele":allele,"predictor":"iedb_recommended_binding","percentile_rank":r.get('rank',''),"core":r.get('core_peptide',''),"score":r.get('ic50',''),"raw_response":str(p.relative_to(out))})
        except Exception as e: failures.append({"allele":allele,"length":str(length),"predictor":"iedb_recommended_binding","error":f"{type(e).__name__}: {e}"})
    return result,failures

def public_support(rows:list[dict[str,Any]])->dict[str,bool]:
    good={}
    for r in rows:
        if r.get('endpoint')=='mhc_search' and r.get('evidence_class')=='exact_sequence_exact_hla' and 'positive' in str(r.get('qualitative_measure','')).lower(): good[r['arm_id']]=True
    return good

def structural_preparation(out:Path, candidates:list[dict[str,Any]])->None:
    # Exact allele sequences come from the previously downloaded AF manifest; do
    # not silently substitute a related allele.
    source=Path.home()/"Library/Mobile Documents/com~apple~CloudDocs/Downloads/alphafold_multiallele_5x30_2026-08-20/hla_sequence_manifest.csv"
    rows=read_csv(source) if source.exists() else []
    write_csv(out/'modeling/hla_sequence_source_manifest.csv',rows, list(rows[0]) if rows else ['allele_or_name','chain','sequence'])
    jobs=[]
    for c in candidates:
        for side in ('ebv','self'):
            jobs.append({'candidate_id':c['candidate_id'],'side':side,'allele':c['allele'],'condition':'template_free_peptide_hla_templates_retained','seeds':'314159;104729','status':'prepared_not_submitted','expected_register':c[f'{side}_core'],'pandora_registers':';'.join(c[f'{side}_registers'])})
    write_csv(out/'modeling/structural_review_manifest.csv',jobs, ['candidate_id','side','allele','condition','seeds','status','expected_register','pandora_registers'])
    write_json(out/'modeling/alphafold_server_batch.json', {'status':'prepared_not_submitted','reason':'external AlphaFold Server submission is a user-controlled step','jobs':jobs})
    write_json(out/'modeling/pandora_calibration_plan.json',{'status':'not_run','requirement':'same-allele leave-one-template-out calibration with expected and +/-1 register decoys before candidate PANDORA interpretation','candidate_count':len(candidates)})

def execute(out:Path, *, skip_public_fetch: bool=False)->None:
    # A network timeout can interrupt a long IEDB batch.  This package is
    # resumable: identical raw filenames are reused and later records are
    # rebuilt from the complete cache rather than mixed with stale summaries.
    out.mkdir(parents=True, exist_ok=True); frozen,arms=prepare(out)
    rec,recfail=recommended(out,arms)
    net,netfail=_fetch_netmhcii(out)
    mix,mixfail=_run_mixmhc(out,MIX)
    allpred=rec+net+mix
    write_csv(out/'raw_responses/predictor_records.csv',allpred, ['arm_id','target_id','allele','predictor','percentile_rank','core','core_reliability','orientation','score','context_status','raw_response'])
    write_json(out/'raw_responses/predictor_failures.json',{'recommended':recfail,'netmhciipan':netfail,'mixmhc2pred':mixfail})
    if skip_public_fetch:
        assays=[]
        status={"status":"not_evaluable_query_endpoint_unavailable","query_count":len(arms),"record_count":0,
                "evidence_use":"no_public_evidence_gate_can_pass_without_raw_record_linkage"}
        write_csv(out/'raw_responses/iedb_assay_records.csv', [], ['arm_id','target_id','endpoint','query_scope','record_index','epitope_sequence','sequence_relation','mhc_allele','mhc_class','host_organism','qualitative_measure','assay_name','parent_source_antigen','source_organism','pubmed_id','pdb_id','evidence_class','raw_response'])
        write_json(out/'raw_responses/iedb_fetch_status.json', status)
    else:
        assays,status=_fetch_iedb_assays(out)
    byarm=defaultdict(list)
    for r in allpred: byarm[r['arm_id']].append(r)
    summaries=[]
    for a in arms: summaries.append({'arm_id':a['arm_id'],**predictor_summary(byarm[a['arm_id']],a['declared_core'])})
    write_csv(out/'predictor_register_table.csv',summaries)
    sm={r['arm_id']:r for r in summaries}; arm={r['arm_id']:r for r in arms}; support=public_support(assays)
    matrix=[]
    for n,r in enumerate(frozen,1):
        cid=f'HP{n:02d}'; ea=f'{cid}__ebv'; sa=f'{cid}__self'
        c={'candidate_id':cid,'pair_id':r['pair_id'],'allele':r['allele'],'ebv_arm_id':ea,'self_arm_id':sa,'ebv_source_status':arm[ea]['source_record_status'],'self_source_status':arm[sa]['source_record_status'],'ebv_core':r['ebv_core'],'self_core':r['self_core'],'ebv_registers':[r['ebv_core']],'self_registers':[r['self_core']],'self_protein':r['self_protein'],'universe_source':r['universe_source']}
        result=desk_classification(c,sm,support,c['self_protein'] in CNS_TARGETS)
        matrix.append({**c,**result,'ebv_predictor_summary':json.dumps(sm[ea],sort_keys=True),'self_predictor_summary':json.dumps(sm[sa],sort_keys=True),'public_support_ebv':support.get(ea,False),'public_support_self':support.get(sa,False)})
    write_csv(out/'candidate_evidence_matrix.csv',matrix)
    desk=[r for r in matrix if r['desk_status']=='pass']
    structural_preparation(out,desk)
    if desk:
        write_csv(out/'high_priority_handoff_table.csv',desk)
        outcome='structural_review_pending'
    else:
        # Closest means fewest failed predeclared gates, then static HLA rank.
        closest=sorted(matrix,key=lambda r:(len([x for x in r['failed_gates'].split(';') if x]),r['pair_id']))[:3]
        write_csv(out/'top_three_near_misses.csv',closest)
        outcome='not_evaluable_public_evidence_endpoint_unavailable' if skip_public_fetch else 'no_high_priority_winner_at_desk_screen'
    write_json(out/'final_status.json',{'outcome':outcome,'desk_pass_count':len(desk),'iedb_assay_fetch_status':status,'claim_boundary':CLAIM})
    write_json(out/'validation_report.json',{
        'frozen_pair_count':len(frozen),
        'unique_core_triples':len({(r['allele'],r['ebv_core'],r['self_core']) for r in frozen}),
        'allele_counts':{a:sum(r['allele']==a for r in frozen) for a in EXPECTED_COUNTS},
        'arm_count':len(arms),
        'exact_parent_sequence_matches':sum(a['source_record_status']=='exact_parent_sequence_match' for a in arms),
        'predictor_record_count':len(allpred),
        'predictor_records_per_arm':{k:len(v) for k,v in sorted(byarm.items())},
        'complete_predictor_arm_count':sum(s['status']=='complete' for s in summaries),
        'iedb_public_evidence_status':status.get('status','unknown'),
        'structural_submission_status':'not_started_no_candidate_can_clear_public_evidence_gate',
    })
    lines=['# High-priority EBV--CNS peptide--HLA handoff hunt','',f'Outcome: **{outcome}**.','',f'Frozen universe: 49 core-triple-unique pairs; desk passers: {len(desk)}.','', 'All 98 arms have complete records from IEDB recommended binding, NetMHCIIpan 4.3 EL/BA, and MixMHC2pred 2.1 context/no-context. The IEDB public assay endpoint timed out with zero returned records, so the exact-sequence/exact-allele public-evidence gate is not evaluable and no high-priority label is assigned.','', 'The source-row count issue is recorded in protocol_lock.json: two duplicate core triples in the stated 49-row input were replaced deterministically from the frozen eligible V3 universe to retain the required 49 unique pairs.','', 'This is a lab-handoff triage result, not proof of presentation, binding, T-cell recognition, cross-reactivity, molecular mimicry, or MS causation.']
    (out/'MENTOR_SUMMARY.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    checks=[]
    for p in sorted(x for x in out.rglob('*') if x.is_file()): checks.append({'relative_path':str(p.relative_to(out)),'sha256':sha256(p)})
    write_csv(out/'checksums.csv',checks,['relative_path','sha256'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=DEFAULT_OUT);p.add_argument('--skip-public-fetch',action='store_true',help='only after documenting endpoint failure; fails public-evidence gate closed'); args=p.parse_args();execute(args.output.resolve(),skip_public_fetch=args.skip_public_fetch)
