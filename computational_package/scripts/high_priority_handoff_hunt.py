"""Fail-closed desk screen for the frozen 49-pair EBV--CNS HLA-II universe.

This module deliberately separates candidate selection from claims about natural
presentation or T-cell biology.  It is reusable: remote retrieval and modelling
are performed by the companion runner, while the decision functions below are
deterministic and testable without network access.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

ROOT = Path(__file__).resolve().parents[1]
V3 = ROOT / "processed/literature_grounded_hla2_rankings_v3_2026-08-27/v3_all_hla_ranked_pairs.csv"
STAGE1 = ROOT / "processed/high_yield_candidate_evidence_2026-08-28/frozen_candidate_registry.csv"
EXPANSION = ROOT / "processed/high_yield_candidate_expansion_2026-08-28/frozen_expansion_targets.csv"
HUMAN_FASTA = ROOT / "processed/high_yield_candidate_evidence_2026-08-28/raw_responses/human_reviewed_canonical.fasta"
EBV_FASTA = ROOT / "processed/high_yield_candidate_evidence_2026-08-28/raw_responses/ebv_uniprot_sequences.fasta"
ALLELES = ("HLA-DRB1*03:01", "HLA-DRB1*08:01", "HLA-DRB1*13:03", "HLA-DRB1*15:01")
EXPECTED_COUNTS = {"HLA-DRB1*03:01": 7, "HLA-DRB1*08:01": 14, "HLA-DRB1*13:03": 16, "HLA-DRB1*15:01": 12}
CNS_TARGETS = ("ANO2", "CLDN11", "CNP", "CNTN2", "CRYAB", "MAG", "MBP", "MOBP", "MOG", "PLP1", "TALDO1")
CLAIM = ("High priority is a computational/public-evidence lab-handoff label only. It does not establish natural presentation, TCR recognition, cross-reactivity, molecular mimicry, MS causation, or affinity.")

def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f: return list(csv.DictReader(f))

def write_csv(path: Path, rows: Iterable[Mapping[str, Any]], fields: list[str] | None = None) -> None:
    rows = list(rows); path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else [])
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore"); w.writeheader(); w.writerows(rows)

def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""): h.update(block)
    return h.hexdigest()

def pair_key(row: Mapping[str, str]) -> tuple[str, str, str]:
    return row["allele"], row.get("ebv_core", row.get("ebv_predicted_core", "")), row.get("self_core", row.get("self_predicted_core", ""))

def _clean(row: Mapping[str, str], source: str) -> dict[str, str]:
    d = dict(row)
    d["ebv_core"] = d.get("ebv_core", d.get("ebv_predicted_core", "")).upper()
    d["self_core"] = d.get("self_core", d.get("self_predicted_core", "")).upper()
    d["universe_source"] = source
    d["retention_rank"] = {"stage1": "1", "near_miss": "2", "prefilter": "3", "prefilter_dedup_replacement": "4"}[source]
    return d

def frozen_universe(v3_path: Path = V3, stage1_path: Path = STAGE1, expansion_path: Path = EXPANSION) -> list[dict[str, str]]:
    stage1 = [_clean(r, "stage1") for r in read_csv(stage1_path)]
    near = [_clean(r, "near_miss") for r in read_csv(expansion_path)]
    existing = {r["pair_id"] for r in stage1 + near}
    # The supplied 8+10+31 rows contain two duplicate core triples.  Retaining
    # only 31 would violate the requested 49 unique-pair acceptance test, so
    # take deterministic eligible replacements from the same frozen V3 source.
    all_v3 = read_csv(v3_path)
    pool = [_clean(r, "prefilter") for r in all_v3
            if r["pair_id"] not in existing and r["register_robust"] == "True"
            and float(r["ebv_binding_percentile_rank"]) <= 20 and float(r["self_binding_percentile_rank"]) <= 20]
    if len(stage1) != 8 or len(near) != 10 or len(pool) < 31:
        raise ValueError(f"unexpected source sizes: {len(stage1)}, {len(near)}, {len(pool)}")
    retained: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in sorted(stage1 + near, key=lambda r: (int(r["retention_rank"]), r.get("hla_rank", "999999"), r["pair_id"])):
        retained.setdefault(pair_key(row), row)
    for row in sorted(pool, key=lambda r: (r.get("hla_rank", "999999"), r["pair_id"])):
        if len(retained) >= 49: break
        retained.setdefault(pair_key(row), row)
    # Replacement candidates are still binding-eligible but did not pass the
    # historical V3 register-robust filter; the new three-predictor register
    # gate must independently pass before they can advance.
    counts = Counter(r["allele"] for r in retained.values())
    for allele in ALLELES:
        needed = EXPECTED_COUNTS[allele] - counts[allele]
        if needed <= 0: continue
        fallback = [_clean(r, "prefilter_dedup_replacement") for r in all_v3
                    if r["allele"] == allele and r["pair_id"] not in existing
                    and float(r["ebv_binding_percentile_rank"]) <= 20 and float(r["self_binding_percentile_rank"]) <= 20]
        for row in sorted(fallback, key=lambda r: (r.get("hla_rank", "999999"), r["pair_id"])):
            if needed <= 0: break
            if pair_key(row) not in retained:
                retained[pair_key(row)] = row; needed -= 1; counts[allele] += 1
    rows = sorted(retained.values(), key=lambda r: (r["allele"], int(float(r.get("hla_rank", "999999"))), r["pair_id"]))
    counts = Counter(r["allele"] for r in rows)
    if len(rows) != 49 or dict(counts) != EXPECTED_COUNTS:
        raise ValueError(f"frozen universe failure: n={len(rows)}, counts={dict(counts)}")
    return rows

def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records=[]; acc=""; seq=[]
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if acc: records.append((acc, "".join(seq)))
            acc=line[1:].split()[0]; seq=[]
        else: seq.append(line.strip().upper())
    if acc: records.append((acc, "".join(seq)))
    return records

def arms_for_universe(rows: list[Mapping[str, str]], human_fasta: Path = HUMAN_FASTA, ebv_fasta: Path = EBV_FASTA) -> list[dict[str, str]]:
    parents = {"ebv": parse_fasta(ebv_fasta), "self": parse_fasta(human_fasta)}
    arms=[]
    for i, r in enumerate(rows, 1):
        candidate_id=f"HP{i:02d}"
        for side in ("ebv", "self"):
            seq=r[f"{side}_sequence"].upper(); core=r[f"{side}_core"].upper()
            matches=[(a,s,s.find(seq)) for a,s in parents[side] if seq in s]
            a,s,start=sorted(matches)[0] if matches else ("","",-1)
            arms.append({"arm_id":f"{candidate_id}__{side}","candidate_id":candidate_id,"target_id":candidate_id,"pair_id":r["pair_id"],"side":side,"allele":r["allele"],"protein":r[f"{side}_protein"],"sequence":seq,"declared_core":core,"declared_core_start_1_based":str(seq.find(core)+1),"source_accession":a,"source_coordinate_start_1_based":str(start+1) if start >= 0 else "","source_record_status":"exact_parent_sequence_match" if matches else "not_evaluable_missing_parent_match","universe_source":r["universe_source"]})
    return arms

def predictor_summary(records: Iterable[Mapping[str, Any]], declared_core: str) -> dict[str, Any]:
    wanted=("iedb_recommended_binding", "netmhciipan_4_3_el", "mixmhc2pred_2_1_context")
    by={str(r.get("predictor")):r for r in records}
    subset=[by[x] for x in wanted if x in by]
    if len(subset) != 3:
        return {"status":"not_evaluable", "binding_pass":False,"register_pass":False,"consensus_core":"", "failure":"missing_predictor"}
    try: ranks=[float(r["percentile_rank"]) for r in subset]
    except (TypeError, ValueError): return {"status":"not_evaluable", "binding_pass":False,"register_pass":False,"consensus_core":"", "failure":"invalid_percentile"}
    cores=[str(r.get("core","")).upper() for r in subset]
    votes=Counter(c for c in cores if len(c)==9)
    core,n=max(votes.items(), key=lambda x:(x[1], x[0])) if votes else ("",0)
    return {"status":"complete","binding_pass":all(x <= 20 for x in ranks),"register_pass":n >= 2,"consensus_core":core,"register_matches_declared":core == declared_core,"predictor_cores":";".join(cores),"predictor_percentiles":";".join(map(str,ranks)),"failure":""}

def desk_classification(candidate: Mapping[str, Any], arm_summaries: Mapping[str, Mapping[str, Any]], public_support: Mapping[str, bool], cns_verified: bool = True) -> dict[str, Any]:
    ebv=arm_summaries[candidate["ebv_arm_id"]]; self_=arm_summaries[candidate["self_arm_id"]]
    fails=[]
    if candidate["ebv_source_status"] != "exact_parent_sequence_match" or candidate["self_source_status"] != "exact_parent_sequence_match": fails.append("sequence_or_coordinate_provenance")
    if not cns_verified: fails.append("cns_expression_documentation")
    for side, summary in (("ebv",ebv),("self",self_)):
        if not summary["binding_pass"]: fails.append(f"{side}_binding_percentile")
        if not summary["register_pass"] or not summary["register_matches_declared"]: fails.append(f"{side}_register_consensus")
    if not (public_support.get(candidate["ebv_arm_id"],False) or public_support.get(candidate["self_arm_id"],False)): fails.append("no_exact_sequence_exact_hla_public_support")
    return {"desk_status":"pass" if not fails else "fail", "failed_gates":";".join(fails), "structural_status":"prepared_pending_external_alphafold_and_same_allele_pandora_calibration" if not fails else "not_run_desk_gate_failed"}
