from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from high_priority_handoff_hunt import EXPECTED_COUNTS, frozen_universe, pair_key, predictor_summary

def test_frozen_universe_has_exact_predeclared_counts():
    rows=frozen_universe()
    assert len(rows)==49
    assert len({pair_key(r) for r in rows})==49
    assert {a:sum(r['allele']==a for r in rows) for a in EXPECTED_COUNTS} == EXPECTED_COUNTS

def test_predictor_gate_boundaries_and_register_vote():
    rows=[{'predictor':'iedb_recommended_binding','percentile_rank':'20','core':'ABCDEFGHI'}, {'predictor':'netmhciipan_4_3_el','percentile_rank':'20','core':'ABCDEFGHI'}, {'predictor':'mixmhc2pred_2_1_context','percentile_rank':'20','core':'XXXXXXXXX'}]
    result=predictor_summary(rows,'ABCDEFGHI')
    assert result['binding_pass'] and result['register_pass'] and result['register_matches_declared']
    rows[2]['percentile_rank']='20.01'
    assert not predictor_summary(rows,'ABCDEFGHI')['binding_pass']
