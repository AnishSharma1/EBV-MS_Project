# Cross-source public-evidence recovery

`iedb_iq_api_arm_census.csv` is a literal exact-peptide/exact-allele retrieval against IEDB's public Query API. The full returned records are preserved in `iedb_iq_api_raw_records.jsonl`.

`literature_recovery_manifest.csv` provides one exact peptide-and-HLA primary-literature query per unique frozen peptide/allele/arm combination. Search hits must be reviewed down to an experimental assay before they can affect the gate.

`pubmed_targeted_discovery.csv` records the completed exact-sequence PubMed title/abstract discovery pass. It returned no candidate citations, which is a retrieval result rather than evidence of absence. `independent_source_retrieval_status.csv` records the current access status for every independent source, including sources that timed out or were rate-limited; these conditions are never converted to biological negatives.

`presentation_source_review_manifest.csv` defines the next targeted review across HLA Ligand Atlas (self peptides only), SysteMHC Atlas, and Human HLA PeptideAtlas. These sources address natural presentation and are deliberately kept separate from the locked exact-allele binding/ligand gate. No candidate rank is changed by this recovery.
