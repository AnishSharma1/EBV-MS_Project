# EBV-MS Verification Logs, Audit Trails, and QA Records

This directory consolidates all execution logs, solver outputs, audit trails, and quality assurance logs across the project.

---

## 1. Directory Structure

```
logs/
├── solver_logs/              # Physics and structure solver runtime logs
│   └── io.mc                 # APBS Poisson-Boltzmann solver output log
│
├── provenance_audits/        # Text provenance and author-verification logs
│   ├── TEXT_PROVENANCE_AUDIT.md
│   ├── TEXT_PROVENANCE_AUDIT_FINAL_PDF.md
│   ├── USER_AUTHORED_GOOGLE_DOCS_MASTER.txt
│   └── PROVENANCE_SHA256.txt
│
├── citation_audits/          # Bibliographic reference audits against external databases
│   ├── CITATION_AUDIT.md
│   └── CITATION_NOTES.md
│
├── quality_assurance/        # Language, grammar, and LaTeX compilation audits
│   ├── GRAMMARLY_REPAIR_LOG.md
│   └── GRAMMARLY_REVIEWED_REPAIRED.md
│
└── package_verification/     # Cryptographic integrity and validation manifests
    └── MANIFEST.sha256
```

---

## 2. Key Audit Roles

- **Solver Logs (`solver_logs/`)**: Retains execution logs from external biophysical tools (APBS, PDB2PQR).
- **Provenance Audits (`provenance_audits/`)**: Documents line-by-line verification that all manuscript prose derives directly from the author's verified text.
- **Citation Audits (`citation_audits/`)**: Cross-verifies all reference keys and assertions against PubMed, Europe PMC, and DOI registries.
- **Quality Assurance (`quality_assurance/`)**: Documents editorial and typographical revisions.
- **Package Verification (`package_verification/`)**: Preserves cryptographic SHA-256 signatures for deterministic reproducibility.
