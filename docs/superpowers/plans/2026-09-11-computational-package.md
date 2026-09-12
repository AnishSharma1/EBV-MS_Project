# EBV-MS Computational Package Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Add a compact, navigable EBV-MS computational package that preserves validated inputs, results, provenance, and claim boundaries without publishing raw model ensembles, cached solvers, credentials, or machine-specific runtimes.

**Architecture:** Keep the historical repository contents intact and add `computational_package/` as the canonical public entry point. Store portable code in `scripts/`, frozen tables and compact summaries in `data/` and `results/`, and make provenance plus intentional omissions explicit in documentation and a checksum manifest.

**Tech Stack:** Python 3, CSV/TSV/JSON/FASTA, Markdown, SHA-256.

**Spec:** User request on 2026-09-11: update the GitHub repository with recent project material and organize it into a full computational package.

## Global Constraints

- Preserve the authoritative 17 GB iCloud archive without moving or editing it.
- Do not publish tokens, credentials, local runtimes, downloaded raw structure ensembles, APBS caches, or machine-specific paths as runnable instructions.
- Preserve the documented boundary: computational prioritization is not evidence of natural presentation, TCR binding/recognition, activation, cross-reactivity, molecular mimicry, or MS mechanism.
- Keep the unsuccessful electrostatics control gate visible; do not reframe a failed scientific gate as a positive candidate result.

---

### Task 1: Establish package contract and navigation

**Files:**
- Create: `README.md`
- Create: `computational_package/README.md`
- Create: `computational_package/docs/SCOPE_AND_LIMITS.md`
- Create: `computational_package/docs/DATA_AND_ARTIFACT_POLICY.md`

**Interfaces:**
- Consumes: frozen project artifacts and original source scripts.
- Produces: a documented public entry point and an explicit artifact policy.

- [x] **Step 1: Create root navigation**
- [x] **Step 2: Define package scope, evidence limits, and run order**
- [x] **Step 3: Record inclusions and omissions**
- [x] **Step 4: Verify every linked package path exists**

### Task 2: Add compact computational artifacts

**Files:**
- Create: `computational_package/scripts/`
- Create: `computational_package/data/`
- Create: `computational_package/results/`
- Create: `computational_package/MANIFEST.sha256`

**Interfaces:**
- Consumes: source scripts; frozen t-cell library; control validation; compact electrostatics records; corrected Stage-1 review; TCRmodel2 reports.
- Produces: a public subset with traceable checksums.

- [x] **Step 1: Copy portable scripts without rewriting their historical behavior**
- [x] **Step 2: Copy compact, versioned tables and reports**
- [x] **Step 3: Exclude raw structures, solver caches, and bundled runtimes**
- [x] **Step 4: Generate SHA-256 manifest**

### Task 3: Add automated package validation

**Files:**
- Create: `computational_package/scripts/verify_package.py`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: package directory and `MANIFEST.sha256`.
- Produces: a nonzero exit code for absent required files, forbidden artifacts, checksum drift, or an altered documented gate status.

- [x] **Step 1: Implement required-file, exclusion, gate-status, and checksum checks**
- [x] **Step 2: Add secret/runtime/cache exclusion rules to `.gitignore`**
- [x] **Step 3: Run the validator**
- [x] **Step 4: Review staged files and publish only the verified package**
