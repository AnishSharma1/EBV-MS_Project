# Data and artifact policy

## Included

The public package includes portable source scripts, small frozen tables, FASTA/JSON/CSV/TSV records, readable Markdown summaries, and SVG figures where available. Each included collection retains its original versioned directory name.

## Deliberately omitted

The authoritative local archive contains large or nonportable artifacts that are intentionally not mirrored here:

- AlphaFold/other raw coordinate ensembles and downloaded server folders.
- APBS raw calculations, sampled surface vectors, model-calculation records, aligned-model collections, and `io.mc` caches.
- Bundled Python runtimes, compiled extensions, temporary files, operating-system metadata, and local inspection output.
- Credentials, tokens, browser sessions, service account data, and local absolute-path configuration.

The omission of raw material does not erase provenance: package tables preserve source identifiers, request manifests, checksums, and method/gate records. Raw artifact access must be requested from the project owner and checked against the corresponding versioned manifest.

## Package policy going forward

New public updates should add a versioned subdirectory, update this package’s top-level documentation, regenerate `MANIFEST.sha256`, and pass `scripts/verify_package.py`. Do not overwrite frozen result directories or relabel exploratory outputs as validated biological findings.
