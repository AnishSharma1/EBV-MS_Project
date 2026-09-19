#!/usr/bin/env python3
"""Write SHA-256 checksums for the manifests/configuration, excluding derived hashes."""
from __future__ import annotations
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
files = sorted(p for p in ROOT.rglob("*") if p.is_file() and p.name != "SHA256SUMS.txt")
lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}" for p in files]
(ROOT / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
print(f"Wrote checksums for {len(files)} files")
