#!/usr/bin/env python3
"""Sequence-aligned component audit for the two TCRModel2 controls.

This intentionally does not alter the existing DockQ chain normalizer.  It
uses its grouped PDB outputs, aligns model and native sequences within pMHC
(R) and TCR (L), and reports only identical-residue C-alpha pairs.  It is a
directional diagnostic: pMHC is the superposition scaffold and TCR placement
is measured after that superposition.

Run from the refresh package directory with its bundled runtime on PYTHONPATH:
  PYTHONPATH="$PWD/runtime" python3 analysis/audit_tcrmodel2_components.py
"""

from collections import OrderedDict
from pathlib import Path

import numpy as np
from Bio import pairwise2
from Bio.SVDSuperimposer import SVDSuperimposer

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}
CONTROLS = (
    ("1YMM", "CAL_1YMM", "CAL_1YMM"),
    ("2WBJ", "CAL_2WBJ_EXCLUDED", "CAL_2WBJ_EXCLUDED"),
)


def residues_by_chain(path: Path):
    """Return ordered [one-letter code, CA coordinate] records by chain."""
    chains, seen = OrderedDict(), set()
    for line in path.read_text().splitlines():
        if line[:6] != "ATOM  ":
            continue
        key = (line[21], line[22:27])
        if key not in seen:
            seen.add(key)
            chains.setdefault(line[21], []).append(
                [AA3_TO_AA1.get(line[17:20].strip(), "X"), None]
            )
        if line[12:16].strip() == "CA":
            chains[line[21]][-1][1] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
    return chains


def identical_ca_pairs(native, model, chain):
    """Globally align sequence then keep only identical, coordinate-present pairs."""
    native_sequence = "".join(item[0] for item in native[chain])
    model_sequence = "".join(item[0] for item in model[chain])
    alignment = pairwise2.align.globalms(
        native_sequence, model_sequence, 2, -1, -5, -0.5, one_alignment_only=True
    )[0]
    native_i = model_i = 0
    pairs = []
    for native_aa, model_aa in zip(alignment.seqA, alignment.seqB):
        ni = native_i if native_aa != "-" else None
        mi = model_i if model_aa != "-" else None
        if native_aa != "-":
            native_i += 1
        if model_aa != "-":
            model_i += 1
        if (
            ni is not None and mi is not None and native_aa == model_aa
            and native[chain][ni][1] is not None and model[chain][mi][1] is not None
        ):
            pairs.append((native[chain][ni][1], model[chain][mi][1]))
    return np.array([x for x, _ in pairs]), np.array([y for _, y in pairs])


def superposition(native_coordinates, model_coordinates):
    superimposer = SVDSuperimposer()
    superimposer.set(native_coordinates, model_coordinates)
    superimposer.run()
    rotation, translation = superimposer.get_rotran()
    return rotation, translation, superimposer.get_rms()


def main():
    package = Path(__file__).resolve().parents[1]
    dockq_root = package / "analysis" / "dockq_calibration"
    print(
        "control,rank,pmhc_identical_ca_pairs,tcr_identical_ca_pairs,"
        "pmhc_ca_rmsd_A,tcr_fold_ca_rmsd_A,"
        "tcr_placement_after_pmhc_alignment_ca_rmsd_A"
    )
    for control, directory, stem in CONTROLS:
        folder = dockq_root / directory
        native = residues_by_chain(folder / f"{stem}.native.canonical.grouped.pdb")
        for model_path in sorted(folder.glob("ranked_*.grouped.pdb")):
            model = residues_by_chain(model_path)
            native_pmhc, model_pmhc = identical_ca_pairs(native, model, "R")
            native_tcr, model_tcr = identical_ca_pairs(native, model, "L")
            rotation, translation, pmhc_rmsd = superposition(native_pmhc, model_pmhc)
            _, _, tcr_fold_rmsd = superposition(native_tcr, model_tcr)
            placed_tcr = model_tcr @ rotation + translation
            placement_rmsd = np.sqrt(np.mean(np.sum((placed_tcr - native_tcr) ** 2, axis=1)))
            print(
                f"{control},{model_path.stem.split('.')[0]},{len(native_pmhc)},"
                f"{len(native_tcr)},{pmhc_rmsd:.3f},{tcr_fold_rmsd:.3f},"
                f"{placement_rmsd:.3f}"
            )


if __name__ == "__main__":
    main()
