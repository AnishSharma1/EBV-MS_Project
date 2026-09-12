# Control peptide-window audit

This audit answers only whether the submitted class-II 11-mers are fair, traceable windows of the deposited controls. It does **not** establish an experimentally correct HLA register or peptide presentation.

| Control | Deposited/raw peptide sequence | Submitted 11-mer | Intended 9-residue middle core | Finding |
|---|---|---|---|---|
| 1YMM | `ENPVVHFFKNIVTPR` | `VHFFKNIVTPR` | `HFFKNIVTP` | The submission is the C-terminal 11-residue window (one V flank, 9-residue middle, one R flank). The deposition has no atom coordinates for terminal R, so coordinate scoring uses ten observed overlapping residues. |
| 2WBJ | source chain begins expression M, then `DFARVHFISAL`, followed by tag/linker residues | `DFARVHFISAL` | `FARVHFISA` | The submission removes the expression-start Met and uses the deposited 11-residue antigen segment (D flank, 9-residue middle, L flank). |

## Consequence

The original 1YMM 15-mer submission was incompatible with TCRModel2's class-II input convention. The repaired 1YMM and canonical 2WBJ submissions are syntactically valid 11-mer windows and sequence-consistent with their controls. However, neither sequence-window check substitutes for experimental peptide-HLA binding or register mapping, and the unresolved 1YMM terminal residue remains an explicit scoring caveat.

