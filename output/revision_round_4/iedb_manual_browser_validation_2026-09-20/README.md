# Manual IEDB browser validation

This is a targeted browser-mediated recovery attempt for the two peptide pairs
that pass the frozen three-predictor 20% binding/register screen: HP36
(BALF5--TALDO1, DRB1*13:03) and HP47 (BALF5--TALDO1, DRB1*15:01).

## Result

The IEDB public web interface was available despite the earlier Query API
failure. Four exact-peptide lookups were evaluated in the browser.

- HP47's BALF5 arm has a positive, exact-sequence/exact-HLA-DRB1*15:01 MHC
  binding record (IEDB assay 1774674). Thus the public-evidence gate is
  evaluable and passes **for that arm**.
- HP47's TALDO1 arm has no exact-sequence IEDB record in the public interface.
- HP36's TALDO1 arm has one exact-sequence positive MHC-ligand record, but
  only at generic `HLA class II`; it is not evidence for DRB1*13:03.
- HP36's BALF5 search returned exact-sequence records, but no
  HLA-DRB1*13:03 row was displayed.

This replaces neither the 98-arm automated retrieval nor its failed API audit.
It establishes that the blanket statement that public evidence is unavailable
or unevaluable is too broad: browser-retrieved evidence is available for the
HP47 viral arm. The remaining 94 arms were not inferred to be negative and are
explicitly `not manually reviewed` rather than treated as missing evidence.

The recovered record supports experimental peptide--HLA binding only. It does
not establish endogenous presentation, recognition of the TALDO1 peptide,
dual-reactive T cells, cross-reactivity, molecular mimicry, or MS causation.
