# Results changelog

## 2026-09-11 00_gate (results/00_gate/20260911T182159Z, commit 2954e1b)

- All four sites pass the Phase 0 gate, so the verdict is go. Three MGB tables only link by patient and order date.
- E1 composite prevalence is 11.1% at Stanford, 13.3% at MGB, 7.9% at UTSW and 9.3% at ECUH. MGB rises to 18.2% under CLSI 2022.
- E2 decision-time prevalence is 1.17%, 1.85%, 0.53% and 0.96%, about a tenfold drop from E1 at every site.
- Mapping UTSW's "Pip/tazobactam" raised UTSW E1 evaluable cultures from 49,615 to 50,327 and resistant cultures from 3,608 to 3,973.
- Flagged: mismatched MGB IDs in three tables, MGB urine and respiratory positivity two to three times the other sites, and 14,397 UTSW positives with no organism or AST.

## 2026-09-11 features (results/features_<site>/20260911T2007*, commit d31b2f2)

- Built the 48 harmonized features for all four sites, plus MGB under CLSI 2022, from the AMIA 2026 definitions linked through the loaders.
- ARMD's "Congestive heart failure" category is a catch-all. Using CCSR "Heart failure" alone cuts heart failure prevalence from 81.5% to 15.6% at Stanford, 41.6% to 9.0% at MGB, 61.4% to 9.6% at UTSW and 72.6% to 14.1% at ECUH.
- Drug capture differs by release. Any prior ceftriaxone is 6.9% at Stanford, 28.8% at MGB, 3.9% at UTSW and 37.8% at ECUH. UTSW has no TMP-SMX or amoxicillin/clavulanate, and ECUH has no piperacillin/tazobactam.
- One resistance phenotype definition gives prior ESBL of 3.3%, 7.3%, 3.3% and 3.9%, and prior MRSA of 1.5%, 5.1%, 2.7% and 2.5%.
- Flagged: the AMIA notebook joined three MGB tables on mismatched culture IDs, so 14 of its 45 features and 4 interactions built on them came from unrelated cultures.
