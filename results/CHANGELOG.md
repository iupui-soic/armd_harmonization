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

## 2026-09-11 01_build_cohorts (results/01_build_cohorts/20260911T210800Z, commit 4f1782d)

- Built E1, E2 and sensitivity cohorts on the derivation partition. E1 has 62,558, 108,572, 40,215 and 69,333 cultures at Stanford, MGB, UTSW and ECUH, and E2 has 595,125, 775,796, 602,062 and 672,137.
- Derivation prevalence is 11.1%, 13.3%, 7.9% and 9.2% under E1 and 1.17%, 1.86%, 0.53% and 0.95% under E2. Negative cultures make up 84.0%, 78.7%, 87.0% and 86.7% of ordered derivation cultures.
- E2 is younger and has fewer prior cultures than E1 at every site. Blood specimens rise (SMD 0.81 to 1.52) and urine specimens fall (SMD -0.83 to -1.41).
- S2 now counts antibacterial records 1 to 3 days before the order. It removes 7.1%, 20.8%, 5.0% and 2.2% of E2 derivation cultures, where counting day 0 had removed 45.1% and 60.9% at Stanford and MGB. Flow steps before the partition no longer count labels.
- Flagged: the plan expected E2 to add low-acuity outpatients, but E2 is more inpatient than E1 at Stanford, MGB and UTSW (SMD 0.14 to 0.78) and only ECUH shifts toward outpatients. Cultures S2 removes have E2 prevalence of 1.03%, 1.05%, 0.52% and 1.13%, against 1.18%, 2.07%, 0.53% and 0.95% for those kept.
