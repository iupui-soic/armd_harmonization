# Results changelog

## 2026-09-11 00_gate (results/00_gate/20260911T182159Z, commit d96f9fc)

- All four sites pass the Phase 0 gate, so the verdict is go. Three MGB tables only link by patient and order date.
- E1 composite prevalence is 11.1% at Stanford, 13.3% at MGB, 7.9% at UTSW and 9.3% at ECUH. MGB rises to 18.2% under CLSI 2022.
- E2 decision-time prevalence is 1.17%, 1.85%, 0.53% and 0.96%, about a tenfold drop from E1 at every site.
- Mapping UTSW's "Pip/tazobactam" raised UTSW E1 evaluable cultures from 49,615 to 50,327 and resistant cultures from 3,608 to 3,973.
- Flagged: mismatched MGB IDs in three tables, MGB urine and respiratory positivity two to three times the other sites, and 14,397 UTSW positives with no organism or AST.

## 2026-09-11 features (results/features_<site>/20260911T2007*, commit 2fefb4e)

- Built the 48 harmonized features for all four sites, plus MGB under CLSI 2022, from the AMIA 2026 definitions linked through the loaders.
- ARMD's "Congestive heart failure" category is a catch-all. Using CCSR "Heart failure" alone cuts heart failure prevalence from 81.5% to 15.6% at Stanford, 41.6% to 9.0% at MGB, 61.4% to 9.6% at UTSW and 72.6% to 14.1% at ECUH.
- Drug capture differs by release. Any prior ceftriaxone is 6.9% at Stanford, 28.8% at MGB, 3.9% at UTSW and 37.8% at ECUH. UTSW has no TMP-SMX or amoxicillin/clavulanate, and ECUH has no piperacillin/tazobactam.
- One resistance phenotype definition gives prior ESBL of 3.3%, 7.3%, 3.3% and 3.9%, and prior MRSA of 1.5%, 5.1%, 2.7% and 2.5%.
- Flagged: the AMIA notebook joined three MGB tables on mismatched culture IDs, so 14 of its 45 features and 4 interactions built on them came from unrelated cultures.

## 2026-09-11 01_build_cohorts (results/01_build_cohorts/20260911T210800Z, commit 060cc8b)

- Built E1, E2 and sensitivity cohorts on the derivation partition. E1 has 62,558, 108,572, 40,215 and 69,333 cultures at Stanford, MGB, UTSW and ECUH, and E2 has 595,125, 775,796, 602,062 and 672,137.
- Derivation prevalence is 11.1%, 13.3%, 7.9% and 9.2% under E1 and 1.17%, 1.86%, 0.53% and 0.95% under E2. Negative cultures make up 84.0%, 78.7%, 87.0% and 86.7% of ordered derivation cultures.
- E2 is younger and has fewer prior cultures than E1 at every site. Blood specimens rise (SMD 0.81 to 1.52) and urine specimens fall (SMD -0.83 to -1.41).
- S2 now counts antibacterial records 1 to 3 days before the order. It removes 7.1%, 20.8%, 5.0% and 2.2% of E2 derivation cultures, where counting day 0 had removed 45.1% and 60.9% at Stanford and MGB. Flow steps before the partition no longer count labels.
- Flagged: the plan expected E2 to add low-acuity outpatients, but E2 is more inpatient than E1 at Stanford, MGB and UTSW (SMD 0.14 to 0.78) and only ECUH shifts toward outpatients. Cultures S2 removes have E2 prevalence of 1.03%, 1.05%, 0.52% and 1.13%, against 1.18%, 2.07%, 0.53% and 0.95% for those kept.

## 2026-09-26 qc (results/qc_<site>/20260926T1204*, commit 5909480, documentation edits uncommitted)

- The QC report now profiles repeat and same-day cultures. Stanford, UTSW and ECUH have no culture 1 to 13 days after another from the same patient, and gaps start at exactly 14 days.
- MGB applies no repeat rule. 35.0% of its cultures follow another culture of the same patient by 1 to 13 days, and 155,950 follow by one day.
- At the Dryad sites every same-day culture of a patient shares one order time. Blood cultures come in same-day sets: 90%, 86%, 91% and 93% of them at Stanford, MGB, UTSW and ECUH.
- The chained 14-day rule changes nothing at the three Dryad sites and keeps 632,625 of 972,798 MGB cultures (65.0%). MGB positivity rises from 21.3% to 26.2% and its blood share falls from 64.5% to 57.4%.
- Flagged: MGB's documented rule reads as the reverse of the other releases, and its data show no rule at all. Asked of the maintainers.

## 2026-09-26 01_build_cohorts (results/01_build_cohorts/20260926T120230Z, commit 5909480, documentation edits uncommitted)

- Added cohorts under the proposed chained 14-day rule with one row per specimen-day (SAP R7, R8), plus order events (S6) and first index order per patient (S7). Per-culture cohorts are unchanged.
- Under the proposed unit, E1 prevalence is 11.1%, 10.7%, 8.1% and 9.2% at Stanford, MGB, UTSW and ECUH, and E2 prevalence is 1.48%, 2.74%, 0.65% and 1.17%. E2 rises because each blood culture set now counts once.
- MGB E1 prevalence falls from 13.3% to 10.7% once repeat cultures go, so repeats after a resistant culture inflated it. Keeping each patient's first index order lowers E1 prevalence to 8.0%, 6.8%, 5.8% and 5.1%.
- The E2 minus E1 inpatient SMD moves from +0.18 to +0.04 at Stanford and from +0.14 to +0.02 at UTSW, so most of the earlier "E2 is more inpatient" surprise came from counting blood sets separately. MGB stays more inpatient (+0.78 to +0.47) and ECUH shifts further to outpatients (-0.32 to -0.45).
- Flow counts within 10 of another published count are now withheld. This hides the S5 counts at MGB and UTSW, which differed from the primary cohorts by fewer than 11.
