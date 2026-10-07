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

## 2026-09-27 04_cascade_sensitivity (results/04_cascade_sensitivity/20260927T100038Z, commit 56e7d84)

- Among derivation Enterobacterales cultures, a sentinel agent is reported for 93.7% at Stanford, 99.9% at MGB, 89.8% at UTSW and 99.97% at ECUH. Cefepime is reported for 31.9% at Stanford and ceftazidime for 11.7% at UTSW.
- Cultures without a sentinel result are 6.3% at Stanford and 10.2% at UTSW, against 130 cultures at MGB and 18 at ECUH. At Stanford they are more often urine (SMD +0.44), inpatient (+0.39) and E. coli (+0.29), with less prior ESBL (-0.29). That fits cascade reporting, which hides broad agents when narrow ones test susceptible.
- The MGB enzyme and CLSI 2022 cross-check cannot inform layer T. Fewer than 11 of the 130 MGB cultures without a sentinel result had an ESBL test or a CLSI 2022 sentinel result.
- Applying every implied-susceptibility rule makes 3,071 of 4,245 not-tested Stanford cultures and 2,531 of 4,584 UTSW cultures evaluable. It adds no resistant culture at Stanford, so Stanford E1 prevalence falls from 11.15% to 10.62%. Strict rules change fewer than 11 cultures at any site, and neither rule set changes ECUH.
- Flagged: the rules only impute susceptible results from susceptible ones, so the lower prevalence is built in. It does not show that the not-tested cultures are susceptible.

## 2026-09-27 02_train_e1 (results/02_train_e1/20260927T113507Z, commit 1b99f33, tag sap-v1.0)

- Trained XGBoost and logistic models under the frozen plan on the derivation partition, with the chained 14-day rule and one row per specimen-day. E1 has 62,468, 93,955, 39,458 and 68,410 units at Stanford, MGB, UTSW and ECUH, and E2 has 471,928, 364,930, 487,157 and 542,116, with the same 6,961, 10,016, 3,179 and 6,323 resistant units.
- The 90%-sensitivity threshold is 0.057, 0.044, 0.041 and 0.046 for the E1-trained XGBoost model and 0.0063, 0.0153, 0.0030 and 0.0049 for the E2-trained model. Both sit below 0.05, the smallest planned net-benefit threshold.
- Platt maps from E1 toward E2 have intercepts of -1.97, -1.75, -2.73 and -1.71 and slopes of 1.10, 0.95, 0.99 and 1.19 for XGBoost. That is close to the shift in log odds of prevalence, with the ranking kept.
- The recalibrated E1 thresholds on E2 are 0.0063, 0.0094, 0.0029 and 0.0049, close to those of the E2-trained models.
- Flagged: slopes near 1 suggest recalibration may recover much of the gap, which is the null judged in SAP section 12. These are cross-fitted derivation predictions and say nothing yet about test performance.

## 2026-09-27 03_evaluate_ladder (results/03_evaluate_ladder/20260927T114044Z, commit dc77b63, tag sap-v1.0)

- First read of the test partition, with 2,000 patient bootstrap draws. On the E2 test cohort, the E1-trained XGBoost model has calibration intercepts of -2.15, -1.67, -2.73 and -2.03 at Stanford, MGB, UTSW and ECUH. Its standardized net benefit is negative at every planned threshold (-1.88, -0.79, -3.99 and -0.89 at 0.10).
- Platt recalibration restores calibration (intercepts -0.02 to 0.04, slopes 0.95 to 1.03) but not ranking. AUROC stays at 0.686, 0.680, 0.653 and 0.753, against 0.812, 0.846, 0.795 and 0.860 for the E2-trained model. The paired differences run from -0.11 to -0.17, with every interval below 0.
- The recalibrated E1 model trails the E2 model in standardized net benefit at 0.10 by 0.10, 0.10, 0.03 and 0.09. Under SAP section 12 the gap is not negligible at any site for either model class, so the result is not null.
- Flagged: the plan expected AUROC to rise from E1 to E2 as an artifact of easy negatives. For the E1-trained model it falls instead (0.724 to 0.686 at Stanford, 0.770 to 0.680 at MGB). Within urine the AUROC gap shrinks to 0.06, so much of it lies between specimen types. The pooled Platt map over-predicts within blood (standardized net benefit at 0.10 of -0.52 at MGB and -0.35 at ECUH), and among first index orders (S7) the E1 model's AUROC falls to 0.50 to 0.58.
- Flagged: a logic bug, fixed in d5099dd, makes this run's summary read the null verdict as undetermined at UTSW and overall. One failed check settles the verdict, and the stored checks give "not negligible" at all four sites.

## 2026-09-27 05_e3_bounds (results/05_e3_bounds/20260927T143411Z, from the 03 run above)

- At the primary threshold of 0.10, net benefit among cultured events is 0.0024, 0.0078, 0.0002 and 0.0021 for the E2-trained XGBoost model at Stanford, MGB, UTSW and ECUH. For the uncalibrated E1-trained model it is -0.028, -0.022, -0.026 and -0.011.
- For the uncalibrated E1 model to break even, uncultured events would need 2.5 to 81 times the cultured risk, across sites and uncultured shares from 0.10 to 0.75.
- The E2-trained model's net benefit turns negative only at uncultured shares of 0.50 or 0.75, and only if uncultured events carry 0.13 to 0.65 times the cultured risk.
- Worst-case PPV bounds are wide at every share. At u = 0.25 they run from 0.03 or less up to 0.93 or more for the E2-trained model, so E3 PPV is not usefully bounded without assumptions.
- Flagged: at UTSW the recalibrated E1 model flags no unit at 0.10, so its PPV bounds span 0 to 1 by construction.

## 2026-09-27 06_figures (results/06_figures/20260927T143412Z)

- Drew figures 1, 2, 3, 5 and 6 and table 4 from the runs above, with the chained 14-day rule and specimen-day units, as PDF, PNG and CSV. Figure and table files stay local.

## 2026-09-27 02_train_e1, urine only (results/02_train_e1/20260927T190025Z, commit 6b431eb, exploratory E-3)

- Trained XGBoost and logistic models on the urine units of the derivation partition only, for the exploratory amendment in SAP section 17 (E-3).
- Urine E1 has 56,176, 87,021, 37,541 and 62,613 units at Stanford, MGB, UTSW and ECUH, with 6,065, 8,849, 2,845 and 5,652 resistant. Urine E2 has 292,939, 199,057, 351,972 and 376,977 units.
- Platt maps from urine E1 toward urine E2 have XGBoost intercepts of -1.14, -0.97, -1.40 and -1.17, about half the pooled shift, because resistance is more common in urine than across all specimens.
- The 90%-sensitivity threshold of the urine E2 XGBoost model is 0.0095, 0.0201, 0.0037 and 0.0076.
- Flagged: this run was chosen after the primary test results were seen, so it and the evaluation it feeds are exploratory.

## 2026-09-27 03_exploratory (results/03_exploratory/20260927T190415Z, commit 23a55e6, SAP section 17 E-1 to E-5)

- E-1: the approximated Tumbarello score ranks better at order time (AUROC 0.70 to 0.75) than on E1 (0.61 to 0.67). At its published cut-off it flags 234 to 397 of every 1,000 orders, with PPV of 1.5% to 4.9% and negative standardized net benefit at 0.10 at every site. At the score's sensitivity the E2 XGBoost model flags 87 to 164 fewer orders per 1,000.
- E-2: prior ESBL alone flags 29 to 56 orders per 1,000 with sensitivity of 25% to 39%, and the E2 model saves only 4 to 18 flags per 1,000 at that sensitivity. Adding recent beta-lactam or fluoroquinolone exposure flags 128 to 233 per 1,000, with negative net benefit at 0.10 at every site. R3, a five-predictor logistic model trained on all orders, is calibrated (intercepts 0.00 to 0.05), has AUROC 0.70 to 0.81, and matches the XGBoost model's flags at 0.10 within 2 per 1,000.
- E-3: within urine, urine-trained E1 models reach AUROC 0.72 to 0.79 at order time against 0.79 to 0.85 for urine-trained E2 models, a paired gap of 0.055 to 0.064. After Platt recalibration the standardized net benefit gap at 0.10 is 0.01 to 0.04.
- E-4: at the 10% threshold, the E1 model starts broad coverage in 165 to 366 of every 1,000 orders to cover 2.9 to 17.1 resistant cultures, and the E2 model in 6 to 38 orders to cover 0.8 to 10.8. To catch 90% of resistant cultures, R3 flags 576 to 744 per 1,000 against 436 to 549 for the E2 XGBoost model. E-5: every conclusion above holds at the three Dryad sites alone.
- Flagged: within urine, the E1 model's AUROC rises at order time at three sites, as the plan originally expected. The pooled fall comes from mixing specimen types. Much of the E2 advantage comes from training on all orders rather than from model complexity. All of E-1 to E-5 is exploratory.

## 2026-10-07 02_train_e1, blood only (results/02_train_e1/20261007T101226Z, commit dc8f0c2, exploratory E-7)

- Trained XGBoost and logistic models on the blood units of the derivation partition only, for the second exploratory amendment in SAP section 17 (E-7).
- Blood E1 has 4,384, 4,614, 1,514 and 4,937 units at Stanford, MGB, UTSW and ECUH, with 649, 674, 228 and 533 resistant. Blood E2 has 127,173, 153,193, 116,971 and 146,883 units.
- Platt maps from blood E1 toward blood E2 have XGBoost intercepts of -2.89, -3.15, -2.82 and -2.71, a larger shift than the pooled maps, because resistance is rarer among blood orders.
- The 90%-sensitivity threshold of the blood E2 XGBoost model is 0.0028, 0.0023, 0.0014 and 0.0022, about half the share of blood orders that are resistant.
- Flagged: the Platt slopes are 1.37, 1.27, 1.96 and 1.38, well above the pooled and urine slopes near 1. With 228 to 674 resistant units the blood models are noisy, UTSW most of all. This run was chosen after test results were seen, so it is exploratory.

## 2026-10-07 03_exploratory_2 (results/03_exploratory_2/20261007T101904Z, commit b5ae1b2, SAP section 17 E-6 to E-8)

- E-6: with ceftriaxone resistance as the outcome, E2 prevalence is 1.32%, 2.56%, 0.44% and 1.10% at Stanford, MGB, UTSW and ECUH. On E2 the E1-trained XGBoost model has calibration intercepts of -2.27, -1.76, -3.09 and -2.13 and AUROC of 0.706, 0.696, 0.686 and 0.760, against 0.821, 0.855, 0.821 and 0.866 for the E2-trained model. After Platt recalibration it trails in standardized net benefit at 0.10 by 0.11, 0.11, 0.03 and 0.09. Every primary conclusion holds.
- E-7, blood models: blood-trained E2 models reach AUROC of 0.701, 0.777, 0.741 and 0.753 on blood test units, no better than the pooled E2 model within blood (0.696, 0.778, 0.739, 0.758). Blood-trained E1 models reach 0.652, 0.747, 0.680 and 0.731, a paired gap of 0.022 to 0.061. At a threshold of 0.01 the blood E2 model flags 60, 59, 0 and 39 of every 1,000 blood orders and catches 32%, 38%, 0% and 35% of resistant blood cultures, with standardized net benefit of 0.21, 0.25, 0.00 and 0.24. Almost nothing is flagged at 0.05.
- E-7, thresholds and specimen type: with the 90%-sensitivity threshold set within each specimen type, the pooled E2 model flags 595, 633, 584 and 443 of every 1,000 urine orders against 700, 787, 709 and 639 for the pooled E1 model, the same as the urine-trained models, and 699, 544 and 671 of every 1,000 blood orders at Stanford, MGB and ECUH against 765, 702 and 763. UTSW blood is withheld for small cells. Specimen type alone reaches AUROC of 0.622, 0.688, 0.585 and 0.606 on E2 when fitted on E2, and 0.378, 0.323, 0.416 and 0.396 when fitted on E1, below chance, because blood carries more resistance than urine among positive cultures and less among orders. The E2 XGBoost model adds 0.16 to 0.25 over specimen type alone.
- E-8: at 6 points or more the approximated Tumbarello score flags 88, 87, 40 and 55 of every 1,000 orders, with sensitivity of 23%, 23%, 15% and 24% and PPV of 3.9%, 7.2%, 2.5% and 5.3%. Standardized net benefit at 0.10 stays negative at every site (-0.41, -0.10, -0.52, -0.24). At that sensitivity the E2 XGBoost model flags 31 to 74 fewer orders per 1,000.
- Flagged: ceftriaxone is not reported for every Enterobacterales culture, so the ceftriaxone label leaves 15%, 15%, 9.5% and 1.9% of E1 test units not evaluable. The blood thresholds set on derivation units reach only 85% to 87% sensitivity on test, where each site has 71 to 180 resistant blood units. Within a specimen type the E1 to E2 AUROC gap is 0.02 to 0.06, about half the pooled gap of 0.11 to 0.17, so specimen mix explains about half of it and not most. All of E-6 to E-8 is exploratory.
