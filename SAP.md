# Statistical analysis plan: Paper A

**Status: Draft v0.1 for review. Not frozen.**

**Freezing.** The project lead approves the review items in section 16. The status line then becomes
"Frozen v1.0" with a date, and the commit is tagged `sap-v1.0`, with an optional OSF registration.
No model is trained and no test partition is read before that tag exists. Any change after freezing
goes into the deviations log with a date and a reason.

Working title: *Clinical prediction models for antimicrobial resistance estimate the wrong quantity:
an estimand analysis across four US health systems.*

## 1. Objective

The paper quantifies the gap between the estimand that published EHR models of ESBL phenotype
target (E1) and the estimand that deployment at culture order needs (E2). It measures what happens
when a model is trained under E1 and used under E2, and tests whether recalibration closes the gap.

## 2. Data

| site | release | version |
|---|---|---|
| Stanford | Dryad doi:10.5061/dryad.jq2bvq8kp | 20 |
| MGB | PhysioNet doi:10.13026/2r5k-b955 | 1.0.0 |
| UTSW | Dryad doi:10.5061/dryad.0rxwdbsd5 | 5 |
| ECUH | Dryad doi:10.5061/dryad.7sqv9s55x | 4 |

- **Population.** Adults with urine, blood or respiratory cultures, as defined by each release. Every
  release excludes repeat cultures from the same patient within 14 days.
- **Harmonization.** armd-common at the freezing commit. `docs/HARMONIZATION.md` documents every
  decision.
- **Unit of analysis.** The culture episode. A patient can contribute several cultures.
- **Time zero.** The culture order time. MGB records only the order date.
- **Splits.** The patient-level manifests in `benchmarks/manifests/`, pinned by `SHA256SUMS`. One of
  five outer folds is the test partition, and derivation patients have five inner folds. Patients
  were stratified by ESBL history.

## 3. Estimands

| | estimand | conditioning | identification |
|---|---|---|---|
| E1 | status quo | ordered, grew Enterobacterales, at least one sentinel agent reported | direct |
| E2 | decision time | ordered | direct, through negative cultures |
| E3 | population | empiric prescribing event | not identified in ARMD; bounds only |

- **Outcome Y.** A culture is positive if any Enterobacterales isolate is resistant (R, laboratory
  interpretation) to at least one of ceftriaxone, ceftazidime, cefepime or piperacillin/tazobactam.
  E1 uses `labels.esbl_phenotype_composite` and E2 uses `labels.esbl_decision_time`.
- **E2 zeros.** E2 sets Y = 0 for negative cultures and for positive cultures without an
  Enterobacterales isolate.
- **Non-evaluable cultures** are excluded from both estimands and counted in the flow. These are
  Enterobacterales cultures with no sentinel result, and cultures whose isolates are unidentified.
- **E3.** ARMD contains cultured patients only, so E3 is analysed as bounds and a tipping point
  (section 7.5).
- **Estimand statement.** Each result names its estimand: population, variable, handling of
  non-evaluable cultures, and summary measure. TRIPOD+AI does not require this; the paper adds it.

## 4. Predictors

- **Features.** Every site and every task uses all 48 harmonized features (`features.FEATURES_48`),
  including cross-site pairs. This was decided on 2026-09-11. Per-site capture of each feature is
  reported alongside the results.
- **Timing.** All features are available at time zero. History features count only events at least
  one day before the order. Nursing home visits keep day 0.
- **MGB.** Prior resistance features use the laboratory interpretation.
- **Missing values.** These follow the feature definitions and are not imputed:
  - absent history is 0;
  - a missing age band takes the middle band;
  - a missing ADI score is 0, with an indicator.

## 5. Cohorts and flow

- **Cohort builder.** Cohorts come from `cohorts.build_cohort`. The flow per site runs from ordered
  cultures, to cultures evaluable under the estimand, to the partition, and it feeds a TRIPOD+AI
  diagram. Labels are counted only from the partition step on, so no step before evaluation counts
  outcomes in the test partition.
- **Descriptive timing.** Descriptive statistics before freezing (`analysis/01_build_cohorts.py`) use
  the derivation partition only. The full-cohort Table 1 is produced after freezing, together with
  the test evaluation.

## 6. Models

**Primary model.** XGBoost with the tuned hyperparameters of the AMIA 2026 ARMD-MGB model, held fixed
with no re-tuning:

| parameter | value |
|---|---|
| n_estimators | 700 |
| max_depth | 6 |
| learning_rate | 0.02001 |
| min_child_weight | 6 |
| subsample | 0.8807 |
| colsample_bytree | 0.5222 |
| gamma | 2.596 |
| reg_alpha | 8.608 |
| reg_lambda | 5.384 |
| objective | binary:logistic |

Two choices need review (item R1):
- **Class weighting.** `scale_pos_weight` is set to 1 instead of the AMIA balanced ratio (5.86).
  Calibration and net benefit are primary here, and reweighting shifts the probability scale. The AMIA
  paper found discrimination unchanged without imbalance correction (AUROC 0.762 against 0.763).
- **Fixed hyperparameters.** They were tuned on the AMIA MGB cohort, where 18 of 45 features had been
  joined to the wrong cultures. They are still kept fixed, to limit analyst degrees of freedom.

**Secondary model.** L2-penalized logistic regression with C = 1.0. Features are standardized with
derivation means and standard deviations, and there is no class weighting.

**Training.**
- One model per site and training estimand, fitted on all derivation patients with seed 20260911.
- Cross-fitted predictions for recalibration and threshold selection use the frozen inner folds.
- Test predictions are cached, and the bootstrap resamples predictions rather than refitting models.

## 7. Analyses

### 7.1 Analysis 1: prevalence and composition shift

For each site:
- number of cultures, patients and label prevalence under E1 and E2;
- the share of ordered cultures that were negative or grew no Enterobacterales;
- covariate distributions with standardized mean differences between E2 and E1 for age 65 or older,
  sex, inpatient, outpatient and ED setting, specimen type, prior culture history, prior ESBL, 90-day
  third-generation cephalosporin and fluoroquinolone exposure, nursing home, ADI missingness,
  comorbidity count and number of prior organisms.

### 7.2 Analysis 2: the train and evaluate grid (headline)

| training estimand | evaluated on E1 test | evaluated on E2 test |
|---|---|---|
| E1 | reproduces published practice | what deployment does |
| E2 | reported for completeness | the corrected model |

The primary comparison is the E1-trained model against the E2-trained model, both evaluated on the
E2 test partition.

### 7.3 Analysis 3: recalibration

- **Maps.** The E1-trained model is recalibrated toward E2 with Platt scaling (primary) and isotonic
  regression (secondary).
- **Fitting.** The maps are fitted on cross-fitted predictions. For each inner fold, an E1 model
  trained on the other folds predicts that fold's E2 derivation cultures, then the map is applied to
  test predictions.
- **Reporting.** Net benefit, calibration and discrimination are reported before and after
  recalibration. If recalibration closes the gap, the paper says so.

### 7.4 Analysis 4: cascade reporting (layer T)

- **Reporting rates.** Per site, the share of Enterobacterales cultures with each sentinel agent
  reported as S, I or R.
- **Non-evaluable cultures.** Compared with evaluable cultures on setting, organism group, specimen
  and prior history, using standardized mean differences.
- **MGB cross-check.** Among Enterobacterales cultures, compare ESBL enzyme test positivity and the
  CLSI 2022 phenotype between cultures with and without a reported sentinel agent.
- **Imputation.** The estimand is reported both on evaluable cultures only and with strict implied
  susceptibility applied.

### 7.5 Analysis 5: E3 bounds

- **Bounds.** Manski worst-case bounds for population prevalence, PPV and NPV at the primary threshold.
  These assume that a share u of empiric prescribing events has no culture, with u in {0.10, 0.25,
  0.50}, and leave resistance among the uncultured unrestricted.
- **Tipping point.** The smallest relative risk of ESBL among uncultured compared with cultured
  patients at which net benefit at the primary threshold changes sign.
- **Reporting.** One figure, framed as the remaining unidentified layer. The parameter ranges need
  justification from the literature (item R4).

### 7.6 Analysis 6: multi-site replication

- **Replication.** Analyses 7.1 to 7.3 run at every site and are reported site by site, without
  pooling.
- **Cross-site transport (secondary).** Train on the derivation partition of site A and evaluate on
  the test partition of site B, for every ordered pair, under E1 and E2.

## 8. Metrics

- **Primary: net benefit.** NB(t) = TP/n - FP/n x t/(1 - t), at thresholds t of 0.05, 0.10, 0.15 and
  0.20, under E2. Standardized net benefit (NB divided by prevalence) is reported alongside, because
  E2 prevalence is near 1%. Decision curves cover t from 0.01 to 0.40.
- **Secondary: calibration.** Calibration intercept (calibration in the large) and slope, from
  logistic recalibration on the logit of the prediction, plus calibration curves and the Brier score.
- **Secondary: discrimination.** AUROC and AUPRC.
- **Operating point.** Sensitivity, specificity, PPV, NPV and net benefit at a 90%-sensitivity
  threshold. The threshold comes from the cross-fitted derivation predictions of the model's training
  cohort, as it would be set before deployment, and is then applied unchanged to the test partition
  (item R3).
- **AUROC artifact.** AUROC can rise from E1 to E2 because E2 adds many easy negatives. That is a
  prevalence and case-mix artifact, not an improvement. The abstract states this, and the results
  lead with calibration, PPV and net benefit.

## 9. Uncertainty

- **Bootstrap.** 2,000 patient-clustered resamples of each site's test partition, keeping all of a
  patient's cultures together, with seed 20260911 and percentile 95% intervals.
- **Differences.** Differences between cells evaluated on the same cohort use paired resamples.
- **Cross-site pairs.** These resample patients of the target site.

## 10. Subgroups

- **Setting:** `setting_primary`. ED and urgent care are only distinguishable at MGB.
- **Specimen type.**
- **Prior culture history:** any prior organism or prior resistant organism at least one day before
  the order.
- **ADI:** within-site quintile, with cutpoints from the derivation partition and missing as its own
  group.

Subgroup results are descriptive, with no subgroup hypothesis tests. Counts from 1 to 10 are
suppressed.

## 11. Sensitivity analyses

| id | variation |
|---|---|
| S1 | Intermediate counts as resistant (non-susceptibility) |
| S2 | Exclude cultures with an antibacterial record 1 to 3 days before the order, because prior antibiotics can suppress growth in negative cultures |
| S3 | CLSI 2022 interpretation at MGB, for labels and prior resistance features |
| S4 | IDSA definition for E1, and the MGB enzyme label for E1 on the MGB subset that has ESBL tests |
| S5 | Strict implied susceptibility applied before labelling |

S2 counts whole days, because every release records exposure times in days (item R6).
- **Day 0 does not count.** A day 0 record cannot be placed before the order. Stanford has day 0
  records for 41% of ordered cultures. At MGB, 37% have a last dose on the order day and 24% a last
  dose after it, and the two groups overlap. These include treatment started after the culture was
  drawn. UTSW and ECUH hold no day 0 records.
- **Effect.** Counting day 0 and later would exclude 45% of ordered cultures at Stanford and 61% at
  MGB, against 5% at UTSW and 2% at ECUH. Records 1 to 3 days before the order exist for 7%, 21%, 5%
  and 2%.
- **Missed exposure.** Doses on the order day before the culture, and MGB courses that began earlier
  and continued past the order.
- **Source.** The table the antibiotic features use, so the drug list per release matches section 4.

## 12. What counts as a null result

This is proposed and needs review (item R2).
- **Per site.** The estimand gap is negligible at a site if two conditions hold on the E2 test
  partition:
  - At every primary threshold, the difference in standardized net benefit between the E1-trained
    and E2-trained XGBoost models is within ±0.05, with its 95% interval within ±0.10.
  - The E1-trained model's calibration slope interval includes 1 and its intercept interval
    includes 0.
- **Overall.** The paper reports a null result if the gap is negligible at every site that passed the
  gate.

A null result is published with the same analyses and a retitled framing.

## 13. Reporting and provenance

- **Reporting guideline.** TRIPOD+AI, with the estimand statement added.
- **Provenance.** Every number traces to a `results/**/manifest.json`, which records the git commit,
  seed, input digests, package versions and row counts.
- **Small cells.** Counts from 1 to 10 are suppressed.

## 14. Limitations to state

- MGB may exclude cultures that grew an organism without AST, which would shrink its E2 denominator.
  This is being raised with the maintainers.
- Three MGB tables link by patient and order date, because their culture IDs do not match the cohort.
  MGB order times are dates only.
- Antibiotic class, procedure and comorbidity-count features depend on what each release records.
- Sex is inferred at Stanford, UTSW and ECUH.
- ARMD's "Congestive heart failure" category is a catch-all. Heart failure uses the CCSR category.
- Comparisons with the AMIA 2026 model must note that 18 of its 45 features were joined to the wrong
  cultures.

## 15. Analysis scripts

| script | role | reads test partition |
|---|---|---|
| analysis/00_gate.py | Phase 0 feasibility gate | no |
| analysis/01_build_cohorts.py | cohorts and Analysis 1 on derivation | no |
| analysis/02_train_e1.py | models for both estimands, cross-fitted predictions, recalibration maps | no |
| analysis/03_evaluate_ladder.py | Analyses 2, 3 and 6, full-cohort Table 1 | yes, the only one |
| analysis/04_cascade_sensitivity.py | Analysis 4 and sensitivity analyses | no |
| analysis/05_e3_bounds.py | Analysis 5 | no |
| analysis/06_figures.py | figures from saved results | no |

S1 to S5 are evaluated on the test partition inside `03_evaluate_ladder.py`. Script 04 prepares their
cohorts and derivation-side summaries only.

## 16. Review items before freezing

1. **R1.** XGBoost with `scale_pos_weight = 1`, and the AMIA hyperparameters kept fixed rather than
   re-tuned.
2. **R2.** The null-result margins in section 12.
3. **R3.** Taking the 90%-sensitivity threshold from the training cohort's cross-fitted derivation
   predictions.
4. **R4.** The uncultured shares and relative-risk ranges for E3, with references.
5. **R5.** Derivation-only descriptive statistics before freezing, and the full-cohort Table 1 at
   evaluation.
6. **R6.** S2 counts antibacterial records 1 to 3 days before the order and ignores day 0 (section 11).

## 17. Deviations log

None yet.
