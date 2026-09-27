# armd-common harmonization decisions

This document records every choice that turns four ARMD releases into one schema, with the reason
for each. It is the substance of the resource preprint. Counts come from aggregate profiling, the
Phase 0 gate (`results/00_gate/`) and the per-site QC reports (`results/qc_<site>/`). Counts from 1
to 10 are never reported.

## Releases

| site | source | version | license | downloaded |
|---|---|---|---|---|
| Stanford | Dryad doi:10.5061/dryad.jq2bvq8kp | 20 (2025-10-22) | CC0-1.0 | 2026-02-26 |
| MGB | PhysioNet doi:10.13026/2r5k-b955 | 1.0.0 (2025-12-05) | PhysioNet Credentialed Health Data License 1.5.0 | 2026-02-23 |
| UTSW | Dryad doi:10.5061/dryad.0rxwdbsd5 | 5 (2025-09-26) | CC0-1.0 | 2026-09-11 |
| ECUH | Dryad doi:10.5061/dryad.7sqv9s55x | 4 (2025-11-10) | CC0-1.0 | 2026-09-11 |

`src/armd_common/data_sources.yaml` pins the size and SHA-256 of every file. Stanford matches all
17 Dryad digests and MGB matches all 12 entries in its SHA256SUMS.txt. UTSW and ECUH were verified
against the Dryad digests during download.

## Staging

Each CSV becomes Parquet once, under the site folder in `/data0`. Staging keeps raw column names and
changes only three things.

1. The tokens `""`, `NA`, `null`, `Null` and `NULL` become NULL. Stanford, UTSW and ECUH write
   missing values as `Null`.
2. `anon_id`, `pat_enc_csn_id_coded` and `order_proc_id_coded` are always VARCHAR, so joins never
   depend on sniffed types.
3. Standard CSV quoting is set explicitly. The releases quote only the rare fields that need it, so
   the sniffer often sees no quotes in its sample. It then fails on a quoted comma (UTSW cohort, 2
   lines; ECUH antibiotic class exposure, 518 lines) or silently keeps the quote characters (Stanford
   microbial resistance, 2,841 lines). ECUH comorbidity has 1.3 million quoted lines.

If a column fails to parse under its sniffed type, only that column is retried as VARCHAR, and the
retry is recorded in `_conversion.json`. No column needed a retry at any site. Error messages from
the CSV reader quote raw lines, so the converter reports only a line number and a column name.

## Keys

`culture_id` is `<site>:<anon_id>:<order_proc_id_coded>` and `patient_id` is `<site>:<anon_id>`.

- **The patient is part of the culture key.** MGB reuses `order_proc_id_coded` across patients:
  2,633 values belong to two patients each.
- **MGB rows repeat once per linked encounter.** 118,408 MGB cultures link to more than one
  encounter, and every cohort and feature row is repeated for each. Loaders collapse these rows
  before counting or labelling.
- **Three MGB tables need a different link.** In `comorbidity`, `nursing_home_visits` and
  `prior_procedures`, encounter and culture IDs come from a different mapping than the cohort. No
  ordering of the three ID columns links any of their rows to cohort cultures. Yet every
  (anon_id, order date) pair in those files exists in the cohort: 407,018 of 407,018, 386,657 of
  386,657 and 78,259 of 78,259. The MGB loader links them by patient and order date. Their contents
  are measured relative to the order date, so a patient's same-day cultures share the same values
  and nothing is lost at MGB's day resolution. This looks like a release defect and should be
  reported to the dataset maintainers.
- **The MGB defect is in the release, not our copy.** It reproduces on the raw CSVs read with
  Python's csv module, outside the staging pipeline, and every row has the header's field count.
  All 12 local files match PhysioNet's SHA256SUMS.txt and the byte sizes in the server's directory
  listing. A fresh copy of SHA256SUMS.txt fetched from PhysioNet on 2026-09-11 is byte-identical to
  the one in the local release, and all 12 files verify against it. The three files also share one alternative ID set: every triple in
  `nursing_home_visits` (446,041) and `prior_procedures` (92,171) appears in `comorbidity`. Random
  corruption would not produce consistent IDs across files.

## Culture positivity and the E2 denominator

Positivity is taken from each release: `was_positive` at the Dryad sites and the absence of `neg_cx`
at MGB. The definitions are not the same, and this matters most for the decision-time estimand.

| site | cultures | negative | urine positive | blood positive | respiratory positive |
|---|---|---|---|---|---|
| Stanford | 751,075 | 630,798 | 23.8% | 4.5% | 21.1% |
| MGB | 972,798 | 765,693 | 48.3% | 5.3% | 68.3% |
| UTSW | 774,149 | 673,063 | 16.3% | 5.4% | 33.9% |
| ECUH | 840,776 | 728,750 | 19.1% | 5.1% | 16.7% |

- The Dryad sites keep a negative culture as one row with no organism. UTSW and ECUH describe
  positivity as "the inclusion of susceptibility results".
- In MGB every non-negative row carries AST, so a culture that grew an organism without AST cannot
  appear as positive. Such cultures are either absent from the release or folded into negatives.
- MGB is the outlier. Blood positivity matches the other sites, but urine and respiratory positivity
  are two to three times higher. That fits a release that drops growth without AST, such as mixed
  urine flora and normal respiratory flora. MGB also has the fewest respiratory cultures (36,821).

**Open question.** If MGB drops growth-without-AST cultures, its E2 denominator is smaller than the
other sites' by construction. Cross-site E2 comparisons then need a matched definition. This should
be resolved with the MGB dataset authors before Analysis 6.

**Positive cultures without an identified organism** number 1,510 at Stanford, 533 at MGB, 14,906 at
UTSW and 620 at ECUH. At UTSW, 14,397 of them have no organism and no susceptibility rows at all,
even though UTSW defines positivity by the presence of susceptibility results. Such cultures get
`organism_unidentified = TRUE` and are non-evaluable under every label. No negative culture carries
an organism name at any site.

## Organisms

Rules live in `organisms.py` and run over each site's distinct organism names (314 at Stanford, 97 at
MGB, 395 at UTSW and 265 at ECUH).

- Names are uppercased. Legacy `ZZ` and `ZZZ` prefixes and lab annotations in parentheses or braces
  are removed, for example `(CARBAPENEM RESISTANT)`, `{MRSA}` and `(MUCOID)`.
- Groups are ecoli, klebsiella, proteus, enterobacter, other_enterobacterales, non_enterobacterales
  and none (negative culture).
- Grouping follows current taxonomy. `ENTEROBACTER AEROGENES` is renamed `KLEBSIELLA AEROGENES` and
  falls in klebsiella. This is a choice. An AmpC-oriented grouping would put it with Enterobacter.
- Other Enterobacterales are matched by genus, including Citrobacter, Serratia, Morganella,
  Providencia, Hafnia, Raoultella, Pantoea, Kluyvera, Cronobacter, Salmonella, Shigella, Yersinia and
  non-coli Escherichia.
- Generic descriptors get a best guess and `organism_ambiguous = TRUE`. "Enteric" and
  "lactose fermenting" gram negative rods count as other_enterobacterales. "Non lactose fermenting"
  gram negative rods, bare "gram negative rods" and "gram negative organism" count as
  non_enterobacterales. Named non-fermenters are non_enterobacterales and not ambiguous.
- Descriptive names such as `MIXED GENITAL FLORA` and UTSW names with a drug suffix such as
  `CANDIDA GLABRATA ANIDULAFUNGIN` fall in non_enterobacterales.
- `idsa_esbl_species` flags E. coli, K. pneumoniae (including subsp. ozaenae), K. oxytoca and
  P. mirabilis.
- A polymicrobial culture takes the first group present in the order above, so E. coli plus
  Pseudomonas is ecoli. Isolate-level groups stay in `isolate_ast`.

## Antibiotics

Canonical codes are the American Society for Microbiology abbreviations that MGB ships in
`AST_code`, including its combination suffixes such as `CAZ_AVI` and `TOL_TAZ`. MGB codes are used
as shipped. Names at the other sites are normalized and mapped, with aliases for release spellings.

- Stanford: `Imipenem/Ebactam` is read as imipenem/relebactam, and `Ticacarcillin/Clavulanic Acid`
  is a misspelling of ticarcillin/clavulanate.
- UTSW: `Pip/tazobactam` (a sentinel agent), `Amox/clavulanate`, `Trimeth/sulf`, `Polymixin e`
  (colistin), `Fluconizole` and others.
- ECUH: `Augmentin`, `Amp/sulbactam`, `Trimeth-sulfamethoxazole` and others. ECUH also reports
  cefotaxime and ceftazidime with clavulanate (`CTX_CLV`, `CAZ_CLV`), which are the CLSI ESBL
  confirmatory tests, and mCIM and eCIM carbapenemase tests.
- Some names stay unmapped on purpose, listed in `vocab.UNMAPPED_BY_DESIGN`. These are gentamicin
  synergy screens, the cefoxitin screen, isoniazid critical-concentration variants, and two UTSW
  names with an undocumented `rm` suffix.

The sentinel agents for the composite label are CRO, CAZ, FEP and TZP. The QC report shows how often
each one is reported for Enterobacterales cultures at each site.

## Susceptibility

| raw | canonical |
|---|---|
| Susceptible | S |
| Intermediate | I |
| Resistant | R |
| Susceptible dose dependent (either spelling) | I |
| Non-susceptible | R |
| Inconclusive, Synergism, anything unrecognized | other |
| NULL | NULL |

- Susceptible dose-dependent sits in the intermediate band for reporting, so it maps to I. Under the
  non-susceptibility rule it counts as a hit.
- Non-susceptible is only reported for drugs without an R breakpoint, so it maps to R.
- Labels count only S, I and R as tested. NT, other and NULL never make a culture evaluable.

**Duplicates and conflicts.** The same culture, organism and agent can carry different results,
mostly from multiple isolates of one species. This happens in 8,882 Stanford groups, and in 19,795
MGB groups under the lab interpretation (13,477 under CLSI 2022). Labels pool results, so any
resistant result makes the culture resistant.

**Breakpoint era.** The plan assumed MGB is uniformly reinterpreted to CLSI 2022. The release ships
both the lab interpretation (`AST_pheno`) and a CLSI 2022 reinterpretation (`CLSI_2022_pheno`), and
both are kept.

- `sir` is the lab interpretation at every site, which is comparable to Stanford's contemporaneous
  results.
- `sir_clsi2022` exists only at MGB and supports a within-site sensitivity analysis.
- The shift is large. Ceftriaxone moves from S to R in 3,467 rows. Piperacillin/tazobactam and
  cefepime move from S to SDD in 32,800 and 17,203 rows.
- MGB E1 composite prevalence is 13.3% under the lab interpretation and 18.2% under CLSI 2022.
  CLSI 2022 also has fewer evaluable cultures, because it is missing where no MIC was available.

**Preliminary results (MGB).** A preliminary result is dropped when a final result exists for the
same culture, organism and agent (77,786 groups). Preliminary-only results are kept with
`preliminary = TRUE` (44,871 groups).

## Enzyme tests (MGB only)

MGB reports ESBL, carbapenemase, beta-lactamase and PBP2a tests as their own rows. These rows also
carry an `AST_code` (`ESBL`, `CARBA`, `BLAC`, `MEC`), so the loader keeps them out of `isolate_ast`.
A few other MGB test rows, such as mCIM and inducible clindamycin resistance, have no enzyme class.
They stay in `isolate_ast`, with their Positive or Negative results mapped to other.
ESBL tests cover 8,181 cultures, almost all E. coli, K. pneumoniae or K. oxytoca. Only 24 cultures
have both a positive and a negative ESBL result, and pooling makes them positive. The
`esbl_enzyme_class` label is therefore evaluable for a small, selected subset of MGB cultures.

## Setting

| site | raw structure |
|---|---|
| Stanford | One row per culture with IP, OP, ER and ICU flags. Emergency visits appear only together with IP. |
| MGB | One flag per linked encounter (IP, OP, ER, UC, day surgery). Flags are OR-ed across a culture's encounters. |
| UTSW, ECUH | Stanford layout. ED-only cultures number 27 at UTSW and fewer than 11 at ECUH. |

- `setting_primary` takes the first setting present in this order: inpatient, ed, urgent_care,
  outpatient (which includes day surgery), then unknown. Raw flags stay in `setting_raw_*`.
- Stanford cultures fall into OP only (275,296), IP with ER (261,060), IP only (152,443), IP with
  ER and ICU (29,164), IP with ICU (21,446) and no flag (11,658), plus two combinations below the
  small-cell limit. An ED visit that ended in discharge is either absent or recorded as outpatient.
- So `ed` is effectively empty at all three Dryad sites, and only MGB separates ED cultures
  (114,704). Setting is unknown for 11,658 Stanford, 229 UTSW and 3,662 ECUH cultures.
- ICU exists only at the Dryad sites, and urgent care and day surgery only at MGB.

## Age

All four sites share the same nine bands: 18-24, 25-34, 35-44, 45-54, 55-64, 65-74, 75-84, 85-89 and
90+. No coarsening is needed and no age value is unmapped.

## Sex

- MGB discloses Female, Male, Unknown and X. X covers 11 cultures and maps to unknown because the
  canonical enum has no other slot.
- Stanford, UTSW and ECUH release a 0 or 1 code without saying which is which. The mapping is
  inferred and recorded in `src/armd_common/mappings/sex_codes.yaml`: code 0 is female and code 1 is
  male at all three sites. `sex` carries the mapped value, `sex_code_masked` keeps the raw code, and
  `sex_inferred` is TRUE for every value that came from the inferred mapping. The data creators have
  not confirmed it.
- `sex_inference.py` tests each masked code against MGB's disclosed sex with two independent lines of
  evidence, and a mapping is recorded only when both agree (run
  `results/sex_inference/20260911T194822Z`, commit 38ca855).

| site | male-only diagnosis, code 0 / code 1 | female-only diagnosis, code 0 / code 1 | pattern correlation with MGB |
|---|---|---|---|
| Stanford | 0.02% / 24.2% | 32.2% / 0.2% | 0.89 |
| UTSW | fewer than 11 patients / 13.5% | 15.2% / 0.08% | 0.92 |
| ECUH | fewer than 11 patients / 14.9% | 20.9% / 0.11% | 0.90 |
| MGB reference, female / male | fewer than 11 patients / 15.3% | 16.3% / fewer than 11 patients | |

- Male-only diagnoses are prostate hyperplasia, prostate, testis and penis cancers and erectile
  dysfunction. Female-only diagnoses are female reproductive cancers, endometriosis, ovarian cysts,
  menstrual disorders and pregnancy. All four releases use the same CCSR names for these.
- The pattern line compares the share of code 0 with MGB's female share across 26 or 27 specimen
  and age strata, under both possible mappings.
- Stanford's code 0 covers 189,864 patients, exactly the female count published in the ARMD-Stanford
  data descriptor. No patient has conflicting codes at any site.

## Area Deprivation Index

- The Dryad sites ship an integer national percentile (1 to 100) and a state decile.
- MGB ships non-integer scores between about 1.3 and 78.7 and non-integer state ranks, which look
  like averages over several block groups.
- Scores are kept as DOUBLE without rescaling, with `adi_missing`.

| site | national score missing | state rank missing |
|---|---|---|
| Stanford | 21.7% | 21.7% |
| MGB | 32.7% | 32.7% |
| UTSW | 4.1% | 57.0% |
| ECUH | 4.9% | 90.0% |

The state rank is not usable at UTSW or ECUH.

## Time

- Stanford order times resolve to the second in `order_time_jittered_utc`. UTSW and ECUH use
  `order_time_jittered`, also to the second. MGB order times are dates.
- Dates are shifted per patient. MGB dates run from 1931 to 2108, and UTSW and ECUH shift by up to
  two months. Calendar time is not usable, so every temporal feature must use the release's
  days-to-culture columns.
- Stanford, UTSW and ECUH exclude repeat cultures within 14 days, and MGB does not. See the next
  section.

## Repeat cultures and the unit of analysis

An earlier version of this document said every release excludes repeat cultures within 14 days. That
is true of three releases and not of MGB. This section records what each release documents, what the
data show, what published rules do, and the options armd-common proposes. Figures come from the QC
runs `results/qc_<site>/20260926T1204*` (`qc.repeat_profile`), over all ordered cultures and without
labels, unless stated otherwise.

### What the releases document

| site | documented rule | source |
|---|---|---|
| Stanford | "repeated cultures from the same patient within a two-week period were excluded". The Dryad page adds: "Patients with prior microbiological cultures within two weeks before the current culture were excluded to avoid overlapping data and ensure distinct clinical events." | Nateghi Haredasht 2025, Methods; Dryad Methods |
| UTSW | "we accounted for patients with an active infection that might have multiple cultures taken in a short time period by excluding patients with prior microbiological cultures within the two weeks before the encounter" | Dryad Methods |
| ECUH | "Patients with possible active infections were identified by having multiple cultures within the previous two weeks from the encounter and were thus excluded from the data set" | Dryad Methods |
| MGB | "Exclusion: Microbiology cultures with no previous culture within the last 14 days." Read literally, this is the reverse of the other three rules. | PhysioNet, Cohort Selection |

- No release says how it handles same-day cultures, different specimen types or blood culture sets.
- UTSW and ECUH anchor the window on the encounter and word it as excluding patients. Stanford
  anchors it on the culture. The data below show the same result at all three.
- MGB says dates "are randomly shifted for each patient, encounter, and microbiology culture". The
  three Dryad releases shift dates per patient.

### What the data show

- **Stanford, UTSW and ECUH.** No patient has a culture 1 to 13 days after another culture. Gaps start
  at exactly 14 days.
- **Same-day cultures at those three sites share one order time.** 321,658 Stanford, 345,365 UTSW and
  359,625 ECUH cultures share an exact order time with another culture of the same patient, and no
  patient has two different order times on one day. The Dryad rule therefore acts on order events and
  keeps every culture in the index order.
- **MGB applies no repeat rule.** 35.0% of MGB cultures have a culture from the same patient 1 to 13
  days earlier, and 35.5% counting 14 days. The gaps peak at one day (155,950 cultures) and fall
  smoothly after that.
- **MGB intervals look preserved.** 66.2% of MGB cultures share an order date with another culture of
  the same patient. Random shifts per culture would erase that clustering and the peak at one day. The
  maintainers have not confirmed this.
- **Blood cultures come in sets.** Among blood cultures, 90% at Stanford, 86% at MGB, 91% at UTSW and
  93% at ECUH share a patient and day with another blood culture.

### What published rules do

The surveillance standards count isolates.
- **CLSI M39.** Duplicates are eliminated "by including only the first isolate of a species, patient,
  and/or analysis period, regardless of specimen source or antimicrobial susceptibility profile"
  (CLSI M39, 5th edition, as summarised in Simner 2022, Table 1).
- **NHSN Repeat Infection Timeframe.** The RIT is "a 14-day timeframe during which no new infections
  of the same type are reported". "Additional pathogens recovered during the RIT from the same type of
  infection are added to the event", and for BSI, UTI and PNEU "the RIT will apply at the major type of
  infection" (NHSN 2026a).
- **NHSN AR Option.** For invasive specimens "there should be 14 days with no positive culture result
  from the laboratory for the patient and specific organism" before another event is entered. The
  manual's worked example chains the window: a culture that was not reported still restarts the 14
  days (NHSN 2026b, Table 1). Same-day duplicates are defined within the invasive or non-invasive
  source type.
- **WHO GLASS.** GLASS reports "only one result for each surveillance period (for example, 12 months)
  ... for each patient per surveyed specimen type and surveyed pathogen". Blood and urine from the
  same patient are both kept. "Repeated negative results for the same specimen type in the same
  patient should also be de-duplicated" (WHO 2023, section 9.5.1).

Blood culture guidance treats sets drawn together as one episode.
- IDSA and ASM define a set as all bottles from one venipuncture and recommend several sets for each
  septic episode (Baron 2013; Miller 2018).
- A blood culture series is "All BC sets collected from the same patient over a period of time",
  usually 24 or 48 hours (Palavecino 2024).

Prediction studies use several units.
- **AMR-UTI, the MGB urinary tract infection dataset (Oberst 2020).** "To mimic the empiric treatment
  setting, we restrict to the first specimen from an infection and exclude any specimens taken within
  a 14 day period from the same body site as duplicates. Specimens taken on the same day are merged or
  kept separate, depending on whether they come from the same or different body sites
  (respectively)."
- **Deelen 2021.** "Multiple sets of blood cultures obtained on the same calendar day were considered
  part of the same infection episode." Its comparison group includes episodes whose blood cultures were
  negative, and episodes are anchored to the start of empiric antibiotics. It is the closest published
  design to the E2 estimand.
- **Yuan 2025.** The first positive blood culture with an Enterobacterales species per 14-day period.
- **Amrollahi 2024, on ARMD-Stanford.** "we excluded culture orders if the patient had a positive
  result within the preceding 14 days".
- **Goodman 2016, Tumbarello 2011, Hebert 2020 and Howard 2024.** One episode or culture per patient.

This evidence has three limits.
- Every standard above counts positive isolates. Only GLASS mentions negative results, and it
  de-duplicates them in the same way. No published rule defines the unit for a cohort of culture orders
  that includes negatives, so the choice below extends the ARMD rule and says so.
- The standards disagree on specimen type. M39 pools specimen sources. The NHSN RIT, the NHSN AR Option
  and GLASS keep infection or specimen types apart.
- Most published models keep one culture per patient or per episode. TRIPOD+AI requires every record
  of one individual to stay in one partition (Collins 2024). The patient-level splits already meet
  that, whatever unit is chosen.

### Options in armd-common

They are the `repeat_rule` and `unit` options of `cohorts.CohortSpec`. On 2026-09-27 the project
adopted `chained_14d` and `specimen_day` for the benchmark tasks and OrderGap. `CohortSpec` itself still
defaults to `release` and `culture`, so a caller reproduces the releases unless it opts in.

**Repeat rule.**

| value | rule | support |
|---|---|---|
| `chained_14d` (proposed default) | An order date is an index if the patient has no culture 1 to 13 days earlier, whether or not that earlier culture was kept. | The Stanford Dryad wording and the NHSN AR Option worked example. A 14-day window is also used by the NHSN RIT, AMR-UTI, Yuan 2025 and Amrollahi 2024. |
| `release` | Keep cultures as released. | MGB then keeps repeat cultures that the other releases drop. |
| `first_per_patient` | Keep each patient's first index order only. | Goodman 2016, Tumbarello 2011, Hebert 2020, Howard 2024. Sensitivity analysis only. |

- **Why a chained window.** It is what the Stanford release describes. It also reproduces all three
  Dryad releases exactly, so it changes only MGB.
- **Why every culture counts.** The rule looks at any earlier culture, not only positive ones. The
  outcome is unknown at order time, and GLASS de-duplicates negatives as well. Only then can E1 and E2
  be drawn from one cohort.
- **Effect at MGB.** 632,625 of 972,798 cultures are kept (65.0%), on 410,763 order dates, and no
  patient drops out. Positivity rises from 21.3% to 26.2% and the blood share falls from 64.5% to
  57.4%. Enterobacterales cultures fall from 135,929 to 120,154.
- **MGB is still the outlier after the rule.** The other sites are 16.0%, 13.1% and 13.3% positive,
  with blood shares of 38.8%, 36.9% and 40.7%. Repeat cultures explain only part of the gap described
  under Culture positivity.

**Unit of analysis.**

| value | one row per | support |
|---|---|---|
| `specimen_day` (proposed default) | patient, specimen type and order date | Sets drawn together belong to one septic episode (Baron 2013; Miller 2018; Deelen 2021). AMR-UTI merges same-day specimens from one body site and keeps different sites apart. The NHSN RIT, the NHSN AR Option and GLASS all keep infection or specimen types apart. |
| `order_event` | patient and order date | One empiric decision covers every culture ordered together. Deelen 2021 anchors episodes on the start of empiric antibiotics. No standard pools specimen types like this, so it is a sensitivity analysis. |
| `culture` | culture, as released | The release grain. Each blood culture set counts separately. |

The order date is the finest grain every site shares. MGB records dates only, and at the Dryad sites
all same-day cultures share one order time. Counts after `chained_14d`:

| site | cultures | specimen-days | order events | order events with two or more specimen types |
|---|---|---|---|---|
| Stanford | 751,075 | 596,787 | 584,403 | 12,072 |
| MGB | 632,625 | 457,123 | 410,763 | 44,744 |
| UTSW | 774,149 | 629,476 | 572,417 | 56,943 |
| ECUH | 840,776 | 678,420 | 645,966 | 30,896 |

**Labels for units with several cultures.**
- A unit is labelled 1 if any of its cultures is labelled 1. This matches how isolates are already
  pooled within a culture.
- A unit is labelled 0 only when every culture in it is evaluable and none is labelled 1. Otherwise
  the unit is non-evaluable.
- Under E1, cultures that are negative or grew no Enterobacterales are left out of the unit instead of
  making it non-evaluable.
- Features are measured from the order date, so they mostly agree within a unit. A unit takes the
  maximum of each flag and count, so it has a history if any of its cultures records one.
  `days_since_prior_org` is 0 when there is no prior organism, so a unit takes its most recent
  non-zero value (`cohorts.unit_feature_sql`).
- Specimen indicators describe the unit's specimen type under `specimen_day`, and any specimen present
  under `order_event`. Ward flags differ within a few hundred order events, where same-day cultures
  came from different encounters, and the unit takes any ward present.
- **Flagged: MGB history differs between cultures ordered on the same date.** Across all MGB cultures,
  276,054 order events hold more than one culture. In 8,246 of them the cultures disagree on
  `days_since_prior_org`, in 737 on `num_prior_orgs` and in 76 on `prior_ESBL`. These counts come from a
  one-off aggregate check on 2026-09-26 that is not yet part of QC. MGB links these tables
  by culture ID, and same-day cultures should share their history at day resolution, so MGB seems to
  compute each culture's history separately. The Dryad sites show no such disagreement.

**Still depends on MGB.** The maintainers are being asked which repeat rule they applied and whether
intervals between a patient's cultures are preserved. If intervals are not preserved, no window rule
can be applied at MGB.

## Feature tables

| table | Stanford | MGB | UTSW | ECUH |
|---|---|---|---|---|
| labs | yes, 14-day window | no | yes, 14-day window | yes, 30-day window |
| vitals | yes | no | no | yes |
| prior antibiotics | prior_med and class and subtype exposure | prior_abx | as Stanford | as Stanford |
| prior resistance | microbial_resistance | prior_micro | as Stanford | as Stanford |
| prior organisms | prior_infecting_organism | prior_org | as Stanford | as Stanford |
| comorbidity | CCSR and Elixhauser components | ICD-10 and Elixhauser category, date-linked | as Stanford | as Stanford |

- The labs window is the median `Period_Day`.
- **Nursing home visits change shape across sites.** Stanford keeps only visits (7,628 rows). MGB,
  UTSW and ECUH keep a row for nearly every culture and leave the value NULL when there was no visit.
  Only 6,509 of 446,141 MGB rows, 3,062 of 761,866 UTSW rows and 35,762 of 806,338 ECUH rows carry a
  value. A feature must test the value, not the presence of a row.
- The gate measures coverage by non-NULL days-to-culture for that reason.
- Stanford's `prior_infecting_organism_days_to_culutre` is misspelled. The loaders resolve column
  names from candidate lists.

## Harmonized features

`features.py` builds 48 features per culture: the 45 features of the AMIA 2026 ARMD-MGB model
(`~/ebsl-amr`), plus urine, blood and respiratory specimen indicators. Every site gets the same names
and definitions. The tables are written to `features/harmonized48_<interpretation>/` inside each
site folder in the data root. The runs are `results/features_<site>/20260911T2007*` and
`20260911T2008*`, on commit 2fefb4e.

### Deliberate departures from the AMIA notebook

- **Links.** The notebook joined comorbidity, prior procedures and nursing home visits on
  `order_proc_id_coded` alone. Those MGB files use a different culture ID mapping, so 14 notebook
  features came from unrelated cultures: six comorbidity flags, the category count, five procedure
  flags, the 30-day procedure flag and nursing home. So did the 4 interactions built on them. Here
  every table links through the loader.
- **Heart failure.** The notebook matched "heart failure|congestive heart". ARMD's "Congestive heart
  failure" category is a catch-all, as shown below. Heart failure now uses the CCSR "Heart failure"
  category, and the catch-all is left out of every comorbidity feature, including the category count.
- **Timing.** Prior antibiotics, resistance, organisms and procedures count only events at least one
  day before the order. The notebook kept MGB antibiotic last doses on or after the culture day,
  which is 9.9% of distinct MGB exposure rows. UTSW and ECUH put the index culture's own resistance
  results at day 0. Nursing home visits keep day 0, which marks residence.
- **Antibiotic classes.** Each class is a fixed set of drug codes taken from MGB's class membership,
  then applied to every site's medication names. Brand names, salts and formulations are normalized,
  and topical products are dropped.

| feature | drug codes |
|---|---|
| T_fq_90d | CIP, DLX, GAT, LVX, MXF, NOR, OFX |
| T_ceph3_90d | CAZ, CPT, CRO, CTX, FDC, FEP |
| T_carb_90d | DOR, ETP, IPM, MEM |
| T_glyco_90d | VAN |
| T_sulfa_90d | SXT, TMP |
| T_esp_90d | AMC, SAM, TZP |
| T_amino_90d | AMK, GEN, KAN, PLZ, STR, TOB |

- **Prior resistant organisms.** Every site uses one phenotype definition based on prior
  susceptibility results:
  - **ESBL:** an Enterobacterales isolate resistant to ceftriaxone, ceftazidime or cefepime.
  - **CRE:** an Enterobacterales isolate resistant to a carbapenem.
  - **MRSA:** *S. aureus* resistant to an anti-staphylococcal beta-lactam, or an organism name that
    says MRSA.
  - **VRE:** *E. faecalis* or *E. faecium* resistant to vancomycin.
- **Checks against MGB's own flags.**
  - MRSA agrees almost exactly: 49,682 cultures are flagged by both, with only 299 and 194
    disagreements.
  - The ESBL definition recovers 71% of MGB's flag, so MGB's own flag is broader.
  - An all-enterococci VRE rule flagged 16,631 cultures that MGB did not. It was then restricted to
    *E. faecalis* and *E. faecium*. Agreement after the restriction was not re-measured.
- **Prior organisms.** Genus counts use the 14 genera recorded at all four sites, because the Dryad
  tables record genus only.
- **Ward.** Flags are combined across a culture's encounters. The notebook kept the first row.

### The heart failure catch-all

At MGB, the only release with ICD-10 codes, "Congestive heart failure" spans 390 ICD-10 prefixes. The
next broadest category spans 59. Its most common codes are hyperlipidemia (E78.5), anxiety (F41.9),
depression (F32.9), hyponatremia, anemia and type 2 diabetes. The CCSR "Heart failure" category holds
I50, I11.0 and I13.0. The same pattern shows at every site:

| site | "Congestive heart failure", ages 18-24 / 45-54 / 75-84 | "Heart failure", same ages | heart failure feature, notebook rule / fixed |
|---|---|---|---|
| Stanford | 61.8% / 82.8% / 91.2% | 3.4% / 11.0% / 25.6% | 81.5% / 15.6% |
| MGB | 27.9% / 38.5% / 48.1% | 0.9% / 4.7% / 15.1% | 41.6% / 9.0% |
| UTSW | 40.5% / 61.7% / 73.3% | 2.6% / 8.7% / 15.1% | 61.4% / 9.6% |
| ECUH | 40.5% / 72.1% / 85.3% | 0.3% / 10.4% / 23.2% | 72.6% / 14.1% |

The other categories the features use hold the right codes: liver disease, lymphoma, metastatic
cancer, obesity and renal failure.

### What each feature measures depends on the release

| antibiotic class, 90 days | Stanford | MGB | UTSW | ECUH |
|---|---|---|---|---|
| T_ceph3_90d | 1.9% | 29.4% | 1.3% | 10.9% |
| T_esp_90d | 7.9% | 13.6% | 5.1% | 0.01% |
| T_sulfa_90d | 6.7% | 10.2% | 0.04% | 0.01% |
| T_fq_90d | 10.7% | 15.7% | 15.1% | 0.9% |
| T_glyco_90d | 4.0% | 22.1% | 9.7% | 4.2% |
| T_carb_90d | 0.6% | 5.7% | 0.3% | 0.3% |
| T_amino_90d | 1.1% | 1.2% | 1.8% | 0.2% |

| drug, any prior exposure | Stanford | MGB | UTSW | ECUH |
|---|---|---|---|---|
| ceftriaxone | 6.9% | 28.8% | 3.9% | 37.8% |
| cefepime | 1.9% | 27.5% | 2.9% | 7.8% |
| piperacillin/tazobactam | 15.1% | 17.2% | 17.6% | not recorded |
| vancomycin | 10.8% | 35.2% | 25.6% | 17.2% |
| meropenem | 0.9% | 8.1% | 0.9% | 0.01% |
| TMP-SMX | 18.8% | 23.4% | not recorded | 0.3% |
| amoxicillin/clavulanate | 15.1% | 16.7% | not recorded | 0.2% |
| ciprofloxacin | 17.5% | 23.6% | 28.4% | 5.5% |

| other features | Stanford | MGB | UTSW | ECUH |
|---|---|---|---|---|
| proc_cvc | 18.5% | 1.6% | 5.3% | 0.2% |
| proc_surgical | 20.9% | 9.5% | 1.9% | 0.2% |
| elixhauser_count, median | 26 | 2 | 7 | 11 |
| prior_ESBL | 3.3% | 7.3% | 3.3% | 3.9% |
| prior_MRSA | 1.5% | 5.1% | 2.7% | 2.5% |
| nursing_home_90d | 1.0% | 1.3% | 0.4% | 4.3% |

- **Each release records a different drug list.**
  - Stanford and UTSW barely record IV cephalosporins and carbapenems.
  - UTSW has no TMP-SMX or amoxicillin/clavulanate.
  - ECUH has no piperacillin/tazobactam and very few oral drugs.
  - MGB records inpatient and outpatient drugs alike.
  - The broader Dryad `prior_med` table holds the same antibiotic rows, so no better source exists.
- **Procedure and comorbidity capture differ by an order of magnitude,** in line with how much history
  each release carries.
- **Definitions match across sites for:** age, sex (inferred at three sites), specimen, the inpatient
  and outpatient flags, the four prior resistance flags, prior organism counts and the individual
  comorbidity flags. Their capture still varies with each release's history depth.
- **Release-limited features:**
  - the seven antibiotic class features and the interactions built on them;
  - the procedure features and `elixhauser_count`;
  - `hosp_ward_ER`, because the Dryad sites flag ED only together with inpatient;
  - `hosp_ward_UC`, which only MGB records;
  - the ADI score scale, which differs at MGB, and `adi_high`, because state rank is mostly missing at
    UTSW and ECUH.
- **MGB interpretation.** Under CLSI 2022, MGB prior MRSA falls from 5.1% to 3.0%, because many
  oxacillin results have no CLSI 2022 interpretation.

## Recent antibiotic exposure

A sensitivity analysis excludes cultures with recent antibacterial exposure
(`CohortSpec(exclude_recent_antibiotics=True)`). Every release
records exposure times in whole days, and day 0 holds different things at each site. Shares are of
ordered cultures in each release's cohort table, measured on 2026-09-11 without labels:

| cultures with an antibacterial record | Stanford | MGB | UTSW | ECUH |
|---|---|---|---|---|
| last dose after the order | none recorded | 23.6% | none recorded | none recorded |
| day 0 | 41.2% | 36.8% | none | none |
| days 1 to 3 | 7.1% | 20.8% | 5.0% | 2.2% |
| day 0 to 3, or last dose up to day 3 at MGB | 45.1% | 61.0% | 5.0% | 2.2% |

- Stanford's `prior_med` holds the day 0 records. Its `antibiotic_subtype_exposure` table starts at
  day 1 and otherwise matches. At UTSW and ECUH the two tables are identical and start at day 1.
- MGB's `prior_abx` records the last dose only, so a course cannot be dated by its first dose.
- A day 0 record cannot be placed before the order, and it includes treatment started after the
  culture was drawn. Counting it would remove 45% and 61% of cultures at Stanford and MGB, and almost
  nothing at UTSW and ECUH.
- The exclusion therefore counts days 1 to 3 at every site, from the table the antibiotic features use. It misses
  doses on the order day before the culture, and MGB courses that began earlier and continued past
  the order. This was adopted on 2026-09-27.

## Implied susceptibility

`implied.py` applies Stanford's rules to any site's `isolate_ast`, and it runs only when a caller asks
for it. The rules file has 50 rows. One is an exact duplicate, so 49 rules remain: 23 intrinsic
resistance, 7 equivalence, 16 "susceptible if susceptible to", 2 "otherwise resistant" and 1
combined rule.

- Reported S, I and R results are never overwritten.
- Rules read reported results only, so imputations never chain.
- Conflicting reports for a source agent use the worst result.
- Every imputed row has `imputed = TRUE` and the rule text.
- `strict=True` applies intrinsic resistance and equivalences only. `strict=False` applies every rule.
- Some Stanford rules deserve review before use. For example, one rule makes Acinetobacter
  ertapenem susceptible when ceftriaxone is susceptible, although Acinetobacter is intrinsically
  resistant to ertapenem.

Stanford also ships a per-row implied table. It uses raw lab strings such as `No Interpretation`,
`See Comment` and `Susceptible - Dose Dependent` instead of the cohort's standardized values, so it
is not used.

## Phase 0 gate

Run `results/00_gate/20260911T182159Z` on commit d96f9fc. All four sites pass, so the verdict is go.

- Demographics, ward and ADI cover every negative culture at every site.
- Every event table records dated events for negative cultures, and none falls below the 10%
  relative-coverage threshold.
- For tables linked by culture ID that carry an order time, the order date matches the cohort for
  every negative culture.
- Every negative culture has an order time.

The estimands below cover whole cohorts, including patients who will land in the test partition.
They describe the cohort and did not inform any modelling decision.

| site | E1 evaluable | resistant | E1 prevalence | E1 not tested | E2 evaluable | E2 prevalence |
|---|---|---|---|---|---|---|
| Stanford | 78,306 | 8,713 | 11.1% | 5,300 | 744,265 | 1.17% |
| MGB, lab interpretation | 135,746 | 18,001 | 13.3% | 183 | 971,521 | 1.85% |
| MGB, CLSI 2022 | 108,838 | 19,791 | 18.2% | 27,091 | 944,613 | 2.10% |
| UTSW | 50,327 | 3,973 | 7.9% | 5,756 | 753,109 | 0.53% |
| ECUH | 86,680 | 8,066 | 9.3% | 19 | 840,085 | 0.96% |

- E1 is the composite label and E2 is the decision-time label. Resistant cultures are the same
  under both, because E2 only adds negative and non-Enterobacterales cultures as Y = 0.
- Prevalence drops roughly tenfold from E1 to E2 at every site.
- Under the IDSA definition, E1 prevalence is 11.5% at Stanford, 13.0% at MGB, 4.9% at UTSW and 7.5%
  at ECUH.
- The MGB enzyme label is evaluable for 8,179 cultures, with prevalence 6.8%.

**Sentinel agent reporting (layer T).** The table gives the share of Enterobacterales cultures with
an S, I or R result for each agent.

| site | Enterobacterales cultures | CRO | CAZ | FEP | TZP | any sentinel |
|---|---|---|---|---|---|---|
| Stanford | 83,606 | 80.0% | 67.4% | 31.7% | 89.2% | 93.7% |
| MGB, lab interpretation | 135,929 | 86.3% | 50.8% | 72.0% | 95.9% | 99.9% |
| MGB, CLSI 2022 | 135,929 | 79.9% | 62.5% | 70.8% | 77.6% | 80.1% |
| UTSW | 56,083 | 81.5% | 11.7% | 58.7% | 87.3% | 89.7% |
| ECUH | 86,699 | 98.1% | 52.3% | 56.6% | 99.7% | 100.0% |

"Resistant to at least one of four agents" rests on a different panel at each site. Cefepime is
rarely reported at Stanford and ceftazidime rarely at UTSW. At MGB, ceftazidime is listed without a
lab interpretation for 12.3% of Enterobacterales cultures. CLSI 2022 interprets many of those MICs
but leaves piperacillin/tazobactam uninterpreted for 18.5%. This is the cascade-reporting layer that
Analysis 4 quantifies.

## Splits

The splits were frozen on 2026-09-11 with seed 20260911 (`python -m armd_common.splits`).

- **Strata.** Patients fall into three strata: no evaluable E1 culture, evaluable but never
  resistant, and ever resistant.
- **Folds.** One of five outer folds is the 20% test partition, and derivation patients get five
  inner folds.
- **Balance.** The ever-resistant share is identical in the test and derivation partitions at every
  site.

| site | patients | strata 0 / 1 / 2 | test patients | ever resistant |
|---|---|---|---|---|
| Stanford | 283,715 | 235,358 / 42,287 / 6,070 | 56,743 | 2.14% |
| MGB | 226,659 | 153,730 / 63,366 / 9,563 | 45,332 | 4.22% |
| UTSW | 237,258 | 206,305 / 27,937 / 3,016 | 47,452 | 1.27% |
| ECUH | 261,217 | 210,603 / 45,111 / 5,503 | 52,244 | 2.11% |

- Manifests live in `benchmarks/manifests/`, and `benchmarks/manifests/SHA256SUMS` pins all four.
- The Stanford, UTSW and ECUH manifests are committed. MGB `anon_id` values are small integers, so
  its patient hashes can be reversed. Its manifest stays local until PhysioNet's terms allow release.
- Regenerating a split gives identical bytes, and `splits.py` refuses to overwrite a frozen manifest.

## Decisions

The project lead decided the following on 2026-09-11:

- **MGB interpretation.** The lab interpretation is the primary definition at MGB, and CLSI 2022 is
  a sensitivity analysis. The AMIA model used CLSI 2022, which raises MGB's E1 prevalence from 13.3%
  to 18.2%.
- **Features in cross-site tasks.** Every task uses all 48 features, cross-site pairs included, and
  per-site capture is reported alongside the results.
- **Masked sex codes.** They are mapped by inference, as described under Sex.
- **Splits.** They are frozen with the stratification described under Splits.

Decided on 2026-09-27:

- **Recent antibiotic window.** Days 1 to 3 at every site, as described under Recent antibiotic
  exposure.
- **Repeat rule.** `chained_14d`, as described under Repeat cultures and the unit of analysis.
  `first_per_patient` and the release as is are sensitivity analyses.
- **Unit of analysis.** `specimen_day`, with `order_event` as a sensitivity analysis.
- **Implied susceptibility.** The sensitivity analysis applies every rule. The rules that can fill a
  sentinel agent for Enterobacterales only impute susceptibility from ceftriaxone, cefotaxime or
  cefazolin, so the flagged Acinetobacter rule does not affect these labels.

Still open:

1. **K. aerogenes grouping.** Klebsiella by taxonomy (current) or Enterobacter by AmpC behaviour.
2. **ESBL confirmatory label at ECUH.** ECUH reports clavulanate combination tests, which could
   support a confirmatory-test label alongside MGB's enzyme label.

## Questions for the ARMD-MGB maintainers

1. Encounter and culture IDs in `comorbidity`, `nursing_home_visits` and `prior_procedures` come
   from a different mapping than the cohort. Is this intended?
2. Were cultures that grew an organism without AST excluded? The answer changes MGB's E2
   denominator.
3. Which repeat-culture exclusion, if any, was applied? The project page says "Exclusion:
   Microbiology cultures with no previous culture within the last 14 days", but 35.5% of cultures have
   a culture from the same patient 1 to 14 days earlier.
4. Are intervals between one patient's cultures preserved, given that dates are shifted "for each
   patient, encounter, and microbiology culture"?

## References

- Amrollahi F, Nateghi Haredasht F, Vansomphone A, et al. Machine learning-based prediction of
  antimicrobial susceptibility: a step towards precision antimicrobial stewardship. AMIA Annu Symp
  Proc. 2024;2024:138-146. PMID 41726396.
- Baron EJ, Miller JM, Weinstein MP, et al. A guide to utilization of the microbiology laboratory for
  diagnosis of infectious diseases: 2013 recommendations by IDSA and ASM. Clin Infect Dis.
  2013;57(4):e22-e121. doi:10.1093/cid/cit278
- CDC NHSN (2026a). Patient Safety Component Manual, chapter 2: Identifying healthcare-associated
  infections for NHSN surveillance. January 2026.
  https://www.cdc.gov/nhsn/pdfs/pscmanual/2psc_identifyinghais_nhsncurrent.pdf
- CDC NHSN (2026b). Antimicrobial Use and Resistance (AUR) Module. July 2026.
  https://www.cdc.gov/nhsn/pdfs/pscmanual/11pscaurcurrent.pdf
- CLSI. M39: Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data. 5th ed.
  Clinical and Laboratory Standards Institute; 2022. Not read directly. Quoted through Simner 2022.
- Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical
  prediction models that use regression or machine learning methods. BMJ. 2024;385:e078378.
  doi:10.1136/bmj-2023-078378
- Cooper LN, et al. Antimicrobial Resistance Microbiological Dataset (ARMD-ECUH). Dryad, version 4.
  doi:10.5061/dryad.7sqv9s55x
- Cooper LN, et al. Antimicrobial Resistance Microbiological Dataset (ARMD-UTSW). Dryad, version 5.
  doi:10.5061/dryad.0rxwdbsd5
- Deelen JWT, et al. An international prospective cohort study to validate 2 prediction rules for
  infections caused by third-generation cephalosporin-resistant Enterobacterales. Clin Infect Dis.
  2021;73(11):e4475-e4483. doi:10.1093/cid/ciaa950
- Goodman KE, et al. A clinical decision tree to predict whether a bacteremic patient is infected with
  an extended-spectrum beta-lactamase-producing organism. Clin Infect Dis. 2016;63(7):896-903.
  doi:10.1093/cid/ciw425
- Hebert C, et al. Prediction of antibiotic susceptibility for urinary tract infection in a hospital
  setting. Antimicrob Agents Chemother. 2020;64(7):e02236-19. doi:10.1128/aac.02236-19
- Howard A, et al. Personalised antimicrobial susceptibility testing with clinical prediction modelling
  informs appropriate antibiotic use. Nat Commun. 2024;15:9924. doi:10.1038/s41467-024-54192-3
- Miller JM, Binnicker MJ, Campbell S, et al. A guide to utilization of the microbiology laboratory for
  diagnosis of infectious diseases: 2018 update by IDSA and ASM. Clin Infect Dis. 2018;67(6):e1-e94.
  doi:10.1093/cid/ciy381
- Nateghi Haredasht F, et al. Antibiotic Resistance Microbiology Dataset (ARMD): a resource for
  antimicrobial resistance from EHRs. Sci Data. 2025;12:1299. doi:10.1038/s41597-025-05649-7. Dryad,
  version 20. doi:10.5061/dryad.jq2bvq8kp
- Oberst M, Boominathan S, Zhou H, Kanjilal S, Sontag D. AMR-UTI: antimicrobial resistance in urinary
  tract infections (version 1.0.0). PhysioNet. 2020. doi:10.13026/se6w-f455
- Palavecino EL, Campodonico VL, She RC. Laboratory approaches to determining blood culture
  contamination rates: an ASM Laboratory Practices Subcommittee report. J Clin Microbiol.
  2024;62(2):e01028-23. doi:10.1128/jcm.01028-23
- Simner PJ, Hindler JA, et al. What's new in antibiograms? Updating CLSI M39 guidance with current
  trends. J Clin Microbiol. 2022;60(10):e02210-21. doi:10.1128/jcm.02210-21
- Tumbarello M, et al. Identifying patients harboring extended-spectrum-beta-lactamase-producing
  Enterobacteriaceae on hospital admission: derivation and validation of a scoring system. Antimicrob
  Agents Chemother. 2011;55(7):3485-3490. doi:10.1128/aac.00009-11
- Wei Z, Kanjilal S. Antibiotic Resistance Microbiology Dataset Mass General Brigham (ARMD-MGB)
  (version 1.0.0). PhysioNet. 2025. doi:10.13026/2r5k-b955
- WHO. GLASS manual for antimicrobial resistance surveillance in common bacteria causing human
  infection. Geneva: World Health Organization; 2023. ISBN 978-92-4-007660-0
- Yuan K, et al. Machine learning and clinician predictions of antibiotic resistance in
  Enterobacterales bloodstream infections. J Infect. 2025;90(2):106388.
  doi:10.1016/j.jinf.2024.106388
