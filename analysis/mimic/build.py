"""MIMIC-IV 3.1 for OrderGap's external validation (SAP amendment E-9). Not part of armd-common.

This module turns MIMIC-IV hosp tables into the canonical culture and susceptibility tables that
armd-common's label constructors expect, builds the 48 harmonized features with the definitions and
windows of SAP section 4, and adds the columns E-9 needs: clinicians' empiric choices from eMAR, the
exact Tumbarello items, and insurance and race for fairness subgroups. It reuses armd-common's
labels, repeat rule, unit logic, organism classifier and drug vocabularies, and changes none of them.

Every function takes a DuckDB connection on which `create_views` has registered the tables as
mimic_<table> views, so tests can register small hand-built tables under the same names.

The test partition is sealed: only analysis/09_mimic_evaluate.py may register it.
"""
from __future__ import annotations

import re
import sys
import uuid
from dataclasses import asdict, dataclass
from pathlib import Path

import duckdb
import yaml

from armd_common import labels, organisms, splits, vocab
from armd_common.cohorts import CEFTRIAXONE_AGENTS, SETTING_PRECEDENCE, unit_feature_sql
from armd_common.episodes import kept_cultures, unit_key_sql, unit_labels
from armd_common.features import (ANTIBIOTIC_CLASSES, CARBAPENEMS, COMMON_PRIOR_GENERA, DOUBLE_FEATURES, ESBL_AGENTS,
                                  FEATURES_48, MRSA_AGENTS, MRSA_NAME_PATTERN, NURSING_HOME_WINDOW_DAYS,
                                  PRIOR_MIN_DAYS, PROCEDURE_WINDOW_DAYS, VRE_SPECIES, ANTIBIOTIC_WINDOW_DAYS)
from armd_common.provenance import sha256_file

SITE = "mimic"
PARQUET = Path("/data0/mimiciv/3.1/parquet")
WORK = Path("/data0/mimiciv/3.1/ordergap")
TMP = Path("/data0/mimiciv/3.1/duckdb_tmp")
SPLIT_PATH = WORK / "split_manifest.csv.gz"
SEALED_READER = "09_mimic_evaluate.py"
CODE_LISTS = Path(__file__).with_name("code_lists.yaml")
TABLES = ("patients", "admissions", "transfers", "microbiologyevents", "diagnoses_icd", "procedures_icd",
          "prescriptions", "emar")
REPEAT_RULE, UNIT = "chained_14d", "specimen_day"

# Specimen map, fixed from the specimen and test vocabulary on 2026-10-29 before any label was built.
# Only culture tests count; Gram stains, screens, serology, viral, fungal and mycobacterial tests do not.
CULTURE_TESTS = {
    "urine": (("URINE", "URINE,SUPRAPUBIC ASPIRATE", "URINE,KIDNEY"), ("URINE CULTURE", "REFLEX URINE CULTURE")),
    "blood": (("BLOOD CULTURE",), ("Blood Culture, Routine",)),
    "respiratory": (("SPUTUM", "BRONCHOALVEOLAR LAVAGE", "Mini-BAL", "BRONCHIAL WASHINGS", "TRACHEAL ASPIRATE",
                     "BRONCHIAL BRUSH", "BRONCHIAL BRUSH - PROTECTED", "ASPIRATE"), ("RESPIRATORY CULTURE",)),
}
# A cancelled culture has no result. It leaves the cohort instead of counting as negative.
NOT_AN_ORGANISM = ("CANCELLED",)
# MIMIC reports single letters. D (susceptible dose-dependent) is never resistant, so it counts as I.
SIR = {"S": "S", "I": "I", "R": "R", "D": "I"}
# MIMIC antibiotic names that armd-common's vocabulary does not map.
ANTIBIOTIC_OVERRIDES = {"PIPERACILLIN/TAZO": "TZP", "TRIMETHOPRIM/SULFA": "SXT", "PENICILLIN G": "PEN",
                        "DELAFLOXICIN": "DLX"}
MEDICATION_OVERRIDES = ((re.compile(r"^sulfameth/trimethoprim", re.I), "SXT"),
                        (re.compile(r"^piperacillin-tazobactam", re.I), "TZP"),
                        (re.compile(r"^cefpodoxime", re.I), "CPD"))
# Prescription routes that deliver an antibacterial systemically.
SYSTEMIC_ROUTES = ("IV", "PO/NG", "PO", "IV DRIP", "IM", "ORAL", "NG", "IV BOLUS", "PB", "IVPCA", "G TUBE", "J TUBE",
                   "OG", "NG/OG", "PO OR NG")
# eMAR events that mean a dose was given.
ADMINISTERED = ("Administered", "Started", "Restarted", "Delayed Administered", "Administered Bolus from IV Drip",
                "Administered in Other Location", "Delayed Started", "Partial Administered")
CARBAPENEM_CODES = CARBAPENEMS
BROAD_CODES = CARBAPENEMS + ("TZP", "FEP")
NURSING_HOME_ADMISSION = ("TRANSFER FROM SKILLED NURSING FACILITY",)
NURSING_HOME_DISCHARGE = ("SKILLED NURSING FACILITY",)
TRANSFER_ADMISSION = ("TRANSFER FROM HOSPITAL", "TRANSFER FROM SKILLED NURSING FACILITY")
EMAR_YEAR_GROUPS = ("2014 - 2016", "2017 - 2019", "2020 - 2022")
CLINICIAN_WINDOW_HOURS = 24
# Genus stems, so 'STAPH AUREUS' and 'BETA STREPTOCOCCUS GROUP B' count under their genus.
GENUS_STEMS = {g: g[:max(5, len(g) - 4)] for g in COMMON_PRIOR_GENERA}
GENUS_STEMS.update({"STAPHYLOCOCCUS": "STAPH", "STREPTOCOCCUS": "STREPTOCOC", "ENTEROCOCCUS": "ENTEROCOC"})


# Connection and views -------------------------------------------------------------------------------

def connect(memory: str = "64GB", threads: int = 10) -> duckdb.DuckDBPyConnection:
    TMP.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(config={"temp_directory": str(TMP), "memory_limit": memory, "threads": threads,
                                  "preserve_insertion_order": False})


def create_views(con, parquet_dir: Path = PARQUET) -> None:
    for table in TABLES:
        con.execute(f"CREATE OR REPLACE VIEW mimic_{table} AS SELECT * FROM read_parquet('{parquet_dir / (table + '.parquet')}')")


def input_digests(parquet_dir: Path = PARQUET) -> dict:
    raw = Path("/data0/mimiciv/3.1/hosp")
    return {t: sha256_file(raw / f"{t}.csv.gz", cache_path=WORK / "digests" / f"{t}.sha256") for t in TABLES}


def _quoted(values) -> str:
    return ", ".join("'" + str(v).replace("'", "''") + "'" for v in values)


def antibiotic_code(name: str | None) -> str | None:
    if name is None:
        return None
    return ANTIBIOTIC_OVERRIDES.get(name.strip().upper()) or vocab.antibiotic_code(name)


def medication_code(name: str | None) -> str | None:
    if name is None:
        return None
    for pattern, code in MEDICATION_OVERRIDES:
        if pattern.search(name.strip()):
            return code
    return vocab.medication_code(name)


def _lookup(con, name: str, values, mapper) -> str:
    rows = sorted({(v, mapper(v)) for v in values if v is not None})
    con.execute(f"CREATE OR REPLACE TEMP TABLE {name} (raw VARCHAR, value VARCHAR)")
    if rows:
        con.executemany(f"INSERT INTO {name} VALUES (?, ?)", rows)
    return name


def _distinct(con, table: str, column: str) -> list:
    return [r[0] for r in con.execute(f"SELECT DISTINCT {column} FROM {table}").fetchall()]


def code_lists(path: Path = CODE_LISTS) -> dict:
    return yaml.safe_load(path.read_text())


# Cultures and susceptibility ---------------------------------------------------------------------------

def culture_rows_sql() -> str:
    parts = []
    for specimen, (specs, tests) in CULTURE_TESTS.items():
        parts.append(f"(spec_type_desc IN ({_quoted(specs)}) AND test_name IN ({_quoted(tests)}) AND '{specimen}' = s.specimen)")
    specimens = " UNION ALL ".join(f"SELECT '{s}' AS specimen" for s in CULTURE_TESTS)
    return f"""
        SELECT m.*, s.specimen AS specimen_type
        FROM mimic_microbiologyevents m
        JOIN ({specimens}) s ON {' OR '.join(parts)}
    """


def build_cultures(con) -> dict:
    """Create mimic_culture_episodes and mimic_isolate_ast. One culture is one specimen's culture test.

    culture_episodes: culture_id, patient_id, subject_id, micro_hadm_id, order_time, specimen_type,
                      culture_positive, organism_unidentified. Cancelled cultures are left out.
    isolate_ast: culture_id, is_enterobacterales, idsa_esbl_species, antibiotic_code, sir. Isolates
                 without susceptibility rows appear once with a NULL antibiotic_code.
    """
    con.execute(f"CREATE OR REPLACE TEMP TABLE _mimic_culture_rows AS {culture_rows_sql()}")
    org = organisms.register_lookup(con, _distinct(con, "_mimic_culture_rows", "org_name"), table="mimic_organism_lookup")
    abx = _lookup(con, "mimic_ab_lookup", _distinct(con, "_mimic_culture_rows", "ab_name"), antibiotic_code)
    not_org = _quoted(NOT_AN_ORGANISM)
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE mimic_culture_all AS
        SELECT CAST(micro_specimen_id AS VARCHAR) AS culture_id,
               'mimic:' || subject_id AS patient_id,
               subject_id,
               MAX(NULLIF(hadm_id, '')) AS micro_hadm_id,
               MIN(COALESCE(TRY_CAST(NULLIF(charttime, '') AS TIMESTAMP), TRY_CAST(NULLIF(chartdate, '') AS TIMESTAMP))) AS order_time,
               MIN(specimen_type) AS specimen_type,
               COALESCE(BOOL_OR(org_name IS NOT NULL AND org_name <> '' AND org_name NOT IN ({not_org})), FALSE) AS culture_positive,
               COALESCE(BOOL_OR(org_name IN ({not_org})), FALSE) AS cancelled,
               FALSE AS organism_unidentified
        FROM _mimic_culture_rows
        GROUP BY micro_specimen_id, subject_id
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE mimic_culture_episodes AS
        SELECT * EXCLUDE (cancelled) FROM mimic_culture_all
        WHERE order_time IS NOT NULL AND NOT (cancelled AND NOT culture_positive)
    """)
    sir_case = "CASE upper(trim(r.interpretation)) " + " ".join(f"WHEN '{k}' THEN '{v}'" for k, v in SIR.items()) + \
               " WHEN '' THEN 'NT' ELSE (CASE WHEN r.interpretation IS NULL THEN 'NT' ELSE 'other' END) END"
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE mimic_isolate_ast AS
        SELECT DISTINCT CAST(r.micro_specimen_id AS VARCHAR) AS culture_id,
               o.is_enterobacterales, o.idsa_esbl_species,
               a.value AS antibiotic_code,
               CASE WHEN a.value IS NULL THEN NULL ELSE {sir_case} END AS sir
        FROM _mimic_culture_rows r
        JOIN mimic_culture_episodes e ON e.culture_id = CAST(r.micro_specimen_id AS VARCHAR)
        JOIN {org} o ON o.organism_raw = r.org_name
        LEFT JOIN {abx} a ON a.raw = r.ab_name
        WHERE r.org_name IS NOT NULL AND r.org_name NOT IN ({not_org})
    """)
    counts = dict(con.execute("""
        SELECT 'cultures', COUNT(*) FROM mimic_culture_all
        UNION ALL SELECT 'cancelled', COUNT(*) FILTER (WHERE cancelled AND NOT culture_positive) FROM mimic_culture_all
        UNION ALL SELECT 'no_time', COUNT(*) FILTER (WHERE order_time IS NULL) FROM mimic_culture_all
        UNION ALL SELECT 'kept', COUNT(*) FROM mimic_culture_episodes
    """).fetchall())
    return {"culture_episodes": "mimic_culture_episodes", "isolate_ast": "mimic_isolate_ast", "counts": counts}


# Encounter context ---------------------------------------------------------------------------------------

def build_context(con, episodes: str = "mimic_culture_episodes") -> str:
    """Per culture: setting flags, matched admission, age band, sex, insurance and race groups, year group,
    and whether eMAR covers the encounter. Writes mimic_context."""
    age_case = " ".join(f"WHEN age BETWEEN {lo} AND {hi} THEN '{band}'" for band, (lo, hi) in {
        "18-24": (18, 24), "25-34": (25, 34), "35-44": (35, 44), "45-54": (45, 54), "55-64": (55, 64),
        "65-74": (65, 74), "75-84": (75, 84), "85-89": (85, 89)}.items())
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _mimic_adm AS
        SELECT hadm_id, subject_id, TRY_CAST(admittime AS TIMESTAMP) AS admittime, TRY_CAST(dischtime AS TIMESTAMP) AS dischtime,
               TRY_CAST(NULLIF(edregtime, '') AS TIMESTAMP) AS edregtime, TRY_CAST(NULLIF(edouttime, '') AS TIMESTAMP) AS edouttime,
               admission_location, discharge_location, insurance, race
        FROM mimic_admissions
    """)
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE mimic_context AS
        WITH enc AS (
            SELECT e.culture_id, arg_min(a.hadm_id, a.admittime) AS hadm_id
            FROM {episodes} e JOIN _mimic_adm a ON a.subject_id = e.subject_id
             AND e.order_time BETWEEN COALESCE(LEAST(a.edregtime, a.admittime), a.admittime) AND a.dischtime
            GROUP BY e.culture_id
        ),
        ed AS (
            SELECT DISTINCT e.culture_id
            FROM {episodes} e JOIN mimic_transfers t ON t.subject_id = e.subject_id AND t.eventtype = 'ED'
             AND e.order_time BETWEEN TRY_CAST(t.intime AS TIMESTAMP) AND TRY_CAST(t.outtime AS TIMESTAMP)
        ),
        last_adm AS (
            SELECT e.culture_id, arg_max(a.insurance, a.admittime) AS insurance, arg_max(a.race, a.admittime) AS race
            FROM {episodes} e JOIN _mimic_adm a ON a.subject_id = e.subject_id AND a.admittime <= e.order_time + INTERVAL 1 DAY
            GROUP BY e.culture_id
        ),
        emar_adm AS (SELECT DISTINCT hadm_id FROM mimic_emar WHERE hadm_id IS NOT NULL AND hadm_id <> '')
        SELECT e.culture_id,
               COALESCE(e.micro_hadm_id, n.hadm_id) AS hadm_id,
               CAST(COALESCE(e.micro_hadm_id, n.hadm_id) IS NOT NULL AS INTEGER) AS hosp_ward_IP,
               CAST(d.culture_id IS NOT NULL AS INTEGER) AS hosp_ward_ER,
               CAST(COALESCE(e.micro_hadm_id, n.hadm_id) IS NULL AND d.culture_id IS NULL AS INTEGER) AS hosp_ward_OP,
               CASE WHEN COALESCE(e.micro_hadm_id, n.hadm_id) IS NOT NULL THEN 'inpatient'
                    WHEN d.culture_id IS NOT NULL THEN 'ed' ELSE 'outpatient' END AS setting_primary,
               p.age, CASE WHEN p.age >= 90 THEN '90+' {age_case} ELSE NULL END AS age_band, p.age_exact,
               CASE p.gender WHEN 'M' THEN 'm' WHEN 'F' THEN 'f' END AS sex,
               p.anchor_year_group,
               CASE WHEN COALESCE(a.insurance, l.insurance) IN ('Medicare', 'Medicaid', 'Private') THEN COALESCE(a.insurance, l.insurance)
                    ELSE 'Other or unknown' END AS insurance_group,
               CASE WHEN COALESCE(a.race, l.race) LIKE 'WHITE%' OR COALESCE(a.race, l.race) = 'PORTUGUESE' THEN 'White'
                    WHEN COALESCE(a.race, l.race) LIKE 'BLACK%' THEN 'Black'
                    WHEN COALESCE(a.race, l.race) LIKE 'HISPANIC%' OR COALESCE(a.race, l.race) = 'SOUTH AMERICAN' THEN 'Hispanic or Latino'
                    WHEN COALESCE(a.race, l.race) LIKE 'ASIAN%' THEN 'Asian'
                    ELSE 'Other or unknown' END AS race_group,
               a.admission_location,
               CAST(COALESCE(e.micro_hadm_id, n.hadm_id) IS NOT NULL AND p.anchor_year_group IN ({_quoted(EMAR_YEAR_GROUPS)})
                    AND m.hadm_id IS NOT NULL AS INTEGER) AS emar_covered
        FROM {episodes} e
        LEFT JOIN enc n USING (culture_id)
        LEFT JOIN ed d USING (culture_id)
        LEFT JOIN last_adm l USING (culture_id)
        LEFT JOIN _mimic_adm a ON a.hadm_id = COALESCE(e.micro_hadm_id, n.hadm_id)
        LEFT JOIN emar_adm m ON m.hadm_id = COALESCE(e.micro_hadm_id, n.hadm_id)
        LEFT JOIN (
            SELECT e2.culture_id, pt.gender, pt.anchor_year_group,
                   TRY_CAST(pt.anchor_age AS INTEGER) + year(e2.order_time) - TRY_CAST(pt.anchor_year AS INTEGER) AS age,
                   TRY_CAST(pt.anchor_age AS INTEGER) + year(e2.order_time) - TRY_CAST(pt.anchor_year AS INTEGER) AS age_exact
            FROM {episodes} e2 JOIN mimic_patients pt ON pt.subject_id = e2.subject_id
        ) p USING (culture_id)
    """)
    return "mimic_context"


# Features --------------------------------------------------------------------------------------------------

def _prefix_match(column: str, codes) -> str:
    codes = [c for c in codes if c]
    if not codes:
        return "FALSE"
    return "(" + " OR ".join(f"starts_with({column}, '{c}')" for c in codes) + ")"


def _icd_match(codes9, codes10, code: str = "d.icd_code", version: str = "d.icd_version") -> str:
    return f"(({version} = '9' AND {_prefix_match(code, codes9)}) OR ({version} = '10' AND {_prefix_match(code, codes10)}))"


def build_drug_events(con) -> str:
    """mimic_drug_events: subject_id, event_time, code, source for every systemic antibacterial
    prescription start and eMAR administration. mimic_drug_days keeps one row per subject, day and code."""
    rx = _lookup(con, "mimic_rx_lookup", _distinct(con, "mimic_prescriptions", "drug"), medication_code)
    em = _lookup(con, "mimic_emar_lookup", _distinct(con, "mimic_emar", "medication"), medication_code)
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE mimic_drug_events AS
        SELECT r.subject_id, TRY_CAST(r.starttime AS TIMESTAMP) AS event_time, l.value AS code, 'prescription' AS source
        FROM mimic_prescriptions r JOIN {rx} l ON l.raw = r.drug
        WHERE l.value IS NOT NULL AND r.route IN ({_quoted(SYSTEMIC_ROUTES)})
        UNION ALL
        SELECT m.subject_id, TRY_CAST(m.charttime AS TIMESTAMP), l.value, 'emar'
        FROM mimic_emar m JOIN {em} l ON l.raw = m.medication
        WHERE l.value IS NOT NULL AND m.event_txt IN ({_quoted(ADMINISTERED)})
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE mimic_drug_days AS
        SELECT DISTINCT subject_id, CAST(event_time AS DATE) AS d, code FROM mimic_drug_events WHERE event_time IS NOT NULL
    """)
    return "mimic_drug_events"


def _genus_case(name: str) -> str:
    return "CASE " + " ".join(f"WHEN position('{stem}' IN {name}) > 0 THEN '{genus}'" for genus, stem in GENUS_STEMS.items()) + " END"


def build_admission_codes(con, lists: dict) -> None:
    """Per admission: comorbidity flags, Elixhauser categories and Charlson categories from its diagnoses,
    so cultures join to admissions rather than to every code."""
    elix, charl, hf = lists["elixhauser_quan2005"], lists["charlson_quan2005"], lists["heart_failure_ccsr"]
    comorbid = {
        "comorbid_heart_failure": (hf["icd9"], hf["icd10"]),
        "comorbid_liver_disease": (elix["liver_disease"]["icd9"], elix["liver_disease"]["icd10"]),
        "comorbid_lymphoma": (elix["lymphoma"]["icd9"], elix["lymphoma"]["icd10"]),
        "comorbid_metastatic_cancer": (elix["metastatic_cancer"]["icd9"], elix["metastatic_cancer"]["icd10"]),
        "comorbid_obesity": (elix["obesity"]["icd9"], elix["obesity"]["icd10"]),
        "comorbid_renal_failure": (elix["renal_failure"]["icd9"], elix["renal_failure"]["icd10"]),
    }
    flags = ", ".join(f"MAX(CAST({_icd_match(*codes)} AS INTEGER)) AS {n}" for n, codes in comorbid.items())
    con.execute(f"CREATE OR REPLACE TEMP TABLE mimic_adm_flags AS SELECT d.hadm_id, {flags} FROM mimic_diagnoses_icd d GROUP BY d.hadm_id")
    elix_rows = " UNION ALL ".join(f"SELECT DISTINCT d.hadm_id, '{k}' AS category FROM mimic_diagnoses_icd d WHERE {_icd_match(v['icd9'], v['icd10'])}"
                                   for k, v in elix.items())
    con.execute(f"CREATE OR REPLACE TEMP TABLE mimic_adm_elix AS {elix_rows}")
    charl_rows = " UNION ALL ".join(f"SELECT DISTINCT d.hadm_id, '{k}' AS category, {int(v['weight'])} AS weight FROM mimic_diagnoses_icd d "
                                    f"WHERE {_icd_match(v['icd9'], v['icd10'])}" for k, v in charl.items())
    con.execute(f"CREATE OR REPLACE TEMP TABLE mimic_adm_charlson AS {charl_rows}")


def _prior_admissions(episodes: str) -> str:
    return f"""
        SELECT e.culture_id, a.hadm_id
        FROM {episodes} e JOIN _mimic_adm a ON a.subject_id = e.subject_id
        WHERE date_diff('day', CAST(a.dischtime AS DATE), CAST(e.order_time AS DATE)) >= {PRIOR_MIN_DAYS}"""


def build_features(con, lists: dict, episodes: str = "mimic_culture_episodes", context: str = "mimic_context") -> str:
    """mimic_features: culture_id and FEATURES_48, with the definitions and windows of SAP section 4."""
    build_drug_events(con)
    build_admission_codes(con, lists)
    days = "date_diff('day', {d}, CAST(e.order_time AS DATE))"
    abx_flags = ", ".join(f"MAX(CAST(x.code IN ({_quoted(codes)}) AS INTEGER)) AS {name}" for name, codes in ANTIBIOTIC_CLASSES.items())
    abx = f"""
        SELECT e.culture_id, {abx_flags}
        FROM {episodes} e JOIN mimic_drug_days x ON x.subject_id = e.subject_id
        WHERE {days.format(d='x.d')} BETWEEN {PRIOR_MIN_DAYS} AND {ANTIBIOTIC_WINDOW_DAYS}
        GROUP BY e.culture_id"""

    # Prior organisms and resistance from every earlier microbiology result, any specimen type, by day.
    org_all = organisms.register_lookup(con, _distinct(con, "mimic_microbiologyevents", "org_name"), table="mimic_organism_lookup_all")
    ab_all = _lookup(con, "mimic_ab_lookup_all", _distinct(con, "mimic_microbiologyevents", "ab_name"), antibiotic_code)
    vre = " OR ".join(f"starts_with(o.organism_std, '{s}')" for s in VRE_SPECIES)
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE mimic_micro_days AS
        SELECT m.subject_id,
               CAST(COALESCE(TRY_CAST(NULLIF(m.charttime, '') AS TIMESTAMP), TRY_CAST(NULLIF(m.chartdate, '') AS TIMESTAMP)) AS DATE) AS d,
               {_genus_case('o.organism_std')} AS genus,
               MAX(CAST(upper(trim(m.interpretation)) = 'R' AND COALESCE(o.is_enterobacterales, FALSE) AND a.value IN ({_quoted(ESBL_AGENTS)}) AS INTEGER)) AS esbl,
               MAX(CAST((starts_with(o.organism_std, 'STAPHYLOCOCCUS AUREUS') AND ((upper(trim(m.interpretation)) = 'R' AND a.value IN ({_quoted(MRSA_AGENTS)}))
                                                                                     OR m.spec_type_desc = 'MRSA SCREEN'))
                        OR regexp_matches(m.org_name, '{MRSA_NAME_PATTERN}') AS INTEGER)) AS mrsa,
               MAX(CAST(upper(trim(m.interpretation)) = 'R' AND COALESCE(o.is_enterobacterales, FALSE) AND a.value IN ({_quoted(CARBAPENEMS)}) AS INTEGER)) AS cre,
               MAX(CAST(upper(trim(m.interpretation)) = 'R' AND ({vre}) AND a.value = 'VAN' AS INTEGER)) AS vre
        FROM mimic_microbiologyevents m
        JOIN {org_all} o ON o.organism_raw = m.org_name
        LEFT JOIN {ab_all} a ON a.raw = m.ab_name
        WHERE m.org_name IS NOT NULL AND m.org_name NOT IN ({_quoted(NOT_AN_ORGANISM)})
        GROUP BY 1, 2, 3
    """)
    res = f"""
        SELECT e.culture_id, MAX(p.esbl) AS prior_ESBL, MAX(p.mrsa) AS prior_MRSA, MAX(p.cre) AS prior_CRE, MAX(p.vre) AS prior_VRE
        FROM {episodes} e JOIN mimic_micro_days p ON p.subject_id = e.subject_id
        WHERE {days.format(d='p.d')} >= {PRIOR_MIN_DAYS}
        GROUP BY e.culture_id"""
    orgs = f"""
        SELECT e.culture_id, COUNT(DISTINCT p.genus) AS num_prior_orgs, MIN({days.format(d='p.d')}) AS days_since_prior_org
        FROM {episodes} e JOIN mimic_micro_days p ON p.subject_id = e.subject_id
        WHERE p.genus IS NOT NULL AND {days.format(d='p.d')} >= {PRIOR_MIN_DAYS}
        GROUP BY e.culture_id"""

    proc_lists = lists["procedures"]
    proc_map = {"proc_cvc": "cvc", "proc_mechvent": "mechvent", "proc_dialysis": "dialysis",
                "proc_urethral": "urethral_catheter", "proc_surgical": "surgical"}
    proc_match = {n: _icd_match(proc_lists[k].get("icd9", []), proc_lists[k].get("icd10pcs", []), "x.icd_code", "x.icd_version")
                  for n, k in proc_map.items()}
    proc_flags = ", ".join(f"MAX(CAST({m} AS INTEGER)) AS {n}" for n, m in proc_match.items())
    any_proc = " OR ".join(proc_match.values())
    proc_day = "CAST(TRY_CAST(x.chartdate AS TIMESTAMP) AS DATE)"
    proc = f"""
        SELECT e.culture_id, {proc_flags},
               MAX(CAST({days.format(d=proc_day)} <= {PROCEDURE_WINDOW_DAYS} AND ({any_proc}) AS INTEGER)) AS any_proc_30d
        FROM {episodes} e JOIN mimic_procedures_icd x ON x.subject_id = e.subject_id
        WHERE {days.format(d=proc_day)} >= {PRIOR_MIN_DAYS}
        GROUP BY e.culture_id"""

    nh = f"""
        SELECT e.culture_id, 1 AS nursing_home_90d
        FROM {episodes} e JOIN _mimic_adm a ON a.subject_id = e.subject_id
        WHERE (a.admission_location IN ({_quoted(NURSING_HOME_ADMISSION)})
               AND {days.format(d='CAST(a.admittime AS DATE)')} BETWEEN 0 AND {NURSING_HOME_WINDOW_DAYS})
           OR (a.discharge_location IN ({_quoted(NURSING_HOME_DISCHARGE)})
               AND {days.format(d='CAST(a.dischtime AS DATE)')} BETWEEN 0 AND {NURSING_HOME_WINDOW_DAYS})
        GROUP BY e.culture_id"""

    con.execute(f"CREATE OR REPLACE TEMP TABLE mimic_prior_adm AS {_prior_admissions(episodes)}")
    comorbid = ("comorbid_heart_failure", "comorbid_liver_disease", "comorbid_lymphoma", "comorbid_metastatic_cancer",
                "comorbid_obesity", "comorbid_renal_failure")
    com = f"""
        SELECT f.culture_id, {', '.join(f'f.{n}' for n in comorbid)}, COALESCE(c.elixhauser_count, 0) AS elixhauser_count
        FROM (SELECT p.culture_id, {', '.join(f'MAX(a.{n}) AS {n}' for n in comorbid)}
              FROM mimic_prior_adm p JOIN mimic_adm_flags a USING (hadm_id) GROUP BY p.culture_id) f
        LEFT JOIN (SELECT p.culture_id, COUNT(DISTINCT x.category) AS elixhauser_count
                   FROM mimic_prior_adm p JOIN mimic_adm_elix x USING (hadm_id) GROUP BY p.culture_id) c USING (culture_id)"""

    parts = {"abx": abx, "res": res, "org": orgs, "proc": proc, "nh": nh, "com": com}
    owners = {
        "abx": tuple(ANTIBIOTIC_CLASSES),
        "res": ("prior_ESBL", "prior_MRSA", "prior_CRE", "prior_VRE"),
        "org": ("num_prior_orgs", "days_since_prior_org"),
        "proc": (*proc_map, "any_proc_30d"),
        "nh": ("nursing_home_90d",),
        "com": ("elixhauser_count", *comorbid),
    }
    joined = "\n".join(f"LEFT JOIN ({sql}) {alias} USING (culture_id)" for alias, sql in parts.items())
    history = ",\n".join(f"COALESCE({alias}.{name}, 0) AS {name}" for alias, names in owners.items() for name in names)
    age_case = " ".join(f"WHEN '{band}' THEN {i}" for i, band in enumerate(vocab.AGE_BANDS))
    base = f"""
        SELECT e.culture_id,
               CASE c.age_band {age_case} ELSE 4 END AS age_encoded,
               COALESCE(CAST(c.sex = 'm' AS INTEGER), 0) AS gender_male,
               0 AS adi_score_clean, 0 AS adi_high, 1 AS adi_missing,
               c.hosp_ward_IP, c.hosp_ward_OP, c.hosp_ward_ER, 0 AS hosp_ward_UC,
               CAST(e.specimen_type = 'urine' AS INTEGER) AS specimen_urine,
               CAST(e.specimen_type = 'blood' AS INTEGER) AS specimen_blood,
               CAST(e.specimen_type = 'respiratory' AS INTEGER) AS specimen_respiratory,
               {history}
        FROM {episodes} e
        JOIN {context} c USING (culture_id)
        {joined}
    """
    derived = {
        "age_65plus": "CAST(age_encoded >= 5 AS INTEGER)",
        "fq_and_ceph3": "CAST(T_fq_90d > 0 AND T_ceph3_90d > 0 AS INTEGER)",
        "multi_prior_resistant": "prior_ESBL + prior_CRE + prior_MRSA + prior_VRE",
        "any_immunocomp": "CAST(comorbid_lymphoma = 1 OR comorbid_metastatic_cancer = 1 AS INTEGER)",
        "any_invasive_proc": "CAST(proc_dialysis = 1 OR proc_mechvent = 1 OR proc_urethral = 1 OR proc_cvc = 1 AS INTEGER)",
        "prior_esbl_and_ceph3": "CAST(prior_ESBL > 0 AND T_ceph3_90d > 0 AS INTEGER)",
        "hospital_multiabx": "CAST(hosp_ward_IP > 0 AND (T_fq_90d + T_ceph3_90d + T_carb_90d + T_glyco_90d) >= 2 AS INTEGER)",
        "renal_urine": "CAST(comorbid_renal_failure > 0 AND specimen_urine = 1 AS INTEGER)",
        "nursing_home_flag": "CAST(nursing_home_90d > 0 AS INTEGER)",
    }
    select = ", ".join(f"CAST({derived.get(n, n)} AS {'DOUBLE' if n in DOUBLE_FEATURES else 'INTEGER'}) AS {n}" for n in FEATURES_48)
    con.execute(f"CREATE OR REPLACE TEMP TABLE mimic_features AS SELECT culture_id, {select} FROM ({base}) b")
    return "mimic_features"


EXTRA_FLAGS = ("emar_covered", "clin_carbapenem", "clin_carbapenem_started_after", "clin_broad", "treated_at_order",
               "tum_hosp_12m", "tum_transfer", "tum_charlson4", "tum_age70", "tum_catheter_30d")
EXTRA_CATEGORIES = ("insurance_group", "race_group", "anchor_year_group")


def build_extras(con, lists: dict, episodes: str = "mimic_culture_episodes", context: str = "mimic_context") -> str:
    """mimic_extras: per culture, eMAR coverage, clinicians' empiric choices, exact Tumbarello items and
    fairness subgroups (SAP amendment E-9). Needs build_features first (drug events, prior admissions)."""
    w = CLINICIAN_WINDOW_HOURS
    # eMAR doses within a day of collection: join by subject and calendar day, then compare times.
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE _mimic_emar_events AS
        SELECT subject_id, event_time, CAST(event_time AS DATE) AS d, code FROM mimic_drug_events
        WHERE source = 'emar' AND event_time IS NOT NULL
    """)
    hours = "date_diff('minute', e.order_time, x.event_time) / 60.0"
    clin = f"""
        SELECT e.culture_id,
               MAX(CAST(x.code IN ({_quoted(CARBAPENEM_CODES)}) AND {hours} BETWEEN -{w} AND {w} AS INTEGER)) AS clin_carbapenem,
               MAX(CAST(x.code IN ({_quoted(CARBAPENEM_CODES)}) AND {hours} > 0 AND {hours} <= {w} AS INTEGER)) AS carb_after,
               MAX(CAST(x.code IN ({_quoted(CARBAPENEM_CODES)}) AND {hours} BETWEEN -{w} AND 0 AS INTEGER)) AS carb_before,
               MAX(CAST(x.code IN ({_quoted(BROAD_CODES)}) AND {hours} BETWEEN -{w} AND {w} AS INTEGER)) AS clin_broad,
               MAX(CAST({hours} BETWEEN -{w} AND {w} AS INTEGER)) AS treated_at_order
        FROM {episodes} e
        JOIN _mimic_emar_events x ON x.subject_id = e.subject_id
         AND x.d BETWEEN CAST(e.order_time AS DATE) - 1 AND CAST(e.order_time AS DATE) + 1
        GROUP BY e.culture_id"""
    hosp = f"""
        SELECT e.culture_id, 1 AS tum_hosp_12m
        FROM {episodes} e JOIN _mimic_adm a ON a.subject_id = e.subject_id
        WHERE date_diff('day', CAST(a.dischtime AS DATE), CAST(e.order_time AS DATE)) BETWEEN {PRIOR_MIN_DAYS} AND 365
        GROUP BY e.culture_id"""
    charl = lists["charlson_quan2005"]
    roles = {r: _charlson_key(charl, r) for r in ("mild_liver", "severe_liver", "diabetes_uncomplicated", "diabetes_complicated",
                                                 "malignancy", "metastatic")}
    # Quan 2005 hierarchy: the severe form replaces the mild one in liver disease, diabetes and cancer.
    charlson = f"""
        WITH cats AS (SELECT DISTINCT p.culture_id, c.category, c.weight FROM mimic_prior_adm p JOIN mimic_adm_charlson c USING (hadm_id))
        SELECT culture_id, SUM(weight) AS charlson FROM cats r
        WHERE NOT EXISTS (SELECT 1 FROM cats s WHERE s.culture_id = r.culture_id AND (
              (r.category = '{roles['mild_liver']}' AND s.category = '{roles['severe_liver']}') OR
              (r.category = '{roles['diabetes_uncomplicated']}' AND s.category = '{roles['diabetes_complicated']}') OR
              (r.category = '{roles['malignancy']}' AND s.category = '{roles['metastatic']}')))
        GROUP BY culture_id"""
    cath = lists["procedures"]["urethral_catheter"]
    catheter = f"""
        SELECT e.culture_id, 1 AS tum_catheter_30d
        FROM {episodes} e JOIN mimic_procedures_icd x ON x.subject_id = e.subject_id
        WHERE date_diff('day', CAST(TRY_CAST(x.chartdate AS TIMESTAMP) AS DATE), CAST(e.order_time AS DATE)) BETWEEN {PRIOR_MIN_DAYS} AND 30
          AND {_icd_match(cath.get('icd9', []), cath.get('icd10pcs', []), 'x.icd_code', 'x.icd_version')}
        GROUP BY e.culture_id"""
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE mimic_extras AS
        SELECT e.culture_id, c.emar_covered, c.insurance_group, c.race_group, c.anchor_year_group,
               COALESCE(k.clin_carbapenem, 0) AS clin_carbapenem,
               CAST(COALESCE(k.carb_after, 0) = 1 AND COALESCE(k.carb_before, 0) = 0 AS INTEGER) AS clin_carbapenem_started_after,
               COALESCE(k.clin_broad, 0) AS clin_broad,
               COALESCE(k.treated_at_order, 0) AS treated_at_order,
               COALESCE(h.tum_hosp_12m, 0) AS tum_hosp_12m,
               CAST(COALESCE(c.admission_location IN ({_quoted(TRANSFER_ADMISSION)}), FALSE) AS INTEGER) AS tum_transfer,
               CAST(COALESCE(ch.charlson, 0) >= 4 AS INTEGER) AS tum_charlson4,
               CAST(COALESCE(c.age_exact >= 70, FALSE) AS INTEGER) AS tum_age70,
               COALESCE(u.tum_catheter_30d, 0) AS tum_catheter_30d
        FROM {episodes} e
        JOIN {context} c USING (culture_id)
        LEFT JOIN ({clin}) k USING (culture_id)
        LEFT JOIN ({hosp}) h USING (culture_id)
        LEFT JOIN ({charlson}) ch USING (culture_id)
        LEFT JOIN ({catheter}) u USING (culture_id)
    """)
    return "mimic_extras"


def _charlson_key(charl: dict, role: str) -> str:
    """Category name in the code list for each role in the Quan 2005 hierarchy."""
    patterns = {"mild_liver": r"mild.*liver|liver.*mild", "severe_liver": r"(moderate|severe).*liver|liver.*(moderate|severe)",
                "diabetes_uncomplicated": r"diabetes.*(without|uncomplicated)", "diabetes_complicated": r"diabetes.*(with_|complicat)(?!.*without)",
                "malignancy": r"^(any_)?malignan|^cancer|^malignancy", "metastatic": r"metasta"}
    for key in charl:
        if re.search(patterns[role], key) and not (role == "diabetes_complicated" and re.search(r"without|uncomplicated", key)):
            return key
    raise KeyError(f"no Charlson category for {role} in {sorted(charl)}")


# Split and cohorts -------------------------------------------------------------------------------------------

def _caller_may_read_test() -> bool:
    return Path(sys.argv[0]).name == SEALED_READER


def write_split(con, seed: int, path: Path = SPLIT_PATH, episodes: str = "mimic_culture_episodes",
                features: str = "mimic_features") -> dict:
    """Patient split stratified by an earlier ESBL phenotype, a feature, so no outcome is read. Frozen once."""
    rows = con.execute(f"""
        SELECT e.patient_id, MAX(f.prior_ESBL) AS stratum
        FROM {episodes} e JOIN {features} f USING (culture_id) GROUP BY e.patient_id
    """).fetchall()
    assignments = splits.assign([r[0] for r in rows], [int(r[1]) for r in rows], seed=seed)
    if path.exists():
        existing = path.read_bytes()
        tmp = path.with_suffix(".check.csv.gz")
        splits.write_manifest(SITE, assignments, path=tmp, force=True)
        same = tmp.read_bytes() == existing
        tmp.unlink()
        if not same:
            raise RuntimeError(f"{path} is frozen and the regenerated split differs")
    else:
        splits.write_manifest(SITE, assignments, path=path)
    return {"path": str(path), "sha256": sha256_file(path), "patients": len(rows),
            "strata": {str(k): sum(1 for r in rows if int(r[1]) == k) for k in (0, 1)}}


def register_split(con, include_test: bool = False, path: Path = SPLIT_PATH) -> str:
    if include_test and not _caller_may_read_test():
        raise splits.SealedPartitionError(f"only analysis/{SEALED_READER} may read the MIMIC test partition")
    name = f"_mimic_split_{uuid.uuid4().hex[:8]}"
    keep = "TRUE" if include_test else "partition <> 'test'"
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE {name} AS
        SELECT patient_id_hash, partition, TRY_CAST(NULLIF(fold, '') AS INTEGER) AS fold
        FROM read_csv('{path}', header = true, all_varchar = true) WHERE {keep}
    """)
    return name


@dataclass(frozen=True)
class MimicSpec:
    estimand: str
    label: str = "composite"

    @property
    def name(self) -> str:
        return f"{self.estimand.lower()}_{self.label}_{REPEAT_RULE}_{UNIT}"


def _label(con, spec: MimicSpec, built: dict) -> labels.LabelResult:
    agents = CEFTRIAXONE_AGENTS if spec.label == "ceftriaxone" else labels.SENTINEL_AGENTS
    build = labels.esbl_phenotype_composite if spec.estimand == "E1" else labels.esbl_decision_time
    return build(con, built["culture_episodes"], built["isolate_ast"], agents=agents)


def build_cohort(con, spec: MimicSpec, built: dict, partitions: tuple[str, ...], features: str = "mimic_features",
                 context: str = "mimic_context", extras: str = "mimic_extras") -> tuple[str, list]:
    """One row per evaluable unit in the partitions, with the columns 02_train_e1 and 03_evaluate_ladder read."""
    episodes = built["culture_episodes"]
    split = register_split(con, include_test="test" in partitions)
    label = _label(con, spec, built)
    tag = uuid.uuid4().hex[:8]
    kept = kept_cultures(con, episodes, REPEAT_RULE, f"_mimic_kept_{tag}")
    units = unit_labels(con, label.frame.alias, episodes, kept, UNIT, spec.estimand, f"_mimic_units_{tag}")
    in_partitions = ", ".join(f"'{p}'" for p in partitions)
    split_join = f"JOIN {split} s ON s.patient_id_hash = sha256(u.patient_id)"
    flow = []

    def record(step: str, where: str, with_split: bool) -> None:
        join = split_join if with_split else ""
        positives = "COUNT(*) FILTER (WHERE u.label = 1)" if with_split else "NULL"
        n, c, p, pos = con.execute(f"SELECT COUNT(*), SUM(u.n_cultures), COUNT(DISTINCT u.patient_id), {positives} "
                                   f"FROM {units} u {join} WHERE {where}").fetchone()
        flow.append({"step": step, "units": n, "cultures": c or 0, "patients": p, "label_positive": pos})

    n, p = con.execute(f"SELECT COUNT(*), COUNT(DISTINCT patient_id) FROM {episodes}").fetchone()
    flow.append({"step": "ordered cultures", "units": None, "cultures": n, "patients": p, "label_positive": None})
    record(f"index orders under {REPEAT_RULE}", "TRUE", with_split=False)
    record(f"evaluable for {spec.estimand}", "u.evaluable", with_split=False)
    conditions = f"u.evaluable AND s.partition IN ({in_partitions})"
    record(f"in partition {', '.join(partitions)}", conditions, with_split=True)

    setting_rank = "CASE m.setting_primary " + " ".join(f"WHEN '{v}' THEN {i}" for i, v in enumerate(SETTING_PRECEDENCE)) + " END"
    table = f"mimic_cohort_{spec.name}_{tag}"
    feature_columns = ", ".join(f"f.{n}" for n in FEATURES_48)
    feature_aggregates = ", ".join(f"{unit_feature_sql(n)} AS {n}" for n in FEATURES_48)
    extra_flags = ", ".join(f"MAX(m.{n}) AS {n}" for n in EXTRA_FLAGS)
    extra_categories = ", ".join(f"MIN(m.{n}) AS {n}" for n in EXTRA_CATEGORIES)
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE {table} AS
        WITH m AS (
            SELECT {unit_key_sql(UNIT)} AS unit_id, e.culture_id, c.setting_primary, e.specimen_type, c.age_band, c.sex,
                   {', '.join(f'x.{n}' for n in EXTRA_FLAGS + EXTRA_CATEGORIES)}, {feature_columns}
            FROM {episodes} e
            JOIN {kept} k USING (culture_id)
            JOIN {features} f USING (culture_id)
            JOIN {context} c USING (culture_id)
            JOIN {extras} x USING (culture_id)
        )
        SELECT '{SITE}' AS site, u.patient_id, u.unit_id, MIN(m.culture_id) AS culture_id, u.n_cultures,
               s.partition, s.fold, u.label, u.reason,
               arg_min(m.setting_primary, {setting_rank}) AS setting_primary,
               CASE WHEN COUNT(DISTINCT m.specimen_type) = 1 THEN MIN(m.specimen_type) ELSE 'multiple' END AS specimen_type,
               MIN(m.age_band) AS age_band, MIN(m.sex) AS sex,
               CAST(MAX(m.num_prior_orgs) > 0 OR MAX(m.multi_prior_resistant) > 0 AS INTEGER) AS prior_culture_history,
               'missing' AS adi_quintile,
               {extra_flags}, {extra_categories}, {feature_aggregates}
        FROM {units} u
        {split_join}
        JOIN m USING (unit_id)
        WHERE {conditions}
        GROUP BY u.patient_id, u.unit_id, u.n_cultures, s.partition, s.fold, u.label, u.reason
    """)
    return table, flow


def provenance(spec: MimicSpec) -> dict:
    return {"spec": asdict(spec), "name": spec.name, "repeat_rule": REPEAT_RULE, "unit": UNIT,
            "split_manifest": {"path": str(SPLIT_PATH), "sha256": sha256_file(SPLIT_PATH) if SPLIT_PATH.exists() else None}}
