"""Hand-built fixtures for the MIMIC-IV mapping into armd-common's label constructors (SAP amendment E-9).

Rule 4 of the project: every label path needs susceptible, intermediate, resistant, not-tested and
organism-absent cases before use. The MIMIC loader adds its own specimen map, interpretation codes,
antibiotic overrides and cancelled-culture rule, so each is covered here.
"""
import sys
from pathlib import Path

import duckdb
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))

from armd_common import labels  # noqa: E402
from armd_common.cohorts import CEFTRIAXONE_AGENTS  # noqa: E402
from armd_common.splits import SealedPartitionError  # noqa: E402
from mimic import build  # noqa: E402

COLUMNS = ("micro_specimen_id", "subject_id", "hadm_id", "chartdate", "charttime", "spec_type_desc", "test_name",
           "org_name", "ab_name", "interpretation")
SENTINELS = ("CEFTRIAXONE", "CEFTAZIDIME", "CEFEPIME", "PIPERACILLIN/TAZO")


def urine(specimen, org=None, results=(), spec="URINE", test="URINE CULTURE"):
    base = (str(specimen), "1", None, "2150-01-01 00:00:00", f"2150-01-{specimen:02d} 10:00:00", spec, test)
    if org is None:
        return [base + (None, None, None)]
    if not results:
        return [base + (org, None, None)]
    return [base + (org, ab, sir) for ab, sir in results]


def all_sentinels(sir, override=None):
    return [(ab, (override or {}).get(ab, sir)) for ab in SENTINELS]


ROWS = (
    urine(1, "ESCHERICHIA COLI", all_sentinels("S"))                                     # susceptible
    + urine(2, "ESCHERICHIA COLI", all_sentinels("S", {"CEFTRIAXONE": "I"}))             # intermediate
    + urine(3, "ESCHERICHIA COLI", all_sentinels("S", {"CEFTRIAXONE": "R"}))             # resistant
    + urine(4, "ESCHERICHIA COLI", [("NITROFURANTOIN", "S"), ("GENTAMICIN", "S")])       # not tested
    + urine(5, "STAPH AUREUS COAG +", [("OXACILLIN", "R")])                              # organism absent
    + urine(6)                                                                           # negative
    + urine(7, "CANCELLED")                                                              # cancelled
    + urine(8, "ESCHERICHIA COLI", all_sentinels("S", {"PIPERACILLIN/TAZO": "R"}))       # piperacillin-tazobactam only
    + urine(9, "KLEBSIELLA PNEUMONIAE", all_sentinels("S", {"CEFTRIAXONE": "D"}))        # dose-dependent
    + urine(10, None, spec="SEROLOGY/BLOOD", test="HIV-1 Viral Load")                    # not a culture
    + urine(11, "ESCHERICHIA COLI", spec="BLOOD CULTURE", test="Aerobic Bottle Gram Stain")  # Gram stain only
    + urine(12, "ESCHERICHIA COLI", all_sentinels("S", {"CEFTRIAXONE": "R"}), spec="BLOOD CULTURE",
            test="Blood Culture, Routine")                                               # resistant blood
    + urine(13, "ESCHERICHIA COLI")                                                      # isolate with no AST rows
)


@pytest.fixture()
def con():
    c = duckdb.connect()
    c.execute(f"CREATE TABLE mimic_microbiologyevents ({', '.join(f'{col} VARCHAR' for col in COLUMNS)})")
    c.executemany(f"INSERT INTO mimic_microbiologyevents VALUES ({', '.join('?' for _ in COLUMNS)})", ROWS)
    yield c
    c.close()


def reasons(con, result):
    return dict(con.execute(f"SELECT culture_id, reason FROM {result.frame.alias}").fetchall())


def label_map(con, result):
    return dict(con.execute(f"SELECT culture_id, label FROM {result.frame.alias}").fetchall())


def test_specimen_map_and_cancelled(con):
    built = build.build_cultures(con)
    kept = {r[0] for r in con.execute("SELECT culture_id FROM mimic_culture_episodes").fetchall()}
    assert kept == {"1", "2", "3", "4", "5", "6", "8", "9", "12", "13"}
    assert built["counts"]["cancelled"] == 1
    specimen = dict(con.execute("SELECT culture_id, specimen_type FROM mimic_culture_episodes").fetchall())
    assert specimen["12"] == "blood" and specimen["1"] == "urine"
    positive = dict(con.execute("SELECT culture_id, culture_positive FROM mimic_culture_episodes").fetchall())
    assert positive["6"] is False and positive["5"] is True


def test_decision_time_labels(con):
    built = build.build_cultures(con)
    result = labels.esbl_decision_time(con, built["culture_episodes"], built["isolate_ast"])
    r, y = reasons(con, result), label_map(con, result)
    assert r["1"] == "not_resistant" and y["1"] == 0
    assert r["2"] == "not_resistant" and y["2"] == 0      # intermediate is not resistant
    assert r["3"] == "resistant" and y["3"] == 1
    assert r["4"] == "not_tested" and y["4"] is None      # no sentinel result: not evaluable
    assert r["5"] == "organism_absent" and y["5"] == 0
    assert r["6"] == "culture_negative" and y["6"] == 0
    assert r["8"] == "resistant" and y["8"] == 1          # PIPERACILLIN/TAZO maps to TZP
    assert r["9"] == "not_resistant" and y["9"] == 0      # D is never resistant
    assert r["12"] == "resistant"
    assert r["13"] == "not_tested" and y["13"] is None


def test_positive_culture_labels(con):
    built = build.build_cultures(con)
    result = labels.esbl_phenotype_composite(con, built["culture_episodes"], built["isolate_ast"])
    r, y = reasons(con, result), label_map(con, result)
    assert y["1"] == 0 and y["2"] == 0 and y["3"] == 1 and y["8"] == 1 and y["9"] == 0
    assert y["4"] is None and r["4"] == "not_tested"
    assert y["5"] is None and r["5"] == "organism_absent"   # outside E1
    assert y["6"] is None and r["6"] == "culture_negative"  # outside E1


def test_ceftriaxone_label(con):
    built = build.build_cultures(con)
    result = labels.esbl_decision_time(con, built["culture_episodes"], built["isolate_ast"], agents=CEFTRIAXONE_AGENTS)
    y = label_map(con, result)
    assert y["3"] == 1 and y["8"] == 0                      # piperacillin-tazobactam alone is not ceftriaxone resistance


def test_overrides():
    assert build.antibiotic_code("PIPERACILLIN/TAZO") == "TZP"
    assert build.antibiotic_code("TRIMETHOPRIM/SULFA") == "SXT"
    assert build.antibiotic_code("CEFTRIAXONE") == "CRO"
    assert build.medication_code("Sulfameth/Trimethoprim DS") == "SXT"
    assert build.medication_code("Piperacillin-Tazobactam Na") == "TZP"
    assert build.medication_code("Meropenem") == "MEM"


def test_test_partition_is_sealed(tmp_path):
    with pytest.raises(SealedPartitionError):
        build.register_split(duckdb.connect(), include_test=True, path=tmp_path / "split.csv.gz")


STUB_LISTS = {
    "elixhauser_quan2005": {k: {"icd9": ["9999"], "icd10": ["Z999"]} for k in
                            ("liver_disease", "lymphoma", "metastatic_cancer", "obesity", "renal_failure", "congestive_heart_failure")},
    "charlson_quan2005": {"mild_liver_disease": {"icd9": ["5712"], "icd10": ["K703"], "weight": 1},
                          "moderate_severe_liver_disease": {"icd9": ["5722"], "icd10": ["K704"], "weight": 3},
                          "diabetes_without_complications": {"icd9": ["2500"], "icd10": ["E119"], "weight": 1},
                          "diabetes_with_complications": {"icd9": ["2504"], "icd10": ["E112"], "weight": 2},
                          "any_malignancy": {"icd9": ["174"], "icd10": ["C50"], "weight": 2},
                          "metastatic_solid_tumor": {"icd9": ["196"], "icd10": ["C78"], "weight": 6}},
    "heart_failure_ccsr": {"icd9": ["428"], "icd10": ["I50"]},
    "procedures": {k: {"icd9": ["9999"], "icd10pcs": ["ZZZZZZZ"]} for k in ("cvc", "mechvent", "dialysis", "urethral_catheter", "surgical")},
}
STUB_LISTS["elixhauser_quan2005"]["renal_failure"] = {"icd9": ["585"], "icd10": ["N18"]}
STUB_LISTS["elixhauser_quan2005"]["metastatic_cancer"] = {"icd9": ["196"], "icd10": ["C77", "C78", "C79", "C80"]}


def _table(con, name, columns, rows):
    con.execute(f"CREATE TABLE {name} ({', '.join(f'{c} VARCHAR' for c in columns)})")
    if rows:
        con.executemany(f"INSERT INTO {name} VALUES ({', '.join('?' for _ in columns)})", rows)


def test_features_and_extras_smoke(con):
    # Patient 1: culture 3 (resistant E. coli, day 3) and culture 1 (day 1). Culture 12 is blood on day 12.
    _table(con, "mimic_patients", ("subject_id", "gender", "anchor_age", "anchor_year", "anchor_year_group", "dod"),
           [("1", "F", "72", "2150", "2017 - 2019", None)])
    _table(con, "mimic_admissions", ("subject_id", "hadm_id", "admittime", "dischtime", "deathtime", "admission_type",
                                     "admit_provider_id", "admission_location", "discharge_location", "insurance", "language",
                                     "marital_status", "race", "edregtime", "edouttime", "hospital_expire_flag"),
           [("1", "100", "2149-12-01 00:00:00", "2149-12-10 00:00:00", None, "URGENT", None, "TRANSFER FROM HOSPITAL",
             "SKILLED NURSING FACILITY", "Medicare", None, None, "WHITE", None, None, "0"),
            ("1", "101", "2150-01-12 11:00:00", "2150-01-20 00:00:00", None, "URGENT", None, "EMERGENCY ROOM", "HOME",
             "Medicaid", None, None, "WHITE", "2150-01-12 09:00:00", "2150-01-12 11:00:00", "0")])
    _table(con, "mimic_transfers", ("subject_id", "hadm_id", "transfer_id", "eventtype", "careunit", "intime", "outtime"),
           [("1", "101", "1", "ED", "Emergency Department", "2150-01-12 09:00:00", "2150-01-12 11:00:00")])
    _table(con, "mimic_diagnoses_icd", ("subject_id", "hadm_id", "seq_num", "icd_code", "icd_version"),
           [("1", "100", "1", "N186", "10"), ("1", "100", "2", "C787", "10"), ("1", "100", "3", "C509", "10")])
    _table(con, "mimic_procedures_icd", ("subject_id", "hadm_id", "seq_num", "chartdate", "icd_code", "icd_version"), [])
    _table(con, "mimic_prescriptions", ("subject_id", "hadm_id", "drug", "route", "starttime"),
           [("1", "100", "Ciprofloxacin HCl", "PO", "2149-12-05 09:00:00"), ("1", "100", "Ciprofloxacin HCl", "TP", "2149-12-05 09:00:00")])
    _table(con, "mimic_emar", ("subject_id", "hadm_id", "medication", "event_txt", "charttime"),
           [("1", "101", "Meropenem", "Administered", "2150-01-12 12:00:00"),
            ("1", "101", "Meropenem", "Not Given", "2150-01-12 14:00:00")])
    built = build.build_cultures(con)
    build.build_context(con)
    build.build_features(con, STUB_LISTS)
    build.build_extras(con, STUB_LISTS)
    f = {r[0]: r[1:] for r in con.execute(
        "SELECT culture_id, T_fq_90d, prior_ESBL, nursing_home_90d, comorbid_renal_failure, comorbid_metastatic_cancer, "
        "hosp_ward_IP, hosp_ward_ER, hosp_ward_OP, age_65plus, gender_male FROM mimic_features").fetchall()}
    assert f["1"][0] == 1                       # ciprofloxacin 27 days before, oral route
    assert f["3"][1] == 0 and f["12"][1] == 1   # culture 3 is resistant; culture 12 nine days later has an earlier ESBL phenotype
    assert f["1"][2] == 1                       # discharged to a nursing facility 22 days before
    assert f["1"][3] == 1 and f["1"][4] == 1    # comorbidities from the earlier admission
    assert f["12"][5] == 1 and f["12"][6] == 1  # blood culture drawn in the ED, then admitted
    assert f["1"][7] == 1                       # urine before any admission: outpatient
    assert f["1"][8] == 1 and f["1"][9] == 0    # age 72, female
    x = {r[0]: r[1:] for r in con.execute(
        "SELECT culture_id, emar_covered, clin_carbapenem, clin_carbapenem_started_after, treated_at_order, tum_hosp_12m, "
        "tum_transfer, tum_charlson4, tum_age70, insurance_group FROM mimic_extras").fetchall()}
    assert x["12"][0] == 1 and x["12"][1] == 1 and x["12"][2] == 1 and x["12"][3] == 1  # meropenem 2 hours after collection
    assert x["1"][4] == 1                       # discharged 22 days before
    assert x["12"][5] == 0                      # current admission came from the ED
    assert x["1"][6] == 1                       # renal 2 + metastatic 6, malignancy dropped by the hierarchy
    assert x["1"][7] == 1 and x["12"][8] == "Medicaid"
