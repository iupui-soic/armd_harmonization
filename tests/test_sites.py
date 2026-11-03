"""Site loaders against tiny synthetic Parquet releases that mimic each site's quirks."""

import duckdb
import pytest

from armd_common import labels
from armd_common.db import connect
from armd_common.schema import validate
from armd_common.sites import MGB2Loader, MGBLoader, StanfordLoader, UTSWLoader


def write_table(root, site_dir, table, ddl, rows, part=0):
    con = duckdb.connect()
    con.execute("SET TimeZone = 'UTC'")
    con.execute(f"CREATE TABLE t ({ddl})")
    if rows:
        placeholders = ", ".join("?" for _ in rows[0])
        con.executemany(f"INSERT INTO t VALUES ({placeholders})", rows)
    out = root / site_dir / "parquet" / table
    # A later part adds a file to a table that already exists.
    out.mkdir(parents=True, exist_ok=part > 0)
    path = out / f"data_{part}.parquet"
    con.execute(f"COPY t TO '{path}' (FORMAT parquet)")
    con.close()


def rows_by_id(con, table, columns):
    return {r[0]: r[1:] for r in con.execute(f"SELECT culture_id, {columns} FROM {table} ORDER BY 1").fetchall()}


def reasons(result):
    return {cid: (label, reason) for cid, label, _, reason in result.frame.order("culture_id").fetchall()}


IDS = "anon_id VARCHAR, pat_enc_csn_id_coded VARCHAR, order_proc_id_coded VARCHAR"


def build_dryad_release(root, site_dir, time_column):
    cohort = (
        f"{IDS}, {time_column} TIMESTAMPTZ, ordering_mode VARCHAR, culture_description VARCHAR, was_positive BIGINT, "
        "organism VARCHAR, antibiotic VARCHAR, susceptibility VARCHAR"
    )
    write_table(root, site_dir, "cohort", cohort, [
        ("p1", "e1", "101", "2020-01-01 10:00:00+00", "Outpatient", "URINE", 0, None, None, None),
        ("p1", "e2", "102", "2020-02-01 08:30:00+00", "Inpatient", "URINE", 1, "ESCHERICHIA COLI", "Ceftriaxone", "Resistant"),
        ("p1", "e2", "102", "2020-02-01 08:30:00+00", "Inpatient", "URINE", 1, "ESCHERICHIA COLI", "Ampicillin", "Resistant"),
        ("p2", "e3", "103", "2019-05-05 00:00:00+00", "Inpatient", "BLOOD", 1, "STAPHYLOCOCCUS AUREUS", "Oxacillin", "Susceptible"),
        ("p2", "e4", "104", "2019-06-06 00:00:00+00", None, "RESPIRATORY", 1, None, None, None),
        ("p3", "e5", "105", "2021-03-03 12:00:00+00", "Outpatient", "URINE", 1, "KLEBSIELLA PNEUMONIAE", "Nitrofurantoin", "Susceptible"),
        ("p3", "e6", "106", "2021-04-04 12:00:00+00", "Inpatient", "URINE", 1, "ZZZENTEROBACTER AEROGENES", "Cefepime", "Intermediate"),
        ("p3", "e6", "106", "2021-04-04 12:00:00+00", "Inpatient", "URINE", 1, "PSEUDOMONAS AERUGINOSA", "Ceftazidime", "Resistant"),
    ])
    write_table(root, site_dir, "ward_info", f"{IDS}, hosp_ward_IP BIGINT, hosp_ward_OP BIGINT, hosp_ward_ER BIGINT, hosp_ward_ICU BIGINT", [
        ("p1", "e1", "101", 0, 1, 0, 0),
        ("p1", "e2", "102", 1, 0, 1, 0),
        ("p2", "e3", "103", 1, 0, 0, 1),
        ("p2", "e4", "104", 0, 0, 0, 0),
        ("p3", "e5", "105", 0, 1, 0, 0),
        ("p3", "e6", "106", 1, 0, 0, 0),
    ])
    write_table(root, site_dir, "demographics", f"{IDS}, age VARCHAR, gender BIGINT", [
        ("p1", "e1", "101", "18-24 years", 0),
        ("p1", "e2", "102", "18-24 years", 0),
        ("p2", "e3", "103", "above 90", 1),
        ("p2", "e4", "104", "above 90", 1),
        ("p3", "e5", "105", "55-64 years", None),
        ("p3", "e6", "106", "55-64 years", None),
    ])
    write_table(root, site_dir, "adi_scores", f"{IDS}, adi_score BIGINT, adi_state_rank BIGINT", [
        ("p1", "e1", "101", 5, 1),
        ("p1", "e2", "102", 5, 1),
        ("p2", "e3", "103", None, None),
        ("p2", "e4", "104", None, None),
        ("p3", "e5", "105", 88, 9),
        ("p3", "e6", "106", 88, 9),
    ])


@pytest.fixture()
def stanford(tmp_path, monkeypatch):
    monkeypatch.setenv("ARMD_DATA_ROOT", str(tmp_path))
    build_dryad_release(tmp_path, "armd-stanford", "order_time_jittered_utc")
    loader = StanfordLoader(connect())
    return loader, loader.build()


def test_stanford_schema(stanford):
    loader, built = stanford
    for kind, table in built.items():
        assert validate(loader.con, table, kind) == [], kind


def test_stanford_culture_episodes(stanford):
    loader, built = stanford
    got = rows_by_id(loader.con, built["culture_episodes"],
                     "culture_positive, organism_group, has_enterobacterales, organism_unidentified, setting_primary, specimen_type, n_organisms")
    assert got == {
        "stanford:p1:101": (False, "none", False, False, "outpatient", "urine", 0),
        "stanford:p1:102": (True, "ecoli", True, False, "inpatient", "urine", 1),
        "stanford:p2:103": (True, "non_enterobacterales", False, False, "inpatient", "blood", 1),
        "stanford:p2:104": (True, None, False, True, "unknown", "respiratory", 0),
        "stanford:p3:105": (True, "klebsiella", True, False, "outpatient", "urine", 1),
        "stanford:p3:106": (True, "klebsiella", True, False, "inpatient", "urine", 2),
    }
    details = rows_by_id(loader.con, built["culture_episodes"],
                         "CAST(order_time AS VARCHAR), setting_raw_ed, setting_raw_icu, age_band, sex, sex_code_masked, adi_score, adi_missing")
    assert details["stanford:p1:102"] == ("2020-02-01 08:30:00", True, False, "18-24", "f", "0", 5.0, False)
    assert details["stanford:p2:103"] == ("2019-05-05 00:00:00", False, True, "90+", "m", "1", None, True)
    assert details["stanford:p3:105"][3:6] == ("55-64", "unknown", None)
    inferred = rows_by_id(loader.con, built["culture_episodes"], "sex_inferred")
    assert inferred["stanford:p1:102"] == (True,) and inferred["stanford:p3:105"] == (False,)


def test_stanford_isolate_ast(stanford):
    loader, built = stanford
    rows = loader.con.execute(
        f"SELECT culture_id, antibiotic_code, sir, organism_group, preliminary FROM {built['isolate_ast']} ORDER BY 1, 2"
    ).fetchall()
    assert rows == [
        ("stanford:p1:102", "AMP", "R", "ecoli", False),
        ("stanford:p1:102", "CRO", "R", "ecoli", False),
        ("stanford:p2:103", "OXA", "S", "non_enterobacterales", False),
        ("stanford:p3:105", "NIT", "S", "klebsiella", False),
        ("stanford:p3:106", "CAZ", "R", "non_enterobacterales", False),
        ("stanford:p3:106", "FEP", "I", "klebsiella", False),
    ]


STANFORD_E2 = {
    "p1:101": (0, "culture_negative"),
    "p1:102": (1, "resistant"),
    "p2:103": (0, "organism_absent"),
    "p2:104": (None, "organism_unidentified"),
    "p3:105": (None, "not_tested"),
    "p3:106": (0, "not_resistant"),
}


def test_stanford_decision_time_labels(stanford):
    loader, built = stanford
    got = reasons(labels.esbl_decision_time(loader.con, built["culture_episodes"], built["isolate_ast"]))
    assert got == {f"stanford:{k}": v for k, v in STANFORD_E2.items()}


def test_utsw_order_time_spelling(tmp_path, monkeypatch):
    monkeypatch.setenv("ARMD_DATA_ROOT", str(tmp_path))
    build_dryad_release(tmp_path, "armd-utsw", "order_time_jittered")
    loader = UTSWLoader(connect())
    built = loader.build()
    for kind, table in built.items():
        assert validate(loader.con, table, kind) == [], kind
    got = reasons(labels.esbl_decision_time(loader.con, built["culture_episodes"], built["isolate_ast"]))
    assert got == {f"utsw:{k}": v for k, v in STANFORD_E2.items()}
    assert loader.con.execute(f"SELECT CAST(order_date AS VARCHAR) FROM {built['culture_keys']} WHERE culture_id = 'utsw:p1:102'").fetchone()[0] == "2020-02-01"


MGB_COHORT = (
    f"{IDS}, culture_description VARCHAR, organism VARCHAR, neg_cx VARCHAR, mult_org_ast VARCHAR, has_AST VARCHAR, "
    "prelim_AST VARCHAR, AST_panel VARCHAR, enzyme_class VARCHAR, enzyme VARCHAR, AST_code VARCHAR, antibiotic VARCHAR, "
    "AST_inequality VARCHAR, AST_val1 DOUBLE, AST_val2 DOUBLE, AST_pheno VARCHAR, CLSI_2022_pheno VARCHAR, "
    "order_time_jittered_utc_shifted DATE"
)


def mgb_row(anon, enc, proc, spec, organism=None, neg=None, prelim=None, enzyme_class=None, code=None, pheno=None, clsi=None):
    antibiotic = {"CRO": "ceftriaxone", "FEP": "cefepime", "CAZ": "ceftazidime"}.get(code)
    has_ast = None if neg else "X"
    return (anon, enc, proc, spec, organism, neg, None, has_ast, prelim, None if neg else "broth_microdilution",
            enzyme_class, None, code, antibiotic, "<=" if code else None, 1.0 if code else None, None, pheno, clsi, "2090-05-01")


@pytest.fixture()
def mgb(tmp_path, monkeypatch):
    monkeypatch.setenv("ARMD_DATA_ROOT", str(tmp_path))
    ecoli = "ESCHERICHIA COLI"
    rows = []
    # Culture 500 of patient 1 is linked to two encounters, so every row appears twice.
    for enc in ("10", "11"):
        rows += [
            mgb_row("1", enc, "500", "URINE", ecoli, code="CRO", pheno="Susceptible", clsi="Resistant"),
            mgb_row("1", enc, "500", "URINE", ecoli, prelim="X", code="CRO", pheno="Susceptible", clsi="Susceptible"),
            mgb_row("1", enc, "500", "URINE", ecoli, prelim="X", code="FEP", pheno="Susceptible", clsi="Susceptible dose dependent"),
            mgb_row("1", enc, "500", "URINE", ecoli, enzyme_class="ESBL", code="ESBL", pheno="Positive"),
        ]
    # Patient 2 reuses order_proc_id_coded 500 for a different, negative culture on the same day as culture 501.
    rows.append(mgb_row("2", "20", "500", "BLOOD", neg="X"))
    rows.append(mgb_row("2", "21", "501", "RESPIRATORY_TRACT", code="CAZ", pheno="Resistant", clsi="Resistant"))
    write_table(tmp_path, "armd-mgb", "cohort", MGB_COHORT, rows)
    write_table(tmp_path, "armd-mgb", "ward_info",
                f"{IDS}, order_time_jittered_utc_shifted DATE, hosp_ward_IP BIGINT, hosp_ward_OP BIGINT, hosp_ward_ER BIGINT, hosp_ward_UC BIGINT, hosp_ward_day_surg BIGINT", [
                    ("1", "10", "500", "2090-05-01", 0, 0, 1, 0, 0),
                    ("1", "11", "500", "2090-05-01", 1, 0, 0, 0, 0),
                    ("2", "20", "500", "2090-05-01", 0, 1, 0, 0, 0),
                    ("2", "21", "501", "2090-05-01", 0, 0, 0, 1, 0),
                ])
    write_table(tmp_path, "armd-mgb", "demographics", f"{IDS}, age VARCHAR, gender VARCHAR", [
        ("1", "10", "500", "65-74 years", "Female"),
        ("1", "11", "500", "65-74 years", "Female"),
        ("2", "20", "500", "25-34 years", "Male"),
        ("2", "21", "501", "above 90", "X"),
    ])
    write_table(tmp_path, "armd-mgb", "adi_scores", f"{IDS}, adi_score DOUBLE, adi_state_rank DOUBLE, order_time_jittered_utc_shifted DATE", [
        ("1", "10", "500", 12.5, 2.0, "2090-05-01"),
        ("1", "11", "500", 12.5, 2.0, "2090-05-01"),
        ("2", "20", "500", None, None, "2090-05-01"),
        ("2", "21", "501", 50.0, 7.0, "2090-05-01"),
    ])
    # Comorbidity uses encounter and culture IDs from another mapping, as in the real release.
    write_table(tmp_path, "armd-mgb", "comorbidity", f"{IDS}, ICD10 VARCHAR, category VARCHAR, order_time_jittered_utc_shifted DATE", [
        ("1", "77", "9001", "N39.0", "renal", "2090-05-01"),
        ("2", "78", "9002", "E11.9", "diabetes", "2090-05-01"),
        ("2", "79", "9003", "E11.9", "diabetes", "2090-05-01"),
    ])
    loader = MGBLoader(connect())
    return loader, loader.build()


def test_mgb_schema(mgb):
    loader, built = mgb
    for kind, table in built.items():
        assert validate(loader.con, table, kind) == [], kind


def test_mgb_culture_episodes(mgb):
    loader, built = mgb
    got = rows_by_id(loader.con, built["culture_episodes"],
                     "n_encounters, culture_positive, organism_group, organism_unidentified, setting_primary, setting_raw_ed, "
                     "sex, sex_code_masked, age_band, adi_score, order_time_resolution, CAST(order_time AS VARCHAR)")
    assert got == {
        "mgb:1:500": (2, True, "ecoli", False, "inpatient", True, "f", None, "65-74", 12.5, "day", "2090-05-01 00:00:00"),
        "mgb:2:500": (1, False, "none", False, "outpatient", False, "m", None, "25-34", None, "day", "2090-05-01 00:00:00"),
        "mgb:2:501": (1, True, None, True, "urgent_care", False, "unknown", None, "90+", 50.0, "day", "2090-05-01 00:00:00"),
    }


def test_mgb_isolate_ast_collapses_encounters_and_preliminary(mgb):
    loader, built = mgb
    rows = loader.con.execute(
        f"SELECT culture_id, antibiotic_code, sir, sir_clsi2022, preliminary, n_source_rows FROM {built['isolate_ast']} ORDER BY 1, 2"
    ).fetchall()
    assert rows == [
        ("mgb:1:500", "CRO", "S", "R", False, 2),
        ("mgb:1:500", "FEP", "S", "I", True, 2),
        ("mgb:2:501", "CAZ", "R", "R", False, 1),
    ]


def test_mgb_labels(mgb):
    loader, built = mgb
    ep, ast, enz = built["culture_episodes"], built["isolate_ast"], built["enzyme_tests"]
    assert reasons(labels.esbl_phenotype_composite(loader.con, ep, ast))["mgb:1:500"] == (0, "not_resistant")
    assert reasons(labels.esbl_phenotype_composite(loader.con, ep, ast, sir_column="sir_clsi2022"))["mgb:1:500"] == (1, "resistant")
    assert reasons(labels.esbl_enzyme_class(loader.con, ep, enz))["mgb:1:500"] == (1, "resistant")
    e2 = reasons(labels.esbl_decision_time(loader.con, ep, ast))
    assert e2["mgb:2:500"] == (0, "culture_negative")
    assert e2["mgb:2:501"] == (None, "organism_unidentified")


def test_mgb_date_linked_feature_table(mgb):
    loader, _ = mgb
    assert loader.feature_link("comorbidity")[2] == "patient_and_order_date"
    assert loader.feature_link("ward_info")[2] == "culture_id"
    rows = loader.con.execute(f"SELECT culture_id, ICD10 FROM ({loader.feature_rows_sql('comorbidity')}) ORDER BY 1, 2").fetchall()
    # Patient 2 had two cultures that day, so both receive the patient's comorbidity rows.
    assert rows == [
        ("mgb:1:500", "N39.0"),
        ("mgb:2:500", "E11.9"),
        ("mgb:2:500", "E11.9"),
        ("mgb:2:501", "E11.9"),
        ("mgb:2:501", "E11.9"),
    ]


# ARMD-MGB 2.0.0 follows the Stanford layout with its own quirks ----------------------------------

def build_mgb2_release(root):
    site = "armd-mgb/v2"
    cohort = (
        f"{IDS}, ordering_mode VARCHAR, culture_description VARCHAR, was_positive BIGINT, organism VARCHAR, "
        "antibiotic VARCHAR, susceptibility VARCHAR, order_time_jittered_utc_shifted TIMESTAMP"
    )
    ecoli = "ESCHERICHIA COLI"
    rows = [
        # Culture 500 spans two encounters, so every row appears twice.
        ("1", "10", "500", "Inpatient", "URINE", 1, ecoli, "Ceftriaxone", "Resistant", "2090-05-01 10:15:00"),
        ("1", "11", "500", "Inpatient", "URINE", 1, ecoli, "Ceftriaxone", "Resistant", "2090-05-01 10:15:00"),
        ("1", "10", "500", "Inpatient", "URINE", 1, ecoli, "Cefepime", "Unknown", "2090-05-01 10:15:00"),
        ("1", "11", "500", "Inpatient", "URINE", 1, ecoli, "Cefepime", "Unknown", "2090-05-01 10:15:00"),
        ("2", "20", "600", "Urgent Care", "BLOOD", 0, None, None, None, "2091-01-02 08:00:00"),
        ("3", "30", "700", "Emergency", "RESPIRATORY_TRACT", 1, "KLEBSIELLA PNEUMONIAE", "Cefepime", "Susceptible", "2092-03-03 23:30:00"),
    ]
    write_table(root, site, "cohort", cohort, rows)
    write_table(root, site, "ward_info",
                f"{IDS}, hosp_ward_IP BIGINT, hosp_ward_OP BIGINT, hosp_ward_ER BIGINT, hosp_ward_ICU BIGINT, order_time_jittered_utc_shifted TIMESTAMP", [
                    ("1", "10", "500", 1, 0, 0, None, "2090-05-01 10:15:00"),
                    ("1", "11", "500", 1, 0, 0, None, "2090-05-01 10:15:00"),
                    ("2", "20", "600", 0, 1, 0, None, "2091-01-02 08:00:00"),
                    ("3", "30", "700", 0, 0, 1, None, "2092-03-03 23:30:00"),
                ])
    write_table(root, site, "demographics", f"{IDS}, age VARCHAR, gender VARCHAR", [
        ("1", "10", "500", "85-89 years", "Female"),
        ("1", "11", "500", "85-89 years", "Female"),
        ("2", "20", "600", "above 90", "Others"),
        ("3", "30", "700", "25-34 years", "Male"),
    ])
    write_table(root, site, "adi_scores", f"{IDS}, adi_score DOUBLE, adi_state_rank DOUBLE, order_time_jittered_utc_shifted TIMESTAMP", [
        ("1", "10", "500", 12.0, 2.0, "2090-05-01 10:15:00"),
        ("1", "11", "500", 12.0, 2.0, "2090-05-01 10:15:00"),
        ("3", "30", "700", 90.0, 9.0, "2092-03-03 23:30:00"),
    ])
    # 2.0.0 fixed the ID mismatch, so comorbidity rows carry the cohort's own encounter and culture IDs.
    write_table(root, site, "comorbidity", f"{IDS}, ICD10 VARCHAR, category VARCHAR, order_time_jittered_utc_shifted TIMESTAMP", [
        ("1", "10", "500", "I50.9", "Heart failure", "2090-05-01 10:15:00"),
        ("1", "11", "500", "I50.9", "Heart failure", "2090-05-01 10:15:00"),
        ("1", "10", "500", "E66.9", "Obesity", "2090-05-01 10:15:00"),
        ("3", "30", "700", "N18.6", "Renal failure", "2092-03-03 23:30:00"),
    ])
    write_table(root, site, "prior_procedures", f"{IDS}, procedure_description VARCHAR, procedure_days_culture DOUBLE, order_time_jittered_utc_shifted TIMESTAMP", [
        ("1", "10", "500", "urethral_catheter", 12.0, "2090-05-01 10:15:00"),
        ("3", "30", "700", "cvc", 0.0, "2092-03-03 23:30:00"),
    ])
    write_table(root, site, "antibiotic_subtype_exposure",
                f"{IDS}, medication_name VARCHAR, medication_time_to_cultureTime DOUBLE, medication_category VARCHAR, order_time_jittered_utc_shifted TIMESTAMP", [
                    ("1", "10", "500", "CEFTRIAXONE", 20.0, "CEF2", "2090-05-01 10:15:00"),
                    ("3", "30", "700", "PIPERACILLIN-TAZOBACTAM", 100.0, "PIP1", "2092-03-03 23:30:00"),
                ])
    write_table(root, site, "microbial_resistance",
                f"{IDS}, organism VARCHAR, antibiotic VARCHAR, resistant_time_to_cultureTime DOUBLE, order_time_jittered_utc_shifted TIMESTAMP", [
                    ("1", "10", "500", ecoli, "Ceftriaxone", 300.0, "2090-05-01 10:15:00"),
                ])
    write_table(root, site, "prior_infecting_organism",
                f"{IDS}, prior_infecting_organism_days_to_culture DOUBLE, prior_organism VARCHAR, order_time_jittered_utc_shifted TIMESTAMP", [
                    ("1", "10", "500", 300.0, "Escherichia", "2090-05-01 10:15:00"),
                    ("1", "10", "500", 45.0, "Enterococcus", "2090-05-01 10:15:00"),
                ])
    write_table(root, site, "nursing_home_visits", f"{IDS}, nursing_home_visit_culture DOUBLE, order_time_jittered_utc_shifted TIMESTAMP", [
        ("1", "10", "500", None, "2090-05-01 10:15:00"),
        ("3", "30", "700", 30.0, "2092-03-03 23:30:00"),
    ])


@pytest.fixture()
def mgb2(tmp_path, monkeypatch):
    monkeypatch.setenv("ARMD_DATA_ROOT", str(tmp_path))
    build_mgb2_release(tmp_path)
    loader = MGB2Loader(connect())
    return loader, loader.build()


def test_mgb2_schema(mgb2):
    loader, built = mgb2
    for kind, table in built.items():
        assert validate(loader.con, table, kind) == [], kind
    assert loader.con.execute(f"SELECT COUNT(*) FROM {built['enzyme_tests']}").fetchone()[0] == 0


def test_mgb2_culture_episodes(mgb2):
    loader, built = mgb2
    got = rows_by_id(loader.con, built["culture_episodes"],
                     "n_encounters, culture_positive, organism_group, setting_primary, setting_raw_urgent_care, setting_raw_icu, "
                     "sex, sex_code_masked, age_band, adi_score, order_time_resolution, CAST(order_time AS VARCHAR), ordering_mode_raw")
    assert got == {
        "mgb2:1:500": (2, True, "ecoli", "inpatient", False, None, "f", None, "85-89", 12.0, "second", "2090-05-01 10:15:00", "Inpatient"),
        "mgb2:2:600": (1, False, "none", "urgent_care", True, None, "unknown", None, "90+", None, "second", "2091-01-02 08:00:00", "Urgent Care"),
        "mgb2:3:700": (1, True, "klebsiella", "ed", False, None, "m", None, "25-34", 90.0, "second", "2092-03-03 23:30:00", "Emergency"),
    }


def test_mgb2_isolate_ast_collapses_encounters(mgb2):
    loader, built = mgb2
    rows = loader.con.execute(
        f"SELECT culture_id, antibiotic_code, sir, sir_raw, n_source_rows FROM {built['isolate_ast']} ORDER BY 1, 2"
    ).fetchall()
    assert rows == [
        ("mgb2:1:500", "CRO", "R", "Resistant", 2),
        ("mgb2:1:500", "FEP", "other", "Unknown", 2),
        ("mgb2:3:700", "FEP", "S", "Susceptible", 1),
    ]


def test_mgb2_labels(mgb2):
    loader, built = mgb2
    ep, ast = built["culture_episodes"], built["isolate_ast"]
    assert reasons(labels.esbl_phenotype_composite(loader.con, ep, ast))["mgb2:1:500"] == (1, "resistant")
    e2 = reasons(labels.esbl_decision_time(loader.con, ep, ast))
    assert e2["mgb2:2:600"] == (0, "culture_negative")
    assert e2["mgb2:3:700"] == (0, "not_resistant")


def test_mgb2_feature_tables_link_by_culture(mgb2):
    loader, _ = mgb2
    for table in ("comorbidity", "nursing_home_visits", "prior_procedures"):
        assert loader.feature_link(table)[2] == "culture_id"
    rows = loader.con.execute(f"SELECT culture_id, ICD10 FROM ({loader.feature_rows_sql('comorbidity')}) ORDER BY 1, 2").fetchall()
    assert rows == [("mgb2:1:500", "E66.9"), ("mgb2:1:500", "I50.9"), ("mgb2:1:500", "I50.9"), ("mgb2:3:700", "N18.6")]
