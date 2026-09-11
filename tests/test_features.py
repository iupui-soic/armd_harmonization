"""Harmonized features against synthetic releases that exercise every inclusion and exclusion rule."""

import pytest

from armd_common import features
from armd_common.db import connect
from armd_common.sites import MGBLoader, StanfordLoader
from test_sites import IDS, MGB_COHORT, build_dryad_release, mgb_row, write_table


def feature_rows(loader, table, names):
    columns = ", ".join(names)
    return {r[0]: dict(zip(names, r[1:])) for r in loader.con.execute(f"SELECT culture_id, {columns} FROM {table} ORDER BY 1").fetchall()}


@pytest.fixture()
def stanford(tmp_path, monkeypatch):
    monkeypatch.setenv("ARMD_DATA_ROOT", str(tmp_path))
    site = "armd-stanford"
    build_dryad_release(tmp_path, site, "order_time_jittered_utc")
    write_table(tmp_path, site, "antibiotic_subtype_exposure", f"{IDS}, medication_name VARCHAR, medication_time_to_cultureTime DOUBLE", [
        ("p1", "e2", "102", "Cefepime In", 30.0),
        ("p1", "e2", "102", "Levaquin", 0.5),
        ("p1", "e2", "102", "Piperacillin-Tazobactam-Dextrs", 100.0),
        ("p2", "e3", "103", "Bactrim Ds", 10.0),
        ("p2", "e3", "103", "Vancomycin In Dextrose", 5.0),
        ("p2", "e3", "103", "Ertapenem", 20.0),
        ("p2", "e3", "103", "Ciprofloxacin-Dexamethasone", 3.0),
    ])
    write_table(tmp_path, site, "microbial_resistance", f"{IDS}, organism VARCHAR, antibiotic VARCHAR, resistant_time_to_culturetime DOUBLE", [
        ("p1", "e2", "102", "ESCHERICHIA COLI", "Ceftriaxone", 200.0),
        ("p1", "e2", "102", "ESCHERICHIA COLI", "Meropenem", 0.0),
        ("p2", "e3", "103", "STAPH AUREUS {MRSA}", "Erythromycin", 50.0),
        ("p2", "e3", "103", "ENTEROCOCCUS FAECIUM -  VANCO RESISTANT", "Vancomycin", 30.0),
        ("p3", "e5", "105", "ENTEROCOCCUS GALLINARUM", "Vancomycin", 30.0),
        ("p3", "e5", "105", "METHICILLIN RESISTANT STAPHYLOCOCCUS LUGDUNENSIS", "Cefoxitin", 10.0),
        ("p3", "e6", "106", "KLEBSIELLA PNEUMONIAE", "Ertapenem", 15.0),
    ])
    write_table(tmp_path, site, "prior_infecting_organism", f"{IDS}, prior_organism VARCHAR, prior_infecting_organism_days_to_culutre DOUBLE", [
        ("p1", "e2", "102", "Escherichia", 200.0),
        ("p1", "e2", "102", "Klebsiella", 400.0),
        ("p1", "e2", "102", "CONS", 10.0),
        ("p1", "e2", "102", "Escherichia", 0.2),
    ])
    write_table(tmp_path, site, "prior_procedures", f"{IDS}, procedure_description VARCHAR, procedure_time_to_culturetime DOUBLE", [
        ("p2", "e3", "103", "cvc", 10.0),
        ("p2", "e3", "103", "urethral_catheter", 45.0),
        ("p3", "e6", "106", "parenteral_nutrition", 5.0),
        ("p3", "e6", "106", "dialysis", 0.5),
    ])
    write_table(tmp_path, site, "nursing_home_visits", f"{IDS}, nursing_home_visit_culture DOUBLE", [
        ("p1", "e1", "101", 0.0),
        ("p1", "e2", "102", None),
        ("p2", "e4", "104", 120.0),
    ])
    write_table(tmp_path, site, "comorbidity", f"{IDS}, comorbidity_component VARCHAR, comorbidity_component_start_days_culture DOUBLE", [
        ("p2", "e3", "103", "Congestive heart failure", 100.0),
        ("p2", "e3", "103", "Heart failure", 100.0),
        ("p2", "e3", "103", "Metastatic cancer", 50.0),
        ("p2", "e3", "103", "Renal failure", 20.0),
        ("p3", "e6", "106", "Obesity", 300.0),
        ("p3", "e6", "106", "Acute and unspecified renal failure", 10.0),
    ])
    loader = StanfordLoader(connect())
    built = loader.build()
    return loader, features.build_features(loader, built)


def test_feature_table_shape(stanford):
    loader, table = stanford
    columns = [r[0] for r in loader.con.execute(f"DESCRIBE {table}").fetchall()]
    assert columns == ["culture_id", *features.FEATURES_48]
    assert len(features.FEATURES_45) == 45 and len(features.FEATURES_48) == 48
    assert loader.con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0] == 6


def test_stanford_features(stanford):
    loader, table = stanford
    got = feature_rows(loader, table, features.FEATURES_48)
    expected_nonzero = {
        "stanford:p1:101": {"hosp_ward_OP": 1, "nursing_home_90d": 1, "nursing_home_flag": 1, "adi_score_clean": 5.0, "specimen_urine": 1},
        "stanford:p1:102": {
            "hosp_ward_IP": 1, "hosp_ward_ER": 1, "T_ceph3_90d": 1, "prior_ESBL": 1, "num_prior_orgs": 2,
            "days_since_prior_org": 200.0, "multi_prior_resistant": 1, "prior_esbl_and_ceph3": 1, "specimen_urine": 1,
            "adi_score_clean": 5.0,
        },
        "stanford:p2:103": {
            "age_encoded": 8, "age_65plus": 1, "gender_male": 1, "adi_missing": 1, "hosp_ward_IP": 1,
            "comorbid_heart_failure": 1, "comorbid_metastatic_cancer": 1, "comorbid_renal_failure": 1, "elixhauser_count": 4,
            "T_carb_90d": 1, "T_glyco_90d": 1, "T_sulfa_90d": 1, "prior_MRSA": 1, "prior_VRE": 1, "multi_prior_resistant": 2,
            "proc_cvc": 1, "proc_urethral": 1, "any_proc_30d": 1, "any_immunocomp": 1, "any_invasive_proc": 1,
            "hospital_multiabx": 1, "specimen_blood": 1,
        },
        "stanford:p2:104": {"age_encoded": 8, "age_65plus": 1, "gender_male": 1, "adi_missing": 1, "specimen_respiratory": 1},
        "stanford:p3:105": {"age_encoded": 4, "adi_score_clean": 88.0, "adi_high": 1, "hosp_ward_OP": 1, "specimen_urine": 1},
        "stanford:p3:106": {
            "age_encoded": 4, "adi_score_clean": 88.0, "adi_high": 1, "hosp_ward_IP": 1, "prior_CRE": 1,
            "multi_prior_resistant": 1, "comorbid_obesity": 1, "comorbid_renal_failure": 1, "elixhauser_count": 2,
            "renal_urine": 1, "specimen_urine": 1,
        },
    }
    mismatches = {}
    for culture_id, nonzero in expected_nonzero.items():
        expected = {name: nonzero.get(name, 0) for name in features.FEATURES_48}
        diff = {name: (got[culture_id][name], value) for name, value in expected.items() if got[culture_id][name] != value}
        if diff:
            mismatches[culture_id] = diff
    assert mismatches == {}, f"(got, expected) per feature: {mismatches}"


@pytest.fixture()
def mgb(tmp_path, monkeypatch):
    monkeypatch.setenv("ARMD_DATA_ROOT", str(tmp_path))
    site = "armd-mgb"
    write_table(tmp_path, site, "cohort", MGB_COHORT, [
        mgb_row("1", "10", "500", "URINE", "ESCHERICHIA COLI", code="CRO", pheno="Susceptible", clsi="Susceptible"),
        mgb_row("2", "20", "600", "BLOOD", neg="X"),
    ])
    write_table(tmp_path, site, "ward_info",
                f"{IDS}, order_time_jittered_utc_shifted DATE, hosp_ward_IP BIGINT, hosp_ward_OP BIGINT, hosp_ward_ER BIGINT, hosp_ward_UC BIGINT, hosp_ward_day_surg BIGINT", [
                    ("1", "10", "500", "2090-05-01", 0, 0, 0, 1, 0),
                    ("2", "20", "600", "2090-05-01", 1, 0, 0, 0, 0),
                ])
    write_table(tmp_path, site, "demographics", f"{IDS}, age VARCHAR, gender VARCHAR", [
        ("1", "10", "500", "65-74 years", "Female"),
        ("2", "20", "600", "25-34 years", "Male"),
    ])
    write_table(tmp_path, site, "adi_scores", f"{IDS}, adi_score DOUBLE, adi_state_rank DOUBLE, order_time_jittered_utc_shifted DATE", [
        ("1", "10", "500", 12.5, 9.0, "2090-05-01"),
        ("2", "20", "600", None, None, "2090-05-01"),
    ])
    write_table(tmp_path, site, "prior_abx", f"{IDS}, last_dose_to_culture DOUBLE, drug_code VARCHAR, drug_class VARCHAR", [
        ("1", "10", "500", -1.0, "CIP", "fluoroquinolone"),
        ("1", "10", "500", 0.0, "VAN", "glycopeptide"),
        ("1", "10", "500", 5.0, "CRO", "extended_spectrum_cephalosporin"),
        ("1", "10", "500", 20.0, "CZA", "beta_lactam_combo"),
        ("2", "20", "600", 60.0, "TZP", "extended_spectrum_penicillin"),
        ("2", "20", "600", 91.0, "GEN", "aminoglycoside"),
    ])
    write_table(tmp_path, site, "prior_micro", f"{IDS}, organism VARCHAR, drug_code VARCHAR, AST_pheno VARCHAR, CLSI_2022_pheno VARCHAR, prior_AST_time_to_culture DOUBLE", [
        ("1", "10", "500", "KLEBSIELLA PNEUMONIAE", "CRO", "Susceptible", "Resistant", 10.0),
        ("1", "10", "500", "STAPHYLOCOCCUS AUREUS", "OXA", "Resistant", "Resistant", 300.0),
        ("2", "20", "600", "ESCHERICHIA COLI", "MEM", "Resistant", "Resistant", 40.0),
    ])
    # Comorbidity IDs come from another mapping, as in the real release, so it links by patient and date.
    write_table(tmp_path, site, "comorbidity", f"{IDS}, ICD10 VARCHAR, category VARCHAR, order_time_jittered_utc_shifted DATE", [
        ("1", "77", "9001", "N18.6", "Chronic kidney disease", "2090-05-01"),
        ("1", "77", "9001", "E66.9", "Obesity", "2090-05-01"),
    ])
    loader = MGBLoader(connect())
    return loader, loader.build()


def test_mgb_features_lab(mgb):
    loader, built = mgb
    table = features.build_features(loader, built, "lab")
    names = ["T_fq_90d", "T_ceph3_90d", "T_glyco_90d", "T_esp_90d", "T_amino_90d", "prior_ESBL", "prior_MRSA", "prior_CRE",
             "comorbid_obesity", "elixhauser_count", "hosp_ward_UC", "hosp_ward_IP", "gender_male", "age_encoded", "age_65plus",
             "adi_high", "adi_missing"]
    got = feature_rows(loader, table, names)
    assert got["mgb:1:500"] == dict(zip(names, [0, 1, 0, 0, 0, 0, 1, 0, 1, 2, 1, 0, 0, 5, 1, 1, 0]))
    assert got["mgb:2:600"] == dict(zip(names, [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1]))


def test_mgb_features_clsi2022(mgb):
    loader, built = mgb
    table = features.build_features(loader, built, "clsi2022")
    assert feature_rows(loader, table, ["prior_ESBL"])["mgb:1:500"] == {"prior_ESBL": 1}


def test_rejects_unknown_interpretation(mgb):
    loader, built = mgb
    with pytest.raises(ValueError):
        features.build_features(loader, built, "eucast")
