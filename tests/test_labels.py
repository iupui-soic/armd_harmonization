"""Label constructors against hand-built fixtures.

Every constructor is checked on susceptible, intermediate, resistant, not-tested and
organism-absent cultures, plus negative cultures, unidentified isolates, pooled isolates and
non-sentinel resistance.
"""

import duckdb
import pytest

from armd_common import labels

# culture_id, culture_positive, organism_unidentified
CULTURES = [
    ("neg", False, False),
    ("staph_cro_r", True, False),
    ("unidentified", True, True),
    ("ecoli_s", True, False),
    ("ecoli_cro_i", True, False),
    ("kpn_caz_r", True, False),
    ("ecoli_nit_only", True, False),
    ("ecoli_other_nt", True, False),
    ("ecoli_two_isolates", True, False),
    ("pae_r_plus_ecoli_s", True, False),
    ("ecoli_amp_r", True, False),
    ("ecl_cro_r", True, False),
    ("pmi_cro_i", True, False),
    ("unidentified_plus_ecoli_r", True, True),
]

# culture_id, is_enterobacterales, idsa_esbl_species, antibiotic_code, sir
AST = [
    ("staph_cro_r", False, False, "CRO", "R"),
    ("ecoli_s", True, True, "CRO", "S"),
    ("ecoli_s", True, True, "CAZ", "S"),
    ("ecoli_cro_i", True, True, "CRO", "I"),
    ("ecoli_cro_i", True, True, "CAZ", "S"),
    ("ecoli_cro_i", True, True, "FEP", "S"),
    ("ecoli_cro_i", True, True, "TZP", "S"),
    ("kpn_caz_r", True, True, "CAZ", "R"),
    ("kpn_caz_r", True, True, "CRO", "S"),
    ("ecoli_nit_only", True, True, "NIT", "S"),
    ("ecoli_other_nt", True, True, "CRO", "other"),
    ("ecoli_other_nt", True, True, "CAZ", "NT"),
    ("ecoli_two_isolates", True, True, "CRO", "S"),
    ("ecoli_two_isolates", True, True, "CRO", "R"),
    ("pae_r_plus_ecoli_s", False, False, "CAZ", "R"),
    ("pae_r_plus_ecoli_s", True, True, "CRO", "S"),
    ("ecoli_amp_r", True, True, "AMP", "R"),
    ("ecoli_amp_r", True, True, "CRO", "S"),
    ("ecl_cro_r", True, False, "CRO", "R"),
    ("pmi_cro_i", True, True, "CRO", "I"),
    ("unidentified_plus_ecoli_r", True, True, "CRO", "R"),
    ("unidentified_plus_ecoli_r", None, None, "CRO", "S"),
]


@pytest.fixture()
def con():
    con = duckdb.connect()
    con.execute("CREATE TABLE cultures (culture_id VARCHAR, culture_positive BOOLEAN, organism_unidentified BOOLEAN)")
    con.executemany("INSERT INTO cultures VALUES (?, ?, ?)", CULTURES)
    con.execute("CREATE TABLE ast (culture_id VARCHAR, is_enterobacterales BOOLEAN, idsa_esbl_species BOOLEAN, antibiotic_code VARCHAR, sir VARCHAR)")
    con.executemany("INSERT INTO ast VALUES (?, ?, ?, ?, ?)", AST)
    yield con
    con.close()


def collect(result):
    rows = result.frame.order("culture_id").fetchall()
    return {culture_id: (label, evaluable, reason) for culture_id, label, evaluable, reason in rows}


def expected(reasons, rule):
    return {cid: (rule.get(reason), reason in rule, reason) for cid, reason in reasons.items()}


COMPOSITE_REASONS = {
    "neg": "culture_negative",
    "staph_cro_r": "organism_absent",
    "unidentified": "organism_unidentified",
    "ecoli_s": "not_resistant",
    "ecoli_cro_i": "not_resistant",
    "kpn_caz_r": "resistant",
    "ecoli_nit_only": "not_tested",
    "ecoli_other_nt": "not_tested",
    "ecoli_two_isolates": "resistant",
    "pae_r_plus_ecoli_s": "not_resistant",
    "ecoli_amp_r": "not_resistant",
    "ecl_cro_r": "resistant",
    "pmi_cro_i": "not_resistant",
    "unidentified_plus_ecoli_r": "resistant",
}


def test_composite(con):
    result = labels.esbl_phenotype_composite(con, "cultures", "ast")
    assert collect(result) == expected(COMPOSITE_REASONS, labels.E1_RULE)
    assert result.provenance["estimand"] == "E1"
    assert result.provenance["reason_counts"]["resistant"] == 4


def test_composite_nonsusceptible_counts_intermediate(con):
    reasons = {**COMPOSITE_REASONS, "ecoli_cro_i": "resistant", "pmi_cro_i": "resistant"}
    result = labels.esbl_phenotype_composite(con, "cultures", "ast", nonsusceptible=True)
    assert collect(result) == expected(reasons, labels.E1_RULE)


def test_idsa(con):
    reasons = {
        **COMPOSITE_REASONS,
        "ecoli_cro_i": "resistant",
        "kpn_caz_r": "not_resistant",
        "ecl_cro_r": "organism_absent",
        "pmi_cro_i": "resistant",
    }
    result = labels.esbl_phenotype_idsa(con, "cultures", "ast")
    assert collect(result) == expected(reasons, labels.E1_RULE)


def test_decision_time(con):
    result = labels.esbl_decision_time(con, "cultures", "ast")
    got = collect(result)
    assert got == expected(COMPOSITE_REASONS, labels.E2_RULE)
    assert got["neg"] == (0, True, "culture_negative")
    assert got["staph_cro_r"] == (0, True, "organism_absent")
    assert got["unidentified"] == (None, False, "organism_unidentified")
    assert got["ecoli_nit_only"] == (None, False, "not_tested")
    assert result.provenance["estimand"] == "E2"


def test_label_null_exactly_when_not_evaluable(con):
    for build in (labels.esbl_phenotype_composite, labels.esbl_phenotype_idsa, labels.esbl_decision_time):
        rows = build(con, "cultures", "ast").frame.fetchall()
        assert all((label is None) == (not evaluable) for _, label, evaluable, _ in rows)


def test_alternate_sir_column(con):
    con.execute("ALTER TABLE ast ADD COLUMN sir_clsi2022 VARCHAR")
    con.execute("UPDATE ast SET sir_clsi2022 = CASE WHEN culture_id = 'ecoli_s' AND antibiotic_code = 'CRO' THEN 'R' ELSE sir END")
    got = collect(labels.esbl_phenotype_composite(con, "cultures", "ast", sir_column="sir_clsi2022"))
    assert got["ecoli_s"] == (1, True, "resistant")


def test_enzyme_class(con):
    con.execute("CREATE TABLE enzymes (culture_id VARCHAR, is_enterobacterales BOOLEAN, enzyme_class VARCHAR, enzyme_result VARCHAR)")
    con.executemany("INSERT INTO enzymes VALUES (?, ?, ?, ?)", [
        ("staph_cro_r", False, "beta_lactamase", "positive"),
        ("ecoli_s", True, "ESBL", "negative"),
        ("kpn_caz_r", True, "ESBL", "positive"),
        ("ecoli_cro_i", True, "ESBL", "other"),
        ("ecoli_nit_only", True, "carbapenemase", "negative"),
        ("ecoli_two_isolates", True, "ESBL", "negative"),
        ("ecoli_two_isolates", True, "ESBL", "positive"),
        ("ecl_cro_r", True, "ESBL", "negative"),
    ])
    reasons = {
        "neg": "culture_negative",
        "staph_cro_r": "organism_absent",
        "unidentified": "organism_unidentified",
        "ecoli_s": "not_resistant",
        # an indeterminate ESBL result is not a test result
        "ecoli_cro_i": "not_tested",
        "kpn_caz_r": "resistant",
        "ecoli_nit_only": "not_tested",
        "ecoli_other_nt": "organism_absent",
        "ecoli_two_isolates": "resistant",
        "pae_r_plus_ecoli_s": "organism_absent",
        "ecoli_amp_r": "organism_absent",
        # phenotypic CRO resistance does not override a negative ESBL test
        "ecl_cro_r": "not_resistant",
        "pmi_cro_i": "organism_absent",
        "unidentified_plus_ecoli_r": "organism_unidentified",
    }
    assert collect(labels.esbl_enzyme_class(con, "cultures", "enzymes")) == expected(reasons, labels.E1_RULE)


def test_rejects_duplicate_cultures(con):
    con.execute("INSERT INTO cultures VALUES ('neg', FALSE, FALSE)")
    with pytest.raises(labels.LabelInputError):
        labels.esbl_decision_time(con, "cultures", "ast")


def test_rejects_null_positivity(con):
    con.execute("INSERT INTO cultures VALUES ('mystery', NULL, FALSE)")
    with pytest.raises(labels.LabelInputError):
        labels.esbl_decision_time(con, "cultures", "ast")


def test_accepts_relations(con):
    result = labels.esbl_decision_time(con, con.table("cultures"), con.table("ast"))
    assert collect(result) == expected(COMPOSITE_REASONS, labels.E2_RULE)


def test_rejects_unsafe_names(con):
    with pytest.raises(labels.LabelInputError):
        labels.esbl_decision_time(con, "cultures; DROP TABLE ast", "ast")
    with pytest.raises(labels.LabelInputError):
        labels.esbl_decision_time(con, "cultures", "ast", sir_column="sir) OR (1=1")
