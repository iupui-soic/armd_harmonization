import duckdb
import pytest

from armd_common import implied
from armd_common.schema import ISOLATE_AST, ddl

RULES = implied.parse_rules([
    ("KLEBSIELLA", "Ampicillin", "Resistant"),
    ("ESCHERICHIA COLI", "Ceftriaxone", "Same as Cefotaxime if reported"),
    ("ESCHERICHIA COLI", "Cefepime", "Susceptible if susceptible to Ceftriaxone, Cefotaxime, or Cefazolin"),
    ("ESCHERICHIA COLI", "Doxycycline", "Susceptible if susceptible to Tetracycline, otherwise Resistant"),
    ("ESCHERICHIA COLI", "Cefpodoxime", "Same as Ceftriaxone if reported; otherwise, Susceptible if susceptible to Cefazolin"),
    ("KLEBSIELLA", "Ampicillin", "Resistant"),
])


def test_parse_rule_kinds():
    by_target = {r.target: r for r in RULES}
    assert len(RULES) == 5
    assert by_target["AMP"].kind == "intrinsic" and by_target["AMP"].organism_pattern == "KLEBSIELLA"
    assert (by_target["CRO"].kind, by_target["CRO"].same_as) == ("same_as", "CTX")
    assert (by_target["FEP"].kind, by_target["FEP"].sources) == ("susceptible_if", ("CRO", "CTX", "CFZ"))
    assert (by_target["DOX"].kind, by_target["DOX"].sources) == ("susceptible_if_else_r", ("TET",))
    assert (by_target["CPD"].kind, by_target["CPD"].same_as, by_target["CPD"].sources) == ("same_as_else_susceptible_if", "CRO", ("CFZ",))
    assert implied.parse_rule(0, "Proteus species", "Tigecycline", "Resistant").organism_pattern == "PROTEUS"


def test_parse_rule_rejects_unknown_text():
    with pytest.raises(ValueError):
        implied.parse_rule(0, "SERRATIA", "Cefazolin", "Ask the lab")
    with pytest.raises(ValueError):
        implied.parse_rule(0, "SERRATIA", "Unobtainium", "Resistant")


def isolate(culture, organism, code, sir):
    group = "klebsiella" if organism.startswith("KLEBSIELLA") else "ecoli"
    row = dict.fromkeys(ISOLATE_AST)
    row.update(site="stanford", patient_id="stanford:p", culture_id=culture, organism_raw=organism, organism_std=organism,
               organism_group=group, is_enterobacterales=True, idsa_esbl_species=True, organism_ambiguous=False,
               antibiotic_raw=code, antibiotic_code=code, sir_raw=sir, sir=sir, preliminary=False, n_source_rows=1)
    return tuple(row.values())


@pytest.fixture()
def con():
    con = duckdb.connect()
    con.execute(f"CREATE TABLE ast ({ddl(ISOLATE_AST)})")
    ecoli, kpn = "ESCHERICHIA COLI", "KLEBSIELLA PNEUMONIAE"
    con.executemany(f"INSERT INTO ast VALUES ({', '.join('?' for _ in ISOLATE_AST)})", [
        isolate("kpn_cro_s", kpn, "CRO", "S"),
        isolate("kpn_amp_reported", kpn, "AMP", "S"),
        isolate("ec_ctx_s", ecoli, "CTX", "S"),
        isolate("ec_ctx_s", ecoli, "CFZ", "S"),
        isolate("ec_tet_r", ecoli, "TET", "R"),
        isolate("ec_tet_s_cro_r", ecoli, "TET", "S"),
        isolate("ec_tet_s_cro_r", ecoli, "CRO", "R"),
        isolate("ec_tet_conflict", ecoli, "TET", "S"),
        isolate("ec_tet_conflict", ecoli, "TET", "R"),
    ])
    return con


def imputed(con, table):
    rows = con.execute(f"SELECT culture_id, antibiotic_code, sir FROM {table} WHERE imputed ORDER BY 1, 2").fetchall()
    reported = con.execute(f"SELECT COUNT(*) FROM {table} WHERE NOT imputed").fetchone()[0]
    assert reported == con.execute("SELECT COUNT(*) FROM ast").fetchone()[0]
    return rows


def test_strict_applies_intrinsic_and_equivalence_only(con):
    table = implied.apply_implied_susceptibility(con, "ast", RULES, strict=True)
    assert imputed(con, table) == [
        ("ec_ctx_s", "CRO", "S"),
        ("ec_tet_s_cro_r", "CPD", "R"),
        ("kpn_cro_s", "AMP", "R"),
    ]


def test_non_strict_applies_every_rule_without_chaining(con):
    table = implied.apply_implied_susceptibility(con, "ast", RULES, strict=False)
    assert imputed(con, table) == [
        ("ec_ctx_s", "CPD", "S"),
        ("ec_ctx_s", "CRO", "S"),
        ("ec_ctx_s", "FEP", "S"),
        ("ec_tet_conflict", "DOX", "R"),
        ("ec_tet_r", "DOX", "R"),
        ("ec_tet_s_cro_r", "CPD", "R"),
        ("ec_tet_s_cro_r", "DOX", "S"),
        ("kpn_cro_s", "AMP", "R"),
    ]
    rule = con.execute(f"SELECT imputation_rule FROM {table} WHERE culture_id = 'kpn_cro_s' AND imputed").fetchone()[0]
    assert rule == "Resistant"
