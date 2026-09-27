"""analysis/04_cascade_sensitivity.py on hand-built tables and the synthetic Stanford release."""

import importlib.util

import duckdb
import pytest

from armd_common.config import REPO_ROOT
from armd_common.qc import SMALL
from test_cohorts import stanford  # noqa: F401  (fixture)

spec = importlib.util.spec_from_file_location("cascade", REPO_ROOT / "analysis" / "04_cascade_sensitivity.py")
cascade = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cascade)

FEATURE_COLUMNS = ("hosp_ward_IP", "hosp_ward_OP", "hosp_ward_ER", "specimen_urine", "specimen_blood", "specimen_respiratory",
                   "num_prior_orgs", "multi_prior_resistant", "prior_ESBL")


def tables(con, rows):
    """rows: (culture_id, organism_group, sentinel sir or None, clsi sir or None, IP, prior_ESBL, enzyme result or None)"""
    con.execute("CREATE TABLE e (culture_id VARCHAR, culture_positive BOOLEAN, organism_unidentified BOOLEAN, organism_group VARCHAR)")
    con.execute("CREATE TABLE a (culture_id VARCHAR, is_enterobacterales BOOLEAN, idsa_esbl_species BOOLEAN, antibiotic_code VARCHAR, sir VARCHAR, sir_clsi2022 VARCHAR)")
    con.execute("CREATE TABLE z (culture_id VARCHAR, is_enterobacterales BOOLEAN, enzyme_class VARCHAR, enzyme_result VARCHAR)")
    con.execute(f"CREATE TABLE f (culture_id VARCHAR, {', '.join(c + ' INTEGER' for c in FEATURE_COLUMNS)})")
    for cid, group, sir, clsi, ip, esbl, enzyme in rows:
        con.execute("INSERT INTO e VALUES (?, TRUE, FALSE, ?)", [cid, group])
        # Without a lab result the sentinel agent is listed uninterpreted, as MGB does for ceftazidime.
        con.execute("INSERT INTO a VALUES (?, TRUE, TRUE, ?, ?, ?)", [cid, "CRO" if sir else "CAZ", sir, clsi])
        con.execute("INSERT INTO f VALUES (?, ?, ?, 0, 1, 0, 0, 0, 0, ?)", [cid, ip, 1 - ip, esbl])
        if enzyme:
            con.execute("INSERT INTO z VALUES (?, TRUE, 'ESBL', ?)", [cid, enzyme])


def test_not_tested_cultures_compared_with_evaluable():
    con = duckdb.connect()
    rows = [(f"ev{i}", "ecoli", "R" if i % 5 == 0 else "S", None, i % 2, int(i < 40), None) for i in range(150)]
    rows += [(f"nt{i}", "klebsiella", None, None, 1, int(i < 3), None) for i in range(50)]
    tables(con, rows)
    out = cascade.evaluable_comparison(con, "e", "a", "f")
    assert out["evaluable_cultures"] == 150 and out["not_tested_cultures"] == 50 and out["not_tested_share"] == 0.25
    ip = out["covariates"]["hosp_ward_IP"]
    assert ip["not_tested_mean"] == 1.0 and ip["smd_not_tested_minus_evaluable"] > 0
    assert out["covariates"]["organism_klebsiella"]["not_tested_mean"] == 1.0
    # Three not-tested cultures with prior ESBL is a small cell, so the row is withheld.
    assert out["covariates"]["prior_ESBL"] is None


def test_mgb_cross_check():
    con = duckdb.connect()
    rows = [(f"r{i}", "ecoli", "S", "R" if i < 20 else "S", 0, 0, "positive" if i < 15 else ("negative" if i < 60 else None)) for i in range(100)]
    rows += [(f"n{i}", "ecoli", None, "R" if i < 30 else "S", 0, 0, "positive" if i < 25 else ("negative" if i < 50 else None)) for i in range(80)]
    tables(con, rows)
    out = cascade.mgb_cross_check(con, "e", "a", "z")
    assert out["sentinel_reported"] == {"cultures": 100, "esbl_tested": 60, "esbl_positive_share": 0.25,
                                        "clsi2022_evaluable": 100, "clsi2022_resistant_share": 0.2}
    assert out["no_sentinel_result"]["esbl_positive_share"] == 0.5
    assert out["no_sentinel_result"]["clsi2022_resistant_share"] == pytest.approx(30 / 80, abs=1e-4)


def test_run_site_on_the_synthetic_release(stanford):
    loader, _, manifest = stanford
    out = cascade.run_site("stanford", loader=loader, manifest=manifest)
    assert set(out) == {"reporting", "not_tested_vs_evaluable", "implied_susceptibility"}
    # The synthetic release is tiny, so every count is a small cell.
    assert out["not_tested_vs_evaluable"]["evaluable_cultures"] == SMALL
    assert set(out["implied_susceptibility"]) == {"none", "strict", "all rules"}
