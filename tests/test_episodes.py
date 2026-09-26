"""Repeat rules and unit-level labels on hand-built cultures."""

import duckdb
import pytest

from armd_common import episodes, labels

# culture_id, patient_id, order_time, specimen_type
REPEATS = [
    ("a1", "a", "2020-01-01 09:00:00", "urine"),
    ("a2", "a", "2020-01-01 09:00:00", "blood"),
    ("a3", "a", "2020-01-01 09:00:00", "blood"),
    # Four days after the index: a repeat.
    ("a4", "a", "2020-01-05 09:00:00", "urine"),
    # Fourteen days after the index but ten after a4. The window chains, so this is still a repeat.
    ("a5", "a", "2020-01-15 09:00:00", "urine"),
    # Exactly fourteen days after a5: a new index.
    ("a6", "a", "2020-01-29 09:00:00", "urine"),
    ("b1", "b", "2020-03-01 09:00:00", "blood"),
    ("b2", "b", "2020-03-14 09:00:00", "blood"),
]

# culture_id, patient_id, specimen_type, reason. Every culture is ordered on the same day.
UNITS = [
    # A blood culture set with one resistant bottle.
    ("x1", "p1", "blood", "resistant"),
    ("x2", "p1", "blood", "culture_negative"),
    # Urine not resistant, with a negative blood culture set the same day.
    ("x3", "p2", "urine", "not_resistant"),
    ("x4", "p2", "blood", "culture_negative"),
    ("x5", "p2", "blood", "culture_negative"),
    # One urine not resistant and one urine with no sentinel result.
    ("x6", "p3", "urine", "not_resistant"),
    ("x7", "p3", "urine", "not_tested"),
    # A resistant culture outweighs a culture with no sentinel result.
    ("x8", "p4", "urine", "not_tested"),
    ("x9", "p4", "urine", "resistant"),
    # An unidentified organism next to a negative culture.
    ("x10", "p5", "blood", "organism_unidentified"),
    ("x11", "p5", "blood", "culture_negative"),
    # No Enterobacterales anywhere in the order.
    ("x12", "p6", "respiratory", "organism_absent"),
    ("x13", "p6", "urine", "culture_negative"),
]


@pytest.fixture()
def con():
    con = duckdb.connect()
    con.execute("CREATE TABLE repeats (culture_id VARCHAR, patient_id VARCHAR, order_time TIMESTAMP, specimen_type VARCHAR)")
    con.executemany("INSERT INTO repeats VALUES (?, ?, ?, ?)", REPEATS)
    con.execute("CREATE TABLE unit_cultures (culture_id VARCHAR, patient_id VARCHAR, order_time TIMESTAMP, specimen_type VARCHAR, reason VARCHAR)")
    con.executemany("INSERT INTO unit_cultures VALUES (?, ?, TIMESTAMP '2021-06-01 10:00:00', ?, ?)", UNITS)
    for estimand, rule in (("E1", labels.E1_RULE), ("E2", labels.E2_RULE)):
        cases = " ".join(f"WHEN '{r}' THEN {v}" for r, v in rule.items())
        evaluable = ", ".join(f"'{r}'" for r in rule)
        con.execute(f"""
            CREATE TABLE labels_{estimand} AS
            SELECT culture_id, CAST(CASE reason {cases} END AS TINYINT) AS label,
                   reason IN ({evaluable}) AS evaluable, reason
            FROM unit_cultures
        """)
    episodes.kept_cultures(con, "unit_cultures", "release", "all_units")
    return con


def kept(con, rule):
    episodes.kept_cultures(con, "repeats", rule, "kept")
    return sorted(r[0] for r in con.execute("SELECT culture_id FROM kept").fetchall())


def test_release_keeps_everything(con):
    assert kept(con, "release") == sorted(r[0] for r in REPEATS)


def test_chained_14d_window(con):
    # The whole index order is kept, including the blood culture set.
    assert kept(con, "chained_14d") == ["a1", "a2", "a3", "a6", "b1"]


def test_first_per_patient(con):
    assert kept(con, "first_per_patient") == ["a1", "a2", "a3", "b1"]


def test_repeat_rule_needs_order_times(con):
    con.execute("INSERT INTO repeats VALUES ('c1', 'c', NULL, 'urine')")
    with pytest.raises(ValueError):
        episodes.kept_cultures(con, "repeats", "chained_14d", "kept")
    assert "c1" in kept(con, "release")


def test_unknown_options():
    with pytest.raises(ValueError):
        episodes.unit_key_sql("encounter")
    with pytest.raises(ValueError):
        episodes.kept_cultures(duckdb.connect(), "repeats", "anchored_14d", "kept")


def units(con, unit, estimand):
    episodes.unit_labels(con, f"labels_{estimand}", "unit_cultures", "all_units", unit, estimand, "units")
    return {r[0]: r[1:] for r in con.execute("SELECT unit_id, n_cultures, label, evaluable, reason FROM units").fetchall()}


@pytest.mark.parametrize("estimand", ["E1", "E2"])
def test_culture_unit_is_identity(con, estimand):
    got = units(con, "culture", estimand)
    expected = con.execute(f"SELECT culture_id, 1, label, evaluable, reason FROM labels_{estimand}").fetchall()
    assert got == {r[0]: r[1:] for r in expected}


def test_specimen_day_e1(con):
    assert units(con, "specimen_day", "E1") == {
        "p1:2021-06-01:blood": (2, 1, True, "resistant"),
        "p2:2021-06-01:urine": (1, 0, True, "not_resistant"),
        # Negative cultures are not E1 cultures, so the unit leaves E1.
        "p2:2021-06-01:blood": (2, None, False, "culture_negative"),
        # One culture had no sentinel result, so the unit cannot be called not resistant.
        "p3:2021-06-01:urine": (2, None, False, "not_tested"),
        "p4:2021-06-01:urine": (2, 1, True, "resistant"),
        "p5:2021-06-01:blood": (2, None, False, "organism_unidentified"),
        "p6:2021-06-01:respiratory": (1, None, False, "organism_absent"),
        "p6:2021-06-01:urine": (1, None, False, "culture_negative"),
    }


def test_specimen_day_e2(con):
    assert units(con, "specimen_day", "E2") == {
        "p1:2021-06-01:blood": (2, 1, True, "resistant"),
        "p2:2021-06-01:urine": (1, 0, True, "not_resistant"),
        "p2:2021-06-01:blood": (2, 0, True, "culture_negative"),
        "p3:2021-06-01:urine": (2, None, False, "not_tested"),
        "p4:2021-06-01:urine": (2, 1, True, "resistant"),
        # The unidentified organism might be a resistant Enterobacterales, so the unit is not evaluable.
        "p5:2021-06-01:blood": (2, None, False, "organism_unidentified"),
        "p6:2021-06-01:respiratory": (1, 0, True, "organism_absent"),
        "p6:2021-06-01:urine": (1, 0, True, "culture_negative"),
    }


def test_order_event(con):
    e1 = units(con, "order_event", "E1")
    e2 = units(con, "order_event", "E2")
    # The negative blood set leaves the E1 unit, and the urine result decides it.
    assert e1["p2:2021-06-01"] == (3, 0, True, "not_resistant")
    assert e2["p2:2021-06-01"] == (3, 0, True, "not_resistant")
    assert e1["p6:2021-06-01"] == (2, None, False, "organism_absent")
    assert e2["p6:2021-06-01"] == (2, 0, True, "organism_absent")
    assert e1["p1:2021-06-01"] == e2["p1:2021-06-01"] == (2, 1, True, "resistant")
    assert len(e1) == len(e2) == 6
