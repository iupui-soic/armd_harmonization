import duckdb

from armd_common import qc


def test_cell_and_rate_suppress_small_counts():
    assert [qc.cell(n) for n in (0, 1, 10, 11, None)] == [0, "<11", "<11", 11, None]
    assert qc.rate(0, 500) == 0.0
    assert qc.rate(5, 500) == "<11"
    assert qc.rate(50, 5) == "<11"
    assert qc.rate(50, 200) == 0.25
    assert qc.rate(3, 0) is None


def test_value_counts_and_null_shares():
    con = duckdb.connect()
    con.execute("CREATE TABLE t AS SELECT CASE WHEN i < 3 THEN 'rare' WHEN i < 50 THEN 'common' ELSE NULL END AS v, i FROM range(60) r(i)")
    assert qc.value_counts(con, "t", "v") == [["common", 47], [None, "<11"], ["rare", "<11"]]
    assert qc.null_shares(con, "t") == {"v": "<11", "i": 0.0}


def test_sentinel_reporting():
    con = duckdb.connect()
    n = 100
    con.execute("CREATE TABLE episodes AS SELECT 'c' || i AS culture_id, i < ? AS has_enterobacterales FROM range(?) r(i)", [n, n + 5])
    rows = []
    for i in range(n):
        if i < 50:
            rows.append((f"c{i}", True, "CRO", "S"))
        elif i < 80:
            rows.append((f"c{i}", True, "CRO", None))
        if i < 5:
            rows.append((f"c{i}", True, "CAZ", "S"))
        rows.append((f"c{i}", True, "TZP", "R" if i % 2 else "S"))
    rows.append(("c101", False, "CRO", "R"))
    con.execute("CREATE TABLE ast (culture_id VARCHAR, is_enterobacterales BOOLEAN, antibiotic_code VARCHAR, sir VARCHAR)")
    con.executemany("INSERT INTO ast VALUES (?, ?, ?, ?)", rows)
    out = qc.sentinel_reporting(con, {"culture_episodes": "episodes", "isolate_ast": "ast"})
    assert out["CRO"] == {"reported_sir": 0.5, "listed_without_sir": 0.3, "not_listed": 0.2}
    assert out["TZP"] == {"reported_sir": 1.0, "listed_without_sir": 0.0, "not_listed": 0.0}
    # Five cultures with CAZ is a small cell, so its share is suppressed.
    assert out["CAZ"] == {"reported_sir": "<11", "listed_without_sir": 0.0, "not_listed": 0.95}
    assert out["FEP"]["not_listed"] == 1.0
    assert out["any_sentinel_reported_sir"] == 1.0
    assert out["enterobacterales_cultures"] == 100


def test_vocab_gaps_lists_values_without_counts():
    con = duckdb.connect()
    con.execute("CREATE TEMP TABLE x_lookup_antibiotic AS SELECT * FROM (VALUES ('Ceftriaxone', 'CRO'), ('Mystery Drug', NULL), ('Gent screen', NULL)) t(raw, value)")
    con.execute("CREATE TEMP TABLE x_lookup_age AS SELECT * FROM (VALUES ('18-24 years', '18-24')) t(raw, value)")
    con.execute("CREATE TEMP TABLE x_lookup_organism AS SELECT * FROM (VALUES ('GRAM NEGATIVE RODS', TRUE), ('ESCHERICHIA COLI', FALSE)) t(organism_raw, organism_ambiguous)")
    con.execute("CREATE TABLE ast AS SELECT * FROM (VALUES ('CRO', 'S', 'Susceptible', NULL, NULL), ('ZZZ', 'other', 'See Comment', 'other', 'Weird')) "
                "t(antibiotic_code, sir, sir_raw, sir_clsi2022, sir_clsi2022_raw)")
    con.execute("CREATE TABLE episodes AS SELECT * FROM (VALUES ('urine', 'URINE'), ('other', 'WOUND')) t(specimen_type, specimen_raw)")
    gaps = qc.vocab_gaps(con, "x", {"isolate_ast": "ast", "culture_episodes": "episodes"})
    assert gaps == {
        "antibiotic_names_without_code": ["Mystery Drug"],
        "antibiotic_names_unmapped_by_design": ["Gent screen"],
        "antibiotic_codes_outside_vocab": ["ZZZ"],
        "susceptibility_values_mapped_to_other": ["See Comment", "Weird"],
        "specimen_values_mapped_to_other": ["WOUND"],
        "age_values_unmapped": [],
        "sex_values_unmapped": None,
        "ambiguous_organism_names": ["GRAM NEGATIVE RODS"],
    }


def test_publish_flows_withholds_counts_exposed_by_subtraction():
    primary = [
        {"step": "ordered cultures", "units": None, "cultures": 5000, "patients": 900, "label_positive": None},
        {"step": "evaluable", "units": 400, "cultures": 400, "patients": 300, "label_positive": None},
        {"step": "in partition", "units": 320, "cultures": 320, "patients": 240, "label_positive": 40},
        # Excludes three units, so publishing it would expose a cell of three.
        {"step": "no recent antibiotics", "units": 317, "cultures": 317, "patients": 238, "label_positive": 40},
    ]
    sensitivity = [
        {"step": "ordered cultures", "units": None, "cultures": 5000, "patients": 900, "label_positive": None},
        {"step": "evaluable", "units": 401, "cultures": 401, "patients": 301, "label_positive": None},
        {"step": "in partition", "units": 250, "cultures": 260, "patients": 200, "label_positive": 25},
    ]
    out = qc.publish_flows({"primary": primary, "sensitivity": sensitivity})
    assert out["primary"][2] == {"step": "in partition", "units": 320, "cultures": 320, "patients": 240, "label_positive": 40, "prevalence": 0.125}
    assert out["primary"][3]["units"] == qc.WITHHELD and out["primary"][3]["prevalence"] == qc.WITHHELD
    # One more evaluable culture than the primary cohort is withheld. Identical counts are not.
    assert out["sensitivity"][0]["cultures"] == 5000
    assert out["sensitivity"][1]["cultures"] == qc.WITHHELD
    assert out["sensitivity"][2]["prevalence"] == 0.1


def test_repeat_profile():
    con = duckdb.connect()
    rows = []
    for p in range(40):
        rows += [(f"c{p}a", f"p{p}", "2020-01-01 08:00:00", "blood", False, False),
                 (f"c{p}b", f"p{p}", "2020-01-01 08:00:00", "blood", False, False),
                 (f"c{p}c", f"p{p}", "2020-01-01 08:00:00", "urine", True, True),
                 (f"c{p}d", f"p{p}", "2020-01-04 08:00:00", "urine", False, False),
                 (f"c{p}e", f"p{p}", "2020-01-20 08:00:00", "urine", False, False)]
    con.execute("CREATE TABLE e (culture_id VARCHAR, patient_id VARCHAR, order_time TIMESTAMP, specimen_type VARCHAR, "
                "culture_positive BOOLEAN, has_enterobacterales BOOLEAN)")
    con.executemany("INSERT INTO e VALUES (?, ?, ?, ?, ?, ?)", rows)
    out = qc.repeat_profile(con, "e")
    assert out["gap_days"][3] == 40 and out["gap_days"][16] == 40 and out["gap_days"][1] == 0
    assert out["share_with_culture_1_to_13_days_earlier"] == 0.2
    assert out["share_sharing_order_date"] == 0.6
    assert out["blood_share_in_same_day_sets"] == 1.0
    assert out["by_repeat_rule"]["release"]["order_events"] == 120
    assert out["by_repeat_rule"]["chained_14d"] == {
        "cultures": 160, "specimen_days": 120, "order_events": 80, "order_events_with_several_specimen_types": 40,
        "patients": 40, "positive_share": 0.25, "blood_share": 0.5, "enterobacterales_cultures": 40,
    }
