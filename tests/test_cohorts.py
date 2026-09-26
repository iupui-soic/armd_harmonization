"""E1 and E2 cohorts on a synthetic Stanford-layout release with a hand-written split."""

import pytest

from armd_common import cohorts, splits
from armd_common.cohorts import CohortSpec
from armd_common.db import connect
from armd_common.sites import StanfordLoader
from test_features import build_stanford_feature_release
from test_sites import IDS, write_table


@pytest.fixture()
def stanford(tmp_path, monkeypatch):
    monkeypatch.setenv("ARMD_DATA_ROOT", str(tmp_path))
    build_stanford_feature_release(tmp_path)
    # Recent antibiotic records, added as a second file of the exposure table.
    write_table(tmp_path, "armd-stanford", "antibiotic_subtype_exposure", f"{IDS}, medication_name VARCHAR, medication_time_to_cultureTime DOUBLE", [
        ("p1", "e1", "101", "Ceftriaxone", 0.0),
        ("p1", "e2", "102", "Ceftriaxone", 2.0),
        ("p3", "e6", "106", "Fluconazole", 1.0),
    ], part=1)
    # A synthetic rule that exercises implied susceptibility plumbing.
    write_table(tmp_path, "armd-stanford", "implied_susceptibility_rules", "Organism VARCHAR, Antibiotic VARCHAR, Rule VARCHAR", [
        ("KLEBSIELLA", "Cefepime", "Resistant"),
    ])
    manifest = tmp_path / "stanford.csv.gz"
    splits.write_manifest("stanford", [("stanford:p1", "derivation", 1), ("stanford:p2", "test", None), ("stanford:p3", "derivation", 2)], path=manifest)
    loader = StanfordLoader(connect())
    return loader, loader.build(), manifest


@pytest.fixture()
def stanford_repeats(stanford, tmp_path):
    """The same release plus a same-day blood culture set and repeat cultures, added as second files."""
    loader, _, manifest = stanford
    write_table(tmp_path, "armd-stanford", "cohort", (
        f"{IDS}, order_time_jittered_utc TIMESTAMPTZ, ordering_mode VARCHAR, culture_description VARCHAR, was_positive BIGINT, "
        "organism VARCHAR, antibiotic VARCHAR, susceptibility VARCHAR"), [
        # Two blood cultures ordered with urine culture 106.
        ("p3", "e6", "107", "2021-04-04 12:00:00+00", "Inpatient", "BLOOD", 0, None, None, None),
        ("p3", "e6", "108", "2021-04-04 12:00:00+00", "Inpatient", "BLOOD", 0, None, None, None),
        # Six days after 106: a repeat.
        ("p3", "e7", "109", "2021-04-10 12:00:00+00", "Outpatient", "URINE", 1, "ESCHERICHIA COLI", "Ceftriaxone", "Resistant"),
        # Nine days after 102, then ten days after that. Both are repeats because the window chains.
        ("p1", "e8", "110", "2020-02-10 08:30:00+00", "Outpatient", "URINE", 0, None, None, None),
        ("p1", "e9", "111", "2020-02-20 08:30:00+00", "Outpatient", "URINE", 0, None, None, None),
    ], part=1)
    extra = [("p3", "e6", "107"), ("p3", "e6", "108"), ("p3", "e7", "109"), ("p1", "e8", "110"), ("p1", "e9", "111")]
    write_table(tmp_path, "armd-stanford", "ward_info", f"{IDS}, hosp_ward_IP BIGINT, hosp_ward_OP BIGINT, hosp_ward_ER BIGINT, hosp_ward_ICU BIGINT",
                [ids + ((1, 0, 0, 0) if ids[2] in ("107", "108") else (0, 1, 0, 0)) for ids in extra], part=1)
    write_table(tmp_path, "armd-stanford", "demographics", f"{IDS}, age VARCHAR, gender BIGINT",
                [ids + (("18-24 years", 0) if ids[0] == "p1" else ("55-64 years", None)) for ids in extra], part=1)
    write_table(tmp_path, "armd-stanford", "adi_scores", f"{IDS}, adi_score BIGINT, adi_state_rank BIGINT",
                [ids + ((5, 1) if ids[0] == "p1" else (88, 9)) for ids in extra], part=1)
    loader = StanfordLoader(connect())
    return loader, loader.build(), manifest


def rows(loader, result, columns="label, reason, partition"):
    return {r[0]: r[1:] for r in loader.con.execute(f"SELECT culture_id, {columns} FROM {result.table} ORDER BY 1").fetchall()}


def steps(result):
    return [(s["step"], s["cultures"]) for s in result.flow]


def test_e1_primary(stanford):
    loader, built, manifest = stanford
    result = cohorts.build_cohort(loader, built, CohortSpec("E1"), manifest=manifest)
    assert rows(loader, result) == {
        "stanford:p1:102": (1, "resistant", "derivation"),
        "stanford:p3:106": (0, "not_resistant", "derivation"),
    }
    assert steps(result) == [("ordered cultures", 6), ("evaluable for E1", 2), ("in partition derivation", 2)]


def test_e2_primary_excludes_sealed_test_patients(stanford):
    loader, built, manifest = stanford
    result = cohorts.build_cohort(loader, built, CohortSpec("E2"), manifest=manifest)
    assert rows(loader, result) == {
        "stanford:p1:101": (0, "culture_negative", "derivation"),
        "stanford:p1:102": (1, "resistant", "derivation"),
        "stanford:p3:106": (0, "not_resistant", "derivation"),
    }
    assert steps(result) == [("ordered cultures", 6), ("evaluable for E2", 4), ("in partition derivation", 3)]
    # Label counts start at the partition step, so the test patient's cultures are never counted.
    assert [s["label_positive"] for s in result.flow] == [None, None, 1]


def test_subgroups_and_features(stanford):
    loader, built, manifest = stanford
    result = cohorts.build_cohort(loader, built, CohortSpec("E2"), manifest=manifest)
    got = rows(loader, result, "prior_culture_history, adi_quintile, setting_primary, specimen_type, sex, T_ceph3_90d")
    assert got["stanford:p1:101"] == (0, "q1", "outpatient", "urine", "f", 0)
    assert got["stanford:p1:102"] == (1, "q1", "inpatient", "urine", "f", 1)
    assert got["stanford:p3:106"] == (1, "q4", "inpatient", "urine", "unknown", 0)


def test_nonsusceptible_counts_intermediate(stanford):
    loader, built, manifest = stanford
    result = cohorts.build_cohort(loader, built, CohortSpec("E2", nonsusceptible=True), manifest=manifest)
    assert rows(loader, result, "label")["stanford:p3:106"] == (1,)


def test_recent_antibiotic_exclusion(stanford):
    loader, built, manifest = stanford
    result = cohorts.build_cohort(loader, built, CohortSpec("E2", exclude_recent_antibiotics=True), manifest=manifest)
    # Ceftriaxone two days before the order removes culture 102. Culture 101 stays, because a day 0
    # record cannot be placed before the order. Fluconazole is an antifungal and does not count.
    assert set(rows(loader, result)) == {"stanford:p1:101", "stanford:p3:106"}
    assert steps(result)[-1] == ("no antibacterial record 1 to 3 days before the order", 2)


def test_implied_susceptibility_strict(stanford):
    loader, built, manifest = stanford
    result = cohorts.build_cohort(loader, built, CohortSpec("E1", implied_susceptibility="strict"), manifest=manifest)
    # The synthetic intrinsic rule imputes cefepime R for the Klebsiella isolate that had no sentinel result.
    assert rows(loader, result, "label, reason") == {
        "stanford:p1:102": (1, "resistant"),
        "stanford:p3:105": (1, "resistant"),
        "stanford:p3:106": (0, "not_resistant"),
    }


def test_test_partition_is_sealed(stanford):
    loader, built, manifest = stanford
    with pytest.raises(splits.SealedPartitionError):
        cohorts.build_cohort(loader, built, CohortSpec("E2"), partitions=("derivation", "test"), manifest=manifest)


@pytest.mark.parametrize(
    "spec, site",
    [
        (CohortSpec("E2", label="idsa"), "stanford"),
        (CohortSpec("E1", interpretation="clsi2022"), "stanford"),
        (CohortSpec("E1", label="enzyme"), "stanford"),
        (CohortSpec("E1", label="idsa", nonsusceptible=True), "mgb"),
        (CohortSpec("E1", interpretation="clsi2022", implied_susceptibility="strict"), "mgb"),
        (CohortSpec("E3"), "stanford"),
        (CohortSpec("E2", repeat_rule="anchored_14d"), "stanford"),
        (CohortSpec("E2", unit="encounter"), "stanford"),
    ],
)
def test_invalid_specs(spec, site):
    with pytest.raises(ValueError):
        spec.validate(site)


def test_release_culture_keeps_repeats(stanford_repeats):
    loader, built, manifest = stanford_repeats
    result = cohorts.build_cohort(loader, built, CohortSpec("E2"), manifest=manifest)
    assert set(rows(loader, result)) == {f"stanford:p{p}:{c}" for p, c in [
        (1, 101), (1, 102), (1, 110), (1, 111), (3, 106), (3, 107), (3, 108), (3, 109)]}


def test_chained_14d_specimen_day(stanford_repeats):
    loader, built, manifest = stanford_repeats
    result = cohorts.build_cohort(loader, built, CohortSpec("E2", repeat_rule="chained_14d", unit="specimen_day"), manifest=manifest)
    got = {r[0]: r[1:] for r in loader.con.execute(
        f"SELECT unit_id, culture_id, n_cultures, label, reason, specimen_type FROM {result.table}").fetchall()}
    assert got == {
        "stanford:p1:2020-01-01:urine": ("stanford:p1:101", 1, 0, "culture_negative", "urine"),
        "stanford:p1:2020-02-01:urine": ("stanford:p1:102", 1, 1, "resistant", "urine"),
        "stanford:p3:2021-04-04:urine": ("stanford:p3:106", 1, 0, "not_resistant", "urine"),
        # The blood culture set counts once.
        "stanford:p3:2021-04-04:blood": ("stanford:p3:107", 2, 0, "culture_negative", "blood"),
    }
    assert [(s["step"], s["units"], s["cultures"], s["label_positive"]) for s in result.flow] == [
        ("ordered cultures", None, 11, None),
        ("index orders under chained_14d", 7, 8, None),
        ("evaluable for E2", 5, 6, None),
        ("in partition derivation", 4, 5, 1),
    ]


def test_order_event_pools_specimens(stanford_repeats):
    loader, built, manifest = stanford_repeats
    result = cohorts.build_cohort(loader, built, CohortSpec("E2", repeat_rule="chained_14d", unit="order_event"), manifest=manifest)
    got = loader.con.execute(
        f"SELECT n_cultures, label, reason, specimen_type, specimen_urine, specimen_blood, setting_primary FROM {result.table} "
        "WHERE unit_id = 'stanford:p3:2021-04-04'").fetchall()
    assert got == [(3, 0, "not_resistant", "multiple", 1, 1, "inpatient")]


def test_e1_units_ignore_negative_cultures(stanford_repeats):
    loader, built, manifest = stanford_repeats
    culture = cohorts.build_cohort(loader, built, CohortSpec("E1"), manifest=manifest)
    unit = cohorts.build_cohort(loader, built, CohortSpec("E1", repeat_rule="chained_14d", unit="order_event"), manifest=manifest)
    # As released, repeat culture 109 enters E1. The rule drops it, and the negative blood set does
    # not stop the 2021-04-04 order from entering E1.
    assert set(rows(loader, culture)) == {"stanford:p1:102", "stanford:p3:106", "stanford:p3:109"}
    assert {r[0]: r[1:] for r in loader.con.execute(f"SELECT unit_id, n_cultures, label FROM {unit.table}").fetchall()} == {
        "stanford:p1:2020-02-01": (1, 1),
        "stanford:p3:2021-04-04": (3, 0),
    }


def test_unit_feature_aggregation():
    con = connect()
    con.execute("CREATE TABLE m AS SELECT * FROM (VALUES ('u1', 0, 1, 0.0), ('u1', 1, 3, 30.0), ('u1', 0, 2, 200.0), ('u2', 0, 0, 0.0)) "
                "t(unit_id, prior_ESBL, num_prior_orgs, days_since_prior_org)")
    exprs = ", ".join(cohorts.unit_feature_sql(c) for c in ("prior_ESBL", "num_prior_orgs", "days_since_prior_org"))
    got = dict((r[0], r[1:]) for r in con.execute(f"SELECT unit_id, {exprs} FROM m GROUP BY unit_id").fetchall())
    # The most recent prior organism wins, and 0 still means none.
    assert got == {"u1": (1, 3, 30.0), "u2": (0, 0, 0.0)}


def test_spec_names():
    assert CohortSpec("E2", nonsusceptible=True, exclude_recent_antibiotics=True).name == "e2_composite_nonsusceptible_no_antibiotics_3d"
    assert CohortSpec("E1", interpretation="clsi2022").name == "e1_composite_clsi2022"
    assert CohortSpec("E1", label="idsa").name == "e1_idsa"
    assert CohortSpec("E2", repeat_rule="chained_14d", unit="specimen_day").name == "e2_composite_chained_14d_specimen_day"
