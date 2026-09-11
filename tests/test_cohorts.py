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
    write_table(tmp_path, "armd-stanford", "prior_med", f"{IDS}, medication_name VARCHAR, medication_time_to_culturetime DOUBLE, medication_category VARCHAR", [
        ("p1", "e1", "101", "Ceftriaxone", 0.0, "CEF"),
        ("p3", "e6", "106", "Fluconazole", 1.0, "FLU"),
        ("p3", "e5", "105", "Cefepime In", 10.0, "CEF"),
    ])
    # A synthetic rule that exercises implied susceptibility plumbing.
    write_table(tmp_path, "armd-stanford", "implied_susceptibility_rules", "Organism VARCHAR, Antibiotic VARCHAR, Rule VARCHAR", [
        ("KLEBSIELLA", "Cefepime", "Resistant"),
    ])
    manifest = tmp_path / "stanford.csv.gz"
    splits.write_manifest("stanford", [("stanford:p1", "derivation", 1), ("stanford:p2", "test", None), ("stanford:p3", "derivation", 2)], path=manifest)
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
    assert result.flow[-1]["label_positive"] == 1


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
    # Ceftriaxone on the order day removes culture 101. Fluconazole is an antifungal and does not count.
    assert set(rows(loader, result)) == {"stanford:p1:102", "stanford:p3:106"}
    assert steps(result)[-1] == ("no antibiotics in the 3 days before the order", 2)


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
    ],
)
def test_invalid_specs(spec, site):
    with pytest.raises(ValueError):
        spec.validate(site)


def test_spec_names():
    assert CohortSpec("E2", nonsusceptible=True, exclude_recent_antibiotics=True).name == "e2_composite_nonsusceptible_no_antibiotics_3d"
    assert CohortSpec("E1", interpretation="clsi2022").name == "e1_composite_clsi2022"
    assert CohortSpec("E1", label="idsa").name == "e1_idsa"
