import sys
from collections import Counter

import duckdb
import pytest

from armd_common import splits


def fake_patients(n=2000):
    ids = [f"stanford:P{i:05d}" for i in range(n)]
    strata = [2 if i % 10 == 0 else (1 if i % 3 == 0 else 0) for i in range(n)]
    return ids, strata


def test_partition_sizes_and_folds():
    ids, strata = fake_patients()
    out = splits.assign(ids, strata, seed=1)
    partitions = Counter(p for _, p, _ in out)
    assert partitions["test"] == len(ids) // 5
    assert {f for _, p, f in out if p == "test"} == {None}
    assert set(Counter(f for _, p, f in out if p == "derivation")) == {1, 2, 3, 4, 5}


def test_stratification_is_balanced():
    ids, strata = fake_patients()
    by_id = dict(zip(ids, strata))
    out = splits.assign(ids, strata, seed=1)
    test_rate = sum(by_id[pid] == 2 for pid, p, _ in out if p == "test") / sum(p == "test" for _, p, _ in out)
    assert abs(test_rate - 0.1) < 0.01


def test_deterministic_and_order_independent():
    ids, strata = fake_patients()
    first = sorted(splits.assign(ids, strata, seed=20260911))
    reversed_input = sorted(splits.assign(ids[::-1], strata[::-1], seed=20260911))
    assert first == reversed_input
    assert first != sorted(splits.assign(ids, strata, seed=7))


def test_rejects_duplicate_patients():
    with pytest.raises(ValueError):
        splits.assign(["a", "a", "b", "c", "d"], [0, 0, 1, 1, 0])


def test_manifest_is_frozen_and_test_is_sealed(tmp_path, monkeypatch):
    ids, strata = fake_patients(100)
    path = tmp_path / "stanford.csv.gz"
    splits.write_manifest("stanford", splits.assign(ids, strata), path=path)
    with pytest.raises(FileExistsError):
        splits.write_manifest("stanford", splits.assign(ids, strata), path=path)

    visible = splits.read_manifest("stanford", path=path)
    assert len(visible) == 80 and all(r["partition"] == "derivation" for r in visible)

    monkeypatch.setattr(sys, "argv", ["02_train_e1.py"])
    with pytest.raises(splits.SealedPartitionError):
        splits.read_manifest("stanford", include_test=True, path=path)
    monkeypatch.setattr(sys, "argv", ["analysis/03_evaluate_ladder.py"])
    assert len(splits.read_manifest("stanford", include_test=True, path=path)) == 100


def test_hash_matches_duckdb(tmp_path):
    ids, strata = fake_patients(50)
    path = tmp_path / "stanford.csv.gz"
    splits.write_manifest("stanford", splits.assign(ids, strata), path=path)
    con = duckdb.connect()
    con.execute("CREATE TABLE episodes AS SELECT unnest(?) AS patient_id", [ids])
    table = splits.register_manifest(con, "stanford", path=path)
    joined = con.execute(f"SELECT COUNT(*) FROM episodes e JOIN {table} s ON s.patient_id_hash = sha256(e.patient_id)").fetchone()[0]
    assert joined == 40
