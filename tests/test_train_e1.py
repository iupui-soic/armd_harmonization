"""analysis/02_train_e1.py on synthetic cohorts. No real data and no frozen plan are involved."""

import importlib.util
import json

import duckdb
import numpy as np
import pyarrow as pa
import pytest

from armd_common.config import REPO_ROOT
from armd_common.features import FEATURES_48
from armd_common.splits import SealedPartitionError

pytest.importorskip("xgboost")
pytest.importorskip("statsmodels")

spec = importlib.util.spec_from_file_location("train_e1", REPO_ROOT / "analysis" / "02_train_e1.py")
train = importlib.util.module_from_spec(spec)
spec.loader.exec_module(train)

FAST = {"n_estimators": 20}


def write_cohort(path, rows: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.register("t", pa.table(rows))
    con.execute(f"COPY t TO '{path}' (FORMAT parquet)")
    con.close()


def synthetic_cohorts(root, n=4000, seed=0, test_rows=False):
    """E2 is every unit. E1 keeps the units that grew Enterobacterales, which holds every positive."""
    rng = np.random.default_rng(seed)
    patients = rng.integers(0, n // 2, size=n)
    X = rng.normal(size=(n, len(FEATURES_48)))
    grew = rng.random(n) < 0.25 + 0.2 * (X[:, 0] > 1)
    risk = 1 / (1 + np.exp(-(-1.2 + 1.0 * X[:, 1] + 0.8 * X[:, 2])))
    y = (grew & (rng.random(n) < risk)).astype(int)
    base = {
        "unit_id": [f"u{i:05d}" for i in range(n)],
        "patient_id": [f"site:p{p}" for p in patients],
        "partition": ["derivation"] * n,
        "fold": (patients % 5 + 1).astype(int),
        "label": y,
        **{c: X[:, j].astype(float) for j, c in enumerate(FEATURES_48)},
    }
    if test_rows:
        base["partition"][0] = "test"
    names = train.cohort_names("release", "culture")
    write_cohort(root / names["E2"] / "derivation.parquet", base)
    keep = np.flatnonzero(grew)
    write_cohort(root / names["E1"] / "derivation.parquet", {k: np.asarray(v)[keep] for k, v in base.items()})
    return y, grew


def run(tmp_path, **kw):
    return train.train_site("site", tmp_path / "cohorts", tmp_path / "models", "release", "culture",
                            ("xgboost", "logistic"), seed=20260911, threads=2, params=FAST, **kw)


def test_trains_both_estimands_and_recalibrates_e1(tmp_path):
    y, grew = synthetic_cohorts(tmp_path / "cohorts")
    out = run(tmp_path)
    assert set(out["models"]) == {"xgboost_E1", "xgboost_E2", "logistic_E1", "logistic_E2"}
    assert out["cohorts"]["E2"]["units"] == len(y) and out["cohorts"]["E1"]["units"] == int(grew.sum())
    for kind in ("xgboost", "logistic"):
        e1, e2 = out["models"][f"{kind}_E1"], out["models"][f"{kind}_E2"]
        assert "recalibration_to_E2" in e1 and "recalibration_to_E2" not in e2
        recal = e1["recalibration_to_E2"]
        # E2 adds many units with Y = 0, so the E1 model over-predicts there. The Platt map is a
        # maximum likelihood logistic fit, so it restores the observed number of events exactly.
        names = train.cohort_names("release", "culture")
        path = tmp_path / "models" / names["E1"] / f"{kind}_crossfit_on_{names['E2']}.parquet"
        on_e2 = duckdb.connect().execute(f"SELECT prediction, label FROM read_parquet('{path}')").fetchnumpy()
        assert on_e2["prediction"].mean() > on_e2["label"].mean() + 0.05
        platt = train.modeling.PlattMap(**recal["platt"])
        assert platt(on_e2["prediction"]).mean() == pytest.approx(on_e2["label"].mean(), rel=1e-4)
        assert 0 < e2["threshold_90_sensitivity"] < 1 and 0 < recal["threshold_90_sensitivity_platt"] < 1
        # A monotone map keeps the ranking, so the raw threshold maps onto the Platt threshold.
        assert recal["threshold_90_sensitivity_platt"] == pytest.approx(
            float(train.modeling.PlattMap(**recal["platt"])(np.array([recal["threshold_90_sensitivity_uncalibrated"]]))[0]))


def test_predictions_are_cross_fitted_and_complete(tmp_path):
    y, _ = synthetic_cohorts(tmp_path / "cohorts")
    run(tmp_path)
    names = train.cohort_names("release", "culture")
    folder = tmp_path / "models" / names["E2"]
    con = duckdb.connect()
    n, distinct, bad = con.execute(
        f"SELECT COUNT(*), COUNT(DISTINCT unit_id), COUNT(*) FILTER (WHERE prediction <= 0 OR prediction >= 1 OR prediction IS NULL) "
        f"FROM read_parquet('{folder / 'xgboost_crossfit.parquet'}')").fetchone()
    assert (n, distinct, bad) == (len(y), len(y), 0)
    assert (folder / "xgboost.ubj").exists() and (folder / "logistic.pkl").exists()
    e1_folder = tmp_path / "models" / names["E1"]
    assert (e1_folder / f"xgboost_crossfit_on_{names['E2']}.parquet").exists()
    assert (e1_folder / "xgboost_isotonic_to_E2.pkl").exists()


def test_training_is_deterministic(tmp_path):
    synthetic_cohorts(tmp_path / "cohorts")
    first, second = run(tmp_path), run(tmp_path)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)


def test_refuses_any_partition_but_derivation(tmp_path):
    synthetic_cohorts(tmp_path / "cohorts", test_rows=True)
    with pytest.raises(SealedPartitionError):
        run(tmp_path)


def test_refuses_to_train_before_the_plan_is_frozen(monkeypatch):
    monkeypatch.setattr(train, "FROZEN_TAG", "sap-never-tagged")
    with pytest.raises(SystemExit):
        train.require_frozen_plan()
