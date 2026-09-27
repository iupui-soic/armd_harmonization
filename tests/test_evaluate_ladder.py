"""analysis/03_evaluate_ladder.py on synthetic cohorts and models. No real data and no frozen plan."""

import importlib.util
import json
import sys

import numpy as np
import pytest

from armd_common import cohorts, metrics, splits
from armd_common.cohorts import CohortSpec
from armd_common.config import REPO_ROOT
from armd_common.features import FEATURES_48
from armd_common.qc import WITHHELD

pytest.importorskip("xgboost")
pytest.importorskip("statsmodels")


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "analysis" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


evaluate = load("evaluate_ladder", "03_evaluate_ladder.py")
train = load("train_e1", "02_train_e1.py")
from test_train_e1 import FAST, synthetic_cohorts  # noqa: E402
from test_cohorts import stanford  # noqa: E402,F401  (fixture)


def synthetic_test(n=3000, seed=5):
    """A test partition drawn like the training data: E2 is every unit, E1 those that grew Enterobacterales."""
    rng = np.random.default_rng(seed)
    patients = rng.integers(0, n // 2, size=n)
    X = rng.normal(size=(n, len(FEATURES_48))).astype(np.float32)
    grew = rng.random(n) < 0.25 + 0.2 * (X[:, 0] > 1)
    risk = 1 / (1 + np.exp(-(-1.2 + 1.0 * X[:, 1] + 0.8 * X[:, 2])))
    y = (grew & (rng.random(n) < risk)).astype(int)
    e2 = {"unit_id": np.array([f"t{i:05d}" for i in range(n)]), "patient_id": np.array([f"site:q{p}" for p in patients]),
          "y": y, "X": X,
          "setting_primary": np.where(X[:, 3] > 0, "inpatient", "outpatient"),
          "specimen_type": np.where(X[:, 4] > 0.5, "blood", "urine"),
          "prior_culture_history": (X[:, 5] > 0).astype(int),
          "adi_quintile": np.array(["q1", "q2", "q3", "q4", "q5"])[patients % 5]}
    e1 = {k: v[grew] for k, v in e2.items()}
    return {"E1": e1, "E2": e2}


@pytest.fixture(scope="module")
def trained(tmp_path_factory):
    root = tmp_path_factory.mktemp("site")
    synthetic_cohorts(root / "cohorts")
    summary = train.train_site("site", root / "cohorts", root / "models", "release", "culture", ("xgboost", "logistic"),
                               seed=20260911, threads=2, params=FAST)
    names = train.cohort_names("release", "culture")
    variants = {k: evaluate.load_variants(root / "models", names, k, summary) for k in ("xgboost", "logistic")}
    return variants


def test_grid_cells_and_paired_comparisons(trained):
    out = evaluate.evaluate_site(synthetic_test(), trained, n_boot=30, seed=1)
    assert len(out["grid"]) == 12
    cell = out["grid"]["xgboost:E1->E2"]
    lo, hi = cell["interval"]["auroc"]
    assert lo <= cell["point"]["auroc"] <= hi
    # The E1 model over-predicts on E2 and Platt corrects calibration in the large.
    assert cell["point"]["calibration_intercept"] < -0.3
    assert abs(out["grid"]["xgboost:E1_platt->E2"]["point"]["calibration_intercept"]) < 0.3
    assert set(out["comparisons"]) == {f"{k}:{v}->E2 minus {k}:E2->E2" for k in ("xgboost", "logistic")
                                       for v in ("E1", "E1_platt", "E1_isotonic")}
    assert set(cell["curves"]) == {"decision_curve", "treat_all", "calibration_curve"}
    assert set(cell["subgroups"]) == set(evaluate.SUBGROUPS)


def test_point_estimates_match_the_metrics_module(trained):
    test = synthetic_test()
    out = evaluate.evaluate_site(test, trained, n_boot=5, seed=1)
    fn, threshold = trained["logistic"]["E2"]
    expected = metrics.evaluate(fn(test["E2"]["X"]), test["E2"]["y"], threshold)
    got = out["grid"]["logistic:E2->E2"]["point"]
    for k in ("auroc", "auprc", "brier", "calibration_slope", "net_benefit_0.10"):
        assert got[k] == pytest.approx(expected[k], rel=1e-9), k


def test_identical_cells_have_zero_paired_difference():
    rng = np.random.default_rng(0)
    p = rng.uniform(0.01, 0.9, 2000)
    y = (rng.random(2000) < p).astype(int)
    patients = rng.integers(0, 800, 2000)
    cells = {name: {"cohort": "c", "p": p, "y": y, "patients": patients, "threshold": 0.2} for name in ("a", "b")}
    _, diff = evaluate.evaluate_cells(cells, n_boot=20, seed=3, comparisons=[("a", "b")])
    assert diff["a minus b"]["interval"]["auroc"] == (0.0, 0.0)


def test_small_cells_are_withheld():
    p = np.linspace(0.01, 0.99, 400)
    y = np.zeros(400, dtype=int)
    y[-5:] = 1  # five events
    assert evaluate.withheld_metrics(p, y, 0.5) == {"*"}
    counts = evaluate.confusion_counts(p, y, 0.5)
    assert counts["0.05"] == {"tp": WITHHELD, "fp": WITHHELD, "tn": WITHHELD, "fn": WITHHELD}


def test_a_single_small_subgroup_takes_another_with_it():
    rng = np.random.default_rng(2)
    n = 3000
    p = rng.uniform(0.01, 0.9, n)
    y = (rng.random(n) < p).astype(int)
    specimen = np.array(["urine"] * 2000 + ["blood"] * 994 + ["respiratory"] * 6)
    arrays = {"setting_primary": np.array(["inpatient"] * n), "specimen_type": specimen,
              "prior_culture_history": np.zeros(n, dtype=int), "adi_quintile": np.array(["q1"] * n)}
    out = evaluate.subgroups(p, y, 0.3, arrays)["specimen_type"]
    assert out["respiratory"]["units"] == WITHHELD and out["blood"]["units"] == WITHHELD
    assert out["urine"]["units"] == 2000


def test_evaluation_is_deterministic(trained):
    a = evaluate.evaluate_site(synthetic_test(), trained, n_boot=10, seed=4)
    b = evaluate.evaluate_site(synthetic_test(), trained, n_boot=10, seed=4)
    assert json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True, default=str)


def test_cross_site_uses_the_source_models_on_the_target_cohorts(trained):
    per_site = {"a": (synthetic_test(seed=5), trained), "b": (synthetic_test(seed=6), trained)}
    out = evaluate.evaluate_cross_site(per_site, n_boot=5, seed=1)
    assert set(out) == {"a to b", "b to a"}
    assert set(out["a to b"]) == {f"b:{k}:{v}->{v}" for k in ("xgboost", "logistic") for v in ("E1", "E2")}


def test_only_this_script_builds_test_cohorts(stanford, monkeypatch):
    loader, built, manifest = stanford
    with pytest.raises(splits.SealedPartitionError):
        cohorts.build_cohort(loader, built, CohortSpec("E2"), partitions=("test",), manifest=manifest)
    monkeypatch.setattr(sys, "argv", [str(REPO_ROOT / "analysis" / "03_evaluate_ladder.py")])
    result = cohorts.build_cohort(loader, built, CohortSpec("E2"), partitions=("test",), manifest=manifest)
    got = loader.con.execute(f"SELECT DISTINCT partition FROM {result.table}").fetchall()
    assert got == [("test",)]


def test_refuses_to_run_before_the_plan_is_frozen(monkeypatch):
    monkeypatch.setattr(evaluate, "FROZEN_TAG", "sap-never-tagged")
    with pytest.raises(SystemExit):
        evaluate.require_frozen_plan()


def test_specimen_grids_partition_the_test_cohort(trained):
    out = evaluate.evaluate_site(synthetic_test(), trained, n_boot=5, seed=1)
    grids = out["specimen_grids"]
    assert set(grids) == {"urine", "blood"}
    assert all(len(g["grid"]) == 12 for g in grids.values())
    # Every synthetic unit is urine or blood, so the two grids split the E2 test cohort.
    total = out["grid"]["xgboost:E2->E2"]["units"]
    assert grids["urine"]["grid"]["urine:xgboost:E2->E2"]["units"] + grids["blood"]["grid"]["blood:xgboost:E2->E2"]["units"] == total
    assert "urine:xgboost:E1_platt->E2 minus urine:xgboost:E2->E2" in grids["urine"]["comparisons"]


def site_with(diff: float, interval, slope_interval):
    keys = [f"standardized_net_benefit_{t:.2f}" for t in metrics.NB_THRESHOLDS]
    return {
        "comparisons": {"xgboost:E1_platt->E2 minus xgboost:E2->E2": {
            "point": {k: diff for k in keys}, "interval": {k: interval for k in keys}}},
        "grid": {"xgboost:E1_platt->E2": {"interval": {"calibration_slope": slope_interval}}},
    }


def test_null_result():
    assert evaluate.null_result(site_with(0.01, [-0.04, 0.06], [0.9, 1.1]))["negligible"] is True
    # A gap beyond the margin, or a slope interval that excludes 1, is not a null result.
    assert evaluate.null_result(site_with(0.08, [0.02, 0.14], [0.9, 1.1]))["negligible"] is False
    assert evaluate.null_result(site_with(0.01, [-0.04, 0.06], [1.2, 1.5]))["negligible"] is False
    # A withheld value leaves the verdict open, unless another check has already failed.
    assert evaluate.null_result(site_with(WITHHELD, WITHHELD, [0.9, 1.1]))["negligible"] is None
    assert evaluate.null_result(site_with(WITHHELD, WITHHELD, [1.2, 1.5]))["negligible"] is False
    assert evaluate.conjunction([True, None, False]) is False and evaluate.conjunction([True, None]) is None
