"""The exploratory amendment (SAP section 17) in analysis/03_evaluate_ladder.py, on synthetic data."""

import json

import numpy as np
import pytest

from armd_common import metrics
from armd_common.qc import WITHHELD

sklearn_metrics = pytest.importorskip("sklearn.metrics")
from test_evaluate_ladder import evaluate, synthetic_test, train, trained  # noqa: E402,F401  (trained is a fixture)
from test_train_e1 import FAST, synthetic_cohorts  # noqa: E402


def test_flag_rate_at_sensitivity_by_hand():
    p = np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4])
    y = np.array([1, 0, 1, 0, 1, 0])
    # Two of three positives lie at or above 0.7, where three of six orders are flagged.
    assert evaluate.flag_rate_at_sensitivity(p, y, 2 / 3) == pytest.approx(500)
    assert evaluate.flag_rate_at_sensitivity(p, y, 1.0) == pytest.approx(1000 * 5 / 6)
    # A weight of 0 removes the first positive, so reaching all positives still needs 0.5.
    assert evaluate.flag_rate_at_sensitivity(p, y, 1.0, w=np.array([0, 1, 1, 1, 1, 1])) == pytest.approx(800)


def rule_data(n=4000, seed=0):
    rng = np.random.default_rng(seed)
    score = rng.integers(0, 10, n).astype(float)
    y = (rng.random(n) < 0.02 + 0.03 * score).astype(int)
    return score, score >= 5, y, rng.integers(0, n // 2, n)


def test_evaluate_rule_matches_its_definitions():
    score, flags, y, patients = rule_data()
    out = evaluate.evaluate_rule(score, flags, y, patients, n_boot=20, seed=1, key="E2")
    pt = out["point"]
    assert pt["auroc"] == pytest.approx(sklearn_metrics.roc_auc_score(y, score))
    assert pt["flagged_per_1000"] == pytest.approx(1000 * flags.mean())
    assert pt["sensitivity"] == pytest.approx(flags[y == 1].mean())
    t, w = 0.10, 0.10 / 0.90
    tp, fp = int((flags & (y == 1)).sum()), int((flags & (y == 0)).sum())
    assert pt["net_benefit_0.10"] == pytest.approx((tp - fp * w) / len(y))
    lo, hi = out["interval"]["sensitivity"]
    assert lo <= pt["sensitivity"] <= hi


def test_evaluate_rule_withholds_small_cells():
    score, flags, y, patients = rule_data()
    flags = np.zeros_like(flags)
    flags[np.flatnonzero(y == 1)[:5]] = True  # five true positives and no other flags
    out = evaluate.evaluate_rule(score, flags, y, patients, n_boot=5, seed=1, key="E2")
    assert out["point"]["sensitivity"] == WITHHELD and isinstance(out["point"]["auroc"], float)
    assert out["confusion"]["tp"] == WITHHELD


def test_matched_comparison():
    score, _, y, patients = rule_data()
    p_model = score / 10 + np.random.default_rng(3).normal(0, 0.01, len(score))
    # A rule that is the model's own flags at 0.5: at the same sensitivity the model's threshold is its
    # lowest flagged positive, which can only drop negatives, so the model needs no more flags.
    flags = p_model >= 0.5
    out = evaluate.matched_comparison(flags, p_model, y, patients, n_boot=10, seed=1, key="E2")
    assert 0 <= out["point"]["difference_per_1000"] < 5
    # A rule that flags at random needs more flags than the model for the same sensitivity.
    rand = np.random.default_rng(4).random(len(y)) < 0.5
    worse = evaluate.matched_comparison(rand, p_model, y, patients, n_boot=10, seed=1, key="E2")
    assert worse["point"]["difference_per_1000"] > 0


def test_decision_impact_matches_per_thousand():
    score, _, y, patients = rule_data()
    p = score / 10
    out = evaluate.decision_impact(p, y, patients, 0.3, n_boot=10, seed=1, key="E2")
    assert set(out) == {"0.05", "0.10", "0.20", "90_sensitivity"}
    assert out["0.20"]["point"] == pytest.approx(metrics.per_thousand(p, y, 0.20))
    assert out["90_sensitivity"]["threshold"] == 0.3


@pytest.fixture(scope="module")
def urine_setup(tmp_path_factory):
    root = tmp_path_factory.mktemp("urine_site")
    synthetic_cohorts(root / "cohorts")
    names = train.cohort_names("release", "culture")
    summary = train.train_site("site", root / "cohorts", root / "urine_models", "release", "culture",
                               ("xgboost", "logistic"), seed=20260911, threads=2, params=FAST, specimen="urine")
    variants = {k: evaluate.load_variants(root / "urine_models", names, k, summary) for k in ("xgboost", "logistic")}
    return variants, root / "cohorts" / names["E2"] / "derivation.parquet"


def test_exploratory_site(trained, urine_setup):
    urine_variants, derivation = urine_setup
    out = evaluate.exploratory_site(synthetic_test(), trained, urine_variants, derivation, n_boot=5, seed=1)
    assert set(out["rules"]) == {f"{r}->{t}" for r in evaluate.RULES for t in ("E1", "E2")}
    assert {"R3_at_0.10->E1", "R3_at_0.10->E2"} <= set(out["matched"])
    assert set(out["minimal_model"]["grid"]) == {"R3->E1", "R3->E2"}
    assert set(out["minimal_model"]["coefficients"]) == {"intercept", "specimen_urine", "specimen_blood",
                                                         "specimen_respiratory", "prior_ESBL", "recent_blfq"}
    assert len(out["urine_models"]["grid"]) == 12
    assert "urine:xgboost:E1_platt->E2 minus urine:xgboost:E2->E2" in out["urine_models"]["comparisons"]
    assert len(out["decision_impact"]) == 12 + 2 + 12
    json.dumps(out, default=str)
    assert "## " not in evaluate.render_exploratory({"sites": {}})


# Second exploratory amendment, SAP section 17, E-6 to E-8 -------------------------------------------

KINDS = ("xgboost", "logistic")


def test_decision_impact_thresholds_and_decision_curve():
    score, _, y, patients = rule_data()
    p = score / 10
    out = evaluate.decision_impact(p, y, patients, 0.3, n_boot=5, seed=1, key="E2", thresholds=(0.15, 0.25))
    assert set(out) == {"0.15", "0.25", "90_sensitivity"}
    assert out["0.25"]["point"] == pytest.approx(metrics.per_thousand(p, y, 0.25))
    assert set(evaluate.decision_impact(p, y, patients, 0.3, n_boot=5, seed=1, key="E2", thresholds=())) == {"90_sensitivity"}
    curve = evaluate.decision_curve(p, y)
    assert set(curve) == {"decision_curve", "treat_all"}
    assert curve["decision_curve"]["0.30"] == pytest.approx(metrics.net_benefit(p, y, 0.30))
    assert set(evaluate.curves(p, y)) == {"decision_curve", "treat_all", "calibration_curve"}


@pytest.fixture(scope="module")
def second_setup(tmp_path_factory):
    root = tmp_path_factory.mktemp("second_site")
    synthetic_cohorts(root / "cohorts")
    names = train.cohort_names("release", "culture")
    kw = dict(seed=20260911, threads=2, params=FAST)
    pooled = train.train_site("site", root / "cohorts", root / "models", "release", "culture", KINDS, **kw)
    blood = train.train_site("site", root / "cohorts", root / "blood_models", "release", "culture", KINDS, specimen="blood", **kw)
    variants = {k: evaluate.load_variants(root / "models", names, k, pooled) for k in KINDS}
    blood_variants = {k: evaluate.load_variants(root / "blood_models", names, k, blood) for k in KINDS}
    thresholds = {k: evaluate.specimen_thresholds(root / "models", root / "cohorts", names, k)[0] for k in KINDS}
    derivation = {e: root / "cohorts" / names[e] / "derivation.parquet" for e in ("E1", "E2")}
    return root, names, variants, blood_variants, thresholds, derivation


def test_specimen_thresholds_by_hand(second_setup):
    import duckdb

    from armd_common import modeling

    root, names, _, _, thresholds, _ = second_setup
    got, files = evaluate.specimen_thresholds(root / "models", root / "cohorts", names, "xgboost")
    assert got == thresholds["xgboost"] and set(got) == {"urine", "blood"} and len(files) == 3
    cohort = root / "cohorts" / names["E2"] / "derivation.parquet"
    sources = {"E2": root / "models" / names["E2"] / "xgboost_crossfit.parquet",
               "E1": root / "models" / names["E1"] / f"xgboost_crossfit_on_{names['E2']}.parquet"}
    con = duckdb.connect()
    for variant, path in sources.items():
        for specimen in ("urine", "blood"):
            p, y = con.execute(f"SELECT p.prediction, p.label FROM read_parquet('{path}') p JOIN read_parquet('{cohort}') c "
                               f"USING (unit_id) WHERE c.specimen_type = '{specimen}'").fetchnumpy().values()
            t = got[specimen][variant]
            assert t == modeling.sensitivity_threshold(p, y)
            assert (p[y == 1] >= t).mean() >= 0.9 > (p[y == 1] > t).mean()
    con.close()
    assert got["urine"]["E2"] != got["blood"]["E2"]


def ceftriaxone_cohorts(base: dict, seed=11) -> dict:
    """The synthetic test cohorts plus ceftriaxone versions: some units lose their result and some positives turn negative."""
    rng = np.random.default_rng(seed)
    out = dict(base)
    for key in ("E1", "E2"):
        a = base[key]
        y = a["y"].copy()
        y[(y == 1) & (rng.random(len(y)) < 0.2)] = 0
        keep = rng.random(len(y)) > 0.05
        out[f"ceftriaxone:{key}"] = {**{k: np.asarray(v)[keep] for k, v in a.items()}, "y": y[keep]}
    return out


def test_second_exploratory_site(second_setup):
    _, _, variants, blood_variants, thresholds, derivation = second_setup
    cohorts = ceftriaxone_cohorts(synthetic_test())
    out = evaluate.second_exploratory_site(cohorts, variants, blood_variants, thresholds, derivation, n_boot=5, seed=1)

    grid = out["ceftriaxone"]["grid"]
    assert set(grid) == {f"{k}:{v}->ceftriaxone:{t}" for k in KINDS for v in ("E1", "E1_platt", "E2") for t in ("E1", "E2")}
    assert set(out["ceftriaxone"]["comparisons"]) == {f"{k}:{v}->ceftriaxone:E2 minus {k}:E2->ceftriaxone:E2"
                                                      for k in KINDS for v in ("E1", "E1_platt")}
    cell = grid["xgboost:E2->ceftriaxone:E2"]
    assert cell["units"] == len(cohorts["ceftriaxone:E2"]["y"]) < len(cohorts["E2"]["y"])
    assert cell["point"]["prevalence"] == pytest.approx(cohorts["ceftriaxone:E2"]["y"].mean())
    # The primary threshold is kept, since nothing is refitted.
    assert cell["threshold_90_sensitivity"] == variants["xgboost"]["E2"][1]

    blood = out["blood_models"]["grid"]
    assert len(blood) == 12 and set(blood["blood:xgboost:E2->E2"]["curves"]) == {"decision_curve", "treat_all"}
    assert blood["blood:xgboost:E2->E2"]["units"] == int((np.asarray(cohorts["E2"]["specimen_type"]) == "blood").sum())
    assert "blood:xgboost:E1_platt->E2 minus blood:xgboost:E2->E2" in out["blood_models"]["comparisons"]
    assert set(out["decision_impact"]["blood:xgboost:E2->E2"]) == {"0.01", "0.02", "0.05", "0.10", "90_sensitivity"}

    within = out["specimen_thresholds"]
    assert set(within) == {"urine", "blood"}
    assert set(within["blood"]) == {f"pooled:blood:{k}:{v}->E2" for k in KINDS for v in ("E1", "E2")}
    assert within["blood"]["pooled:blood:xgboost:E2->E2"]["threshold_90_sensitivity"] == thresholds["xgboost"]["blood"]["E2"]
    assert set(out["decision_impact"]["pooled:urine:logistic:E1->E2"]) == {"90_sensitivity"}

    only = out["specimen_only"]
    assert set(only["grid"]) == {f"specimen_only:{a}->{b}" for a in ("E1", "E2") for b in ("E1", "E2")}
    assert set(only["comparisons"]) == {"xgboost:E1->E2 minus specimen_only:E1->E2", "xgboost:E2->E2 minus specimen_only:E2->E2"}
    assert set(only["coefficients"]["E2"]) == {"intercept", "specimen_urine", "specimen_blood", "specimen_respiratory"}
    assert set(out["decision_impact"]["specimen_only:E2->E2"]) == {"0.05", "0.10", "0.20", "90_sensitivity"}

    assert set(out["rules"]) == set(out["matched"]) == {"tumbarello_6->E1", "tumbarello_6->E2"}
    assert out["rules"]["tumbarello_6->E2"]["point"]["flagged_per_1000"] <= 1000
    json.dumps(out, default=str)
    report = evaluate.render_second_exploratory({"sites": {"site": out}})
    assert "E-8" in report and "specimen_only:E2->E2" in report and "pooled:blood:xgboost:E2->E2" in report


def test_hide_near_counts():
    grid = {"xgboost:E2->ceftriaxone:E2": {"units": 995, "events": 80}, "xgboost:E2->ceftriaxone:E1": {"units": 300, "events": 80}}
    evaluate._hide_near_counts(grid, {"E2": (1000, 100), "E1": (300, 85)})
    assert grid["xgboost:E2->ceftriaxone:E2"] == {"units": WITHHELD, "events": 80}
    assert grid["xgboost:E2->ceftriaxone:E1"] == {"units": 300, "events": WITHHELD}
