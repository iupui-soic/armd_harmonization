"""Metrics against scikit-learn and statsmodels, on synthetic predictions with ties and weights."""

import numpy as np
import pytest

from armd_common import metrics

sklearn_metrics = pytest.importorskip("sklearn.metrics")
sm = pytest.importorskip("statsmodels.api")


def synthetic(n=5000, seed=0):
    rng = np.random.default_rng(seed)
    p = np.round(rng.beta(1, 8, size=n), 3)  # rounding makes ties
    y = (rng.random(n) < np.clip(1.5 * p, 0, 1)).astype(int)
    w = rng.integers(0, 4, size=n).astype(float)
    return p, y, w


@pytest.mark.parametrize("weighted", [False, True])
def test_ranked_matches_sklearn(weighted):
    p, y, w = synthetic()
    w = w if weighted else None
    ranked = metrics.Ranked(p, y)
    assert ranked.auroc(w) == pytest.approx(sklearn_metrics.roc_auc_score(y, p, sample_weight=w), abs=1e-12)
    assert ranked.auprc(w) == pytest.approx(sklearn_metrics.average_precision_score(y, p, sample_weight=w), abs=1e-12)


def test_weights_equal_repeated_rows():
    p, y, w = synthetic(n=800)
    rep = w.astype(int)
    expanded = metrics.evaluate(np.repeat(p, rep), np.repeat(y, rep), 0.1)
    weighted = metrics.evaluate(p, y, 0.1, w=w)
    assert expanded.keys() == weighted.keys()
    for k in expanded:
        assert weighted[k] == pytest.approx(expanded[k], rel=1e-8, abs=1e-10), k


def test_calibration_matches_statsmodels():
    p, y, _ = synthetic()
    lp = metrics.logit(p)
    slope = sm.GLM(y, sm.add_constant(lp), family=sm.families.Binomial()).fit().params[1]
    intercept = sm.GLM(y, np.ones((len(y), 1)), family=sm.families.Binomial(), offset=lp).fit().params[0]
    got = metrics.calibration(p, y)
    assert got == pytest.approx((intercept, slope), abs=1e-8)


def test_perfect_calibration():
    rng = np.random.default_rng(3)
    p = rng.uniform(0.01, 0.5, size=200_000)
    y = (rng.random(len(p)) < p).astype(int)
    intercept, slope = metrics.calibration(p, y)
    assert intercept == pytest.approx(0, abs=0.02) and slope == pytest.approx(1, abs=0.03)


def test_net_benefit_by_hand():
    p = np.array([0.9, 0.6, 0.3, 0.2, 0.05])
    y = np.array([1, 0, 1, 0, 0])
    # At t = 0.25 three are flagged: two true positives out of 5 and one false positive.
    assert metrics.net_benefit(p, y, 0.25) == pytest.approx((2 - 1 * 0.25 / 0.75) / 5)
    assert metrics.treat_all_net_benefit(y, 0.25) == pytest.approx(0.4 - 0.6 / 3)
    op = metrics.operating_point(p, y, 0.25)
    assert op["sensitivity"] == 1.0 and op["ppv"] == pytest.approx(2 / 3) and op["npv"] == 1.0 and op["specificity"] == pytest.approx(2 / 3)


def test_patient_bootstrap_keeps_patients_together():
    patients = np.array(["a", "a", "b", "c", "c", "c"])
    for w in metrics.patient_weights(patients, 50, seed=1):
        # A patient's units always share one weight, and the draw has as many patients as the data.
        assert w[0] == w[1] and w[3] == w[4] == w[5]
        assert w[0] + w[2] + w[3] == 3


def test_bootstrap_is_seeded():
    patients = np.arange(100) % 37
    a = np.array(list(metrics.patient_weights(patients, 5, seed=9)))
    b = np.array(list(metrics.patient_weights(patients, 5, seed=9)))
    assert np.array_equal(a, b)


def test_percentile_interval_ignores_nan():
    assert metrics.percentile_interval([np.nan, 1, 2, 3, 4, 5], alpha=0.5) == pytest.approx((2, 4))


def test_warm_started_calibration_matches_cold():
    p, y, w = synthetic()
    lp = metrics.logit(p)
    start = metrics.calibration_start(p, y, lp)
    assert start[:2] == pytest.approx(metrics.calibration(p, y), abs=1e-10)
    assert metrics.calibration(p, y, w, lp=lp, start=start) == pytest.approx(metrics.calibration(p, y, w), abs=1e-9)


def test_degenerate_calibration_is_nan_not_an_error():
    # Predictions of exactly 0 and 1 that match the outcome leave no information about the slope.
    p = np.array([0.0, 0.0, 1.0, 1.0])
    y = np.array([0, 0, 1, 1])
    intercept, slope = metrics.calibration(p, y)
    assert np.isnan(slope) or abs(slope) > 1


def test_per_thousand():
    p = np.array([0.9, 0.6, 0.3, 0.2, 0.05])
    y = np.array([1, 0, 1, 0, 0])
    assert metrics.per_thousand(p, y, 0.5) == {"flagged": 400.0, "covered": 200.0, "missed": 200.0}
    # Weights count a unit that many times.
    got = metrics.per_thousand(p, y, 0.5, w=np.array([2, 0, 1, 1, 1]))
    assert got == {"flagged": 400.0, "covered": 400.0, "missed": 200.0}
