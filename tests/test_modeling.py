"""Cross-fitting, recalibration and thresholds on synthetic data."""

import numpy as np
import pytest

from armd_common import modeling

pytest.importorskip("xgboost")
pytest.importorskip("statsmodels")


class MeanModel:
    """Predicts the training prevalence and records which folds it trained on."""

    def fit(self, X, y):
        self.p = y.mean()
        self.seen = set(np.unique(X[:, 0]))
        return self

    def predict_proba(self, X):
        return np.column_stack([1 - np.full(len(X), self.p), np.full(len(X), self.p)])


def synthetic(n=4000, k=48, seed=0, prevalence_shift=0.0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, k))
    logit = -2.0 + 1.2 * X[:, 0] - 0.8 * X[:, 1] + 0.5 * X[:, 2] + prevalence_shift
    y = (rng.random(n) < 1 / (1 + np.exp(-logit))).astype(int)
    folds = rng.integers(1, 6, size=n)
    return X, y, folds


def test_cross_fit_predicts_each_fold_from_the_other_folds():
    folds = np.repeat([1, 2, 3], 4)
    X = folds.reshape(-1, 1).astype(float)
    y = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0])
    fit = modeling.cross_fit(MeanModel, X, y, folds)
    # Fold 1 is predicted by a model trained on folds 2 and 3, whose prevalence is 2 of 8.
    assert np.allclose(fit.predictions[folds == 1], 0.25)
    assert np.allclose(fit.predictions[folds == 2], 0.75)
    assert np.allclose(fit.predictions[folds == 3], 0.5)
    for fold, model in fit.fold_models.items():
        assert fold not in model.seen


def test_predict_by_fold_uses_the_held_out_model():
    folds = np.repeat([1, 2, 3], 4)
    X = folds.reshape(-1, 1).astype(float)
    y = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0])
    fit = modeling.cross_fit(MeanModel, X, y, folds)
    other = np.array([3, 1, 2])
    assert np.allclose(fit.predict_by_fold(other.reshape(-1, 1).astype(float), other), [0.5, 0.25, 0.75])
    with pytest.raises(ValueError):
        fit.predict_by_fold(np.zeros((1, 1)), np.array([4]))


def test_cross_fit_rejects_missing_folds_and_single_class():
    X = np.zeros((4, 1))
    with pytest.raises(ValueError):
        modeling.cross_fit(MeanModel, X, np.array([0, 1, 0, 1]), np.array([1, 2, np.nan, 1]))
    with pytest.raises(ValueError):
        modeling.cross_fit(MeanModel, X, np.array([1, 1, 0, 0]), np.array([1, 1, 2, 2]))


def test_platt_recovers_a_known_map():
    rng = np.random.default_rng(1)
    p = rng.uniform(0.01, 0.6, size=200_000)
    logit_p = np.log(p / (1 - p))
    q = 1 / (1 + np.exp(-(-2.3 + 0.8 * logit_p)))
    y = (rng.random(len(p)) < q).astype(int)
    m = modeling.fit_platt(p, y)
    assert m.intercept == pytest.approx(-2.3, abs=0.05)
    assert m.slope == pytest.approx(0.8, abs=0.03)
    assert np.allclose(m(p[:5]), q[:5], atol=0.01)


def test_isotonic_is_monotone_and_bounded():
    X, y, _ = synthetic()
    p = 1 / (1 + np.exp(-(X[:, 0] - 2)))
    iso = modeling.fit_isotonic(p, y)
    grid = np.linspace(-0.5, 1.5, 50)
    out = iso.predict(grid)
    assert np.all(np.diff(out) >= 0) and out.min() >= 0 and out.max() <= 1


def test_sensitivity_threshold():
    p = np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.99, 0.01])
    y = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0])
    # Nine of ten positives lie at or above 0.1, and 0.1 is the highest such threshold.
    assert modeling.sensitivity_threshold(p, y) == 0.1
    assert (p[y == 1] >= 0.1).mean() >= 0.9
    with pytest.raises(ValueError):
        modeling.sensitivity_threshold(p, np.zeros_like(y))


@pytest.mark.parametrize("kind", modeling.MODEL_KINDS)
def test_models_are_deterministic_for_a_seed(kind):
    X, y, _ = synthetic(n=1500)
    params = {"n_estimators": 25}
    a = modeling.predict(modeling.make_model(kind, 7, params, threads=2).fit(X, y), X)
    b = modeling.predict(modeling.make_model(kind, 7, params, threads=2).fit(X, y), X)
    assert np.array_equal(a, b)
    assert 0 < a.min() and a.max() < 1


def test_xgboost_keeps_fixed_hyperparameters():
    m = modeling.make_model("xgboost", 1)
    got = m.get_params()
    assert {k: got[k] for k in ("n_estimators", "max_depth", "scale_pos_weight")} == {"n_estimators": 700, "max_depth": 6, "scale_pos_weight": 1}
    with pytest.raises(ValueError):
        modeling.make_model("forest", 1)
