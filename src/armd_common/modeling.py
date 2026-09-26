"""Models, cross-fitting, recalibration maps and operating thresholds for the benchmark tasks.

Two model classes:
  xgboost   gradient-boosted trees with fixed hyperparameters (XGB_PARAMS), no re-tuning
  logistic  L2-penalized logistic regression, C = 1.0, on features standardized with the training
            means and standard deviations

Neither uses class weighting, because calibration and net benefit depend on the probability scale.

Cross-fitting uses the frozen inner folds of the split manifest. For each fold, a model trained on
the other folds predicts that fold, so every derivation row gets a prediction from a model that never
saw its patient. The same fold models can predict another cohort's rows fold by fold, which is how an
E1-trained model is recalibrated toward E2 without touching the test partition.

Needs the `model` extra: pip install "armd-common[model]".
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

MODEL_KINDS = ("xgboost", "logistic")

# Tuned hyperparameters of the AMIA 2026 ARMD-MGB model, held fixed.
XGB_PARAMS = {
    "n_estimators": 700,
    "max_depth": 6,
    "learning_rate": 0.02001,
    "min_child_weight": 6,
    "subsample": 0.8807,
    "colsample_bytree": 0.5222,
    "gamma": 2.596,
    "reg_alpha": 8.608,
    "reg_lambda": 5.384,
    "objective": "binary:logistic",
    "scale_pos_weight": 1,
    "tree_method": "hist",
    "eval_metric": "logloss",
}
LOGISTIC_C = 1.0
SENSITIVITY_TARGET = 0.90
_EPS = 1e-7


def make_model(kind: str, seed: int, params: dict | None = None, threads: int = 1):
    """An unfitted estimator. `params` overrides XGB_PARAMS, for tests only.

    XGBoost results repeat for a given seed and thread count, so callers record `threads`.
    """
    if kind == "xgboost":
        from xgboost import XGBClassifier

        return XGBClassifier(**{**XGB_PARAMS, **(params or {})}, random_state=seed, n_jobs=threads)
    if kind == "logistic":
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        return make_pipeline(StandardScaler(), LogisticRegression(C=LOGISTIC_C, max_iter=5000))
    raise ValueError(f"model kind must be one of {MODEL_KINDS}")


def predict(model, X: np.ndarray) -> np.ndarray:
    return model.predict_proba(X)[:, 1]


@dataclass
class CrossFit:
    """Out-of-fold predictions and the fold models that made them."""

    predictions: np.ndarray
    fold_models: dict

    def predict_by_fold(self, X: np.ndarray, folds: np.ndarray) -> np.ndarray:
        """Predict rows of another cohort, each with the model that did not train on its fold."""
        folds = np.asarray(folds)
        unknown = set(np.unique(folds)) - set(self.fold_models)
        if unknown:
            raise ValueError(f"no fold model for folds {sorted(unknown)}")
        out = np.full(len(folds), np.nan)
        for fold, model in self.fold_models.items():
            rows = folds == fold
            if rows.any():
                out[rows] = predict(model, X[rows])
        return out


def cross_fit(factory: Callable[[], object], X: np.ndarray, y: np.ndarray, folds: np.ndarray) -> CrossFit:
    """Fit one model per fold on the other folds and predict the held-out fold."""
    folds = np.asarray(folds)
    if np.isnan(folds.astype(float)).any():
        raise ValueError("every training row needs an inner fold")
    out = np.full(len(y), np.nan)
    fold_models = {}
    for fold in np.unique(folds):
        held = folds == fold
        if y[~held].min() == y[~held].max():
            raise ValueError(f"training rows outside fold {fold} hold one class only")
        model = factory()
        model.fit(X[~held], y[~held])
        out[held] = predict(model, X[held])
        fold_models[fold.item() if hasattr(fold, "item") else fold] = model
    return CrossFit(out, fold_models)


def _logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), _EPS, 1 - _EPS)
    return np.log(p / (1 - p))


@dataclass(frozen=True)
class PlattMap:
    """Logistic recalibration: logit(q) = intercept + slope * logit(p)."""

    intercept: float
    slope: float

    def __call__(self, p: np.ndarray) -> np.ndarray:
        return 1 / (1 + np.exp(-(self.intercept + self.slope * _logit(p))))


def fit_platt(p: np.ndarray, y: np.ndarray) -> PlattMap:
    import statsmodels.api as sm

    X = sm.add_constant(_logit(p), has_constant="add")
    fit = sm.GLM(np.asarray(y, dtype=float), X, family=sm.families.Binomial()).fit()
    intercept, slope = (float(v) for v in fit.params)
    return PlattMap(intercept, slope)


def fit_isotonic(p: np.ndarray, y: np.ndarray):
    """Isotonic recalibration. Its breakpoints come from row-level outcomes, so store it with the
    data, never in a shared results file."""
    from sklearn.isotonic import IsotonicRegression

    return IsotonicRegression(y_min=0.0, y_max=1.0, out_of_bounds="clip").fit(np.asarray(p, dtype=float), np.asarray(y, dtype=float))


def sensitivity_threshold(p: np.ndarray, y: np.ndarray, target: float = SENSITIVITY_TARGET) -> float:
    """The highest threshold t for which predicting positive at p >= t reaches the target sensitivity."""
    positives = np.sort(np.asarray(p, dtype=float)[np.asarray(y) == 1])[::-1]
    if len(positives) == 0:
        raise ValueError("no positive rows to set a sensitivity threshold")
    needed = int(np.ceil(target * len(positives)))
    return float(positives[needed - 1])
