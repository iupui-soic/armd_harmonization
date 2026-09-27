"""Performance metrics for risk predictions, with patient-clustered bootstrap intervals.

Every metric takes optional frequency weights. The bootstrap resamples patients by drawing a count
for each patient and weighting that patient's units by it, so no rows are copied and every cell
evaluated on one cohort can share the same draw (paired differences).

Metrics:
  auroc, auprc         discrimination (auprc is average precision, as in scikit-learn)
  brier                mean squared error of the predicted risk
  calibration          calibration intercept (calibration in the large: logistic model with the
                       logit of the prediction as offset) and slope (logistic model on the logit)
  net_benefit          NB(t) = TP/n - FP/n * t / (1 - t), and standardized NB = NB / prevalence
  operating_point      sensitivity, specificity, PPV, NPV and net benefit at a threshold
"""

from __future__ import annotations

import numpy as np

_EPS = 1e-7
NB_THRESHOLDS = (0.05, 0.10, 0.15, 0.20)
DECISION_CURVE = tuple(round(t, 2) for t in np.arange(0.01, 0.401, 0.01))


def _weights(w, n: int) -> np.ndarray:
    return np.ones(n) if w is None else np.asarray(w, dtype=float)


def logit(p) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), _EPS, 1 - _EPS)
    return np.log(p / (1 - p))


class Ranked:
    """Predictions sorted once, so weighted AUROC and AUPRC cost one pass per resample."""

    def __init__(self, p, y):
        p, y = np.asarray(p, dtype=float), np.asarray(y, dtype=float)
        scores, inverse = np.unique(-p, return_inverse=True)
        self.group = inverse.ravel()
        self.n_groups = len(scores)
        self.y = y

    def _curve(self, w):
        w = _weights(w, len(self.y))
        tp = np.cumsum(np.bincount(self.group, w * self.y, self.n_groups))
        fp = np.cumsum(np.bincount(self.group, w * (1 - self.y), self.n_groups))
        return tp, fp

    def auroc(self, w=None) -> float:
        tp, fp = self._curve(w)
        if tp[-1] == 0 or fp[-1] == 0:
            return float("nan")
        tpr = np.concatenate([[0.0], tp / tp[-1]])
        fpr = np.concatenate([[0.0], fp / fp[-1]])
        return float(np.sum(np.diff(fpr) * (tpr[1:] + tpr[:-1]) / 2))

    def auprc(self, w=None) -> float:
        tp, fp = self._curve(w)
        if tp[-1] == 0:
            return float("nan")
        recall = tp / tp[-1]
        precision = tp / np.maximum(tp + fp, _EPS)
        return float(np.sum(np.diff(np.concatenate([[0.0], recall])) * precision))


def brier(p, y, w=None) -> float:
    w = _weights(w, len(y))
    return float(np.sum(w * (np.asarray(p, dtype=float) - y) ** 2) / w.sum())


def _logistic(X: np.ndarray, y: np.ndarray, w: np.ndarray, offset: np.ndarray, iterations: int = 50) -> np.ndarray:
    beta = np.zeros(X.shape[1])
    for _ in range(iterations):
        mu = 1 / (1 + np.exp(-(X @ beta + offset)))
        grad = X.T @ (w * (y - mu))
        hess = (X * (w * mu * (1 - mu))[:, None]).T @ X
        step = np.linalg.solve(hess, grad)
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return beta


def calibration(p, y, w=None) -> tuple[float, float]:
    """(intercept, slope). The intercept is estimated with the slope fixed at 1."""
    y = np.asarray(y, dtype=float)
    w = _weights(w, len(y))
    lp = logit(p)
    ones = np.ones((len(y), 1))
    intercept = _logistic(ones, y, w, lp)[0]
    slope = _logistic(np.column_stack([ones, lp]), y, w, np.zeros(len(y)))[1]
    return float(intercept), float(slope)


def confusion(p, y, t: float, w=None) -> dict:
    w = _weights(w, len(y))
    flag = np.asarray(p) >= t
    y = np.asarray(y).astype(bool)
    return {"tp": float(w[flag & y].sum()), "fp": float(w[flag & ~y].sum()),
            "tn": float(w[~flag & ~y].sum()), "fn": float(w[~flag & y].sum())}


def net_benefit(p, y, t: float, w=None) -> float:
    c = confusion(p, y, t, w)
    n = sum(c.values())
    return (c["tp"] - c["fp"] * t / (1 - t)) / n


def treat_all_net_benefit(y, t: float, w=None) -> float:
    w = _weights(w, len(y))
    prevalence = float(np.sum(w * y) / w.sum())
    return prevalence - (1 - prevalence) * t / (1 - t)


def _ratio(a: float, b: float) -> float:
    return a / b if b else float("nan")


def operating_point(p, y, t: float, w=None) -> dict:
    c = confusion(p, y, t, w)
    n = sum(c.values())
    return {
        "sensitivity": _ratio(c["tp"], c["tp"] + c["fn"]),
        "specificity": _ratio(c["tn"], c["tn"] + c["fp"]),
        "ppv": _ratio(c["tp"], c["tp"] + c["fp"]),
        "npv": _ratio(c["tn"], c["tn"] + c["fn"]),
        "net_benefit": (c["tp"] - c["fp"] * t / (1 - t)) / n,
    }


def evaluate(p, y, operating_threshold: float, w=None, ranked: Ranked | None = None,
             thresholds=NB_THRESHOLDS) -> dict:
    """Every scalar metric of the SAP for one set of predictions."""
    y = np.asarray(y, dtype=float)
    w = _weights(w, len(y))
    ranked = ranked or Ranked(p, y)
    prevalence = float(np.sum(w * y) / w.sum())
    intercept, slope = calibration(p, y, w)
    out = {
        "prevalence": prevalence,
        "auroc": ranked.auroc(w),
        "auprc": ranked.auprc(w),
        "brier": brier(p, y, w),
        "calibration_intercept": intercept,
        "calibration_slope": slope,
    }
    for t in thresholds:
        nb = net_benefit(p, y, t, w)
        out[f"net_benefit_{t:.2f}"] = nb
        out[f"standardized_net_benefit_{t:.2f}"] = nb / prevalence if prevalence else float("nan")
    for k, v in operating_point(p, y, operating_threshold, w).items():
        out[f"{k}_at_90_sensitivity"] = v
    return out


def patient_weights(patient_ids, n_boot: int, seed: int):
    """Yield one weight vector per bootstrap draw, resampling whole patients with replacement."""
    _, codes = np.unique(np.asarray(patient_ids), return_inverse=True)
    codes = codes.ravel()
    n_patients = int(codes.max()) + 1
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        counts = np.bincount(rng.integers(0, n_patients, n_patients), minlength=n_patients)
        yield counts[codes].astype(float)


def percentile_interval(samples, alpha: float = 0.05) -> tuple[float, float]:
    samples = np.asarray(samples, dtype=float)
    samples = samples[~np.isnan(samples)]
    if len(samples) == 0:
        return float("nan"), float("nan")
    lo, hi = np.quantile(samples, [alpha / 2, 1 - alpha / 2])
    return float(lo), float(hi)
