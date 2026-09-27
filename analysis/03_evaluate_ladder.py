"""Evaluate the estimand ladder on the sealed test partition (plan sections 5.3, 5.4 and 5.7).

This is the only script that may read the test partition (splits.SEALED_READER). It refuses to run
until the analysis plan is frozen, that is until the repository has the git tag sap-v1.0, and the
manifest records the SHA-256 of the local SAP.md.

For each site it builds the E1 and E2 test cohorts with the repeat rule and unit used in training
and predicts them with the models from analysis/02_train_e1.py. Each model class gives six cells:

  E1 model on E1 test              reproduces published practice
  E1 model on E2 test              what deployment does
  E1 model, Platt map, on E2 test  recalibrated toward E2 (primary map)
  E1 model, isotonic, on E2 test   recalibrated toward E2 (secondary map)
  E2 model on E1 test              reported for completeness
  E2 model on E2 test              the corrected model

Each cell reports the planned metrics with 95% percentile intervals from patient-clustered bootstrap
draws. Cells on one test cohort share each draw, so the COMPARISONS are paired differences. The same
grid runs within urine and within blood (SPECIMEN_GRIDS). Subgroups are reported without intervals.
The sensitivity cohorts (S1 to S8) are evaluated with the primary models. Every ordered pair of sites
is evaluated with the source site's models on the target site's test cohorts. Table 1 describes the
full cohorts, derivation and test together. null_result applies the plan's definition of a null
result to the recalibrated comparison.

Predictions are row-level, so they are cached under models/<rule>_<unit>/test/ in the site folder of
the data root. results/03_evaluate_ladder/<timestamp>/ receives manifest.json, summary.json and
report.md. A metric is withheld when a confusion count behind it is between 1 and 10, and the four
counts at a threshold are withheld together.

With --exploratory it runs the SAP's exploratory amendment (section 17, E-1 to E-4) instead, and
writes results/03_exploratory/<timestamp>/. That mode needs the urine-only training run (E-3).

Usage: python analysis/03_evaluate_ladder.py [--train-run results/02_train_e1/<timestamp>]
       [--sites stanford,mgb] [--boot 2000] [--seed 20260911]
       python analysis/03_evaluate_ladder.py --exploratory [--urine-train-run results/02_train_e1/<timestamp>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pickle
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa

from armd_common import features, metrics, modeling, scores
from armd_common.cohorts import TABLE1_COVARIATES, CohortSpec, build_cohort
from armd_common.config import DEFAULT_SEED, MIN_CELL, REPO_ROOT, RESULTS_DIR, get_site
from armd_common.db import connect
from armd_common.provenance import RunManifest, require_git_tag, sha256_file
from armd_common.qc import WITHHELD, cell, publish_flows
from armd_common.sites import get_loader

FROZEN_TAG = "sap-v1.0"
N_BOOT = 2000
SUBGROUPS = ("setting_primary", "specimen_type", "prior_culture_history", "adi_quintile")
# (model variant, test cohort). A variant is a training estimand, with a recalibration map for E1.
GRID = (("E1", "E1"), ("E1", "E2"), ("E1_platt", "E2"), ("E1_isotonic", "E2"), ("E2", "E1"), ("E2", "E2"))
# Paired differences, each variant against the E2 model on the E2 test cohort.
COMPARISONS = (("E1", "E2", "E2"), ("E1_platt", "E2", "E2"), ("E1_isotonic", "E2", "E2"))
SENSITIVITY = {
    "S1": {"nonsusceptible": True},
    "S2": {"exclude_recent_antibiotics": True},
    "S3": {"interpretation": "clsi2022"},
    "S4": {"label": "idsa"},
    "S4_enzyme": {"label": "enzyme"},
    "S5": {"implied_susceptibility": "full"},
    "S6": {"unit": "order_event"},
    "S7": {"repeat_rule": "first_per_patient"},
    "S8": {"repeat_rule": "release", "unit": "culture"},
}
SPECIMEN_GRIDS = ("urine", "blood")
# A null result: every primary threshold within these margins of zero, for the recalibrated comparison.
NULL_MARGIN, NULL_INTERVAL_MARGIN = 0.05, 0.10
NULL_COMPARISON = ("E1_platt", "E2")
SENSITIVITY_VARIANTS = ("E1", "E1_platt", "E2")
CROSS_SITE = (("E1", "E1"), ("E2", "E2"))
OPERATING = ("sensitivity", "specificity", "ppv", "npv", "net_benefit")


def require_frozen_plan() -> dict:
    return require_git_tag(FROZEN_TAG, REPO_ROOT / "SAP.md")


# Small-cell rules ------------------------------------------------------------------------------

def small(n) -> bool:
    return 0 < n < MIN_CELL


def withheld_metrics(p, y, threshold: float) -> set[str]:
    """Names of the metrics that would expose a count from 1 to 10. '*' withholds all of them."""
    events = int(np.sum(y))
    if small(len(y)) or small(events) or small(len(y) - events):
        return {"*"}
    if events == 0 or events == len(y):
        return {"*"}  # nothing to estimate, which is not a disclosure but withheld the same way
    out = set()
    for t in metrics.NB_THRESHOLDS:
        if any(small(v) for v in metrics.confusion(p, y, t).values()):
            out |= {f"net_benefit_{t:.2f}", f"standardized_net_benefit_{t:.2f}"}
    if any(small(v) for v in metrics.confusion(p, y, threshold).values()):
        out |= {f"{k}_at_90_sensitivity" for k in OPERATING}
    return out


def publish(values: dict, withheld: set) -> dict:
    if "*" in withheld:
        return {k: WITHHELD for k in values}
    return {k: WITHHELD if k in withheld else v for k, v in values.items()}


def confusion_counts(p, y, threshold: float) -> dict:
    """Counts at each threshold for the bounds in analysis/05_e3_bounds.py, all four withheld together."""
    out = {}
    for name, t in [(f"{t:.2f}", t) for t in metrics.NB_THRESHOLDS] + [("90_sensitivity", threshold)]:
        c = {k: int(v) for k, v in metrics.confusion(p, y, t).items()}
        out[name] = {k: WITHHELD for k in c} if any(small(v) for v in c.values()) else c
    return out


# Evaluation --------------------------------------------------------------------------------------

def cohort_seed(seed: int, key: str) -> int:
    return int.from_bytes(hashlib.sha256(f"{seed}:{key}".encode()).digest()[:4], "little")


def evaluate_cells(cells: dict, n_boot: int, seed: int, comparisons=()) -> tuple[dict, dict]:
    """Point estimates, bootstrap intervals and paired differences.

    cells maps a name to {"cohort", "p", "y", "patients", "threshold"}. Cells with the same cohort key
    must be row-aligned, and they share every bootstrap draw.
    """
    prepared, by_cohort = {}, {}
    for name, c in cells.items():
        p, y = np.asarray(c["p"], dtype=float), np.asarray(c["y"], dtype=float)
        lp = metrics.logit(p)
        prepared[name] = {"p": p, "y": y, "lp": lp, "ranked": metrics.Ranked(p, y),
                          "start": metrics.calibration_start(p, y, lp), "threshold": float(c["threshold"])}
        by_cohort.setdefault(c["cohort"], []).append(name)
    samples = {name: [] for name in cells}
    for key, names in by_cohort.items():
        for w in metrics.patient_weights(cells[names[0]]["patients"], n_boot, cohort_seed(seed, key)):
            for name in names:
                q = prepared[name]
                samples[name].append(metrics.evaluate(q["p"], q["y"], q["threshold"], w, q["ranked"], lp=q["lp"], start=q["start"]))

    raw, out = {}, {}
    for name, q in prepared.items():
        point = metrics.evaluate(q["p"], q["y"], q["threshold"], ranked=q["ranked"], lp=q["lp"], start=q["start"])
        keys = list(point)
        draws = np.array([[s[k] for k in keys] for s in samples[name]]) if samples[name] else np.empty((0, len(keys)))
        withheld = withheld_metrics(q["p"], q["y"], q["threshold"])
        raw[name] = {"point": point, "draws": draws, "keys": keys, "withheld": withheld}
        intervals = {k: metrics.percentile_interval(draws[:, i]) for i, k in enumerate(keys)} if len(draws) else {}
        out[name] = {
            "cohort": cells[name]["cohort"],
            "units": cell(len(q["y"])),
            "events": cell(int(q["y"].sum())),
            "threshold_90_sensitivity": q["threshold"],
            "point": publish(point, withheld),
            "interval": publish(intervals, withheld),
            "confusion": confusion_counts(q["p"], q["y"], q["threshold"]),
        }
    differences = {}
    for a, b in comparisons:
        if a not in raw or b not in raw or cells[a]["cohort"] != cells[b]["cohort"]:
            continue
        ra, rb = raw[a], raw[b]
        withheld = ra["withheld"] | rb["withheld"]
        point = {k: ra["point"][k] - rb["point"][k] for k in ra["keys"]}
        diff = ra["draws"] - rb["draws"]
        intervals = {k: metrics.percentile_interval(diff[:, i]) for i, k in enumerate(ra["keys"])} if len(diff) else {}
        differences[f"{a} minus {b}"] = {"point": publish(point, withheld), "interval": publish(intervals, withheld)}
    return out, differences


def curves(p, y) -> dict:
    """Decision curve and a smoothed calibration curve for figures. Values only, never counts."""
    import statsmodels.api as sm

    p, y = np.asarray(p, dtype=float), np.asarray(y, dtype=float)
    decision = {}
    for t in metrics.DECISION_CURVE:
        c = metrics.confusion(p, y, t)
        decision[f"{t:.2f}"] = WITHHELD if any(small(v) for v in c.values()) else metrics.net_benefit(p, y, t)
    grid = np.quantile(p, np.linspace(0.05, 0.95, 19))
    observed = sm.nonparametric.lowess(y, p, frac=0.3, it=0, xvals=grid)
    return {
        "decision_curve": decision,
        "treat_all": {f"{t:.2f}": metrics.treat_all_net_benefit(y, t) for t in metrics.DECISION_CURVE},
        "calibration_curve": {"predicted": grid.tolist(), "observed": np.asarray(observed).tolist()},
    }


def subgroups(p, y, threshold: float, arrays: dict) -> dict:
    """Point estimates by subgroup, without intervals and without tests.

    A subgroup with 1 to 10 units, events or non-events is withheld. If only one value of a column is
    withheld, the next smallest is withheld too, so it cannot be recovered from the column total.
    """
    out = {}
    for column in SUBGROUPS:
        values = np.asarray(arrays[column]).astype(str)
        stats = {}
        for value in sorted(set(values)):
            rows = values == value
            stats[value] = (rows, int(rows.sum()), int(y[rows].sum()))
        hidden = {v for v, (_, n, e) in stats.items() if small(n) or small(e) or small(n - e)}
        if len(hidden) == 1:
            rest = sorted((n, v) for v, (_, n, _) in stats.items() if v not in hidden)
            if rest:
                hidden.add(rest[0][1])
        out[column] = {}
        for value, (rows, n, e) in stats.items():
            if value in hidden:
                out[column][value] = {"units": WITHHELD, "events": WITHHELD, "point": {}}
                continue
            withheld = withheld_metrics(p[rows], y[rows], threshold)
            point = metrics.evaluate(p[rows], y[rows], threshold) if "*" not in withheld else {}
            out[column][value] = {"units": cell(n), "events": cell(e), "point": publish(point, withheld)}
    return out


# Models and cohorts ------------------------------------------------------------------------------

def load_model(folder: Path, kind: str, filename: str):
    if kind == "xgboost":
        from xgboost import XGBClassifier

        model = XGBClassifier()
        model.load_model(folder / filename)
        return model
    return pickle.loads((folder / filename).read_bytes())


def load_variants(model_dir: Path, names: dict, kind: str, trained: dict) -> dict:
    """variant -> (function from features to risk, 90%-sensitivity threshold)."""
    m1, m2 = trained["models"][f"{kind}_E1"], trained["models"][f"{kind}_E2"]
    e1 = load_model(model_dir / names["E1"], kind, m1["model_file"])
    e2 = load_model(model_dir / names["E2"], kind, m2["model_file"])
    recal = m1["recalibration_to_E2"]
    platt = modeling.PlattMap(**recal["platt"])
    isotonic = pickle.loads((model_dir / names["E1"] / recal["isotonic_file"]).read_bytes())
    return {
        "E1": (lambda X: modeling.predict(e1, X), m1["threshold_90_sensitivity"]),
        "E1_platt": (lambda X: platt(modeling.predict(e1, X)), recal["threshold_90_sensitivity_platt"]),
        "E1_isotonic": (lambda X: isotonic.predict(modeling.predict(e1, X)), recal["threshold_90_sensitivity_isotonic"]),
        "E2": (lambda X: modeling.predict(e2, X), m2["threshold_90_sensitivity"]),
    }


def cohort_arrays(con, table: str) -> dict:
    columns = ", ".join(("unit_id", "patient_id", "label") + SUBGROUPS + features.FEATURES_48)
    a = con.execute(f"SELECT {columns} FROM {table} ORDER BY unit_id").fetchnumpy()
    return {
        "unit_id": a["unit_id"], "patient_id": a["patient_id"], "y": np.asarray(a["label"]).astype(int),
        "X": np.column_stack([np.asarray(a[c], dtype=np.float32) for c in features.FEATURES_48]),
        **{c: a[c] for c in SUBGROUPS},
    }


def build_test_cohorts(site: str, specs: dict) -> tuple[dict, dict]:
    """Test cohorts for every spec that applies at the site, and their flows."""
    cfg = get_site(site)
    loader = get_loader(site, connect(cfg.tmp_dir))
    built = loader.build()
    feature_tables, cohorts, flows = {}, {}, {}
    for key, spec in specs.items():
        try:
            spec.validate(site)
        except ValueError:
            continue
        if spec.interpretation not in feature_tables:
            feature_tables[spec.interpretation] = features.build_features(loader, built, spec.interpretation)
        result = build_cohort(loader, built, spec, partitions=("test",), features_table=feature_tables[spec.interpretation])
        cohorts[key] = cohort_arrays(loader.con, result.table)
        flows[key] = result.flow
    loader.con.close()
    return cohorts, flows


def table1(site: str, specs: dict) -> dict:
    """Covariate means for the full E1 and E2 cohorts, derivation and test together."""
    cfg = get_site(site)
    loader = get_loader(site, connect(cfg.tmp_dir))
    built = loader.build()
    feats = features.build_features(loader, built, "lab")
    moments = {}
    for estimand in ("E1", "E2"):
        result = build_cohort(loader, built, specs[estimand], partitions=("derivation", "test"), features_table=feats)
        exprs = ", ".join(f"AVG({c}), VAR_SAMP({c})" for c in TABLE1_COVARIATES)
        values = loader.con.execute(f"SELECT COUNT(*), {exprs} FROM {result.table}").fetchone()
        moments[estimand] = {"units": cell(values[0]),
                             **{c: (values[1 + 2 * i], values[2 + 2 * i]) for i, c in enumerate(TABLE1_COVARIATES)}}
    loader.con.close()
    out = {"units": {e: moments[e]["units"] for e in moments}, "covariates": {}}
    for c in TABLE1_COVARIATES:
        (m1, v1), (m2, v2) = moments["E1"][c], moments["E2"][c]
        pooled = math.sqrt((v1 + v2) / 2) if None not in (v1, v2) else 0
        out["covariates"][c] = {"E1_mean": round(m1, 4), "E2_mean": round(m2, 4),
                                "smd_E2_minus_E1": round((m2 - m1) / pooled, 3) if pooled else None}
    return out


def cache_predictions(folder: Path, name: str, arrays: dict, p: np.ndarray) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.register("p", pa.table({"unit_id": arrays["unit_id"], "patient_id": arrays["patient_id"], "label": arrays["y"], "prediction": p}))
    safe = name.replace(":", "_").replace(">", "").replace("-", "_to_")
    con.execute(f"COPY (SELECT * FROM p ORDER BY unit_id) TO '{folder / (safe + '.parquet')}' (FORMAT parquet, COMPRESSION zstd)")
    con.close()


def grid_cells(cohorts: dict, variants: dict, pairs, prefix: str = "") -> dict:
    cells = {}
    for kind, vs in variants.items():
        for variant, test in pairs:
            if test not in cohorts:
                continue
            fn, threshold = vs[variant]
            a = cohorts[test]
            cells[f"{prefix}{kind}:{variant}->{test}"] = {"cohort": f"{prefix}{test}", "p": fn(a["X"]), "y": a["y"],
                                                          "patients": a["patient_id"], "threshold": threshold, "arrays": a}
    return cells


def evaluate_site(cohorts: dict, variants: dict, n_boot: int, seed: int, cache_dir: Path | None = None) -> dict:
    """The grid, paired comparisons, curves, subgroups and sensitivity cohorts for one site."""
    cells = grid_cells(cohorts, variants, GRID)
    comparisons = [(f"{k}:{a}->{t}", f"{k}:{b}->{t}") for k in variants for a, b, t in COMPARISONS]
    grid, differences = evaluate_cells(cells, n_boot, seed, comparisons)
    for name, c in cells.items():
        grid[name]["curves"] = curves(c["p"], c["y"])
        grid[name]["subgroups"] = subgroups(np.asarray(c["p"], dtype=float), np.asarray(c["y"], dtype=float), c["threshold"], c["arrays"])
        if cache_dir is not None:
            cache_predictions(cache_dir, name, c["arrays"], c["p"])
    sensitivity = {}
    for key in sorted(k for k in cohorts if ":" in k):
        cells = grid_cells(cohorts, variants, [(v, key) for v in SENSITIVITY_VARIANTS])
        sensitivity[key], _ = evaluate_cells(cells, n_boot, seed)
    out = {"grid": grid, "comparisons": differences, "sensitivity": sensitivity,
           "specimen_grids": specimen_grids(cohorts, variants, n_boot, seed)}
    out["null_result"] = {k: null_result(out, k) for k in variants}
    return out


def restrict(cohorts: dict, column: str, value: str) -> dict:
    """The E1 and E2 test cohorts limited to units whose `column` equals `value`."""
    out = {}
    for key in ("E1", "E2"):
        if key not in cohorts:
            continue
        a = cohorts[key]
        rows = np.asarray(a[column]).astype(str) == value
        out[key] = {k: np.asarray(v)[rows] for k, v in a.items()}
    return out


def specimen_grids(cohorts: dict, variants: dict, n_boot: int, seed: int) -> dict:
    """The grid and its paired comparisons within each specimen type in SPECIMEN_GRIDS."""
    out = {}
    for specimen in SPECIMEN_GRIDS:
        cells = grid_cells(restrict(cohorts, "specimen_type", specimen), variants, GRID, prefix=f"{specimen}:")
        comparisons = [(f"{specimen}:{k}:{a}->{t}", f"{specimen}:{k}:{b}->{t}") for k in variants for a, b, t in COMPARISONS]
        grid, differences = evaluate_cells(cells, n_boot, seed, comparisons)
        out[specimen] = {"grid": grid, "comparisons": differences}
    return out


def conjunction(values) -> bool | None:
    """False if any value is False, None if none is False but some is unknown, else True."""
    values = list(values)
    if any(v is False for v in values):
        return False
    return None if any(v is None for v in values) else True


def _interval(v):
    return tuple(v) if isinstance(v, (list, tuple)) and len(v) == 2 else None


def null_result(site: dict, model: str = "xgboost") -> dict:
    """The plan's null result at one site, on the E2 test cohort.

    Negligible when, for the Platt-recalibrated E1 model against the E2 model, the difference in
    standardized net benefit is within NULL_MARGIN at every primary threshold with its 95% interval
    inside NULL_INTERVAL_MARGIN, and the recalibrated model's calibration slope interval includes 1.
    Every check must pass, so one failed check settles the verdict as not negligible. A withheld
    value leaves it undetermined only when no check has failed.
    """
    a, b = NULL_COMPARISON
    comparison = site["comparisons"].get(f"{model}:{a}->E2 minus {model}:{b}->E2")
    recalibrated = site["grid"].get(f"{model}:{a}->E2")
    if comparison is None or recalibrated is None:
        return {"model": model, "negligible": None, "checks": {}}
    checks = {}
    for t in metrics.NB_THRESHOLDS:
        key = f"standardized_net_benefit_{t:.2f}"
        point, interval = comparison["point"].get(key), _interval(comparison["interval"].get(key))
        checks[key] = None if not isinstance(point, float) or interval is None else bool(
            abs(point) <= NULL_MARGIN and -NULL_INTERVAL_MARGIN <= interval[0] and interval[1] <= NULL_INTERVAL_MARGIN)
    slope = _interval(recalibrated["interval"].get("calibration_slope"))
    checks["calibration_slope_includes_1"] = None if slope is None else bool(slope[0] <= 1 <= slope[1])
    return {"model": model, "comparison": f"{a} against {b} on E2", "negligible": conjunction(checks.values()), "checks": checks}


def evaluate_cross_site(per_site: dict, n_boot: int, seed: int) -> dict:
    """Site A's models on site B's test cohorts, for every ordered pair, resampling B's patients."""
    out = {}
    for source, (_, variants) in per_site.items():
        for target, (cohorts, _) in per_site.items():
            if source == target:
                continue
            cells = grid_cells({k: cohorts[k] for k in ("E1", "E2") if k in cohorts}, variants, CROSS_SITE, prefix=f"{target}:")
            out[f"{source} to {target}"], _ = evaluate_cells(cells, n_boot, seed)
    return out


# Report ------------------------------------------------------------------------------------------

def fmt(point, interval) -> str:
    if point == WITHHELD or point is None:
        return "withheld"
    lo, hi = interval if isinstance(interval, (list, tuple)) else (float("nan"), float("nan"))
    return f"{point:.3f} ({lo:.3f} to {hi:.3f})"


def render(summary: dict) -> str:
    lines = ["# 03_evaluate_ladder", "", "Test partition. Point estimates with 95% patient-clustered bootstrap intervals.", ""]
    shown = ["auroc", "calibration_intercept", "calibration_slope", "standardized_net_benefit_0.05",
             "standardized_net_benefit_0.10", "ppv_at_90_sensitivity"]
    for site, s in summary["sites"].items():
        lines += [f"## {site}", "", "| cell | units | " + " | ".join(shown) + " |", "|---|---|" + "---|" * len(shown)]
        for name, c in s["grid"].items():
            lines.append(f"| {name} | {c['units']} | " + " | ".join(fmt(c["point"].get(m), c["interval"].get(m)) for m in shown) + " |")
        lines += ["", "Paired differences on the E2 test cohort:", "", "| comparison | " + " | ".join(shown) + " |",
                  "|---|" + "---|" * len(shown)]
        for name, d in s["comparisons"].items():
            lines.append(f"| {name} | " + " | ".join(fmt(d["point"].get(m), d["interval"].get(m)) for m in shown) + " |")
        lines += ["", f"Null result (plan section 12): {s['null_result']}", ""]
    return "\n".join(lines) + "\n"


# Exploratory amendment, SAP section 17 --------------------------------------------------------------

DECISION_THRESHOLDS = (0.05, 0.10, 0.20)
R3_FLAG_THRESHOLD = 0.10  # R3 is a probability model, so it flags at the primary decision threshold
RULES = ("tumbarello", "R1_prior_esbl", "R2_prior_esbl_or_recent_blfq")


def flag_rate_at_sensitivity(p, y, sensitivity: float, w=None) -> float:
    """Orders flagged per 1,000 at the highest threshold whose sensitivity reaches `sensitivity`."""
    p, y = np.asarray(p, dtype=float), np.asarray(y).astype(bool)
    w = np.ones(len(p)) if w is None else np.asarray(w, dtype=float)
    positives = np.flatnonzero(y & (w > 0))
    if len(positives) == 0 or not np.isfinite(sensitivity):
        return float("nan")
    order = positives[np.argsort(-p[positives], kind="stable")]
    reached = np.cumsum(w[order]) / w[order].sum()
    k = min(int(np.searchsorted(reached, sensitivity - 1e-12)), len(order) - 1)
    return float(1000 * w[p >= p[order[k]]].sum() / w.sum())


def _rule_values(score, flags, y, ranked, w=None) -> dict:
    c = metrics.confusion(flags.astype(float), y, 0.5, w)
    n, positives, flagged = sum(c.values()), c["tp"] + c["fn"], c["tp"] + c["fp"]
    nan = float("nan")
    out = {
        "auroc": ranked.auroc(w),
        "flagged_per_1000": 1000 * flagged / n,
        "covered_per_1000": 1000 * c["tp"] / n,
        "missed_per_1000": 1000 * c["fn"] / n,
        "sensitivity": c["tp"] / positives if positives else nan,
        "specificity": c["tn"] / (c["tn"] + c["fp"]) if c["tn"] + c["fp"] else nan,
        "ppv": c["tp"] / flagged if flagged else nan,
        "npv": c["tn"] / (c["tn"] + c["fn"]) if c["tn"] + c["fn"] else nan,
    }
    for t in DECISION_THRESHOLDS:
        nb = (c["tp"] - c["fp"] * t / (1 - t)) / n
        out[f"net_benefit_{t:.2f}"] = nb
        out[f"standardized_net_benefit_{t:.2f}"] = nb / (positives / n) if positives else nan
    return out


def _draws(fn, patients, n_boot: int, seed: int, key: str):
    """Point value of fn(None), and 95% intervals over patient bootstrap draws shared by cohort key."""
    point = fn(None)
    keys = list(point)
    rows = [[d[k] for k in keys] for d in (fn(w) for w in metrics.patient_weights(patients, n_boot, cohort_seed(seed, key)))]
    draws = np.array(rows) if rows else np.empty((0, len(keys)))
    intervals = {k: metrics.percentile_interval(draws[:, i]) for i, k in enumerate(keys)} if len(draws) else {}
    return point, intervals


def _flag_withheld(flags, y) -> tuple[set, dict]:
    """Small-cell rule for a set of flags: everything, or every value but AUROC, or nothing."""
    events = int(np.sum(y))
    counts = {k: int(v) for k, v in metrics.confusion(np.asarray(flags, dtype=float), y, 0.5).items()}
    if small(len(y)) or small(events) or small(len(y) - events) or events in (0, len(y)):
        return {"*"}, {k: WITHHELD for k in counts}
    if any(small(v) for v in counts.values()):
        return {"withhold_all_but_auroc"}, {k: WITHHELD for k in counts}
    return set(), counts


def evaluate_rule(score, flags, y, patients, n_boot: int, seed: int, key: str) -> dict:
    """E-1 and E-2: AUROC of the score, and the rule's flags as a decision, with bootstrap intervals."""
    score, flags, y = np.asarray(score, dtype=float), np.asarray(flags, dtype=bool), np.asarray(y, dtype=float)
    ranked = metrics.Ranked(score, y)
    point, intervals = _draws(lambda w: _rule_values(score, flags, y, ranked, w), patients, n_boot, seed, key)
    withheld, counts = _flag_withheld(flags, y)
    if "withhold_all_but_auroc" in withheld:
        withheld = set(point) - {"auroc"}
    return {"units": cell(len(y)), "events": cell(int(y.sum())), "confusion": counts,
            "point": publish(point, withheld), "interval": publish(intervals, withheld)}


def matched_comparison(flags, p_model, y, patients, n_boot: int, seed: int, key: str) -> dict:
    """Orders flagged per 1,000 by a rule, and by a model at the rule's own sensitivity."""
    flags, y = np.asarray(flags, dtype=bool), np.asarray(y, dtype=float)

    def values(w):
        c = metrics.confusion(flags.astype(float), y, 0.5, w)
        n, positives = sum(c.values()), c["tp"] + c["fn"]
        sensitivity = c["tp"] / positives if positives else float("nan")
        rule_rate = 1000 * (c["tp"] + c["fp"]) / n
        model_rate = flag_rate_at_sensitivity(p_model, y, sensitivity, w)
        return {"rule_sensitivity": sensitivity, "rule_flagged_per_1000": rule_rate,
                "model_flagged_per_1000": model_rate, "difference_per_1000": rule_rate - model_rate}

    point, intervals = _draws(values, patients, n_boot, seed, key)
    withheld, _ = _flag_withheld(flags, y)
    withheld = set(point) if withheld else set()
    return {"point": publish(point, withheld), "interval": publish(intervals, withheld)}


def decision_impact(p, y, patients, threshold90, n_boot: int, seed: int, key: str) -> dict:
    """E-4: per 1,000 orders, broad coverage started and resistant cultures covered and missed."""
    p, y = np.asarray(p, dtype=float), np.asarray(y, dtype=float)
    labels = [(f"{t:.2f}", t) for t in DECISION_THRESHOLDS]
    if threshold90 is not None:
        labels.append(("90_sensitivity", float(threshold90)))

    def values(w):
        out = {}
        for label, t in labels:
            for k, v in metrics.per_thousand(p, y, t, w).items():
                out[f"{label}:{k}"] = v
        return out

    point, intervals = _draws(values, patients, n_boot, seed, key)
    out = {}
    for label, t in labels:
        hidden = any(small(v) for v in metrics.confusion(p, y, t).values()) or small(int(y.sum()))
        out[label] = {"threshold": t, "point": {}, "interval": {}}
        for k in ("flagged", "covered", "missed"):
            out[label]["point"][k] = WITHHELD if hidden else point[f"{label}:{k}"]
            out[label]["interval"][k] = WITHHELD if hidden else intervals.get(f"{label}:{k}")
    return out


def fit_minimal(path: Path):
    """R3: logistic regression on scores.MINIMAL_PREDICTORS, fitted on a derivation cohort file."""
    from sklearn.linear_model import LogisticRegression

    con = duckdb.connect()
    src = f"read_parquet('{str(path).replace(chr(39), chr(39) * 2)}')"
    partitions = [r[0] for r in con.execute(f"SELECT DISTINCT partition FROM {src}").fetchall()]
    if partitions != ["derivation"]:
        raise ValueError(f"{path} must hold the derivation partition only")
    a = con.execute(f"SELECT label, {', '.join(features.FEATURES_48)} FROM {src}").fetchnumpy()
    con.close()
    X = np.column_stack([np.asarray(a[c], dtype=np.float32) for c in features.FEATURES_48])
    y = np.asarray(a["label"]).astype(int)
    D = scores.minimal_design(scores.columns(X))
    model = LogisticRegression(C=1e6, max_iter=5000).fit(D, y)  # effectively unpenalized
    threshold = modeling.sensitivity_threshold(model.predict_proba(D)[:, 1], y)
    coefficients = dict(zip(("intercept",) + scores.MINIMAL_PREDICTORS,
                            [float(model.intercept_[0])] + [float(c) for c in model.coef_[0]]))
    return model, threshold, coefficients


def exploratory_site(cohorts: dict, variants: dict, urine_variants: dict, derivation_e2: Path,
                     n_boot: int, seed: int) -> dict:
    """E-1 to E-4 for one site. Draws are shared by cohort key, so rules and models on one cohort pair up."""
    out = {"rules": {}, "matched": {}, "minimal_model": {}, "urine_models": {}, "decision_impact": {}}
    tests = [t for t in ("E1", "E2") if t in cohorts]
    xgb_e2 = variants["xgboost"]["E2"][0]
    for test in tests:
        a = cohorts[test]
        f = scores.columns(a["X"])
        points = scores.tumbarello_points(f)
        rules = {
            "tumbarello": (points, points >= scores.TUMBARELLO_CUTOFF),
            "R1_prior_esbl": (scores.prior_esbl(f).astype(float), scores.prior_esbl(f)),
            "R2_prior_esbl_or_recent_blfq": (scores.prior_esbl_or_recent_blfq(f).astype(float), scores.prior_esbl_or_recent_blfq(f)),
        }
        p_model = xgb_e2(a["X"])
        for name, (score, flags) in rules.items():
            out["rules"][f"{name}->{test}"] = evaluate_rule(score, flags, a["y"], a["patient_id"], n_boot, seed, test)
            out["matched"][f"{name}->{test}"] = matched_comparison(flags, p_model, a["y"], a["patient_id"], n_boot, seed, test)

    model, threshold, coefficients = fit_minimal(derivation_e2)
    cells = {}
    for test in tests:
        a = cohorts[test]
        p = model.predict_proba(scores.minimal_design(scores.columns(a["X"])))[:, 1]
        cells[f"R3->{test}"] = {"cohort": test, "p": p, "y": a["y"], "patients": a["patient_id"], "threshold": threshold}
        flags = p >= R3_FLAG_THRESHOLD
        out["matched"][f"R3_at_{R3_FLAG_THRESHOLD:.2f}->{test}"] = matched_comparison(
            flags, variants["xgboost"]["E2"][0](a["X"]), a["y"], a["patient_id"], n_boot, seed, test)
    grid, _ = evaluate_cells(cells, n_boot, seed)
    out["minimal_model"] = {"coefficients": coefficients, "threshold_90_sensitivity": threshold, "grid": grid}

    urine = restrict(cohorts, "specimen_type", "urine")
    ucells = grid_cells(urine, urine_variants, GRID, prefix="urine:")
    comparisons = [(f"urine:{k}:{a}->{t}", f"urine:{k}:{b}->{t}") for k in urine_variants for a, b, t in COMPARISONS]
    ugrid, udiff = evaluate_cells(ucells, n_boot, seed, comparisons)
    out["urine_models"] = {"grid": ugrid, "comparisons": udiff}

    impact_cells = {**grid_cells(cohorts, variants, GRID), **cells, **ucells}
    for name, c in impact_cells.items():
        out["decision_impact"][name] = decision_impact(c["p"], c["y"], c["patients"], c["threshold"], n_boot, seed, c["cohort"])
    return out


def render_exploratory(summary: dict) -> str:
    lines = ["# 03_exploratory", "", "Exploratory amendment E-1 to E-4 (SAP section 17). Test partition, not pre-specified.", ""]

    def v(d, k, digits=3):
        x = d.get(k)
        return f"{x:.{digits}f}" if isinstance(x, float) else str(x)

    for site, s in summary["sites"].items():
        lines += [f"## {site}", "", "Rules and scores on the E2 test cohort:", "",
                  "| rule | AUROC | flagged per 1,000 | sensitivity | PPV | standardized NB at 0.10 |", "|---|---|---|---|---|---|"]
        for name, r in s["rules"].items():
            if name.endswith("->E2"):
                pt = r["point"]
                lines.append(f"| {name} | {v(pt, 'auroc')} | {v(pt, 'flagged_per_1000', 0)} | {v(pt, 'sensitivity', 2)} | "
                             f"{v(pt, 'ppv', 3)} | {v(pt, 'standardized_net_benefit_0.10', 2)} |")
        lines += ["", "The E2-trained XGBoost model at each rule's own sensitivity:", "",
                  "| rule | rule sensitivity | rule flags per 1,000 | model flags per 1,000 | difference |", "|---|---|---|---|---|"]
        for name, m in s["matched"].items():
            if name.endswith("->E2"):
                pt = m["point"]
                lines.append(f"| {name} | {v(pt, 'rule_sensitivity', 2)} | {v(pt, 'rule_flagged_per_1000', 0)} | "
                             f"{v(pt, 'model_flagged_per_1000', 0)} | {v(pt, 'difference_per_1000', 0)} |")
        lines += ["", "Urine-specific XGBoost models on urine test units:", "",
                  "| cell | AUROC | calibration intercept | standardized NB at 0.10 |", "|---|---|---|---|"]
        for name, c in s["urine_models"]["grid"].items():
            if ":xgboost:" in name:
                pt = c["point"]
                lines.append(f"| {name} | {v(pt, 'auroc')} | {v(pt, 'calibration_intercept', 2)} | {v(pt, 'standardized_net_benefit_0.10', 2)} |")
        lines.append("")
    return "\n".join(lines) + "\n"


def run_exploratory(args, frozen: dict, train_dir: Path) -> int:
    trained = json.loads((train_dir / "summary.json").read_text())
    train_args = json.loads((train_dir / "manifest.json").read_text())["args"]
    urine_dir = Path(args.urine_train_run) if args.urine_train_run else latest_run("02_train_e1", "urine")
    urine_trained = json.loads((urine_dir / "summary.json").read_text())
    urine_args = json.loads((urine_dir / "manifest.json").read_text())["args"]
    rule, unit, kinds = train_args["repeat_rule"], train_args["unit"], train_args["models"].split(",")
    if (urine_args["repeat_rule"], urine_args["unit"]) != (rule, unit):
        raise SystemExit("the urine training run must use the primary repeat rule and unit")
    sites = args.sites.split(",") if args.sites else list(trained["sites"])

    run = RunManifest("03_exploratory", seed=args.seed, args=vars(args))
    run.set("analysis_plan", frozen)
    run.set("amendment", "SAP section 17, exploratory analyses E-1 to E-5")
    for label, d in (("02_train_e1.summary", train_dir), ("02_train_e1_urine.summary", urine_dir)):
        run.add_input(label, d / "summary.json", sha256_file(d / "summary.json"))
    base = {"repeat_rule": rule, "unit": unit}
    specs = {"E1": CohortSpec("E1", **base), "E2": CohortSpec("E2", **base)}
    names = {e: specs[e].name for e in specs}
    summary = {"partition": "test", "exploratory": True, "amendment": "SAP section 17, E-1 to E-5", "analysis_plan": frozen,
               "train_run": str(train_dir), "urine_train_run": str(urine_dir),
               "tumbarello": {"cutoff": scores.TUMBARELLO_CUTOFF, "max": scores.TUMBARELLO_MAX},
               "r3_flag_threshold": R3_FLAG_THRESHOLD, "sites": {}}
    for site in sites:
        print(f"exploratory: {site}", flush=True)
        root = get_site(site).parquet_dir.parent
        variants = {k: load_variants(root / "models" / f"{rule}_{unit}", names, k, trained["sites"][site]) for k in kinds}
        urine_variants = {k: load_variants(root / "models" / f"{rule}_{unit}_urine", names, k, urine_trained["sites"][site])
                          for k in kinds}
        derivation = root / "cohorts" / names["E2"] / "derivation.parquet"
        run.add_input(f"{site}.{names['E2']}.derivation", derivation, sha256_file(derivation))
        cohorts, flows = build_test_cohorts(site, specs)
        for key, flow in flows.items():
            for step in flow:
                run.count(f"{site}.{key}.{step['step']}", step["cultures"])
        summary["sites"][site] = exploratory_site(cohorts, variants, urine_variants, derivation, args.boot, args.seed)
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (run.dir / "report.md").write_text(render_exploratory(summary))
    print(f"exploratory: wrote {run.dir}")
    return 0


def latest_run(script: str, specimen: str | None = None) -> Path:
    """The latest run of `script` with a summary. For training runs, the one trained on `specimen`
    (None is every specimen type), read from its manifest."""
    def matches(p: Path) -> bool:
        if not (p / "summary.json").exists():
            return False
        manifest = p / "manifest.json"
        return not manifest.exists() or json.loads(manifest.read_text())["args"].get("specimen") == specimen
    runs = sorted(p for p in (RESULTS_DIR / script).iterdir() if matches(p))
    if not runs:
        raise SystemExit(f"no {script} run with a summary.json" + (f" for {specimen}" if specimen else ""))
    return runs[-1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--train-run", help="results/02_train_e1/<timestamp> (default: the latest)")
    parser.add_argument("--sites", help="comma-separated sites (default: every site in the training run)")
    parser.add_argument("--boot", type=int, default=N_BOOT)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--exploratory", action="store_true", help="run the SAP's exploratory amendment E-1 to E-4")
    parser.add_argument("--urine-train-run", help="the urine-only training run for E-3 (default: the latest)")
    args = parser.parse_args()
    frozen = require_frozen_plan()
    train_dir = Path(args.train_run) if args.train_run else latest_run("02_train_e1")
    if args.exploratory:
        return run_exploratory(args, frozen, train_dir)
    trained = json.loads((train_dir / "summary.json").read_text())
    train_args = json.loads((train_dir / "manifest.json").read_text())["args"]
    rule, unit, kinds = train_args["repeat_rule"], train_args["unit"], train_args["models"].split(",")
    sites = args.sites.split(",") if args.sites else list(trained["sites"])

    run = RunManifest("03_evaluate_ladder", seed=args.seed, args=vars(args))
    run.set("analysis_plan", frozen)
    run.add_input("02_train_e1.summary", train_dir / "summary.json", sha256_file(train_dir / "summary.json"))
    base = {"repeat_rule": rule, "unit": unit}
    specs = {"E1": CohortSpec("E1", **base), "E2": CohortSpec("E2", **base)}
    for key, overrides in SENSITIVITY.items():
        for estimand in ("E1", "E2"):
            specs[f"{key}:{estimand}"] = CohortSpec(estimand, **{**base, **overrides})

    summary = {"partition": "test", "analysis_plan": frozen, "train_run": str(train_dir), "sites": {}}
    per_site = {}
    for site in sites:
        print(f"evaluate: {site}", flush=True)
        model_dir = get_site(site).parquet_dir.parent / "models" / f"{rule}_{unit}"
        names = {e: specs[e].name for e in ("E1", "E2")}
        variants = {k: load_variants(model_dir, names, k, trained["sites"][site]) for k in kinds}
        cohorts, flows = build_test_cohorts(site, specs)
        for key, flow in flows.items():
            for step in flow:
                run.count(f"{site}.{key}.{step['step']}", step["cultures"])
        result = evaluate_site(cohorts, variants, args.boot, args.seed, cache_dir=model_dir / "test")
        result["flows"] = publish_flows(flows)
        result["table1"] = table1(site, specs)
        summary["sites"][site] = result
        per_site[site] = ({k: cohorts[k] for k in ("E1", "E2") if k in cohorts}, variants)
    summary["cross_site"] = evaluate_cross_site(per_site, args.boot, args.seed)
    verdicts = [s["null_result"]["xgboost"]["negligible"] for s in summary["sites"].values() if "xgboost" in s["null_result"]]
    summary["null_result_overall"] = conjunction(verdicts) if verdicts else None
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (run.dir / "report.md").write_text(render(summary))
    print(f"evaluate: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
