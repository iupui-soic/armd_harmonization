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

Usage: python analysis/03_evaluate_ladder.py [--train-run results/02_train_e1/<timestamp>]
       [--sites stanford,mgb] [--boot 2000] [--seed 20260911]
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

from armd_common import features, metrics, modeling
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


def _interval(v):
    return tuple(v) if isinstance(v, (list, tuple)) and len(v) == 2 else None


def null_result(site: dict, model: str = "xgboost") -> dict:
    """The plan's null result at one site, on the E2 test cohort.

    Negligible when, for the Platt-recalibrated E1 model against the E2 model, the difference in
    standardized net benefit is within NULL_MARGIN at every primary threshold with its 95% interval
    inside NULL_INTERVAL_MARGIN, and the recalibrated model's calibration slope interval includes 1.
    A withheld value leaves the verdict undetermined.
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
    verdict = None if any(v is None for v in checks.values()) else all(checks.values())
    return {"model": model, "comparison": f"{a} against {b} on E2", "negligible": verdict, "checks": checks}


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


def latest_run(script: str) -> Path:
    runs = sorted(p for p in (RESULTS_DIR / script).iterdir() if (p / "summary.json").exists())
    if not runs:
        raise SystemExit(f"no {script} run with a summary.json")
    return runs[-1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--train-run", help="results/02_train_e1/<timestamp> (default: the latest)")
    parser.add_argument("--sites", help="comma-separated sites (default: every site in the training run)")
    parser.add_argument("--boot", type=int, default=N_BOOT)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    frozen = require_frozen_plan()
    train_dir = Path(args.train_run) if args.train_run else latest_run("02_train_e1")
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
    summary["null_result_overall"] = None if not verdicts or None in verdicts else all(verdicts)
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (run.dir / "report.md").write_text(render(summary))
    print(f"evaluate: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
