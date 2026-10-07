"""Train the OrderGap models on the derivation partition (plan section 5.3, SAP sections 6 and 7.3).

For each site, each model class and each training estimand (E1 and E2) this fits:
- cross-fitted predictions on the frozen inner folds, so every derivation unit is predicted by a
  model that never saw its patient;
- one final model on every derivation unit;
- the 90%-sensitivity threshold from the training cohort's cross-fitted predictions.

For E1-trained models it also predicts the E2 derivation cohort fold by fold with the same fold
models, then fits Platt (primary) and isotonic (secondary) maps toward E2 on those predictions. The
thresholds on the recalibrated scales come from the same predictions.

Inputs are the cohort Parquet files from analysis/01_build_cohorts.py. Models, predictions and
isotonic maps are row-level, so they go under models/<cohort spec>/ in the site folder of the data
root and never into the repository. results/02_train_e1/<timestamp>/ receives manifest.json and
summary.json with counts, thresholds and Platt parameters. No performance metric is computed here.

Two guards:
- The script refuses to train until the analysis plan is frozen, that is until the repository has
  the git tag sap-v1.0. The manifest records the SHA-256 of the local SAP.md it ran under.
- It reads derivation rows only and fails if a cohort file holds any other partition.

With --specimen, every model is trained and recalibrated on that specimen type's units only, and the
outputs go under models/<rule>_<unit>_<specimen>/ (the SAP's exploratory analyses E-3 for urine and
E-7 for blood).

Usage: python analysis/02_train_e1.py [--sites stanford,mgb] [--models xgboost,logistic]
       [--repeat-rule release] [--unit culture] [--specimen urine] [--threads 8] [--seed 20260911]
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import duckdb
import numpy as np
import pyarrow as pa

from armd_common import modeling
from armd_common.cohorts import CohortSpec
from armd_common.config import DEFAULT_SEED, REPO_ROOT, SITES, get_site
from armd_common.episodes import REPEAT_RULES, UNITS
from armd_common.features import FEATURES_48
from armd_common.provenance import RunManifest, require_git_tag, sha256_file
from armd_common.qc import cell
from armd_common.splits import SealedPartitionError

FROZEN_TAG = "sap-v1.0"
ESTIMANDS = ("E1", "E2")
SPECIMENS = ("urine", "blood", "respiratory")


def require_frozen_plan() -> dict:
    return require_git_tag(FROZEN_TAG, REPO_ROOT / "SAP.md")


def cohort_names(repeat_rule: str, unit: str) -> dict:
    return {e: CohortSpec(e, repeat_rule=repeat_rule, unit=unit).name for e in ESTIMANDS}


def load_cohort(path: Path, specimen: str | None = None) -> dict:
    """Derivation units of one cohort as arrays, optionally one specimen type's. Fails on any other partition."""
    con = duckdb.connect()
    src = f"read_parquet('{str(path).replace(chr(39), chr(39) * 2)}')"
    partitions = [r[0] for r in con.execute(f"SELECT DISTINCT partition FROM {src}").fetchall()]
    if partitions != ["derivation"]:
        raise SealedPartitionError(f"{path} holds partitions {partitions}; training reads derivation only")
    missing = con.execute(f"SELECT COUNT(*) - COUNT(label), COUNT(*) - COUNT(fold) FROM {src}").fetchone()
    if any(missing):
        raise ValueError(f"{path} has units without a label or an inner fold")
    columns = ", ".join(FEATURES_48)
    where = "" if specimen is None else f"WHERE specimen_type = '{specimen}'"
    arrays = con.execute(f"SELECT unit_id, patient_id, fold, label, {columns} FROM {src} {where} ORDER BY unit_id").fetchnumpy()
    con.close()
    return {
        "unit_id": arrays["unit_id"],
        "patient_id": arrays["patient_id"],
        "fold": np.asarray(arrays["fold"]).astype(int),
        "y": np.asarray(arrays["label"]).astype(int),
        "X": np.column_stack([np.asarray(arrays[c], dtype=np.float32) for c in FEATURES_48]),
    }


def write_predictions(path: Path, data: dict, prediction: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.register("p", pa.table({"unit_id": data["unit_id"], "patient_id": data["patient_id"], "fold": data["fold"],
                                "label": data["y"], "prediction": prediction}))
    con.execute(f"COPY (SELECT * FROM p ORDER BY unit_id) TO '{path}' (FORMAT parquet, COMPRESSION zstd)")
    con.close()


def save_model(model, kind: str, path_stem: Path) -> Path:
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    if kind == "xgboost":
        path = path_stem.with_suffix(".ubj")
        model.save_model(path)
    else:
        path = path_stem.with_suffix(".pkl")
        path.write_bytes(pickle.dumps(model))
    return path


def train_site(site: str, cohort_dir: Path, out_dir: Path, repeat_rule: str, unit: str, models, seed: int,
               threads: int, run: RunManifest | None = None, params: dict | None = None,
               specimen: str | None = None) -> dict:
    """Train every model for one site, on one specimen type's units if `specimen` is given.
    `params` overrides the XGBoost hyperparameters, for tests only."""
    if specimen is not None and specimen not in SPECIMENS:
        raise ValueError(f"specimen must be one of {SPECIMENS}")
    names = cohort_names(repeat_rule, unit)
    data = {}
    for estimand, name in names.items():
        path = cohort_dir / name / "derivation.parquet"
        if run is not None:
            run.add_input(f"{site}.{name}", path, sha256_file(path))
        data[estimand] = load_cohort(path, specimen)
        if run is not None:
            run.count(f"{site}.{name}.units" + (f".{specimen}" if specimen else ""), len(data[estimand]["y"]))

    out = {"specimen": specimen,
           "cohorts": {e: {"name": n, "units": cell(len(data[e]["y"])), "label_positive": cell(int(data[e]["y"].sum()))}
                       for e, n in names.items()}, "models": {}}
    for kind in models:
        def factory(kind=kind):
            return modeling.make_model(kind, seed, params, threads)

        for estimand in ESTIMANDS:
            d, name = data[estimand], names[estimand]
            fit = modeling.cross_fit(factory, d["X"], d["y"], d["fold"])
            final = factory().fit(d["X"], d["y"])
            stem = out_dir / name / kind
            write_predictions(stem.with_name(f"{kind}_crossfit.parquet"), d, fit.predictions)
            model_path = save_model(final, kind, stem)
            entry = {
                "trained_on": name,
                "model_file": model_path.name,
                "threshold_90_sensitivity": modeling.sensitivity_threshold(fit.predictions, d["y"]),
            }
            if estimand == "E1":
                e2 = data["E2"]
                on_e2 = fit.predict_by_fold(e2["X"], e2["fold"])
                write_predictions(stem.with_name(f"{kind}_crossfit_on_{names['E2']}.parquet"), e2, on_e2)
                platt = modeling.fit_platt(on_e2, e2["y"])
                isotonic = modeling.fit_isotonic(on_e2, e2["y"])
                (stem.with_name(f"{kind}_isotonic_to_E2.pkl")).write_bytes(pickle.dumps(isotonic))
                entry["recalibration_to_E2"] = {
                    "fitted_on": names["E2"],
                    "platt": {"intercept": platt.intercept, "slope": platt.slope},
                    "isotonic_file": f"{kind}_isotonic_to_E2.pkl",
                    "threshold_90_sensitivity_uncalibrated": modeling.sensitivity_threshold(on_e2, e2["y"]),
                    "threshold_90_sensitivity_platt": modeling.sensitivity_threshold(platt(on_e2), e2["y"]),
                    "threshold_90_sensitivity_isotonic": modeling.sensitivity_threshold(isotonic.predict(on_e2), e2["y"]),
                }
            out["models"][f"{kind}_{estimand}"] = entry
            print(f"  {site}.{kind}.{estimand}: trained on {cell(len(d['y']))} units", flush=True)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", help="comma-separated sites (default: every site with built cohorts)")
    parser.add_argument("--models", default=",".join(modeling.MODEL_KINDS))
    parser.add_argument("--repeat-rule", default="release", choices=REPEAT_RULES)
    parser.add_argument("--unit", default="culture", choices=UNITS)
    parser.add_argument("--specimen", choices=SPECIMENS, help="train on one specimen type only (exploratory E-3 and E-7)")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    frozen = require_frozen_plan()
    models = args.models.split(",")
    unknown = set(models) - set(modeling.MODEL_KINDS)
    if unknown:
        parser.error(f"unknown models: {sorted(unknown)}")

    names = cohort_names(args.repeat_rule, args.unit)
    root = {s: get_site(s).parquet_dir.parent for s in (args.sites.split(",") if args.sites else SITES)}
    sites = [s for s, r in root.items() if all((r / "cohorts" / n / "derivation.parquet").exists() for n in names.values())]
    run = RunManifest("02_train_e1", seed=args.seed, args=vars(args))
    run.set("analysis_plan", frozen)
    run.set("xgboost_params", modeling.XGB_PARAMS)
    run.set("logistic_C", modeling.LOGISTIC_C)
    summary = {"partition": "derivation", "analysis_plan": frozen, "specimen": args.specimen, "sites": {}}
    folder = f"{args.repeat_rule}_{args.unit}" + (f"_{args.specimen}" if args.specimen else "")
    for site in sites:
        print(f"train: {site}", flush=True)
        summary["sites"][site] = train_site(site, root[site] / "cohorts", root[site] / "models" / folder,
                                            args.repeat_rule, args.unit, models, args.seed, args.threads, run,
                                            specimen=args.specimen)
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"train: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
