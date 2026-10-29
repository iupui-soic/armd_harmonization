"""Train the OrderGap models on MIMIC-IV derivation units (SAP amendment E-9, retrained grid).

This calls analysis/02_train_e1.py's train_site unchanged, so MIMIC models use the frozen settings of
SAP section 6: XGBoost and logistic models under E1 and E2, cross-fitted predictions on the frozen
inner folds, Platt and isotonic maps from E1 toward E2, and 90%-sensitivity thresholds. With
--specimen the models are trained on urine or blood units only, as in E-3 and E-7.

Models and predictions stay under /data0/mimiciv/3.1/ordergap/models. results/08_mimic_train/<timestamp>/
gets manifest.json and summary.json.

Usage: python analysis/08_mimic_train.py [--specimen urine|blood] [--threads 10] [--seed 20260911]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from armd_common import modeling  # noqa: E402
from armd_common.config import DEFAULT_SEED  # noqa: E402
from armd_common.provenance import RunManifest  # noqa: E402
from mimic import build  # noqa: E402


def load_train_module():
    spec = importlib.util.spec_from_file_location("train_e1", Path(__file__).with_name("02_train_e1.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--specimen", choices=("urine", "blood"))
    parser.add_argument("--models", default=",".join(modeling.MODEL_KINDS))
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    train = load_train_module()
    frozen = train.require_frozen_plan()
    run = RunManifest("08_mimic_train", seed=args.seed, args=vars(args))
    run.set("analysis_plan", frozen)
    run.set("xgboost_params", modeling.XGB_PARAMS)
    run.set("logistic_C", modeling.LOGISTIC_C)
    folder = f"{build.REPEAT_RULE}_{build.UNIT}" + (f"_{args.specimen}" if args.specimen else "")
    print(f"train mimic {folder}", flush=True)
    out = train.train_site(build.SITE, build.WORK / "cohorts", build.WORK / "models" / folder, build.REPEAT_RULE, build.UNIT,
                           args.models.split(","), args.seed, args.threads, run, specimen=args.specimen)
    summary = {"partition": "derivation", "analysis_plan": frozen, "specimen": args.specimen, "model_folder": folder,
               "sites": {build.SITE: out}}
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"train mimic: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
