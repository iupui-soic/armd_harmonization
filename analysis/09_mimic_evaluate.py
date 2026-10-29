"""Evaluate OrderGap at Beth Israel Deaconess with MIMIC-IV (SAP amendment E-9 and its addendum).

This is the only script that reads the MIMIC-IV test partition, and it reads it once. It reuses the
evaluation functions of analysis/03_evaluate_ladder.py on MIMIC arrays:

- retrained grid: the primary grid, urine and blood grids, null criterion (section 7.2, 12);
- E-1 to E-4 and E-6 to E-8 at BIDMC: the ARMD approximation of the Tumbarello score, simple rules,
  the five-input model, urine-only and blood-only models, thresholds within specimen type, the
  specimen-only model, the ceftriaxone outcome and decision counts;
- the Tumbarello score from exact MIMIC items, at 3 and 6 points;
- portability: frozen ARMD models applied unchanged, and the same models recalibrated on MIMIC
  derivation units, against the retrained MIMIC model on the same test units;
- clinicians' empiric choices from eMAR, scored like a model, with paired differences against models;
- subgroups: orders treated at the time of the culture, urine by setting, insurance and race;
- Table 1 for the full MIMIC cohorts.

Results go to results/09_mimic_evaluate/<timestamp>/ as manifest.json, summary.json and report.md,
aggregate values only, with counts from 1 to 10 withheld.

Usage: python analysis/09_mimic_evaluate.py [--boot 2000] [--seed 20260911]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import duckdb
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from armd_common import metrics, modeling, scores  # noqa: E402
from armd_common.config import DEFAULT_SEED, REPO_ROOT, get_site  # noqa: E402
from armd_common.features import FEATURES_48  # noqa: E402
from armd_common.provenance import RunManifest, require_git_tag, sha256_file  # noqa: E402
from armd_common.qc import cell  # noqa: E402
from mimic import build  # noqa: E402

ARMD_SITES = ("stanford", "utsw", "ecuh", "mgb")
ARMD_RUNS = {None: "results/02_train_e1/20260927T113507Z", "urine": "results/02_train_e1/20260927T190025Z",
             "blood": "results/02_train_e1/20261007T101226Z"}
POLICIES = ("clin_carbapenem", "clin_carbapenem_started_after", "clin_broad")
SPECIMEN_THRESHOLDS = {"all": (0.10, 0.20), "urine": (0.10, 0.20), "blood": (0.01, 0.02)}
FAIRNESS = ("insurance_group", "race_group")


def load_evaluate_module():
    spec = importlib.util.spec_from_file_location("evaluate_ladder", Path(__file__).with_name("03_evaluate_ladder.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ev = load_evaluate_module()
COLUMNS = ("unit_id", "patient_id", "label", "specimen_type", "setting_primary", "prior_culture_history", "adi_quintile") \
    + build.EXTRA_FLAGS + build.EXTRA_CATEGORIES


def arrays_from(con, source: str) -> dict:
    a = con.execute(f"SELECT {', '.join(COLUMNS + FEATURES_48)} FROM {source} ORDER BY unit_id").fetchnumpy()
    out = {"unit_id": a["unit_id"], "patient_id": a["patient_id"], "y": np.asarray(a["label"]).astype(int),
           "X": np.column_stack([np.asarray(a[c], dtype=np.float32) for c in FEATURES_48])}
    for c in COLUMNS[3:]:
        out[c] = np.asarray(a[c])
    return out


def latest(script: str, specimen: str | None) -> Path:
    runs = sorted((REPO_ROOT / "results" / script).glob("*/summary.json"))
    for path in reversed(runs):
        if json.loads(path.read_text()).get("specimen") == specimen:
            return path.parent
    raise SystemExit(f"no {script} run for specimen {specimen}")


def subset(a: dict, rows) -> dict:
    return {k: (v[rows] if isinstance(v, np.ndarray) and len(v) == len(rows) else v) for k, v in a.items()}


def mask(a: dict, column: str, values) -> np.ndarray:
    return np.isin(np.asarray(a[column]).astype(str), [str(v) for v in values])


# Portability --------------------------------------------------------------------------------------

def local_platt(fn, derivation: dict):
    """Logistic recalibration (intercept and slope) of a transported model on MIMIC derivation units."""
    p = fn(derivation["X"])
    platt = modeling.fit_platt(p, derivation["y"])
    threshold = modeling.sensitivity_threshold(platt(p), derivation["y"])
    return (lambda X, f=fn, m=platt: m(f(X))), threshold, {"intercept": platt.intercept, "slope": platt.slope}


def transport(test: dict, derivation: dict, retrained: dict, armd: dict, n_boot: int, seed: int, key: str,
              thresholds=ev.DECISION_THRESHOLDS) -> dict:
    """ARMD models unchanged and recalibrated on MIMIC, against the retrained MIMIC model, same test units."""
    cells, maps = {}, {}
    fn, t = retrained["E2"]
    cells["mimic:E2"] = {"cohort": key, "p": fn(test["X"]), "y": test["y"], "patients": test["patient_id"], "threshold": t}
    for site, variants in armd.items():
        for variant, (fn, t) in variants.items():
            if variant not in ("E1", "E1_platt", "E2"):
                continue
            cells[f"{site}:{variant}"] = {"cohort": key, "p": fn(test["X"]), "y": test["y"], "patients": test["patient_id"], "threshold": t}
        for variant in ("E1", "E2"):
            if variant in variants:
                lfn, lt, maps[f"{site}:{variant}_local"] = local_platt(variants[variant][0], derivation)
                cells[f"{site}:{variant}_local"] = {"cohort": key, "p": lfn(test["X"]), "y": test["y"], "patients": test["patient_id"],
                                                    "threshold": lt}
    comparisons = [(name, "mimic:E2") for name in cells if name != "mimic:E2"]
    grid, differences = ev.evaluate_cells(cells, n_boot, seed, comparisons)
    impact = {name: ev.decision_impact(c["p"], c["y"], c["patients"], c["threshold"], n_boot, seed, c["cohort"], thresholds=thresholds)
              for name, c in cells.items()}
    return {"grid": grid, "comparisons": differences, "local_maps": maps, "decision_impact": impact}


# Clinicians' choices --------------------------------------------------------------------------------

def _values(flags: np.ndarray, y: np.ndarray, t: float, w) -> dict:
    c = metrics.confusion(flags.astype(float), y, 0.5, w)
    n, pos = sum(c.values()), c["tp"] + c["fn"]
    nb = (c["tp"] - c["fp"] * t / (1 - t)) / n
    return {"flagged": 1000 * (c["tp"] + c["fp"]) / n, "covered": 1000 * c["tp"] / n, "missed": 1000 * c["fn"] / n,
            "sensitivity": c["tp"] / pos if pos else float("nan"), "standardized_net_benefit": nb / (pos / n) if pos else float("nan")}


def clinician_comparison(a: dict, models: dict, thresholds, n_boot: int, seed: int, key: str) -> dict:
    """Policies from eMAR and models at each threshold on the same units, with paired differences.

    models maps a name to predicted risks. A policy is a fixed set of flags, so its net benefit is
    computed at each threshold's exchange rate. Differences are model minus policy, per draw.
    """
    y = np.asarray(a["y"], dtype=float)
    deciders = {p: np.asarray(a[p]).astype(int) > 0 for p in POLICIES}
    deciders.update({m: None for m in models})

    def values(w):
        out = {}
        for t in thresholds:
            per = {}
            for name in deciders:
                flags = deciders[name] if deciders[name] is not None else np.asarray(models[name]) >= t
                per[name] = _values(flags, y, t, w)
                for k, v in per[name].items():
                    out[f"{t:.2f}:{name}:{k}"] = v
            for m in models:
                for p in POLICIES:
                    out[f"{t:.2f}:{m} minus {p}:standardized_net_benefit"] = per[m]["standardized_net_benefit"] - per[p]["standardized_net_benefit"]
                    out[f"{t:.2f}:{m} minus {p}:flagged"] = per[m]["flagged"] - per[p]["flagged"]
        return out

    point, intervals = ev._draws(values, a["patient_id"], n_boot, seed, key)
    hidden = set()
    for t in thresholds:
        for name in deciders:
            flags = deciders[name] if deciders[name] is not None else np.asarray(models[name]) >= t
            if any(ev.small(v) for v in metrics.confusion(flags.astype(float), y, 0.5).values()):
                hidden.add((t, name))
    out = {"units": cell(len(y)), "events": cell(int(y.sum())), "thresholds": {}}
    for t in thresholds:
        block = {}
        for k, v in point.items():
            if not k.startswith(f"{t:.2f}:"):
                continue
            name, metric = k.split(":", 2)[1:]
            involved = [n for n in deciders if name == n or name.startswith(n + " minus ") or name.endswith(" minus " + n)]
            withheld = any((t, n) in hidden for n in involved)
            block.setdefault(name, {})[metric] = {"point": ev.WITHHELD if withheld else v,
                                                  "interval": ev.WITHHELD if withheld else intervals.get(k)}
        out["thresholds"][f"{t:.2f}"] = block
    return out


# Subgroups -----------------------------------------------------------------------------------------

def grid_on(a: dict, variants: dict, key: str, n_boot: int, seed: int, thresholds=ev.DECISION_THRESHOLDS) -> dict:
    cells = {f"xgboost:{v}": {"cohort": key, "p": variants[v][0](a["X"]), "y": a["y"], "patients": a["patient_id"],
                              "threshold": variants[v][1]} for v in ("E1", "E1_platt", "E2") if v in variants}
    comparisons = [("xgboost:E1_platt", "xgboost:E2"), ("xgboost:E1", "xgboost:E2")]
    grid, differences = ev.evaluate_cells(cells, n_boot, seed, comparisons)
    impact = {n: ev.decision_impact(c["p"], c["y"], c["patients"], c["threshold"], n_boot, seed, key, thresholds=thresholds)
              for n, c in cells.items()}
    return {"grid": grid, "comparisons": differences, "decision_impact": impact}


def fairness(a: dict, variants: dict) -> dict:
    saved = ev.SUBGROUPS
    ev.SUBGROUPS = FAIRNESS
    try:
        return {f"xgboost:{v}": ev.subgroups(np.asarray(variants[v][0](a["X"]), dtype=float), np.asarray(a["y"], dtype=float),
                                             variants[v][1], a) for v in ("E1", "E1_platt", "E2")}
    finally:
        ev.SUBGROUPS = saved


def exact_tumbarello(a: dict) -> np.ndarray:
    f = scores.columns(a["X"])
    return (3 * (a["tum_hosp_12m"] > 0) + 3 * (a["tum_transfer"] > 0) + 2 * scores.recent_blfq(f)
            + 2 * (a["tum_charlson4"] > 0) + 2 * (a["tum_age70"] > 0) + 2 * (a["tum_catheter_30d"] > 0)).astype(int)


def table1(con, table_e1: str, table_e2: str) -> dict:
    moments = {}
    for estimand, table in (("E1", table_e1), ("E2", table_e2)):
        exprs = ", ".join(f"AVG({c}), VAR_SAMP({c})" for c in ev.TABLE1_COVARIATES)
        values = con.execute(f"SELECT COUNT(*), SUM(label), {exprs} FROM {table}").fetchone()
        moments[estimand] = {"units": values[0], "label_positive": values[1],
                             **{c: (values[2 + 2 * i], values[3 + 2 * i]) for i, c in enumerate(ev.TABLE1_COVARIATES)}}
    out = {"units": {e: cell(m["units"]) for e, m in moments.items()},
           "label_positive": {e: cell(int(m["label_positive"] or 0)) for e, m in moments.items()}, "covariates": {}}
    for c in ev.TABLE1_COVARIATES:
        (m1, v1), (m2, v2) = moments["E1"][c], moments["E2"][c]
        pooled = math.sqrt((v1 + v2) / 2) if None not in (v1, v2) else 0
        out["covariates"][c] = {"E1_mean": round(m1, 4), "E2_mean": round(m2, 4),
                                "smd_E2_minus_E1": round((m2 - m1) / pooled, 3) if pooled else None}
    return out


def derivation_arrays(path: Path, specimen: str | None = None) -> dict:
    con = duckdb.connect()
    where = "" if specimen is None else f"WHERE specimen_type = '{specimen}'"
    partitions = [r[0] for r in con.execute(f"SELECT DISTINCT partition FROM read_parquet('{path}')").fetchall()]
    if partitions != ["derivation"]:
        raise ValueError(f"{path} must hold derivation units only")
    a = con.execute(f"SELECT label, {', '.join(FEATURES_48)} FROM read_parquet('{path}') {where} ORDER BY unit_id").fetchnumpy()
    con.close()
    return {"y": np.asarray(a["label"]).astype(int), "X": np.column_stack([np.asarray(a[c], dtype=np.float32) for c in FEATURES_48])}


def render(summary: dict) -> str:
    lines = ["# 09_mimic_evaluate", "", "MIMIC-IV test partition, SAP amendment E-9. Point estimates with 95% patient-bootstrap intervals.", ""]
    g = summary["retrained"]["grid"]
    lines += ["| cell | units | AUROC | intercept | slope | sNB 0.10 |", "|---|---|---|---|---|---|"]
    for name, c in g.items():
        p, i = c["point"], c["interval"]
        lines.append(f"| {name} | {c['units']} | {ev.fmt(p.get('auroc'), i.get('auroc'))} | {ev.fmt(p.get('calibration_intercept'), i.get('calibration_intercept'))} | "
                     f"{ev.fmt(p.get('calibration_slope'), i.get('calibration_slope'))} | {ev.fmt(p.get('standardized_net_benefit_0.10'), i.get('standardized_net_benefit_0.10'))} |")
    lines += ["", "Clinicians' empiric choices (eMAR-covered orders), all specimens:", ""]
    for t, block in summary["clinicians"]["all"]["thresholds"].items():
        lines.append(f"Threshold {t}:")
        for name, m in block.items():
            if " minus " in name:
                continue
            lines.append(f"- {name}: flagged {m['flagged']['point']}, sensitivity {m['sensitivity']['point']}, sNB {m['standardized_net_benefit']['point']}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--boot", type=int, default=ev.N_BOOT)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    frozen = require_git_tag("sap-v1.0", REPO_ROOT / "SAP.md")
    run = RunManifest("09_mimic_evaluate", seed=args.seed, args=vars(args))
    run.set("analysis_plan", frozen)
    run.set("amendment", "SAP section 17, E-9 and its addendum (MIMIC-IV external validation)")
    n_boot, seed = args.boot, args.seed

    # MIMIC training runs and models.
    runs = {sp: latest("08_mimic_train", sp) for sp in (None, "urine", "blood")}
    trained = {sp: json.loads((d / "summary.json").read_text())["sites"][build.SITE] for sp, d in runs.items()}
    for sp, d in runs.items():
        run.add_input(f"08_mimic_train.{sp or 'all'}", d / "summary.json", sha256_file(d / "summary.json"))
    names = {e: build.MimicSpec(e).name for e in ("E1", "E2")}
    folder = {sp: build.WORK / "models" / (f"{build.REPEAT_RULE}_{build.UNIT}" + (f"_{sp}" if sp else "")) for sp in runs}
    kinds = ("xgboost", "logistic")
    variants = {k: ev.load_variants(folder[None], names, k, trained[None]) for k in kinds}
    urine_variants = {k: ev.load_variants(folder["urine"], names, k, trained["urine"]) for k in kinds}
    blood_variants = {k: ev.load_variants(folder["blood"], names, k, trained["blood"]) for k in kinds}
    cohort_dir = build.WORK / "cohorts"
    derivation_paths = {e: cohort_dir / names[e] / "derivation.parquet" for e in names}
    for e, p in derivation_paths.items():
        run.add_input(f"mimic.{names[e]}.derivation", p, sha256_file(p))

    # Test cohorts: the single read of the sealed partition.
    con = build.connect()
    build.create_views(con)
    built = build.build_cultures(con)
    build.build_context(con)
    lists = build.code_lists()
    build.build_features(con, lists)
    build.build_extras(con, lists)
    run.add_input("mimic.split_manifest", build.SPLIT_PATH, sha256_file(build.SPLIT_PATH))
    specs = {"E1": build.MimicSpec("E1"), "E2": build.MimicSpec("E2"),
             "ceftriaxone:E1": build.MimicSpec("E1", "ceftriaxone"), "ceftriaxone:E2": build.MimicSpec("E2", "ceftriaxone")}
    cohorts, flows = {}, {}
    for key, spec in specs.items():
        table, flows[key] = build.build_cohort(con, spec, built, ("test",))
        cohorts[key] = arrays_from(con, table)
        run.count(f"test.{key}.units", len(cohorts[key]["y"]))
    full = {e: build.build_cohort(con, specs[e], built, ("derivation", "test")) for e in ("E1", "E2")}
    summary = {"partition": "test", "exploratory": True, "amendment": "E-9", "analysis_plan": frozen,
               "train_runs": {sp or "all": str(d) for sp, d in runs.items()},
               "flows": {k: [{kk: (cell(vv) if isinstance(vv, int) else vv) for kk, vv in row.items()} for row in f] for k, f in flows.items()},
               "full_flows": {e: [{kk: (cell(vv) if isinstance(vv, int) else vv) for kk, vv in row.items()} for row in full[e][1]] for e in full},
               "table1": table1(con, full["E1"][0], full["E2"][0])}
    con.close()
    primary = {k: cohorts[k] for k in ("E1", "E2")}
    print("mimic evaluate: cohorts built", flush=True)

    # Retrained grid and the ARMD exploratory analyses, at BIDMC.
    summary["retrained"] = ev.evaluate_site(primary, variants, n_boot, seed)
    print("retrained grid done", flush=True)
    summary["exploratory"] = ev.exploratory_site(primary, variants, urine_variants, derivation_paths["E2"], n_boot, seed)
    print("exploratory done", flush=True)
    thresholds = {k: ev.specimen_thresholds(folder[None], cohort_dir, names, k)[0] for k in kinds}
    summary["exploratory_2"] = ev.second_exploratory_site(cohorts, variants, blood_variants, thresholds, derivation_paths, n_boot, seed)
    print("second exploratory done", flush=True)

    # The Tumbarello score from exact MIMIC items.
    summary["tumbarello_exact"] = {"rules": {}, "matched": {}}
    for test in ("E1", "E2"):
        a = cohorts[test]
        points = exact_tumbarello(a)
        for cutoff in (scores.TUMBARELLO_CUTOFF, scores.TUMBARELLO_HIGH_CUTOFF):
            flags = points >= cutoff
            name = f"tumbarello_exact_{cutoff}->{test}"
            summary["tumbarello_exact"]["rules"][name] = ev.evaluate_rule(points, flags, a["y"], a["patient_id"], n_boot, seed, test)
            summary["tumbarello_exact"]["matched"][name] = ev.matched_comparison(flags, variants["xgboost"]["E2"][0](a["X"]), a["y"],
                                                                                 a["patient_id"], n_boot, seed, test)

    # Portability of frozen ARMD models.
    armd = {}
    for sp in (None, "urine", "blood"):
        armd_run = REPO_ROOT / ARMD_RUNS[sp]
        armd_summary = json.loads((armd_run / "summary.json").read_text())
        run.add_input(f"02_train_e1.{sp or 'all'}", armd_run / "summary.json", sha256_file(armd_run / "summary.json"))
        armd[sp] = {s: ev.load_variants(get_site(s).parquet_dir.parent / "models" / (f"chained_14d_specimen_day" + (f"_{sp}" if sp else "")),
                                        names, "xgboost", armd_summary["sites"][s]) for s in ARMD_SITES}
    deriv_all = derivation_arrays(derivation_paths["E2"])
    summary["portability"] = {"all": transport(primary["E2"], deriv_all, variants["xgboost"], armd[None], n_boot, seed, "E2")}
    for sp, th in (("urine", ev.DECISION_THRESHOLDS), ("blood", ev.BLOOD_THRESHOLDS)):
        rows = mask(primary["E2"], "specimen_type", [sp])
        summary["portability"][sp] = transport(subset(primary["E2"], rows), derivation_arrays(derivation_paths["E2"], sp),
                                               (urine_variants if sp == "urine" else blood_variants)["xgboost"],
                                               {s: {"E2": v["E2"]} for s, v in armd[sp].items()}, n_boot, seed, f"{sp}:E2", thresholds=th)
    print("portability done", flush=True)

    # Clinicians' empiric choices on eMAR-covered orders, with the models on the same units.
    e2 = primary["E2"]
    covered = np.asarray(e2["emar_covered"]).astype(int) > 0
    summary["clinicians"] = {}
    groups = {"all": (covered, variants["xgboost"]),
              "urine": (covered & mask(e2, "specimen_type", ["urine"]), urine_variants["xgboost"]),
              "blood": (covered & mask(e2, "specimen_type", ["blood"]), blood_variants["xgboost"])}
    for name, (rows, vs) in groups.items():
        a = subset(e2, rows)
        models = {"all_orders_model": vs["E2"][0](a["X"]), "positive_culture_model_recalibrated": vs["E1_platt"][0](a["X"]),
                  "positive_culture_model": vs["E1"][0](a["X"])}
        summary["clinicians"][name] = clinician_comparison(a, models, SPECIMEN_THRESHOLDS[name], n_boot, seed, f"clin:{name}")
    summary["clinicians"]["eligible_share_of_test_orders"] = round(float(covered.mean()), 4)
    print("clinicians done", flush=True)

    # Subgroups: treated at order; urine by setting; fairness.
    treated = covered & (np.asarray(e2["treated_at_order"]).astype(int) > 0)
    summary["subgroups"] = {"treated_at_order": grid_on(subset(e2, treated), variants["xgboost"], "treated", n_boot, seed)}
    urine = mask(e2, "specimen_type", ["urine"])
    for label, settings in (("urine_inpatient_or_ed", ["inpatient", "ed"]), ("urine_outpatient", ["outpatient"])):
        rows = urine & mask(e2, "setting_primary", settings)
        summary["subgroups"][label] = grid_on(subset(e2, rows), urine_variants["xgboost"], label, n_boot, seed)
    summary["subgroups"]["fairness"] = fairness(e2, variants["xgboost"])
    print("subgroups done", flush=True)

    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (run.dir / "report.md").write_text(render(summary))
    print(f"mimic evaluate: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
