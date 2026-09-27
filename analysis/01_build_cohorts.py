"""Build the OrderGap cohorts at every site and describe them (plan section 5.2, Analysis 1).

For each staged site this builds the primary E1 and E2 cohorts and the sensitivity cohorts named in
SAP.md, on the derivation partition only because the test partition is sealed. Each cohort is written
to cohorts/<cohort name>/derivation.parquet inside the site folder in the data root.

results/01_build_cohorts/<timestamp>/ receives manifest.json, summary.json and report.md with the
TRIPOD+AI flow counts, label prevalence, the share of ordered cultures that were negative or grew no
Enterobacterales, and covariate means with standardized mean differences between E2 and E1, per
culture as released and under the proposed repeat rule and unit (SAP items R7 and R8). Counts from 1
to 10 are suppressed, and so is any flow count within 10 of another published count (qc.publish_flows).

Usage: python analysis/01_build_cohorts.py [--sites stanford,mgb] [--seed 20260911]
"""

from __future__ import annotations

import argparse
import json
import math
import shutil

from armd_common import features, splits
from armd_common.cohorts import PRIMARY_SPECS, TABLE1_COVARIATES, CohortSpec, build_cohort
from armd_common.config import DEFAULT_SEED, SITES, get_site
from armd_common.db import connect
from armd_common.provenance import RunManifest, sha256_file
from armd_common.qc import cell, publish_flows, rate
from armd_common.sites import get_loader

SENSITIVITY_SPECS = (
    CohortSpec("E1", nonsusceptible=True),
    CohortSpec("E2", nonsusceptible=True),
    CohortSpec("E1", label="idsa"),
    CohortSpec("E1", implied_susceptibility="strict"),
    CohortSpec("E2", implied_susceptibility="strict"),
    CohortSpec("E2", exclude_recent_antibiotics=True),
)
# Proposed repeat rule and unit (R7, R8), then S6 (specimen types pooled) and S7 (first index order).
UNIT_SPECS = tuple(
    CohortSpec(estimand, repeat_rule=rule, unit=unit)
    for rule, unit in (("chained_14d", "specimen_day"), ("chained_14d", "order_event"), ("first_per_patient", "specimen_day"))
    for estimand in ("E1", "E2")
)
PROPOSED = ("e1_composite_chained_14d_specimen_day", "e2_composite_chained_14d_specimen_day")
MGB_SPECS = (
    CohortSpec("E1", interpretation="clsi2022"),
    CohortSpec("E2", interpretation="clsi2022"),
    CohortSpec("E1", label="enzyme"),
)

COVARIATES = TABLE1_COVARIATES


def staged(site: str) -> bool:
    return (get_site(site).parquet_dir / "cohort").is_dir()


def moments(con, table: str) -> dict:
    exprs = ", ".join(f"AVG({c}), VAR_SAMP({c})" for c in COVARIATES)
    values = con.execute(f"SELECT {exprs} FROM {table}").fetchone()
    return {c: (values[2 * i], values[2 * i + 1]) for i, c in enumerate(COVARIATES)}


def smd(first: tuple, second: tuple) -> float | None:
    (m1, v1), (m2, v2) = first, second
    if None in (m1, m2, v1, v2):
        return None
    pooled = math.sqrt((v1 + v2) / 2)
    return round((m2 - m1) / pooled, 3) if pooled > 0 else None


def write_parquet(con, table: str, out_dir) -> None:
    tmp = out_dir.with_name(out_dir.name + ".tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    con.execute(f"COPY (SELECT * FROM {table} ORDER BY unit_id) TO '{tmp / 'derivation.parquet'}' (FORMAT parquet, COMPRESSION zstd)")
    shutil.rmtree(out_dir, ignore_errors=True)
    tmp.rename(out_dir)


def run_site(site: str, run: RunManifest) -> dict:
    cfg = get_site(site)
    conversion = json.loads((cfg.parquet_dir / "_conversion.json").read_text())
    for table, record in sorted(conversion["tables"].items()):
        run.add_input(f"{site}.{table}", cfg.csv_path(table), record["sha256"])
    run.add_input(f"{site}.split_manifest", splits.manifest_path(site), sha256_file(splits.manifest_path(site)))

    loader = get_loader(site, connect(cfg.tmp_dir))
    built = loader.build()
    con = loader.con
    feature_tables = {"lab": features.build_features(loader, built, "lab")}
    specs = PRIMARY_SPECS + SENSITIVITY_SPECS + UNIT_SPECS
    if site == "mgb":
        feature_tables["clsi2022"] = features.build_features(loader, built, "clsi2022")
        specs = specs + MGB_SPECS

    out = {"cohorts": {}}
    results = {}
    flows = {}
    for spec in specs:
        result = build_cohort(loader, built, spec, features_table=feature_tables[spec.interpretation])
        results[spec.name] = result
        write_parquet(con, result.table, cfg.parquet_dir.parent / "cohorts" / spec.name)
        for step in result.flow:
            run.count(f"{site}.{spec.name}.{step['step']}", step["cultures"])
            if step["units"] is not None and spec.unit != "culture":
                run.count(f"{site}.{spec.name}.{step['step']}.units", step["units"])
        flows[spec.name] = result.flow
        out["cohorts"][spec.name] = {"provenance": result.provenance}
        print(f"  {site}.{spec.name}: {cell(result.flow[-1]['units'])} units", flush=True)
    for name, flow in publish_flows(flows).items():
        out["cohorts"][name]["flow"] = flow

    e1, e2 = results["e1_composite"], results["e2_composite"]
    split = splits.register_manifest(con, site, table=f"{site}_split_descriptive")
    n_all, negative, no_entero = con.execute(f"""
        SELECT COUNT(*), COUNT(*) FILTER (WHERE NOT e.culture_positive),
               COUNT(*) FILTER (WHERE e.culture_positive AND NOT e.has_enterobacterales AND NOT e.organism_unidentified)
        FROM {built['culture_episodes']} e JOIN {split} s ON s.patient_id_hash = sha256(e.patient_id)
        WHERE s.partition = 'derivation'
    """).fetchone()
    out["analysis_1"] = {
        "ordered_derivation_cultures": cell(n_all),
        "negative_share": rate(negative, n_all),
        "positive_without_enterobacterales_share": rate(no_entero, n_all),
        "covariates": covariate_table(con, e1.table, e2.table),
        "covariates_proposed_unit": covariate_table(con, results[PROPOSED[0]].table, results[PROPOSED[1]].table),
    }
    con.close()
    return out


def covariate_table(con, e1_table: str, e2_table: str) -> dict:
    m1, m2 = moments(con, e1_table), moments(con, e2_table)
    return {c: {"E1_mean": round(m1[c][0], 4) if m1[c][0] is not None else None,
                "E2_mean": round(m2[c][0], 4) if m2[c][0] is not None else None,
                "smd_E2_minus_E1": smd(m1[c], m2[c])} for c in COVARIATES}


def render(summary: dict) -> str:
    lines = ["# 01_build_cohorts", "", "Derivation partition only. The test partition is sealed.", ""]
    for site, s in summary["sites"].items():
        a1 = s["analysis_1"]
        lines += [
            f"## {site}",
            "",
            f"Ordered derivation cultures: {a1['ordered_derivation_cultures']}. Negative: {a1['negative_share']}. "
            f"Positive without Enterobacterales: {a1['positive_without_enterobacterales_share']}.",
            "",
            "| cohort | step | units | cultures | patients | label = 1 | prevalence |",
            "|---|---|---|---|---|---|---|",
        ]
        for name, c in s["cohorts"].items():
            for step in c["flow"]:
                positive, prevalence = ("n/a", "n/a") if step["label_positive"] is None else (step["label_positive"], step["prevalence"])
                units = "n/a" if step["units"] is None else step["units"]
                lines.append(f"| {name} | {step['step']} | {units} | {step['cultures']} | {step['patients']} | {positive} | {prevalence} |")
        lines += ["", "Labels are counted from the partition step on, because earlier steps include the sealed test partition.",
                  "Prevalence is per unit. A withheld step is within 10 of another published count."]
        for key, title in (("covariates", "one row per culture as released"),
                           ("covariates_proposed_unit", "chained_14d repeat rule, one row per specimen-day")):
            lines += ["", f"Covariates, {title}:", "", "| covariate | E1 mean | E2 mean | SMD (E2 minus E1) |", "|---|---|---|---|"]
            for cov, v in a1[key].items():
                lines.append(f"| {cov} | {v['E1_mean']} | {v['E2_mean']} | {v['smd_E2_minus_E1']} |")
        lines.append("")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", help="comma-separated sites (default: every staged site)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    sites = args.sites.split(",") if args.sites else [s for s in SITES if staged(s)]
    run = RunManifest("01_build_cohorts", seed=args.seed, args=vars(args))
    summary = {"partition": "derivation", "sites": {}}
    for site in sites:
        print(f"cohorts: {site}", flush=True)
        summary["sites"][site] = run_site(site, run)
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (run.dir / "report.md").write_text(render(summary))
    print(f"cohorts: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
