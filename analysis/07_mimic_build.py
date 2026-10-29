"""Build the MIMIC-IV cohorts for OrderGap's external validation (SAP amendment E-9). Derivation only.

Steps: convert the hosp tables to Parquet if needed, build cultures and susceptibility results in
armd-common's canonical form, the encounter context, the 48 features and the E-9 extras, freeze the
patient split (stratified by an earlier ESBL phenotype, a feature, so no outcome is read), and write
the E1 and E2 derivation cohorts where analysis/02_train_e1.py's train_site can read them.

Row-level outputs stay under /data0/mimiciv/3.1/ordergap. results/07_mimic_build/<timestamp>/ gets
manifest.json and summary.json with aggregate counts only. No test-partition label is computed here.

Usage: python analysis/07_mimic_build.py [--seed 20260911]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from armd_common.config import DEFAULT_SEED, REPO_ROOT  # noqa: E402
from armd_common.features import FEATURES_48  # noqa: E402
from armd_common.provenance import RunManifest, require_git_tag, sha256_file  # noqa: E402
from armd_common.qc import cell  # noqa: E402
from mimic import build, convert  # noqa: E402

SPECS = (build.MimicSpec("E1"), build.MimicSpec("E2"))


def write_parquet(con, table: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    con.execute(f"COPY (SELECT * FROM {table} ORDER BY unit_id) TO '{path}' (FORMAT parquet, COMPRESSION zstd)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    frozen = require_git_tag("sap-v1.0", REPO_ROOT / "SAP.md")
    run = RunManifest("07_mimic_build", seed=args.seed, args=vars(args))
    run.set("analysis_plan", frozen)

    converted = convert.convert()
    for table, digest in build.input_digests().items():
        run.add_input(f"mimic.{table}", Path(f"/data0/mimiciv/3.1/hosp/{table}.csv.gz"), digest)
    run.add_input("mimic.code_lists", build.CODE_LISTS, sha256_file(build.CODE_LISTS))
    for table, n in converted.items():
        run.count(f"parquet.{table}", n)

    con = build.connect()
    build.create_views(con)
    built = build.build_cultures(con)
    for k, v in built["counts"].items():
        run.count(f"cultures.{k}", v)
    print("cultures built", flush=True)
    build.build_context(con)
    lists = build.code_lists()
    build.build_features(con, lists)
    print("features built", flush=True)
    build.build_extras(con, lists)
    print("extras built", flush=True)

    summary = {"partition": "derivation", "analysis_plan": frozen, "cultures": {k: cell(v) for k, v in built["counts"].items()}}
    # Label-free description of every culture: specimen, growth, setting.
    summary["by_specimen"] = {
        r[0]: {"cultures": cell(r[1]), "grew_any_organism": round(r[2], 4), "inpatient": round(r[3], 4), "ed": round(r[4], 4),
               "outpatient": round(r[5], 4), "emar_covered": round(r[6], 4)}
        for r in con.execute("""
            SELECT e.specimen_type, COUNT(*), AVG(CAST(e.culture_positive AS INTEGER)), AVG(c.hosp_ward_IP), AVG(c.hosp_ward_ER),
                   AVG(c.hosp_ward_OP), AVG(c.emar_covered)
            FROM mimic_culture_episodes e JOIN mimic_context c USING (culture_id) GROUP BY 1 ORDER BY 1""").fetchall()
    }
    split = build.write_split(con, args.seed)
    run.add_input("mimic.split_manifest", Path(split["path"]), split["sha256"])
    summary["split"] = {"sha256": split["sha256"], "patients": cell(split["patients"]),
                        "strata": {k: cell(v) for k, v in split["strata"].items()}}

    summary["cohorts"] = {}
    for spec in SPECS:
        table, flow = build.build_cohort(con, spec, built, ("derivation",))
        path = build.WORK / "cohorts" / spec.name / "derivation.parquet"
        write_parquet(con, table, path)
        run.add_input(f"mimic.{spec.name}.derivation", path, sha256_file(path))
        units, positive, patients = con.execute(f"SELECT COUNT(*), SUM(label), COUNT(DISTINCT patient_id) FROM {table}").fetchone()
        run.count(f"{spec.name}.derivation.units", units)
        means = con.execute(f"SELECT {', '.join(f'AVG({c})' for c in FEATURES_48)} FROM {table}").fetchone()
        extras = con.execute(f"SELECT {', '.join(f'AVG({c})' for c in build.EXTRA_FLAGS)} FROM {table}").fetchone()
        by_specimen = {r[0]: {"units": cell(r[1]), "label_positive": cell(r[2])} for r in con.execute(
            f"SELECT specimen_type, COUNT(*), SUM(label) FROM {table} GROUP BY 1 ORDER BY 1").fetchall()}
        summary["cohorts"][spec.estimand] = {
            "name": spec.name, "flow": [{k: (cell(v) if isinstance(v, int) else v) for k, v in row.items()} for row in flow],
            "units": cell(units), "patients": cell(patients), "label_positive": cell(int(positive or 0)),
            "prevalence": round(positive / units, 5) if units else None, "by_specimen": by_specimen,
            "feature_means": {c: round(v, 4) if v is not None else None for c, v in zip(FEATURES_48, means)},
            "extra_means": {c: round(v, 4) if v is not None else None for c, v in zip(build.EXTRA_FLAGS, extras)},
            "provenance": build.provenance(spec),
        }
        print(f"{spec.name}: {cell(units)} derivation units", flush=True)
    con.close()
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"mimic build: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
