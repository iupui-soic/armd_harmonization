"""Phase 0 feasibility gate (plan section 1).

G1  cultures and patients per site, split by culture positivity
G2  coverage of every staged feature table among positive and negative cultures
G3  whether feature-table timing is anchored to the culture order for negative cultures too
G4  order time availability for negative cultures
G5  Enterobacterales-positive fraction of all ordered cultures, plus E1 and E2 label prevalence
G6  UTSW profile, which runs as a site once the release is staged
G7  MGB access, which holds when the credentialed release is staged
G8  release provenance from data_sources.yaml and the conversion records

Writes results/00_gate/<timestamp>/manifest.json, summary.json and report.md. Only aggregate counts
leave this script, and counts from 1 to 10 are suppressed.

Usage: python analysis/00_gate.py [--sites stanford,mgb] [--seed 20260911]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from armd_common import labels
from armd_common.config import DEFAULT_SEED, SITES, get_site, load_sources
from armd_common.db import connect
from armd_common.provenance import RunManifest
from armd_common.qc import cell, rate
from armd_common.schema import validate
from armd_common.sites import get_loader

# Plan go criterion, applied to tables that should hold one row for every culture.
ONE_TO_ONE_MIN_NEGATIVE_COVERAGE = 0.90
# Event tables only hold rows when an event happened, so coverage below 100% is expected. A table
# keyed to positive cultures alone would show almost no negatives, which this threshold catches.
EVENT_MIN_RELATIVE_COVERAGE = 0.10
MAX_NULL_ORDER_TIME = 0.01
QUANTILES = (0.01, 0.25, 0.5, 0.75, 0.99)


def staged(site: str) -> bool:
    return (get_site(site).parquet_dir / "cohort").is_dir()


def g1_counts(con, ep: str) -> dict:
    total = con.execute(f"SELECT COUNT(*), COUNT(DISTINCT patient_id) FROM {ep}").fetchone()
    split = {
        positive: {"cultures": cell(n), "patients": cell(p)}
        for positive, n, p in con.execute(f"SELECT culture_positive, COUNT(*), COUNT(DISTINCT patient_id) FROM {ep} GROUP BY 1").fetchall()
    }
    return {"cultures": cell(total[0]), "patients": cell(total[1]), "positive": split.get(True), "negative": split.get(False)}


def g2_g3_tables(con, loader, ep: str, keys: str) -> dict:
    out = {}
    staged_tables = set(loader.staged_tables())
    for table, time_column in loader.feature_time_columns.items():
        if table not in staged_tables:
            out[table] = {"status": "not in release"}
            continue
        raw = loader.raw(table)
        columns = loader.columns(table)
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE _gate_keys AS
            SELECT DISTINCT k.culture_id FROM {raw} f JOIN {keys} k USING (anon_id, order_proc_id_coded)
        """)
        n_pos, n_neg, cov_pos, cov_neg = con.execute(f"""
            SELECT COUNT(*) FILTER (WHERE e.culture_positive), COUNT(*) FILTER (WHERE NOT e.culture_positive),
                   COUNT(g.culture_id) FILTER (WHERE e.culture_positive), COUNT(g.culture_id) FILTER (WHERE NOT e.culture_positive)
            FROM {ep} e LEFT JOIN _gate_keys g USING (culture_id)
        """).fetchone()
        orphans = con.execute(f"""
            SELECT COUNT(*) FROM (SELECT DISTINCT anon_id, order_proc_id_coded FROM {raw}) f
            ANTI JOIN {keys} k USING (anon_id, order_proc_id_coded)
        """).fetchone()[0]
        one_to_one = table in loader.one_row_per_culture
        pos_rate = cov_pos / n_pos if n_pos else 0.0
        neg_rate = cov_neg / n_neg if n_neg else 0.0
        relative = neg_rate / pos_rate if pos_rate else None
        if one_to_one:
            keyed = neg_rate >= ONE_TO_ONE_MIN_NEGATIVE_COVERAGE
        else:
            keyed = relative is not None and relative >= EVENT_MIN_RELATIVE_COVERAGE
        entry = {
            "status": "staged",
            "kind": "one_row_per_culture" if one_to_one else "event",
            "coverage_positive": rate(cov_pos, n_pos),
            "coverage_negative": rate(cov_neg, n_neg),
            "negative_to_positive_coverage": round(relative, 3) if relative is not None else None,
            "keys_absent_from_cohort": cell(orphans),
            "keyed_to_negatives": keyed,
        }
        if loader.order_time_column in columns:
            anchor = {}
            for positive, n, same_ts, same_day in con.execute(f"""
                WITH f AS (SELECT DISTINCT anon_id, order_proc_id_coded, CAST({loader.order_time_column} AS TIMESTAMP) AS t FROM {raw})
                SELECT e.culture_positive, COUNT(*),
                       COUNT(*) FILTER (WHERE f.t = e.order_time),
                       COUNT(*) FILTER (WHERE CAST(f.t AS DATE) = CAST(e.order_time AS DATE))
                FROM f JOIN {keys} k USING (anon_id, order_proc_id_coded) JOIN {ep} e USING (culture_id)
                GROUP BY 1
            """).fetchall():
                anchor["positive" if positive else "negative"] = {
                    "culture_time_pairs": cell(n), "same_timestamp": rate(same_ts, n), "same_day": rate(same_day, n)
                }
            entry["order_time_anchor"] = anchor
        if time_column and time_column in columns:
            timing = {}
            q_sql = "[" + ", ".join(str(q) for q in QUANTILES) + "]"
            for positive, n, quantiles, gt0, eq0, lt0, nulls in con.execute(f"""
                SELECT e.culture_positive, COUNT(*), approx_quantile(v, {q_sql}),
                       AVG(CAST(v > 0 AS DOUBLE)), AVG(CAST(v = 0 AS DOUBLE)), AVG(CAST(v < 0 AS DOUBLE)),
                       COUNT(*) - COUNT(v)
                FROM (SELECT anon_id, order_proc_id_coded, TRY_CAST({time_column} AS DOUBLE) AS v FROM {raw}) f
                JOIN {keys} k USING (anon_id, order_proc_id_coded) JOIN {ep} e USING (culture_id)
                GROUP BY 1
            """).fetchall():
                timing["positive" if positive else "negative"] = {
                    "rows": cell(n),
                    "quantiles": dict(zip((f"p{int(q * 100):02d}" for q in QUANTILES), quantiles)) if quantiles else None,
                    "share_gt0": round(gt0, 4) if gt0 is not None else None,
                    "share_eq0": round(eq0, 4) if eq0 is not None else None,
                    "share_lt0": round(lt0, 4) if lt0 is not None else None,
                    "null_rows": cell(nulls),
                }
            entry["days_to_culture"] = {"column": time_column, **timing}
        elif time_column:
            entry["days_to_culture"] = {"column": time_column, "status": "column not in release"}
        out[table] = entry
    return out


def g4_order_time(con, ep: str) -> dict:
    out = {}
    for positive, n, nulls in con.execute(f"SELECT culture_positive, COUNT(*), COUNT(*) - COUNT(order_time) FROM {ep} GROUP BY 1").fetchall():
        out["positive" if positive else "negative"] = {"cultures": cell(n), "null_order_time": cell(nulls), "null_share": rate(nulls, n)}
    return out


def g5_prevalence(con, loader, built: dict) -> dict:
    ep, ast, enz = built["culture_episodes"], built["isolate_ast"], built["enzyme_tests"]
    by_specimen = {}
    for specimen, n, pos, ent, unident in con.execute(f"""
        SELECT COALESCE(specimen_type, 'all'), COUNT(*),
               COUNT(*) FILTER (WHERE culture_positive),
               COUNT(*) FILTER (WHERE has_enterobacterales),
               COUNT(*) FILTER (WHERE organism_unidentified)
        FROM {ep} GROUP BY ROLLUP (specimen_type) ORDER BY 1
    """).fetchall():
        by_specimen[specimen] = {
            "cultures": cell(n),
            "positive_share": rate(pos, n),
            "enterobacterales_share": rate(ent, n),
            "unidentified_share": rate(unident, n),
        }
    variants = {
        "E1_composite_lab": labels.esbl_phenotype_composite(con, ep, ast),
        "E1_idsa_lab": labels.esbl_phenotype_idsa(con, ep, ast),
        "E2_decision_time_lab": labels.esbl_decision_time(con, ep, ast),
    }
    if loader.site == "mgb":
        variants["E1_composite_clsi2022"] = labels.esbl_phenotype_composite(con, ep, ast, sir_column="sir_clsi2022")
        variants["E2_decision_time_clsi2022"] = labels.esbl_decision_time(con, ep, ast, sir_column="sir_clsi2022")
        variants["E1_enzyme_class"] = labels.esbl_enzyme_class(con, ep, enz)
    estimands = {}
    for name, result in variants.items():
        n_eval, n_pos = result.frame.aggregate("COUNT(*) FILTER (WHERE evaluable), COUNT(*) FILTER (WHERE label = 1)").fetchone()
        estimands[name] = {
            "evaluable_cultures": cell(n_eval),
            "label_positive": cell(n_pos),
            "prevalence": rate(n_pos, n_eval),
            "reasons": {k: cell(v) for k, v in result.provenance["reason_counts"].items()},
        }
    return {"by_specimen": by_specimen, "estimands": estimands}


def site_verdict(site_summary: dict) -> dict:
    tables = site_summary["g2_g3"]
    staged_tables = {k: v for k, v in tables.items() if v.get("status") == "staged"}
    failing = sorted(k for k, v in staged_tables.items() if not v["keyed_to_negatives"])
    negative_order = site_summary["g4"].get("negative", {})
    null_share = negative_order.get("null_share")
    order_ok = isinstance(null_share, float) and null_share <= MAX_NULL_ORDER_TIME or null_share == 0
    return {
        "feature_tables_keyed_to_negatives": not failing,
        "tables_failing": failing,
        "order_time_available_for_negatives": bool(order_ok),
        "passes": not failing and bool(order_ok),
    }


def run_site(site: str, manifest: RunManifest) -> dict:
    cfg = get_site(site)
    loader = get_loader(site, connect(cfg.tmp_dir))
    built = loader.build()
    con = loader.con
    conversion = json.loads((cfg.parquet_dir / "_conversion.json").read_text())
    for table, record in conversion["tables"].items():
        manifest.add_input(f"{site}.{table}", cfg.csv_path(table), record["sha256"])
    for kind, table in built.items():
        manifest.count(f"{site}.{kind}.rows", con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
    summary = {
        "schema_problems": {kind: validate(con, table, kind) for kind, table in built.items()},
        "g1": g1_counts(con, built["culture_episodes"]),
        "g2_g3": g2_g3_tables(con, loader, built["culture_episodes"], built["culture_keys"]),
        "g4": g4_order_time(con, built["culture_episodes"]),
        "g5": g5_prevalence(con, loader, built),
    }
    summary["verdict"] = site_verdict(summary)
    con.close()
    return summary


def provenance_summary(sites: list[str]) -> dict:
    sources = load_sources()
    out = {}
    for site in SITES:
        src = {k: v for k, v in sources[site].items() if k != "files"}
        record_path = get_site(site).parquet_dir / "_conversion.json"
        tables = json.loads(record_path.read_text())["tables"] if record_path.exists() else {}
        out[site] = {**src, "files_listed": len(sources[site]["files"]), "tables_staged": len(tables), "gated": site in sites}
    return out


def fmt(value) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def render_report(summary: dict) -> str:
    lines = ["# Phase 0 gate report", "", f"Overall go: **{summary['overall']['go']}**. {summary['overall']['rule']}", ""]
    for site, s in summary["sites"].items():
        v = s["verdict"]
        lines += [f"## {site}", "", f"Passes: **{v['passes']}**. Tables failing G2: {v['tables_failing'] or 'none'}.", ""]
        g1 = s["g1"]
        lines += [
            f"G1: {g1['cultures']} cultures from {g1['patients']} patients. "
            f"Positive {g1['positive']['cultures']}, negative {g1['negative']['cultures']}.",
            "",
            "| table | kind | coverage pos | coverage neg | neg/pos | keyed to negatives | same order day (neg) | days-to-culture p50 pos / neg |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for table, t in s["g2_g3"].items():
            if t.get("status") != "staged":
                lines.append(f"| {table} | not in release | | | | | | |")
                continue
            anchor = t.get("order_time_anchor", {}).get("negative", {}).get("same_day")
            timing = t.get("days_to_culture", {})
            p50 = lambda key: (timing.get(key) or {}).get("quantiles", {}) or {}
            lines.append(
                f"| {table} | {t['kind']} | {fmt(t['coverage_positive'])} | {fmt(t['coverage_negative'])} | "
                f"{fmt(t['negative_to_positive_coverage'])} | {t['keyed_to_negatives']} | {fmt(anchor)} | "
                f"{fmt(p50('positive').get('p50'))} / {fmt(p50('negative').get('p50'))} |"
            )
        lines += ["", f"G4 order time: {s['g4']}", "", "G5 by specimen:", "", "| specimen | cultures | positive | Enterobacterales | unidentified |", "|---|---|---|---|---|"]
        for spec, row in s["g5"]["by_specimen"].items():
            lines.append(f"| {spec} | {row['cultures']} | {fmt(row['positive_share'])} | {fmt(row['enterobacterales_share'])} | {fmt(row['unidentified_share'])} |")
        lines += ["", "| estimand | evaluable | label = 1 | prevalence |", "|---|---|---|---|"]
        for name, row in s["g5"]["estimands"].items():
            lines.append(f"| {name} | {row['evaluable_cultures']} | {row['label_positive']} | {fmt(row['prevalence'])} |")
        problems = {k: p for k, p in s["schema_problems"].items() if p}
        lines += ["", f"Schema problems: {problems or 'none'}", ""]
    lines += ["## Provenance (G7, G8)", ""]
    for site, p in summary["provenance"].items():
        lines.append(f"- {site}: {p['repository']} {p['doi']} version {p['version']} ({p['version_date']}), downloaded {p['downloaded_on']}, {p['tables_staged']} tables staged.")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", help="comma-separated sites (default: every staged site)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    sites = args.sites.split(",") if args.sites else [s for s in SITES if staged(s)]
    manifest = RunManifest("00_gate", seed=args.seed, args=vars(args))
    summary = {
        "thresholds": {
            "one_to_one_min_negative_coverage": ONE_TO_ONE_MIN_NEGATIVE_COVERAGE,
            "event_min_relative_coverage": EVENT_MIN_RELATIVE_COVERAGE,
            "max_null_order_time": MAX_NULL_ORDER_TIME,
        },
        "sites": {},
    }
    for site in sites:
        print(f"gate: {site}", flush=True)
        summary["sites"][site] = run_site(site, manifest)
    passing = [s for s, v in summary["sites"].items() if v["verdict"]["passes"]]
    summary["overall"] = {
        "passing_sites": passing,
        "go": len(passing) >= 2 and bool({"stanford", "mgb"} & set(passing)),
        "rule": "Go needs at least two passing sites, one of them Stanford or MGB.",
    }
    summary["provenance"] = provenance_summary(sites)
    manifest.set("overall", summary["overall"])
    manifest.write()
    (manifest.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (manifest.dir / "report.md").write_text(render_report(summary))
    print(f"gate: wrote {manifest.dir}")
    print(json.dumps(summary["overall"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
