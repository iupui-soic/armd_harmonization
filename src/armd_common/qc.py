"""QC and profiling.

Every count that leaves the package for a report, log or shared file passes through cell(), so
counts from 1 to MIN_CELL - 1 never appear. Zero is not a small cell and is shown.

python -m armd_common.qc --site utsw builds a site's canonical tables and writes
results/qc_<site>/<timestamp>/ with manifest.json, summary.json and report.md. The report covers
staged row counts, schema problems, missingness, category frequencies and vocabulary gaps.
Vocabulary gaps list raw values without counts, so rare codes can be mapped without exposing
small cells.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import vocab
from .config import DEFAULT_SEED, MIN_CELL, SITES

SMALL = f"<{MIN_CELL}"

CATEGORICAL = {
    "culture_episodes": (
        "specimen_type", "culture_positive", "organism_group", "has_enterobacterales", "organism_unidentified",
        "setting_primary", "age_band", "sex", "sex_code_masked", "sex_inferred", "order_time_resolution",
    ),
    "isolate_ast": ("organism_group", "antibiotic_code", "sir", "sir_clsi2022", "preliminary", "ast_method"),
    "enzyme_tests": ("enzyme_class", "enzyme_result", "organism_group"),
}


def cell(n):
    if n is None:
        return None
    n = int(n)
    return SMALL if 0 < n < MIN_CELL else n


def rate(numerator, denominator, digits: int = 4):
    """A proportion, or a marker when the numerator or denominator is a small cell."""
    if not denominator:
        return None
    if 0 < numerator < MIN_CELL or 0 < denominator < MIN_CELL:
        return SMALL
    return round(numerator / denominator, digits)


def value_counts(con, table: str, column: str) -> list[list]:
    rows = con.execute(f"SELECT {column}, COUNT(*) FROM {table} GROUP BY 1 ORDER BY 2 DESC, 1 NULLS LAST").fetchall()
    return [[value, cell(n)] for value, n in rows]


def null_shares(con, table: str) -> dict:
    columns = [row[0] for row in con.execute(f"DESCRIBE {table}").fetchall()]
    total = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    if not total:
        return {c: None for c in columns}
    nulls = con.execute("SELECT " + ", ".join(f"COUNT(*) - COUNT({c})" for c in columns) + f" FROM {table}").fetchone()
    return {c: rate(n, total) for c, n in zip(columns, nulls)}


def _table_exists(con, name: str) -> bool:
    return bool(con.execute("SELECT COUNT(*) FROM duckdb_tables() WHERE table_name = ?", [name]).fetchone()[0])


def vocab_gaps(con, site: str, built: dict) -> dict:
    """Raw values that the vocabularies could not place. Values only, never counts."""

    def values(sql: str) -> list:
        return [row[0] for row in con.execute(sql).fetchall()]

    def unmapped(lookup: str) -> list | None:
        name = f"{site}_lookup_{lookup}"
        return values(f"SELECT raw FROM {name} WHERE value IS NULL ORDER BY 1") if _table_exists(con, name) else None

    ast, episodes = built["isolate_ast"], built["culture_episodes"]
    organism_lookup = f"{site}_lookup_organism"
    antibiotic_gaps = unmapped("antibiotic")
    by_design = None
    if antibiotic_gaps is not None:
        by_design = [n for n in antibiotic_gaps if vocab.normalize_antibiotic(n) in vocab.UNMAPPED_BY_DESIGN]
        antibiotic_gaps = [n for n in antibiotic_gaps if n not in by_design]
    return {
        "antibiotic_names_without_code": antibiotic_gaps,
        "antibiotic_names_unmapped_by_design": by_design,
        "antibiotic_codes_outside_vocab": [
            c for c in values(f"SELECT DISTINCT antibiotic_code FROM {ast} WHERE antibiotic_code IS NOT NULL ORDER BY 1") if c not in vocab.CODES
        ],
        "susceptibility_values_mapped_to_other": values(
            f"SELECT sir_raw FROM {ast} WHERE sir = 'other' UNION SELECT sir_clsi2022_raw FROM {ast} WHERE sir_clsi2022 = 'other' ORDER BY 1"
        ),
        "specimen_values_mapped_to_other": values(f"SELECT DISTINCT specimen_raw FROM {episodes} WHERE specimen_type = 'other' ORDER BY 1"),
        "age_values_unmapped": unmapped("age"),
        "sex_values_unmapped": unmapped("sex"),
        "ambiguous_organism_names": values(f"SELECT organism_raw FROM {organism_lookup} WHERE organism_ambiguous ORDER BY 1")
        if _table_exists(con, organism_lookup)
        else None,
    }


def sentinel_reporting(con, built: dict, sir_column: str = "sir") -> dict:
    """How often each sentinel agent is reported for cultures that grew Enterobacterales (layer T preview)."""
    episodes, ast = built["culture_episodes"], built["isolate_ast"]
    agents = ", ".join(f"'{a}'" for a in vocab.SENTINEL_AGENTS)
    rows = con.execute(f"""
        WITH ent AS (SELECT culture_id FROM {episodes} WHERE has_enterobacterales),
        agents AS (SELECT unnest([{agents}]) AS code),
        res AS (
            SELECT culture_id, antibiotic_code,
                   -- A NULL result makes the IN test NULL, so default to FALSE or NULL-only agents vanish.
                   COALESCE(BOOL_OR({sir_column} IN ('S', 'I', 'R')), FALSE) AS interpretable
            FROM {ast}
            WHERE is_enterobacterales AND antibiotic_code IN ({agents})
            GROUP BY ALL
        )
        SELECT a.code, COUNT(*),
               COUNT(*) FILTER (WHERE r.interpretable),
               COUNT(*) FILTER (WHERE r.culture_id IS NOT NULL AND NOT r.interpretable),
               COUNT(*) FILTER (WHERE r.culture_id IS NULL)
        FROM ent CROSS JOIN agents a
        LEFT JOIN res r ON r.culture_id = ent.culture_id AND r.antibiotic_code = a.code
        GROUP BY a.code ORDER BY a.code
    """).fetchall()
    any_row = con.execute(f"""
        SELECT COUNT(*), COUNT(*) FILTER (WHERE EXISTS (
            SELECT 1 FROM {ast} x WHERE x.culture_id = e.culture_id AND x.is_enterobacterales
              AND x.antibiotic_code IN ({agents}) AND x.{sir_column} IN ('S', 'I', 'R')))
        FROM {episodes} e WHERE e.has_enterobacterales
    """).fetchone()
    out = {
        code: {"reported_sir": rate(sir, n), "listed_without_sir": rate(listed, n), "not_listed": rate(absent, n)}
        for code, n, sir, listed, absent in rows
    }
    out["any_sentinel_reported_sir"] = rate(any_row[1], any_row[0])
    out["enterobacterales_cultures"] = cell(any_row[0])
    return out


def render(summary: dict) -> str:
    lines = [f"# QC report: {summary['site']}", "", "## Staged tables", "", "| table | rows | columns retried as VARCHAR |", "|---|---|---|"]
    for table, t in summary["staged_tables"].items():
        lines.append(f"| {table} | {t['rows']} | {', '.join(t['retried_as_varchar']) or 'none'} |")
    for kind, entry in summary["canonical"].items():
        lines += ["", f"## {kind}", "", f"Rows: {entry['rows']}. Schema problems: {entry['schema_problems'] or 'none'}.", ""]
        missing = {c: s for c, s in entry["null_share"].items() if s not in (0, 0.0, None)}
        lines.append(f"Columns with missing values: {missing or 'none'}")
        for column, counts in entry.get("frequencies", {}).items():
            shown = ", ".join(f"{v}: {n}" for v, n in counts[:25])
            more = f" (+{len(counts) - 25} more)" if len(counts) > 25 else ""
            lines += ["", f"- **{column}**: {shown}{more}"]
    lines += ["", "## Sentinel agent reporting among Enterobacterales cultures", ""]
    for interpretation, table in summary["sentinel_reporting"].items():
        lines.append(f"- {interpretation}: {table}")
    lines += ["", "## Vocabulary gaps", ""]
    lines += [f"- {k}: {v if v is not None else 'not applicable'}" for k, v in summary["vocab_gaps"].items()]
    return "\n".join(lines) + "\n"


def site_report(site: str, seed: int = DEFAULT_SEED) -> Path:
    from .config import get_site
    from .db import connect
    from .provenance import RunManifest
    from .schema import validate
    from .sites import get_loader

    cfg = get_site(site)
    loader = get_loader(site, connect(cfg.tmp_dir))
    built = loader.build()
    con = loader.con
    manifest = RunManifest(f"qc_{site}", seed=seed, args={"site": site})
    conversion = json.loads((cfg.parquet_dir / "_conversion.json").read_text())
    summary = {"site": site, "source": cfg.source, "staged_tables": {}, "canonical": {}}
    for table, record in sorted(conversion["tables"].items()):
        manifest.add_input(f"{site}.{table}", cfg.csv_path(table), record["sha256"])
        summary["staged_tables"][table] = {"rows": cell(record["rows"]), "retried_as_varchar": record["retried_as_varchar"]}
    for kind, table in built.items():
        n = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        manifest.count(f"{kind}.rows", n)
        entry = {"rows": cell(n), "schema_problems": validate(con, table, kind), "null_share": null_shares(con, table)}
        if kind in CATEGORICAL:
            entry["frequencies"] = {c: value_counts(con, table, c) for c in CATEGORICAL[kind]}
        summary["canonical"][kind] = entry
    summary["sentinel_reporting"] = {"sir": sentinel_reporting(con, built)}
    if site == "mgb":
        summary["sentinel_reporting"]["sir_clsi2022"] = sentinel_reporting(con, built, "sir_clsi2022")
    summary["vocab_gaps"] = vocab_gaps(con, site, built)
    manifest.write()
    (manifest.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (manifest.dir / "report.md").write_text(render(summary))
    con.close()
    return manifest.dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Write a QC report for one site's canonical tables.")
    parser.add_argument("--site", required=True, choices=SITES)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args(argv)
    print(f"qc: wrote {site_report(args.site, args.seed)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
