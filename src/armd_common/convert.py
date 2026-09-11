"""Stage a site's raw CSV release as Parquet, once.

Staging keeps raw column names and sniffed types. It maps the release's null tokens to NULL and
forces the shared ID columns to VARCHAR, so joins never depend on what the sniffer guessed. If a
column fails to parse, only that column is retried as VARCHAR. Every CSV is checked against
data_sources.yaml first. Harmonization happens later in the site loaders.

Usage: python -m armd_common.convert --site stanford [--tables cohort,demographics] [--force]
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sys
import time
from pathlib import Path

import duckdb

from .config import DEFAULT_SEED, ID_COLUMNS, NULL_TOKENS, SITES, SiteConfig, get_site
from .db import connect, sql_list
from .provenance import RunManifest, sha256_file

SAMPLE_SIZE = 500_000
ROW_GROUP_SIZE = 500_000
FILE_SIZE = "512MB"
MAX_COLUMN_RETRIES = 8

# DuckDB CSV errors quote the offending raw line. Only a line number and column name may leave
# this module, so raw values never reach logs.
_COLUMN_RE = re.compile(r'converting column "([^"]+)"')
_LINE_RE = re.compile(r"CSV Error on Line: (\d+)")


class ConversionError(RuntimeError):
    pass


def _header(path: Path) -> list[str]:
    with path.open(newline="") as fh:
        return next(csv.reader(fh))


def _safe_message(exc: duckdb.Error) -> str:
    text = str(exc)
    parts = [type(exc).__name__]
    if line := _LINE_RE.search(text):
        parts.append(f"line {line.group(1)}")
    if column := _COLUMN_RE.search(text):
        parts.append(f"column {column.group(1)}")
    return ", ".join(parts)


def _read_csv_sql(path: Path, varchar_columns: set[str]) -> str:
    options = [
        "header = true",
        f"nullstr = {sql_list(NULL_TOKENS)}",
        f"sample_size = {SAMPLE_SIZE}",
    ]
    if varchar_columns:
        types = ", ".join(f"'{c}': 'VARCHAR'" for c in sorted(varchar_columns))
        options.append(f"types = {{{types}}}")
    return f"read_csv('{path}', {', '.join(options)})"


def verify_source(cfg: SiteConfig, table: str) -> str:
    src = cfg.csv_path(table)
    expected = cfg.expected_file(src.name)
    if expected is None:
        raise ConversionError(f"{src.name} is not listed for {cfg.site} in data_sources.yaml")
    if expected["bytes"] != src.stat().st_size:
        raise ConversionError(f"{src.name} size differs from data_sources.yaml, so the release may have changed")
    digest = sha256_file(src, cache_path=cfg.parquet_dir / "_sha256_cache.json")
    if expected["sha256"] != digest:
        raise ConversionError(f"{src.name} digest differs from data_sources.yaml, so the release may have changed")
    return digest


def convert_table(con: duckdb.DuckDBPyConnection, cfg: SiteConfig, table: str, force: bool = False) -> dict:
    src = cfg.csv_path(table)
    out = cfg.parquet_dir / table
    if out.exists() and not force:
        return {"table": table, "status": "exists"}
    tmp = out.with_name(f"{table}.tmp")
    header = _header(src)
    varchar = {c for c in ID_COLUMNS if c in header}
    retried: list[str] = []
    started = time.time()
    for attempt in range(MAX_COLUMN_RETRIES + 1):
        shutil.rmtree(tmp, ignore_errors=True)
        try:
            con.execute(
                f"COPY (SELECT * FROM {_read_csv_sql(src, varchar)}) TO '{tmp}' "
                f"(FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE {ROW_GROUP_SIZE}, FILE_SIZE_BYTES '{FILE_SIZE}')"
            )
            break
        except duckdb.Error as exc:
            shutil.rmtree(tmp, ignore_errors=True)
            column = _COLUMN_RE.search(str(exc))
            if column is None or column.group(1) in varchar or attempt == MAX_COLUMN_RETRIES:
                raise ConversionError(f"{cfg.site}.{table}: {_safe_message(exc)}") from None
            varchar.add(column.group(1))
            retried.append(column.group(1))
    shutil.rmtree(out, ignore_errors=True)
    tmp.rename(out)
    glob = cfg.parquet_glob(table)
    rows = con.execute(f"SELECT COUNT(*) FROM read_parquet('{glob}')").fetchone()[0]
    schema = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{glob}')").fetchall()
    return {
        "table": table,
        "status": "converted",
        "csv": src.name,
        "rows": rows,
        "columns": [[name, typ] for name, typ, *_ in schema],
        "retried_as_varchar": retried,
        "seconds": round(time.time() - started, 1),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage a site's raw CSV release as Parquet.")
    parser.add_argument("--site", required=True, choices=sorted(SITES))
    parser.add_argument("--tables", help="comma-separated staging tables (default: all)")
    parser.add_argument("--force", action="store_true", help="rebuild tables that already exist")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="unused here, kept for a uniform CLI")
    args = parser.parse_args(argv)

    cfg = get_site(args.site)
    tables = args.tables.split(",") if args.tables else list(cfg.tables)
    if unknown := sorted(set(tables) - set(cfg.tables)):
        parser.error(f"unknown tables for {cfg.site}: {unknown}")
    if missing := [t for t in tables if not cfg.csv_path(t).exists()]:
        parser.error(f"raw files missing under {cfg.raw_dir}: {missing}")
    tables.sort(key=lambda t: cfg.csv_path(t).stat().st_size)

    cfg.parquet_dir.mkdir(parents=True, exist_ok=True)
    record_path = cfg.parquet_dir / "_conversion.json"
    record = json.loads(record_path.read_text()) if record_path.exists() else {"site": cfg.site, "tables": {}}
    record["source"] = cfg.source
    manifest = RunManifest(f"convert_{cfg.site}", seed=args.seed, args=vars(args))
    con = connect(cfg.tmp_dir)

    for table in tables:
        digest = verify_source(cfg, table)
        manifest.add_input(table, cfg.csv_path(table), digest)
        result = convert_table(con, cfg, table, force=args.force)
        if result["status"] == "converted":
            result["sha256"] = digest
            record["tables"][table] = result
            record_path.write_text(json.dumps(record, indent=2))
            manifest.count(f"{table}.rows", result["rows"])
        print(
            f"{cfg.site}.{table}: {result['status']} rows={result.get('rows')} "
            f"varchar_retries={result.get('retried_as_varchar')} seconds={result.get('seconds')}",
            flush=True,
        )

    manifest.set("parquet_dir", str(cfg.parquet_dir))
    print("manifest:", manifest.write())
    return 0


if __name__ == "__main__":
    sys.exit(main())
