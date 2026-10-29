"""Converts the MIMIC-IV 3.1 hosp tables OrderGap needs from gzip CSV to Parquet with DuckDB.

Every column is kept as text, so no value depends on type sniffing; later SQL casts what it uses.
Output goes to /data0/mimiciv/3.1/parquet, outside the repository. Nothing is printed but row counts.
"""
from __future__ import annotations

import sys
from pathlib import Path

import duckdb

RAW = Path("/data0/mimiciv/3.1/hosp")
OUT = Path("/data0/mimiciv/3.1/parquet")
TMP = Path("/data0/mimiciv/3.1/duckdb_tmp")
TABLES = ("patients", "admissions", "transfers", "microbiologyevents", "diagnoses_icd", "procedures_icd",
          "prescriptions", "emar")


def convert(tables=TABLES, force: bool = False) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(config={"temp_directory": str(TMP), "memory_limit": "48GB", "threads": 10})
    counts = {}
    for table in tables:
        target = OUT / f"{table}.parquet"
        if target.exists() and not force:
            counts[table] = con.execute(f"SELECT COUNT(*) FROM read_parquet('{target}')").fetchone()[0]
            continue
        tmp = target.with_suffix(".parquet.tmp")
        con.execute(f"""COPY (SELECT * FROM read_csv('{RAW / (table + '.csv.gz')}', header = true, all_varchar = true,
                        max_line_size = 10000000)) TO '{tmp}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 1000000)""")
        tmp.rename(target)
        counts[table] = con.execute(f"SELECT COUNT(*) FROM read_parquet('{target}')").fetchone()[0]
        print(f"converted {table}: {counts[table]:,} rows", flush=True)
    con.close()
    return counts


if __name__ == "__main__":
    print(convert(force="--force" in sys.argv))
