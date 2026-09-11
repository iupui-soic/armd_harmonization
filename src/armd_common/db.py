"""DuckDB connections tuned for a shared server.

Spill files go to the site folder under the data root so nothing large lands on the system disk.
"""

from __future__ import annotations

import os
from pathlib import Path

import duckdb


def connect(tmp_dir: Path | None = None, threads: int | None = None, memory_limit: str | None = None) -> duckdb.DuckDBPyConnection:
    threads = threads or int(os.environ.get("ARMD_DUCKDB_THREADS", "8"))
    memory_limit = memory_limit or os.environ.get("ARMD_DUCKDB_MEMORY", "48GB")
    con = duckdb.connect()
    con.execute(f"SET threads = {int(threads)}")
    con.execute(f"SET memory_limit = '{memory_limit}'")
    con.execute("SET preserve_insertion_order = false")
    if tmp_dir is not None:
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        con.execute(f"SET temp_directory = '{tmp_dir}'")
    return con


def sql_list(values) -> str:
    """Render a Python sequence of strings as a DuckDB list literal."""
    return "[" + ", ".join("'" + str(v).replace("'", "''") + "'" for v in values) + "]"
