"""Site loader interface.

A loader reads one site's staged Parquet and builds the canonical tables from schema.py as TEMP
tables named <site>_<table>. Staged tables are exposed as views named <site>_raw_<table>.

Mapping logic lives in Python (vocab.py and organisms.py). It is applied by joining small lookup
tables built from each column's distinct values, so the row-level work stays in SQL.

Subclasses provide three queries:
  culture_stage_sql   one row per (anon_id, order_proc_id_coded) with the STAGE_COLUMNS below
  isolate_ast_sql     the ISOLATE_AST columns
  enzyme_tests_sql    the ENZYME_TESTS columns, or None when the release has no enzyme tests
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Callable

import duckdb

from .. import organisms, vocab
from ..config import get_site
from ..db import connect
from ..schema import CULTURE_EPISODES, CULTURE_KEYS, ENZYME_TESTS, ISOLATE_AST, ddl, projection

STAGE_COLUMNS = (
    "anon_id", "order_proc_id_coded", "n_encounters", "order_time", "specimen_raw", "culture_positive",
    "any_unidentified", "ordering_mode_raw", "inpatient", "outpatient", "ed", "icu", "urgent_care",
    "day_surgery", "age_raw", "gender_raw", "adi_score", "adi_state_rank",
)


def register_value_lookup(con: duckdb.DuckDBPyConnection, table: str, raw_values, fn: Callable, value_type: str = "VARCHAR") -> str:
    """Create a TEMP table (raw, value) from the distinct non-NULL raw values mapped through fn."""
    rows = [(raw, fn(raw)) for raw in sorted({str(v) for v in raw_values if v is not None})]
    con.execute(f"CREATE OR REPLACE TEMP TABLE {table} (raw VARCHAR, value {value_type})")
    if rows:
        con.executemany(f"INSERT INTO {table} VALUES (?, ?)", rows)
    return table


def specimen_type(raw: str) -> str:
    return vocab.SPECIMEN.get(raw.strip().upper(), "other")


class SiteLoader(ABC):
    site: str = ""
    # Candidate names for the culture order time. Releases spell it differently.
    order_time_columns: tuple[str, ...] = ()
    order_time_resolution: str = "second"
    sex_masked: bool = True
    # Predicate on raw cohort rows that marks a positive culture. Written without a table alias.
    positive_row_sql: str = ""
    # Staged feature tables and candidate names for the column holding days between the event and
    # the culture order. None means the table has no such column.
    feature_time_columns: dict[str, tuple[str, ...] | None] = {}
    one_row_per_culture: tuple[str, ...] = ("demographics", "ward_info", "adi_scores")
    # Tables whose encounter and culture IDs do not match the cohort, so they link by patient and order date.
    date_linked_tables: tuple[str, ...] = ()

    def __init__(self, con: duckdb.DuckDBPyConnection | None = None):
        self.cfg = get_site(self.site)
        self.con = con if con is not None else connect(self.cfg.tmp_dir)
        self._registered = False

    def raw(self, table: str) -> str:
        return f"{self.site}_raw_{table}"

    def table(self, kind: str) -> str:
        return f"{self.site}_{kind}"

    def culture_id_sql(self, alias: str = "") -> str:
        p = f"{alias}." if alias else ""
        return f"('{self.site}:' || {p}anon_id || ':' || {p}order_proc_id_coded)"

    def patient_id_sql(self, alias: str = "") -> str:
        p = f"{alias}." if alias else ""
        return f"('{self.site}:' || {p}anon_id)"

    def staged_tables(self) -> list[str]:
        return [t for t in self.cfg.tables if (self.cfg.parquet_dir / t).is_dir()]

    def register_raw(self) -> None:
        if self._registered:
            return
        needed = ("cohort", *self.one_row_per_culture)
        if missing := [t for t in needed if not (self.cfg.parquet_dir / t).is_dir()]:
            raise FileNotFoundError(
                f"{self.site}: staged Parquet missing for {missing}; run python -m armd_common.convert --site {self.site}"
            )
        for t in self.staged_tables():
            self.con.execute(f"CREATE OR REPLACE VIEW {self.raw(t)} AS SELECT * FROM read_parquet('{self.cfg.parquet_glob(t)}')")
        self._registered = True

    def columns(self, table: str) -> list[str]:
        return [row[0] for row in self.con.execute(f"DESCRIBE {self.raw(table)}").fetchall()]

    def resolve_column(self, table: str, candidates) -> str | None:
        """Return the first candidate present in a staged table, matching names case-insensitively."""
        if not candidates:
            return None
        actual = {c.lower(): c for c in self.columns(table)}
        return next((actual[c.lower()] for c in candidates if c.lower() in actual), None)

    def feature_link(self, table: str) -> tuple[str, str, str]:
        """How a staged feature table (alias f) links to culture_keys (alias k).

        Returns a select list over f that exposes anon_id and link, the join condition, and the
        method name. Date linkage fans a row out to every culture the patient had that day.
        """
        if table in self.date_linked_tables:
            order_time = self.resolve_column(table, self.order_time_columns)
            if order_time is None:
                raise ValueError(f"{self.site}.{table} is date-linked but has no order time column")
            return (
                f"f.anon_id, CAST(f.{order_time} AS DATE) AS link",
                "f.anon_id = k.anon_id AND f.link = k.order_date",
                "patient_and_order_date",
            )
        return ("f.anon_id, f.order_proc_id_coded AS link", "f.anon_id = k.anon_id AND f.link = k.order_proc_id_coded", "culture_id")

    def feature_rows_sql(self, table: str) -> str:
        """Rows of a staged feature table with the canonical culture_id attached."""
        select, on, _ = self.feature_link(table)
        return (
            f"SELECT k.culture_id, f.* EXCLUDE (link) "
            f"FROM (SELECT {select}, * EXCLUDE (anon_id) FROM {self.raw(table)} f) f "
            f"JOIN {self.table('culture_keys')} k ON {on}"
        )

    def distinct(self, table: str, column: str) -> list:
        return [row[0] for row in self.con.execute(f"SELECT DISTINCT {column} FROM {self.raw(table)}").fetchall()]

    def lookup(self, name: str, table: str, column: str, fn: Callable, value_type: str = "VARCHAR") -> str:
        return register_value_lookup(self.con, f"{self.site}_lookup_{name}", self.distinct(table, column), fn, value_type)

    @abstractmethod
    def culture_stage_sql(self) -> str: ...

    @abstractmethod
    def isolate_ast_sql(self, organism_lookup: str) -> str: ...

    def enzyme_tests_sql(self, organism_lookup: str) -> str | None:
        return None

    def build(self) -> dict[str, str]:
        """Build every canonical table and return a map from kind to TEMP table name."""
        self.register_raw()
        org = organisms.register_lookup(self.con, self.distinct("cohort", "organism"), table=f"{self.site}_lookup_organism")
        built = {"isolate_ast": self._materialize("isolate_ast", ISOLATE_AST, self.isolate_ast_sql(org))}
        enzyme_sql = self.enzyme_tests_sql(org)
        if enzyme_sql is None:
            name = self.table("enzyme_tests")
            self.con.execute(f"CREATE OR REPLACE TEMP TABLE {name} ({ddl(ENZYME_TESTS)})")
            built["enzyme_tests"] = name
        else:
            built["enzyme_tests"] = self._materialize("enzyme_tests", ENZYME_TESTS, enzyme_sql)
        built["culture_episodes"] = self._materialize("culture_episodes", CULTURE_EPISODES, self._culture_episodes_sql(org))
        order_time = self.resolve_column("cohort", self.order_time_columns)
        built["culture_keys"] = self._materialize(
            "culture_keys",
            CULTURE_KEYS,
            f"SELECT {self.culture_id_sql()} AS culture_id, anon_id, order_proc_id_coded, CAST(MIN({order_time}) AS DATE) AS order_date "
            f"FROM {self.raw('cohort')} GROUP BY anon_id, order_proc_id_coded",
        )
        return built

    def _materialize(self, kind: str, spec: dict[str, str], query: str) -> str:
        name = self.table(kind)
        self.con.execute(f"CREATE OR REPLACE TEMP TABLE {name} AS SELECT {projection(spec, 'q')} FROM ({query}) q")
        return name

    def _culture_episodes_sql(self, org: str) -> str:
        specimen = self.lookup("specimen", "cohort", "culture_description", specimen_type)
        age = self.lookup("age", "demographics", "age", vocab.age_band)
        if self.sex_masked:
            sex_join, masked_sql = "", "s.gender_raw"
            mapping = vocab.masked_sex_mapping(self.site)
            if mapping:
                whens = " ".join(f"WHEN '{code}' THEN '{value}'" for code, value in mapping.items())
                codes = ", ".join(f"'{code}'" for code in mapping)
                sex_sql = f"CASE s.gender_raw {whens} ELSE 'unknown' END"
                inferred_sql = f"COALESCE(s.gender_raw IN ({codes}), FALSE)"
            else:
                sex_sql, inferred_sql = "NULL", "FALSE"
        else:
            sex = self.lookup("sex", "demographics", "gender", vocab.sex)
            sex_join, sex_sql, masked_sql = f"LEFT JOIN {sex} sx ON sx.raw = s.gender_raw", "COALESCE(sx.value, 'unknown')", "NULL"
            inferred_sql = "FALSE"
        precedence = " ".join(f"WHEN list_contains(o.groups, '{g}') THEN '{g}'" for g in organisms.GROUP_PRECEDENCE)
        return f"""
            WITH s AS ({self.culture_stage_sql()}),
            o AS (
                SELECT c.anon_id, c.order_proc_id_coded,
                       TRUE AS has_identified,
                       COUNT(DISTINCT l.organism_std) AS n_organisms,
                       BOOL_OR(l.is_enterobacterales) AS has_enterobacterales,
                       BOOL_OR(l.idsa_esbl_species) AS has_idsa_esbl_species,
                       BOOL_OR(l.organism_ambiguous) AS organism_ambiguous,
                       LIST(DISTINCT l.organism_group) AS groups
                FROM {self.raw('cohort')} c
                JOIN {org} l ON l.organism_raw = c.organism
                WHERE {self.positive_row_sql}
                GROUP BY c.anon_id, c.order_proc_id_coded
            )
            SELECT '{self.site}' AS site,
                   {self.patient_id_sql('s')} AS patient_id,
                   {self.culture_id_sql('s')} AS culture_id,
                   s.n_encounters,
                   s.order_time,
                   '{self.order_time_resolution}' AS order_time_resolution,
                   COALESCE(sp.value, 'other') AS specimen_type,
                   s.specimen_raw,
                   s.culture_positive,
                   CASE WHEN NOT s.culture_positive THEN 'none'
                        WHEN o.groups IS NULL THEN NULL
                        ELSE CASE {precedence} END END AS organism_group,
                   COALESCE(o.n_organisms, 0) AS n_organisms,
                   COALESCE(o.has_enterobacterales, FALSE) AS has_enterobacterales,
                   COALESCE(o.has_idsa_esbl_species, FALSE) AS has_idsa_esbl_species,
                   s.culture_positive AND (o.has_identified IS NULL OR s.any_unidentified) AS organism_unidentified,
                   COALESCE(o.organism_ambiguous, FALSE) AS organism_ambiguous,
                   CASE WHEN s.inpatient THEN 'inpatient'
                        WHEN s.ed THEN 'ed'
                        WHEN s.urgent_care THEN 'urgent_care'
                        WHEN s.outpatient OR s.day_surgery THEN 'outpatient'
                        ELSE 'unknown' END AS setting_primary,
                   s.inpatient AS setting_raw_inpatient,
                   s.outpatient AS setting_raw_outpatient,
                   s.ed AS setting_raw_ed,
                   s.icu AS setting_raw_icu,
                   s.urgent_care AS setting_raw_urgent_care,
                   s.day_surgery AS setting_raw_day_surgery,
                   s.ordering_mode_raw,
                   ag.value AS age_band,
                   {sex_sql} AS sex,
                   {masked_sql} AS sex_code_masked,
                   {inferred_sql} AS sex_inferred,
                   s.adi_score,
                   s.adi_state_rank,
                   s.adi_score IS NULL AS adi_missing
            FROM s
            LEFT JOIN o USING (anon_id, order_proc_id_coded)
            LEFT JOIN {specimen} sp ON sp.raw = s.specimen_raw
            LEFT JOIN {age} ag ON ag.raw = s.age_raw
            {sex_join}
        """
