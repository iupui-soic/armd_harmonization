"""ESBL label constructors for the estimand ladder.

Each constructor sorts every culture in `cultures` into one reason, then maps reasons to a label
and an evaluable flag. It returns LabelResult(frame, provenance). The frame has one row per input
culture with columns culture_id, label (TINYINT), evaluable (BOOLEAN) and reason. label is NULL
exactly when evaluable is false, and evaluable is the mask that encodes layer T.

Reasons, checked in this order:
  culture_negative       the culture was reported negative
  resistant              a target isolate has a qualifying result that meets the resistance rule
  not_resistant          a target isolate has a qualifying result and none meets the rule
  not_tested             a target isolate exists but no qualifying agent has an S, I or R result
  organism_unidentified  positive, no target isolate, and at least one isolate lacks an organism name
  organism_absent        positive, and no isolate belongs to the target organisms

Results from every target isolate in a culture are pooled, so one resistant isolate or one
resistant duplicate result is enough. NT and other never count as tested.

Expected input columns (see schema.py):
  cultures  culture_id, culture_positive, organism_unidentified
  ast       culture_id, is_enterobacterales, idsa_esbl_species, antibiotic_code, and a SIR column
  enzymes   culture_id, is_enterobacterales, enzyme_class, enzyme_result
"""

from __future__ import annotations

import re
import uuid
from dataclasses import dataclass

import duckdb

from .vocab import SENTINEL_AGENTS

REASONS = ("culture_negative", "resistant", "not_resistant", "not_tested", "organism_unidentified", "organism_absent")

# Status quo (E1): only cultures that grew a target organism and had a qualifying result.
E1_RULE = {"resistant": 1, "not_resistant": 0}

# Decision time (E2): negative cultures and cultures without a target organism contribute Y = 0.
E2_RULE = {"resistant": 1, "not_resistant": 0, "culture_negative": 0, "organism_absent": 0}

_IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class LabelInputError(ValueError):
    pass


@dataclass(frozen=True)
class LabelResult:
    frame: duckdb.DuckDBPyRelation
    provenance: dict

    @property
    def label(self) -> duckdb.DuckDBPyRelation:
        return self.frame.select("culture_id, label")

    @property
    def evaluable_mask(self) -> duckdb.DuckDBPyRelation:
        return self.frame.select("culture_id, evaluable")


def _quoted(values) -> str:
    return ", ".join("'" + str(v).replace("'", "''") + "'" for v in values)


def _source(con: duckdb.DuckDBPyConnection, obj, prefix: str) -> str:
    if isinstance(obj, str):
        if not _IDENTIFIER.match(obj):
            raise LabelInputError(f"not a table or view name: {obj!r}")
        return obj
    name = f"_{prefix}_{uuid.uuid4().hex[:12]}"
    obj.create_view(name)
    return name


def _check_cultures(con: duckdb.DuckDBPyConnection, cultures: str) -> None:
    n, distinct, null_positive = con.execute(
        f"SELECT COUNT(*), COUNT(DISTINCT culture_id), COUNT(*) - COUNT(culture_positive) FROM {cultures}"
    ).fetchone()
    if n != distinct:
        raise LabelInputError("cultures must have one row per culture_id")
    if null_positive:
        raise LabelInputError("culture_positive must not be NULL")


def _build(con, cultures, results, *, constructor, estimand, target, tested_sql, hit_sql, rule, details) -> LabelResult:
    cultures = _source(con, cultures, "cultures")
    results = _source(con, results, "results")
    _check_cultures(con, cultures)
    table = f"_labels_{uuid.uuid4().hex[:12]}"
    label_cases = " ".join(f"WHEN '{reason}' THEN {value}" for reason, value in rule.items())
    con.execute(f"""
        CREATE TEMP TABLE {table} AS
        WITH per_culture AS (
            SELECT culture_id,
                   COALESCE(BOOL_OR({target}), FALSE) AS has_target,
                   COALESCE(BOOL_OR({target} AND ({tested_sql})), FALSE) AS tested,
                   COALESCE(BOOL_OR({target} AND ({hit_sql})), FALSE) AS hit
            FROM {results}
            GROUP BY culture_id
        ),
        reasons AS (
            SELECT c.culture_id,
                   CASE
                       WHEN NOT c.culture_positive THEN 'culture_negative'
                       WHEN p.hit THEN 'resistant'
                       WHEN p.tested THEN 'not_resistant'
                       WHEN p.has_target THEN 'not_tested'
                       WHEN COALESCE(c.organism_unidentified, FALSE) THEN 'organism_unidentified'
                       ELSE 'organism_absent'
                   END AS reason
            FROM {cultures} c
            LEFT JOIN per_culture p USING (culture_id)
        )
        SELECT culture_id,
               CAST(CASE reason {label_cases} END AS TINYINT) AS label,
               list_contains([{_quoted(rule)}], reason) AS evaluable,
               reason
        FROM reasons
    """)
    frame = con.table(table)
    counts = dict(con.execute(f"SELECT reason, COUNT(*) FROM {table} GROUP BY 1").fetchall())
    provenance = {
        "constructor": constructor,
        "estimand": estimand,
        "target": target,
        "rule": rule,
        **details,
        "reason_counts": {r: counts.get(r, 0) for r in REASONS},
    }
    return LabelResult(frame, provenance)


def _sir_rule(agents, nonsusceptible: bool, sir_column: str) -> tuple[str, str, dict]:
    if not _IDENTIFIER.match(sir_column):
        raise LabelInputError(f"not a column name: {sir_column!r}")
    agents = tuple(agents)
    hit_values = ("I", "R") if nonsusceptible else ("R",)
    tested = f"antibiotic_code IN ({_quoted(agents)}) AND {sir_column} IN ('S', 'I', 'R')"
    hit = f"antibiotic_code IN ({_quoted(agents)}) AND {sir_column} IN ({_quoted(hit_values)})"
    details = {"agents": list(agents), "hit_values": list(hit_values), "sir_column": sir_column}
    return tested, hit, details


def esbl_phenotype_composite(con, cultures, ast, agents=SENTINEL_AGENTS, nonsusceptible: bool = False, sir_column: str = "sir") -> LabelResult:
    """E1 label: resistance to at least one sentinel agent in an Enterobacterales isolate.

    Evaluable only for positive cultures with an Enterobacterales isolate and at least one sentinel
    agent reported as S, I or R. With nonsusceptible=True an I result also counts as a hit.
    """
    tested, hit, details = _sir_rule(agents, nonsusceptible, sir_column)
    return _build(con, cultures, ast, constructor="esbl_phenotype_composite", estimand="E1",
                  target="is_enterobacterales", tested_sql=tested, hit_sql=hit, rule=E1_RULE, details=details)


def esbl_phenotype_idsa(con, cultures, ast, nonsusceptible: bool = True, sir_column: str = "sir") -> LabelResult:
    """E1 label under the IDSA definition: ceftriaxone non-susceptible in E. coli, K. pneumoniae,
    K. oxytoca or P. mirabilis. Other Enterobacterales count as organism_absent."""
    tested, hit, details = _sir_rule(("CRO",), nonsusceptible, sir_column)
    return _build(con, cultures, ast, constructor="esbl_phenotype_idsa", estimand="E1",
                  target="idsa_esbl_species", tested_sql=tested, hit_sql=hit, rule=E1_RULE, details=details)


def esbl_enzyme_class(con, cultures, enzymes) -> LabelResult:
    """E1 label from a direct ESBL test on an Enterobacterales isolate. Only ARMD-MGB reports these.

    A culture is resistant if any ESBL test is positive and not_resistant if every ESBL test is
    negative. Phenotypic S, I and R results are ignored by design.
    """
    tested = "enzyme_class = 'ESBL' AND enzyme_result IN ('positive', 'negative')"
    hit = "enzyme_class = 'ESBL' AND enzyme_result = 'positive'"
    return _build(con, cultures, enzymes, constructor="esbl_enzyme_class", estimand="E1",
                  target="is_enterobacterales", tested_sql=tested, hit_sql=hit, rule=E1_RULE,
                  details={"enzyme_class": "ESBL"})


def esbl_decision_time(con, cultures, ast, agents=SENTINEL_AGENTS, nonsusceptible: bool = False, sir_column: str = "sir") -> LabelResult:
    """E2 label: the composite rule over all ordered cultures.

    Negative cultures and positive cultures without an Enterobacterales isolate are Y = 0, because
    there is no resistant organism to cover. Cultures with a target isolate but no sentinel result,
    and cultures whose isolates are unidentified, stay non-evaluable.
    """
    tested, hit, details = _sir_rule(agents, nonsusceptible, sir_column)
    return _build(con, cultures, ast, constructor="esbl_decision_time", estimand="E2",
                  target="is_enterobacterales", tested_sql=tested, hit_sql=hit, rule=E2_RULE, details=details)
