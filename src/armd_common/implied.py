"""Implied susceptibility (plan section 4.2, decision 6). Off unless a caller applies it.

Rules come from Stanford's implied_susceptibility_rules.csv, the only release that ships them, and
are applied to any site's isolate_ast. A rule fills a missing result for one organism and agent
from reported results on other agents for the same culture and organism.

  * Reported S, I or R results are never overwritten.
  * Rules read reported results only, so imputations never chain.
  * When an isolate has conflicting reports for a source agent, the worst result (R, then I, then S)
    is used, so conflicting evidence never produces an imputed S.
  * When several rules target the same organism and agent, the first rule in file order that yields
    a value wins.

Rule kinds:
  intrinsic                    "Resistant"
  same_as                      "Same as X if reported"
  susceptible_if               "Susceptible if susceptible to X, Y, or Z"
  susceptible_if_else_r        "Susceptible if susceptible to X, otherwise Resistant"
  same_as_else_susceptible_if  "Same as X if reported; otherwise, Susceptible if susceptible to Y"

strict=True applies intrinsic resistance and same_as equivalences only, which never carry
information from one drug to a different drug's result. For the combined kind it applies only the
same_as clause. strict=False applies every rule.
"""

from __future__ import annotations

import re
import uuid
from dataclasses import dataclass

import duckdb

from .schema import ISOLATE_AST
from .vocab import antibiotic_code

KINDS = ("intrinsic", "same_as", "susceptible_if", "susceptible_if_else_r", "same_as_else_susceptible_if")
STRICT_KINDS = ("intrinsic", "same_as", "same_as_else_susceptible_if")

_SAME_ELSE = re.compile(r"^same as (?P<same>.+?) if reported;\s*otherwise,?\s*susceptible if susceptible to (?P<srcs>.+)$", re.I)
_SAME = re.compile(r"^same as (?P<same>.+?) if reported$", re.I)
_SUSC_ELSE_R = re.compile(r"^susceptible if susceptible to (?P<srcs>.+?),\s*otherwise resistant$", re.I)
_SUSC = re.compile(r"^susceptible if susceptible to (?P<srcs>.+)$", re.I)
_LIST_SPLIT = re.compile(r",\s*(?:or\s+)?|\s+or\s+", re.I)


@dataclass(frozen=True)
class Rule:
    order: int
    organism_pattern: str
    target: str
    kind: str
    same_as: str | None
    sources: tuple[str, ...]
    text: str


def _code(name: str) -> str:
    code = antibiotic_code(name.strip())
    if code is None:
        raise ValueError(f"no antibiotic code for {name.strip()!r}")
    return code


def _codes(text: str) -> tuple[str, ...]:
    return tuple(_code(part) for part in _LIST_SPLIT.split(text.strip()) if part.strip())


def organism_pattern(raw: str) -> str:
    name = raw.strip().upper()
    return re.sub(r"\s+SPECIES$", "", name)


def parse_rule(order: int, organism: str, antibiotic: str, text: str) -> Rule:
    text = text.strip()
    base = dict(order=order, organism_pattern=organism_pattern(organism), target=_code(antibiotic), text=text)
    if text.lower() == "resistant":
        return Rule(kind="intrinsic", same_as=None, sources=(), **base)
    if m := _SAME_ELSE.match(text):
        return Rule(kind="same_as_else_susceptible_if", same_as=_code(m["same"]), sources=_codes(m["srcs"]), **base)
    if m := _SAME.match(text):
        return Rule(kind="same_as", same_as=_code(m["same"]), sources=(), **base)
    if m := _SUSC_ELSE_R.match(text):
        return Rule(kind="susceptible_if_else_r", same_as=None, sources=_codes(m["srcs"]), **base)
    if m := _SUSC.match(text):
        return Rule(kind="susceptible_if", same_as=None, sources=_codes(m["srcs"]), **base)
    raise ValueError(f"unrecognized implied susceptibility rule: {text!r}")


def parse_rules(rows) -> list[Rule]:
    """Parse (organism, antibiotic, rule text) rows in file order, dropping exact duplicates."""
    seen, rules = set(), []
    for organism, antibiotic, text in rows:
        key = (organism_pattern(organism), antibiotic.strip().lower(), text.strip().lower())
        if key in seen:
            continue
        seen.add(key)
        rules.append(parse_rule(len(rules), organism, antibiotic, text))
    return rules


def stanford_rules(con: duckdb.DuckDBPyConnection) -> list[Rule]:
    from .config import get_site

    glob = get_site("stanford").parquet_glob("implied_susceptibility_rules")
    return parse_rules(con.execute(f'SELECT "Organism", "Antibiotic", "Rule" FROM read_parquet(\'{glob}\')').fetchall())


def apply_implied_susceptibility(con: duckdb.DuckDBPyConnection, ast: str, rules: list[Rule], strict: bool = True) -> str:
    """Return a TEMP table with every reported row plus imputed rows, flagged by imputed and imputation_rule."""
    kinds = STRICT_KINDS if strict else KINDS
    active = [r for r in rules if r.kind in kinds]
    rules_table = f"_implied_rules_{uuid.uuid4().hex[:12]}"
    con.execute(
        f"CREATE TEMP TABLE {rules_table} (rule_order INTEGER, organism_pattern VARCHAR, target VARCHAR, kind VARCHAR, "
        "same_as VARCHAR, sources VARCHAR[], rule_text VARCHAR)"
    )
    if active:
        con.executemany(
            f"INSERT INTO {rules_table} VALUES (?, ?, ?, ?, ?, ?, ?)",
            [(r.order, r.organism_pattern, r.target, r.kind, r.same_as, list(r.sources), r.text) for r in active],
        )
    use_sources = "FALSE" if strict else "TRUE"
    organism_columns = "site, patient_id, culture_id, organism_raw, organism_std, organism_group, is_enterobacterales, idsa_esbl_species, organism_ambiguous"
    out = f"_implied_ast_{uuid.uuid4().hex[:12]}"
    passthrough = ", ".join(ISOLATE_AST)
    con.execute(f"""
        CREATE TEMP TABLE {out} AS
        WITH reported AS (
            SELECT culture_id, organism_raw, antibiotic_code,
                   CASE WHEN BOOL_OR(sir = 'R') THEN 'R' WHEN BOOL_OR(sir = 'I') THEN 'I' ELSE 'S' END AS sir
            FROM {ast}
            WHERE sir IN ('S', 'I', 'R') AND organism_raw IS NOT NULL
            GROUP BY ALL
        ),
        isolates AS (
            SELECT DISTINCT {organism_columns} FROM {ast} WHERE organism_std IS NOT NULL
        ),
        candidates AS (
            SELECT i.*, r.*
            FROM isolates i
            JOIN {rules_table} r ON i.organism_std = r.organism_pattern OR starts_with(i.organism_std, r.organism_pattern || ' ')
            WHERE NOT EXISTS (
                SELECT 1 FROM reported p
                WHERE p.culture_id = i.culture_id AND p.organism_raw = i.organism_raw AND p.antibiotic_code = r.target
            )
        ),
        evidence AS (
            SELECT c.culture_id, c.organism_raw, c.rule_order,
                   ANY_VALUE(p.sir) FILTER (WHERE p.antibiotic_code = c.same_as) AS same_as_sir,
                   BOOL_OR(p.sir = 'S' AND list_contains(c.sources, p.antibiotic_code)) AS any_source_s,
                   BOOL_OR(list_contains(c.sources, p.antibiotic_code)) AS any_source_reported
            FROM candidates c
            LEFT JOIN reported p ON p.culture_id = c.culture_id AND p.organism_raw = c.organism_raw
            GROUP BY c.culture_id, c.organism_raw, c.rule_order
        ),
        evaluated AS (
            SELECT c.*,
                   CASE c.kind
                       WHEN 'intrinsic' THEN 'R'
                       WHEN 'same_as' THEN e.same_as_sir
                       WHEN 'susceptible_if' THEN CASE WHEN e.any_source_s THEN 'S' END
                       WHEN 'susceptible_if_else_r' THEN CASE WHEN e.any_source_s THEN 'S' WHEN e.any_source_reported THEN 'R' END
                       WHEN 'same_as_else_susceptible_if' THEN COALESCE(e.same_as_sir, CASE WHEN {use_sources} AND e.any_source_s THEN 'S' END)
                   END AS imputed_sir
            FROM candidates c
            JOIN evidence e USING (culture_id, organism_raw, rule_order)
        ),
        chosen AS (
            SELECT * FROM evaluated
            WHERE imputed_sir IS NOT NULL
            QUALIFY ROW_NUMBER() OVER (PARTITION BY culture_id, organism_raw, target ORDER BY rule_order) = 1
        )
        SELECT {passthrough}, FALSE AS imputed, CAST(NULL AS VARCHAR) AS imputation_rule FROM {ast}
        UNION ALL
        SELECT {organism_columns},
               CAST(NULL AS VARCHAR) AS antibiotic_raw, target AS antibiotic_code,
               CAST(NULL AS VARCHAR) AS sir_raw, imputed_sir AS sir,
               CAST(NULL AS VARCHAR) AS sir_clsi2022_raw, CAST(NULL AS VARCHAR) AS sir_clsi2022,
               CAST(NULL AS VARCHAR) AS mic_inequality, CAST(NULL AS DOUBLE) AS mic_value, CAST(NULL AS DOUBLE) AS mic_value2,
               CAST(NULL AS VARCHAR) AS ast_method, FALSE AS preliminary, CAST(0 AS BIGINT) AS n_source_rows,
               TRUE AS imputed, rule_text AS imputation_rule
        FROM chosen
    """)
    con.execute(f"DROP TABLE {rules_table}")
    return out
