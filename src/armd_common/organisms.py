"""Organism name normalization and grouping.

Rules run over the distinct raw organism names of a site, which number a few hundred at most, and
the lookup is joined back in SQL. Before matching, names are uppercased, legacy ZZ or ZZZ prefixes
are removed and lab annotations in parentheses or braces are dropped.

Grouping follows current taxonomy, so Klebsiella aerogenes (formerly Enterobacter aerogenes) is in
the klebsiella group. Generic descriptors such as "lactose fermenting gram negative rods" get a
best-guess group and ambiguous=True so QC can count them.
"""

from __future__ import annotations

import re
from dataclasses import astuple, dataclass

import duckdb

ORGANISM_GROUPS = ("ecoli", "klebsiella", "proteus", "enterobacter", "other_enterobacterales", "non_enterobacterales", "none")
ENTEROBACTERALES_GROUPS = ("ecoli", "klebsiella", "proteus", "enterobacter", "other_enterobacterales")

# Culture-level summary uses the first group present in this order.
GROUP_PRECEDENCE = ENTEROBACTERALES_GROUPS + ("non_enterobacterales",)

IDSA_ESBL_SPECIES = ("ESCHERICHIA COLI", "KLEBSIELLA PNEUMONIAE", "KLEBSIELLA OXYTOCA", "PROTEUS MIRABILIS")

OTHER_ENTEROBACTERALES_GENERA = frozenset({
    "BUTTIAUXELLA", "CEDECEA", "CITROBACTER", "CRONOBACTER", "EDWARDSIELLA", "ESCHERICHIA", "EWINGELLA",
    "HAFNIA", "KLUYVERA", "LECLERCIA", "LELLIOTTIA", "MORGANELLA", "PANTOEA", "PLESIOMONAS", "PLURALIBACTER",
    "PROVIDENCIA", "RAHNELLA", "RAOULTELLA", "SALMONELLA", "SERRATIA", "SHIGELLA", "TATUMELLA", "YERSINIA",
})

RENAMED = {
    "ENTEROBACTER AEROGENES": "KLEBSIELLA AEROGENES",
    "STAPH AUREUS": "STAPHYLOCOCCUS AUREUS",
}

_ZZ_PREFIX = re.compile(r"^Z{2,}")
_ANNOTATION = re.compile(r"\s*[({\[][^)}\]]*[)}\]]?")
_SPACES = re.compile(r"\s+")
_LIKELY_ENTERIC = re.compile(r"\bENTERIC GRAM NEGATIVE|^LACTOSE FERMENTING GRAM NEGATIVE")
_GENERIC_GNR = re.compile(r"GRAM NEGATIVE (ROD|ORGANISM|BACILL)")
_NON_FERMENTER = re.compile(r"NON[ -]?FERMENTING|NONFERM")


@dataclass(frozen=True)
class OrganismInfo:
    name_std: str | None
    group: str | None
    is_enterobacterales: bool | None
    idsa_esbl_species: bool | None
    ambiguous: bool


def normalize_name(raw: str | None) -> str | None:
    if raw is None:
        return None
    name = _ZZ_PREFIX.sub("", raw.strip().upper())
    name = _ANNOTATION.sub("", name)
    name = _SPACES.sub(" ", name).strip(" -,")
    for old, new in RENAMED.items():
        if name == old or name.startswith(old + " "):
            name = new + name[len(old):]
    return name or None


def _is_species(name: str, species: str) -> bool:
    return name == species or name.startswith(species + " ")


def classify(raw: str | None) -> OrganismInfo:
    name = normalize_name(raw)
    if name is None:
        return OrganismInfo(None, None, None, None, ambiguous=False)
    genus = name.split(" ", 1)[0]
    ambiguous = False
    if _is_species(name, "ESCHERICHIA COLI"):
        group = "ecoli"
    elif genus == "KLEBSIELLA":
        group = "klebsiella"
    elif genus == "PROTEUS":
        group = "proteus"
    elif genus == "ENTEROBACTER":
        group = "enterobacter"
    elif genus in OTHER_ENTEROBACTERALES_GENERA:
        group = "other_enterobacterales"
    elif _LIKELY_ENTERIC.search(name):
        group, ambiguous = "other_enterobacterales", True
    elif _GENERIC_GNR.search(name) and not _NON_FERMENTER.search(name):
        group, ambiguous = "non_enterobacterales", True
    else:
        group = "non_enterobacterales"
    idsa = any(_is_species(name, s) for s in IDSA_ESBL_SPECIES)
    return OrganismInfo(name, group, group in ENTEROBACTERALES_GROUPS, idsa, ambiguous)


def register_lookup(con: duckdb.DuckDBPyConnection, raw_names, table: str = "organism_lookup") -> str:
    """Create a temp table mapping each distinct raw name to its normalized name and flags."""
    rows = [(raw, *astuple(classify(raw))) for raw in sorted({n for n in raw_names if n is not None})]
    con.execute(
        f"CREATE OR REPLACE TEMP TABLE {table} (organism_raw VARCHAR, organism_std VARCHAR, organism_group VARCHAR, "
        "is_enterobacterales BOOLEAN, idsa_esbl_species BOOLEAN, organism_ambiguous BOOLEAN)"
    )
    if rows:
        con.executemany(f"INSERT INTO {table} VALUES (?, ?, ?, ?, ?, ?)", rows)
    return table
