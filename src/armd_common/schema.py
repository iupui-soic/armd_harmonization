"""Canonical tables that every site loader produces, with validation.

culture_episodes  one row per culture, keyed by culture_id
isolate_ast       one row per distinct culture, organism, antibiotic and result
enzyme_tests      one row per distinct culture, organism and enzyme test (ARMD-MGB only)
culture_keys      culture_id with the raw anon_id, order_proc_id_coded and order date, for joining staged feature tables

IDs are strings prefixed with the site. culture_id is '<site>:<anon_id>:<order_proc_id_coded>'
because ARMD-MGB reuses order_proc_id_coded across patients.
"""

from __future__ import annotations

from .organisms import ORGANISM_GROUPS
from .vocab import AGE_BANDS, SETTINGS, SEX_VALUES, SIR_VALUES, SPECIMEN_TYPES

SITES = ("stanford", "mgb", "utsw", "ecuh")

CULTURE_EPISODES = {
    "site": "VARCHAR",
    "patient_id": "VARCHAR",
    "culture_id": "VARCHAR",
    "n_encounters": "BIGINT",
    "order_time": "TIMESTAMP",
    "order_time_resolution": "VARCHAR",
    "specimen_type": "VARCHAR",
    "specimen_raw": "VARCHAR",
    "culture_positive": "BOOLEAN",
    "organism_group": "VARCHAR",
    "n_organisms": "BIGINT",
    "has_enterobacterales": "BOOLEAN",
    "has_idsa_esbl_species": "BOOLEAN",
    "organism_unidentified": "BOOLEAN",
    "organism_ambiguous": "BOOLEAN",
    "setting_primary": "VARCHAR",
    "setting_raw_inpatient": "BOOLEAN",
    "setting_raw_outpatient": "BOOLEAN",
    "setting_raw_ed": "BOOLEAN",
    "setting_raw_icu": "BOOLEAN",
    "setting_raw_urgent_care": "BOOLEAN",
    "setting_raw_day_surgery": "BOOLEAN",
    "ordering_mode_raw": "VARCHAR",
    "age_band": "VARCHAR",
    "sex": "VARCHAR",
    "sex_code_masked": "VARCHAR",
    "adi_score": "DOUBLE",
    "adi_state_rank": "DOUBLE",
    "adi_missing": "BOOLEAN",
}

_ORGANISM_COLUMNS = {
    "organism_raw": "VARCHAR",
    "organism_std": "VARCHAR",
    "organism_group": "VARCHAR",
    "is_enterobacterales": "BOOLEAN",
    "idsa_esbl_species": "BOOLEAN",
    "organism_ambiguous": "BOOLEAN",
}

ISOLATE_AST = {
    "site": "VARCHAR",
    "patient_id": "VARCHAR",
    "culture_id": "VARCHAR",
    **_ORGANISM_COLUMNS,
    "antibiotic_raw": "VARCHAR",
    "antibiotic_code": "VARCHAR",
    "sir_raw": "VARCHAR",
    "sir": "VARCHAR",
    "sir_clsi2022_raw": "VARCHAR",
    "sir_clsi2022": "VARCHAR",
    "mic_inequality": "VARCHAR",
    "mic_value": "DOUBLE",
    "mic_value2": "DOUBLE",
    "ast_method": "VARCHAR",
    "preliminary": "BOOLEAN",
    "n_source_rows": "BIGINT",
}

ENZYME_TESTS = {
    "site": "VARCHAR",
    "patient_id": "VARCHAR",
    "culture_id": "VARCHAR",
    **_ORGANISM_COLUMNS,
    "enzyme_class": "VARCHAR",
    "enzyme": "VARCHAR",
    "enzyme_result_raw": "VARCHAR",
    "enzyme_result": "VARCHAR",
    "ast_method": "VARCHAR",
    "preliminary": "BOOLEAN",
    "n_source_rows": "BIGINT",
}

CULTURE_KEYS = {"culture_id": "VARCHAR", "anon_id": "VARCHAR", "order_proc_id_coded": "VARCHAR", "order_date": "DATE"}

SPECS = {
    "culture_episodes": CULTURE_EPISODES,
    "isolate_ast": ISOLATE_AST,
    "enzyme_tests": ENZYME_TESTS,
    "culture_keys": CULTURE_KEYS,
}

ENUMS = {
    "culture_episodes": {
        "site": SITES,
        "order_time_resolution": ("second", "day"),
        "specimen_type": SPECIMEN_TYPES,
        "organism_group": ORGANISM_GROUPS,
        "setting_primary": SETTINGS,
        "age_band": AGE_BANDS,
        "sex": SEX_VALUES,
    },
    "isolate_ast": {"site": SITES, "organism_group": ORGANISM_GROUPS, "sir": SIR_VALUES, "sir_clsi2022": SIR_VALUES},
    "enzyme_tests": {"site": SITES, "organism_group": ORGANISM_GROUPS, "enzyme_result": ("positive", "negative", "other")},
    "culture_keys": {},
}

NOT_NULL = {
    "culture_episodes": ("site", "patient_id", "culture_id", "culture_positive", "specimen_type", "setting_primary", "order_time_resolution"),
    "isolate_ast": ("site", "patient_id", "culture_id", "preliminary"),
    "enzyme_tests": ("site", "patient_id", "culture_id", "enzyme_class", "enzyme_result"),
    "culture_keys": ("culture_id", "anon_id", "order_proc_id_coded", "order_date"),
}

UNIQUE = {"culture_episodes": "culture_id", "culture_keys": "culture_id"}


def ddl(spec: dict[str, str]) -> str:
    return ", ".join(f"{col} {typ}" for col, typ in spec.items())


def projection(spec: dict[str, str], alias: str) -> str:
    """SELECT list casting every canonical column to its declared type, in declared order."""
    return ", ".join(f"CAST({alias}.{col} AS {typ}) AS {col}" for col, typ in spec.items())


def validate(con, table: str, kind: str) -> list[str]:
    """Return schema problems in a canonical table. Messages never carry counts."""
    spec = SPECS[kind]
    problems: list[str] = []
    actual = {name: typ for name, typ, *_ in con.execute(f"DESCRIBE {table}").fetchall()}
    missing = [c for c in spec if c not in actual]
    extra = [c for c in actual if c not in spec]
    if missing or extra:
        problems.append(f"columns differ from spec: missing {missing}, extra {extra}")
    elif list(actual) != list(spec):
        problems.append("column order differs from spec")
    problems += [f"{c} is {actual[c]}, expected {t}" for c, t in spec.items() if c in actual and actual[c] != t]
    for col, allowed in ENUMS[kind].items():
        if col not in actual:
            continue
        allowed_sql = ", ".join(f"'{v}'" for v in allowed)
        if con.execute(f"SELECT COUNT(*) FROM {table} WHERE {col} IS NOT NULL AND {col} NOT IN ({allowed_sql})").fetchone()[0]:
            problems.append(f"{col} has values outside {list(allowed)}")
    for col in NOT_NULL[kind]:
        if col in actual and con.execute(f"SELECT COUNT(*) FROM {table} WHERE {col} IS NULL").fetchone()[0]:
            problems.append(f"{col} has NULLs")
    if key := UNIQUE.get(kind):
        n, distinct = con.execute(f"SELECT COUNT(*), COUNT(DISTINCT {key}) FROM {table}").fetchone()
        if n != distinct:
            problems.append(f"{key} is not unique")
    return problems
