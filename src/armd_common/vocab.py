"""Controlled vocabularies shared by every site loader.

Covers specimen type, age band, sex, susceptibility categories and antibiotic codes. Canonical
antibiotic codes are the American Society for Microbiology abbreviations that ARMD-MGB ships in
AST_code, including its suffixes for combinations such as CAZ_AVI. Sites that only ship names are
mapped through normalized names and a short alias list.
"""

from __future__ import annotations

import re

SPECIMEN_TYPES = ("urine", "blood", "respiratory", "other")
SPECIMEN = {"URINE": "urine", "BLOOD": "blood", "RESPIRATORY": "respiratory", "RESPIRATORY_TRACT": "respiratory"}

AGE_BANDS = ("18-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85-89", "90+")

SEX_VALUES = ("f", "m", "unknown")
# MGB discloses sex. "X" covers 11 cultures and has no slot in the canonical enum.
SEX = {"female": "f", "male": "m", "unknown": "unknown", "x": "unknown"}

SETTINGS = ("inpatient", "ed", "urgent_care", "outpatient", "unknown")

SIR_VALUES = ("S", "I", "R", "NT", "other")
SIR = {
    "susceptible": "S",
    "intermediate": "I",
    "resistant": "R",
    # Susceptible dose-dependent sits in the intermediate band for reporting purposes.
    "susceptible dose dependent": "I",
    "susceptible dose-dependent": "I",
    # Only drugs without an R breakpoint report non-susceptible, so it is treated as R.
    "non-susceptible": "R",
    "not tested": "NT",
    "inconclusive": "other",
    "synergism": "other",
    "positive": "other",
    "negative": "other",
}

SENTINEL_AGENTS = ("CRO", "CAZ", "FEP", "TZP")

CODES = {
    "AMB": "amphotericin_b",
    "AMC": "amoxicillin_clavulanate",
    "AMK": "amikacin",
    "AMP": "ampicillin",
    "AMX": "amoxicillin",
    "ATM": "aztreonam",
    "AZM": "azithromycin",
    "CAZ": "ceftazidime",
    "CAZ_AVI": "ceftazidime_avibactam",
    "CEC": "cefaclor",
    "CEF": "cephalothin",
    "CFM": "cefixime",
    "CFZ": "cefazolin",
    "CHL": "chloramphenicol",
    "CIP": "ciprofloxacin",
    "CLI": "clindamycin",
    "CLR": "clarithromycin",
    "CPD": "cefpodoxime",
    "CPT": "ceftaroline",
    "CRO": "ceftriaxone",
    "CST": "colistin",
    "CTT": "cefotetan",
    "CTX": "cefotaxime",
    "CXM": "cefuroxime",
    "DAP": "daptomycin",
    "DOR": "doripenem",
    "DOX": "doxycycline",
    "ERV": "eravacycline",
    "ERY": "erythromycin",
    "ETP": "ertapenem",
    "FDC": "cefiderocol",
    "FEP": "cefepime",
    "FOF": "fosfomycin",
    "FOX": "cefoxitin",
    "GEN": "gentamicin",
    "IPM": "imipenem",
    "IPM_REL": "imipenem_relebactam",
    "KAN": "kanamycin",
    "LVX": "levofloxacin",
    "LZD": "linezolid",
    "MEM": "meropenem",
    "MEM_VAB": "meropenem_vaborbactam",
    "MIN": "minocycline",
    "MTZ": "metronidazole",
    "MXF": "moxifloxacin",
    "NIT": "nitrofurantoin",
    "OXA": "oxacillin",
    "PEN": "penicillin",
    "PIP": "piperacillin",
    "PMB": "polymyxin_b",
    "QUD": "quinupristin_dalfopristin",
    "RIF": "rifampin",
    "STR": "streptomycin",
    "SXM": "sulfamethoxazole",
    "SXT": "trimethoprim_sulfamethoxazole",
    "TET": "tetracycline",
    "TGC": "tigecycline",
    "TIC": "ticarcillin",
    "TIM": "ticarcillin_clavulanate",
    "TLV": "telavancin",
    "TOB": "tobramycin",
    "TOL_TAZ": "ceftolozane_tazobactam",
    "TZP": "piperacillin_tazobactam",
    "VAN": "vancomycin",
}

# Release-specific spellings, keyed by normalized name.
ALIASES = {
    "amoxicillin_clavulanic_acid": "AMC",
    "cephalexin_cephalothin": "CEF",
    # Stanford ships "Imipenem/Ebactam". Read as imipenem/relebactam, the only imipenem combination on US panels.
    "imipenem_ebactam": "IPM_REL",
    # Stanford ships the misspelling "Ticacarcillin/Clavulanic Acid".
    "ticacarcillin_clavulanic_acid": "TIM",
    "ticarcillin_clavulanic_acid": "TIM",
    "rifampicin": "RIF",
}

_BY_NAME = {name: code for code, name in CODES.items()}
_NON_ALNUM = re.compile(r"[^a-z0-9]+")
_AGE_RANGE = re.compile(r"^\s*(\d+)\s*-\s*(\d+)")
_AGE_TOP = re.compile(r"(above|over|>=?)\s*(\d+)|(\d+)\s*(\+|and above|or older)")


def normalize_antibiotic(name: str) -> str:
    return _NON_ALNUM.sub("_", name.strip().lower()).strip("_")


def antibiotic_code(name: str | None) -> str | None:
    if name is None:
        return None
    key = normalize_antibiotic(name)
    return ALIASES.get(key) or _BY_NAME.get(key)


def sir(raw: str | None) -> str | None:
    """Harmonize a susceptibility string. None stays None and unknown strings become 'other'."""
    if raw is None:
        return None
    return SIR.get(raw.strip().lower(), "other")


def age_band(raw: str | None) -> str | None:
    if raw is None:
        return None
    if match := _AGE_RANGE.match(raw):
        band = f"{int(match.group(1))}-{int(match.group(2))}"
        return band if band in AGE_BANDS else None
    if match := _AGE_TOP.search(raw.lower()):
        top = int(match.group(2) or match.group(3))
        return "90+" if top in (89, 90) else None
    return None


def sex(raw: str | None) -> str | None:
    if raw is None:
        return "unknown"
    return SEX.get(raw.strip().lower())


def sql_case(expr: str, mapping: dict[str, str], default: str = "NULL", lower: bool = True) -> str:
    """Render a mapping as a SQL CASE over expr, matching keys case-insensitively by default."""
    target = f"lower(trim({expr}))" if lower else expr
    whens = " ".join(
        "WHEN '{}' THEN '{}'".format(k.replace("'", "''"), v.replace("'", "''")) for k, v in mapping.items()
    )
    return f"CASE {target} {whens} ELSE {default} END"
