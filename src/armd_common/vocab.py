"""Controlled vocabularies shared by every site loader.

Covers specimen type, age band, sex, susceptibility categories and antibiotic codes. Canonical
antibiotic codes are the American Society for Microbiology abbreviations that ARMD-MGB ships in
AST_code, including its suffixes for combinations such as CAZ_AVI. Sites that only ship names are
mapped through normalized names and a short alias list.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import yaml

MAPPINGS_DIR = Path(__file__).with_name("mappings")

SPECIMEN_TYPES = ("urine", "blood", "respiratory", "other")
SPECIMEN = {"URINE": "urine", "BLOOD": "blood", "RESPIRATORY": "respiratory", "RESPIRATORY_TRACT": "respiratory"}

AGE_BANDS = ("18-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85-89", "90+")

SEX_VALUES = ("f", "m", "unknown")
# MGB discloses sex. "X" covers 11 cultures in 1.0.0 and has no slot in the canonical enum. 2.0.0
# collapses everything but Female and Male into "Others".
SEX = {"female": "f", "male": "m", "unknown": "unknown", "x": "unknown", "others": "unknown"}

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
    # ARMD-MGB 2.0.0 writes "Unknown" for a lab phenotype its build did not recognize.
    "unknown": "other",
    "inconclusive": "other",
    "synergism": "other",
    "positive": "other",
    "negative": "other",
}

SENTINEL_AGENTS = ("CRO", "CAZ", "FEP", "TZP")

CODES = {
    "5FC": "flucytosine",
    "CDR": "cefdinir",
    "CFR": "cefadroxil",
    "DAL": "dalbavancin",
    "DCX": "dicloxacillin",
    "DLX": "delafloxacin",
    "DPS": "dapsone",
    "FDX": "fidaxomicin",
    "LEX": "cephalexin",
    "MTH": "methenamine",
    "NAF": "nafcillin",
    "RFX": "rifaximin",
    "TZD": "tedizolid",
    "AFG": "anidulafungin",
    "ATM_AVI": "aztreonam_avibactam",
    "BDQ": "bedaquiline",
    "CAP": "capreomycin",
    "CAS": "caspofungin",
    # Clavulanate combinations are the CLSI ESBL confirmatory tests.
    "CAZ_CLV": "ceftazidime_clavulanate",
    "CLO": "clofazimine",
    "CTX_CLV": "cefotaxime_clavulanate",
    "DCS": "cycloserine",
    # EDTA-modified and modified carbapenem inactivation methods.
    "ECIM": "ecim",
    "EMB": "ethambutol",
    "ETO": "ethionamide",
    "FLC": "fluconazole",
    "GAT": "gatifloxacin",
    "ICR": "inducible_clindamycin_resistance",
    "INH": "isoniazid",
    "ISA": "isavuconazole",
    "ITC": "itraconazole",
    "MCIM": "mcim",
    "MCZ": "miconazole",
    "MFG": "micafungin",
    "NAL": "nalidixic_acid",
    "NOR": "norfloxacin",
    "OFX": "ofloxacin",
    "OMC": "omadacycline",
    "PLZ": "plazomicin",
    "POS": "posaconazole",
    "PZA": "pyrazinamide",
    "RFB": "rifabutin",
    "SAM": "ampicillin_sulbactam",
    "SPT": "spectinomycin",
    "TMP": "trimethoprim",
    "TRB": "terbinafine",
    "VRC": "voriconazole",
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
    # ARMD-MGB 2.0.0 reports one combined carbapenem result under this name in its culture and prior
    # resistance tables. It is read as meropenem: every use of carbapenem results in armd-common is
    # by class (prior carbapenem resistance, prior carbapenem exposure), so the agent does not matter.
    "meropenem_imipenem": "MEM",
    # UTSW spellings.
    "pip_tazobactam": "TZP",
    "amox_clavulanate": "AMC",
    "amoxicillin_k_clavulanate_aug": "AMC",
    "trimeth_sulf": "SXT",
    "ticarcil_clavulanate": "TIM",
    "amikacin_iv": "AMK",
    "polymixin_e": "CST",
    "5_fluorocytosine": "5FC",
    "caspofungin_acetate": "CAS",
    "fluconizole": "FLC",
    "bedaquilline": "BDQ",
    # ECUH spellings.
    "amp_sulbactam": "SAM",
    "augmentin": "AMC",
    "trimeth_sulfamethoxazole": "SXT",
    "ticarillin_clavulanic_acid": "TIM",
    "cefotaxime_clavulanic_acid": "CTX_CLV",
    "ceftazidime_clavulanic_acid": "CAZ_CLV",
    "inducible_clindamycin": "ICR",
    # Both.
    "synercid": "QUD",
    "sulfamethoxazole_trimethoprim": "SXT",
}

# ARMD-MGB prior_abx spells some combinations differently from its AST_code column.
PRIOR_ABX_CODE_ALIASES = {"CZA": "CAZ_AVI", "C/T": "TOL_TAZ", "I-R": "IPM_REL", "MVB": "MEM_VAB"}

# Brand names seen in prior medication tables.
MEDICATION_BRANDS = {
    "augmentin": "AMC",
    "avelox": "MXF",
    "bactrim": "SXT",
    "cipro": "CIP",
    "flagyl": "MTZ",
    "hiprex": "MTH",
    "keflex": "LEX",
    "levaquin": "LVX",
    "macrobid": "NIT",
    "macrodantin": "NIT",
    "monurol": "FOF",
    "rocephin": "CRO",
    "septra": "SXT",
    "vancocin": "VAN",
    "xifaxan": "RFX",
    "zithromax": "AZM",
    "zosyn": "TZP",
    "zyvox": "LZD",
}

# Salt, strength and formulation words dropped from medication names before lookup.
_FORMULATION_TOKENS = frozenset({
    "acetate", "aug", "axetil", "calcium", "chloride", "citrate", "combo", "dextrose", "dextrs", "diluent", "ds",
    "er", "ethylsuccinate", "hcl", "hippurate", "hyclate", "in", "iso", "iv", "k", "macrocrystal", "macrocrystals",
    "mandelate", "monohydrate", "nacl", "phosphate", "po", "pot", "potassium", "sod", "sodium", "sulfate",
    "tromethamine", "xr",
})

# Topical, otic and ophthalmic products do not count as systemic antibiotic exposure.
_TOPICAL_MARKERS = frozenset({"dexamethasone", "hydrocortisone", "ophthalmic", "otic", "silver", "topical"})

# Names left without a code on purpose, keyed by normalized name. They are screens or method
# variants whose result is not the drug's own susceptibility.
UNMAPPED_BY_DESIGN = {
    "gent_screen": "high-level gentamicin synergy screen for enterococci",
    "gentamicin_screen": "high-level gentamicin synergy screen for enterococci",
    "cefoxitin_screen": "staphylococcal mecA surrogate screen",
    "isoniazid_high": "isoniazid at the high critical concentration",
    "isoniazid_low": "isoniazid at the low critical concentration",
    "tobramycin_rm": "UTSW method variant with an undocumented rm suffix",
    "vancomycin_rm": "UTSW method variant with an undocumented rm suffix",
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


def medication_code(name: str | None) -> str | None:
    """Map a prior medication name (brand, salt or formulation included) to a canonical code.

    Topical, otic and ophthalmic products return None.
    """
    if name is None:
        return None
    tokens = _NON_ALNUM.sub(" ", name.lower()).split()
    if not tokens or _TOPICAL_MARKERS.intersection(tokens):
        return None
    for token in tokens:
        if token in MEDICATION_BRANDS:
            return MEDICATION_BRANDS[token]
    core = [t for t in tokens if t not in _FORMULATION_TOKENS]
    return antibiotic_code("_".join(core)) if core else None


def prior_abx_code(code: str | None) -> str | None:
    """Canonical code for an ARMD-MGB prior_abx drug_code."""
    if code is None:
        return None
    return PRIOR_ABX_CODE_ALIASES.get(code, code)


@lru_cache(maxsize=1)
def _sex_code_file() -> dict:
    return yaml.safe_load((MAPPINGS_DIR / "sex_codes.yaml").read_text())


def masked_sex_mapping(site: str) -> dict[str, str] | None:
    """Inferred code-to-sex mapping for a site that masks sex, or None when there is none."""
    entry = (_sex_code_file().get("sites") or {}).get(site)
    return {str(code): sex for code, sex in entry["mapping"].items()} if entry else None


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
