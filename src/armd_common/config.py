"""Paths, table maps and shared constants for every ARMD site.

All data lives under ARMD_DATA_ROOT (default /data0). Nothing here points inside the repo.
Release metadata (DOI, version, sizes, digests) lives in data_sources.yaml next to this file.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results"
SOURCES_PATH = Path(__file__).with_name("data_sources.yaml")
DEFAULT_SEED = 20260911

# Smallest count that may be printed, logged or written to a shared file.
MIN_CELL = 11

# Tokens the raw releases use for missing values. Stanford writes "Null" and "null".
NULL_TOKENS = ("", "NA", "null", "Null", "NULL")

# Identifier columns shared by every table. Staging forces them to VARCHAR.
ID_COLUMNS = ("anon_id", "pat_enc_csn_id_coded", "order_proc_id_coded")


def load_dotenv(path: Path = REPO_ROOT / ".env") -> None:
    """Read KEY=VALUE lines from .env into os.environ. Values already set win."""
    if not path.exists():
        return
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        os.environ.setdefault(key.strip(), value.strip().strip("'\""))


def data_root() -> Path:
    load_dotenv()
    return Path(os.environ.get("ARMD_DATA_ROOT", "/data0"))


@lru_cache(maxsize=1)
def load_sources() -> dict:
    return yaml.safe_load(SOURCES_PATH.read_text())


@dataclass(frozen=True)
class SiteConfig:
    site: str
    raw_subdir: str
    parquet_subdir: str
    tables: dict[str, str]

    @property
    def raw_dir(self) -> Path:
        return data_root() / self.raw_subdir

    @property
    def parquet_dir(self) -> Path:
        return data_root() / self.parquet_subdir

    @property
    def tmp_dir(self) -> Path:
        return self.parquet_dir.parent / ".duckdb_tmp"

    @property
    def source(self) -> dict:
        return {k: v for k, v in load_sources()[self.site].items() if k != "files"}

    def expected_file(self, filename: str) -> dict | None:
        return load_sources()[self.site]["files"].get(filename)

    def csv_path(self, table: str) -> Path:
        return self.raw_dir / self.tables[table]

    def parquet_glob(self, table: str) -> str:
        return str(self.parquet_dir / table / "*.parquet")


_DRYAD_STANDARD_TABLES = {
    "cohort": "microbiology_cultures_cohort.csv",
    "adi_scores": "microbiology_cultures_adi_scores.csv",
    "antibiotic_class_exposure": "microbiology_cultures_antibiotic_class_exposure.csv",
    "antibiotic_subtype_exposure": "microbiology_cultures_antibiotic_subtype_exposure.csv",
    "comorbidity": "microbiology_cultures_comorbidity.csv",
    "demographics": "microbiology_cultures_demographics.csv",
    "labs": "microbiology_cultures_labs.csv",
    "microbial_resistance": "microbiology_cultures_microbial_resistance.csv",
    "nursing_home_visits": "microbiology_cultures_nursing_home_visits.csv",
    "prior_med": "microbiology_cultures_prior_med.csv",
    "ward_info": "microbiology_cultures_ward_info.csv",
}

SITES: dict[str, SiteConfig] = {
    "stanford": SiteConfig(
        site="stanford",
        raw_subdir="armd-stanford",
        parquet_subdir="armd-stanford/parquet",
        tables={
            **_DRYAD_STANDARD_TABLES,
            "prior_infecting_organism": "microbiology_culture_prior_infecting_organism.csv",
            "prior_procedures": "microbiology_cultures_priorprocedures.csv",
            "vitals": "microbiology_cultures_vitals.csv",
            "implied_susceptibility": "microbiology_cultures_implied_susceptibility.csv",
            "implied_susceptibility_rules": "implied_susceptibility_rules.csv",
        },
    ),
    "mgb": SiteConfig(
        site="mgb",
        raw_subdir="armd-mgb/physionet.org/files/armd-mgb/1.0.0",
        parquet_subdir="armd-mgb/parquet",
        tables={
            "cohort": "microbiology_cohort_deid_tj_updated.csv",
            "adi_scores": "ADI_deid_tj.csv",
            "comorbidity": "comorbidity_deid_tj.csv",
            "demographics": "demographics_deid_tj.csv",
            "nursing_home_visits": "nursing_home_visits_deid_tj.csv",
            "prior_abx": "prior_abx_deid_tj.csv",
            "prior_micro": "prior_micro_deid_tj.csv",
            "prior_org": "prior_org_deid_tj.csv",
            "prior_procedures": "prior_procedures_deid_tj.csv",
            "ward_info": "ward_type_deid_tj.csv",
        },
    ),
    "utsw": SiteConfig(
        site="utsw",
        raw_subdir="armd-utsw",
        parquet_subdir="armd-utsw/parquet",
        tables={
            **_DRYAD_STANDARD_TABLES,
            "prior_infecting_organism": "microbiology_cultures_prior_infecting_organism.csv",
            "prior_procedures": "microbiology_cultures_prior_procedures.csv",
        },
    ),
    "ecuh": SiteConfig(
        site="ecuh",
        raw_subdir="armd-ecuh",
        parquet_subdir="armd-ecuh/parquet",
        tables={
            "cohort": "00_microbiology_cultures_cohort.csv",
            "adi_scores": "01_microbiology_cultures_adi_scores.csv",
            "antibiotic_class_exposure": "02_microbiology_cultures_antibiotic_class_exposure.csv",
            "antibiotic_subtype_exposure": "03_microbiology_cultures_antibiotic_subtype_exposure.csv",
            "comorbidity": "04_microbiology_cultures_comorbidity.csv",
            "demographics": "05_microbiology_cultures_demographics.csv",
            "labs": "06_microbiology_cultures_labs.csv",
            "microbial_resistance": "07_microbiology_cultures_microbial_resistance.csv",
            "nursing_home_visits": "08_microbiology_cultures_nursing_home_visits.csv",
            "prior_infecting_organism": "09_microbiology_cultures_prior_infecting_organism.csv",
            "prior_med": "010_microbiology_cultures_prior_med.csv",
            "prior_procedures": "011_microbiology_cultures_prior_procedures.csv",
            "ward_info": "012_microbiology_cultures_ward_info.csv",
            "vitals": "013_microbiology_cultures_vitals.csv",
        },
    ),
}


def get_site(site: str) -> SiteConfig:
    try:
        return SITES[site]
    except KeyError:
        raise ValueError(f"unknown site {site!r}; expected one of {sorted(SITES)}") from None
