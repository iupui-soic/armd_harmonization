"""Harmonized 48-feature set for ESBL prediction at culture order time.

The first 45 features follow the AMIA 2026 ARMD-MGB model (~/ebsl-amr), rebuilt from each release so
a feature means the same thing at every site. Three specimen indicators complete the 48-feature
extended model. Differences from the AMIA notebook are deliberate:

  * Feature tables link to cultures through the loader: patient plus culture ID, or patient plus
    order date for the three MGB tables whose IDs do not match the cohort. The notebook joined on
    order_proc_id_coded alone, which attached other patients' comorbidities, procedures and nursing
    home visits.
  * Antibiotic classes are fixed sets of drug codes, taken from MGB's class membership and applied to
    every site's medication names. The Dryad class table puts all beta-lactams in one class.
  * Prior resistant organisms come from each site's prior susceptibility table with one phenotype
    definition, because the Dryad releases record prior organisms at genus level only.
  * History tables drop events on the order day (days >= 1), because UTSW and ECUH place the index
    culture's own resistance results at day 0. Nursing home visits keep day 0, which marks residence.
  * Ward flags are combined across a culture's encounters. The notebook kept the first row.

Missing values follow the notebook. Absent history is 0, a missing age band is the middle band (4),
and a missing ADI score is 0 with an indicator.

Usage: python -m armd_common.features --site stanford [--interpretation lab|clsi2022]
Writes Parquet under the site folder in the data root and results/features_<site>/<timestamp>/.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys

from . import organisms, vocab
from .config import DEFAULT_SEED, get_site
from .db import connect
from .provenance import RunManifest
from .qc import cell, rate
from .sites import get_loader
from .sites.base import register_value_lookup

FEATURES_45 = (
    "age_encoded", "gender_male", "age_65plus",
    "adi_score_clean", "adi_high", "adi_missing",
    "comorbid_heart_failure", "comorbid_liver_disease", "comorbid_lymphoma",
    "comorbid_metastatic_cancer", "comorbid_obesity", "comorbid_renal_failure", "elixhauser_count",
    "hosp_ward_IP", "hosp_ward_OP", "hosp_ward_ER", "hosp_ward_UC",
    "T_fq_90d", "T_ceph3_90d", "T_carb_90d", "T_glyco_90d", "T_sulfa_90d", "T_esp_90d", "T_amino_90d",
    "prior_ESBL", "prior_MRSA", "prior_CRE", "prior_VRE", "num_prior_orgs", "days_since_prior_org",
    "proc_cvc", "proc_mechvent", "proc_dialysis", "proc_urethral", "proc_surgical", "any_proc_30d",
    "nursing_home_90d",
    "fq_and_ceph3", "multi_prior_resistant", "any_immunocomp", "any_invasive_proc",
    "prior_esbl_and_ceph3", "hospital_multiabx", "renal_urine", "nursing_home_flag",
)
SPECIMEN_FEATURES = ("specimen_urine", "specimen_blood", "specimen_respiratory")
FEATURES_48 = FEATURES_45 + SPECIMEN_FEATURES
DOUBLE_FEATURES = frozenset({"adi_score_clean", "days_since_prior_org"})
NON_BINARY_FEATURES = frozenset({"age_encoded", "adi_score_clean", "elixhauser_count", "num_prior_orgs", "days_since_prior_org", "multi_prior_resistant"})

PRIOR_MIN_DAYS = 1
ANTIBIOTIC_WINDOW_DAYS = 90
PROCEDURE_WINDOW_DAYS = 30
NURSING_HOME_WINDOW_DAYS = 90

# MGB class membership, extended with drugs of the same class that MGB never recorded.
ANTIBIOTIC_CLASSES = {
    "T_fq_90d": ("CIP", "DLX", "GAT", "LVX", "MXF", "NOR", "OFX"),
    "T_ceph3_90d": ("CAZ", "CPT", "CRO", "CTX", "FDC", "FEP"),
    "T_carb_90d": ("DOR", "ETP", "IPM", "MEM"),
    "T_glyco_90d": ("VAN",),
    "T_sulfa_90d": ("SXT", "TMP"),
    "T_esp_90d": ("AMC", "SAM", "TZP"),
    "T_amino_90d": ("AMK", "GEN", "KAN", "PLZ", "STR", "TOB"),
}

COMORBIDITY_RULES = {
    "comorbid_heart_failure": "heart failure|congestive heart",
    "comorbid_liver_disease": "liver disease",
    "comorbid_lymphoma": "lymphoma",
    "comorbid_metastatic_cancer": "metastatic cancer",
    "comorbid_obesity": "obesity",
    "comorbid_renal_failure": "renal failure|renal sclerosis",
}
COMORBIDITY_COLUMN = {"mgb": "category"}

PROCEDURES = {
    "proc_cvc": "cvc",
    "proc_mechvent": "mechvent",
    "proc_dialysis": "dialysis",
    "proc_urethral": "urethral_catheter",
    "proc_surgical": "surgical_procedure",
}

ESBL_AGENTS = ("CRO", "CAZ", "FEP")
CARBAPENEMS = ("DOR", "ETP", "IPM", "MEM")
# Any of these beta-lactams reported resistant in S. aureus implies methicillin resistance.
MRSA_AGENTS = ("AMC", "CAZ", "CFZ", "CRO", "ETP", "FEP", "FOX", "IPM", "MEM", "OXA", "SAM", "TZP")
MRSA_NAME_PATTERN = r"(?i)\bmrsa\b|methicillin resistant staph[a-z.]* aureus"
VRE_SPECIES = ("ENTEROCOCCUS FAECALIS", "ENTEROCOCCUS FAECIUM")

# Genera recorded in the prior organism tables of all four releases.
COMMON_PRIOR_GENERA = (
    "ACINETOBACTER", "CITROBACTER", "ENTEROBACTER", "ENTEROCOCCUS", "ESCHERICHIA", "KLEBSIELLA", "MORGANELLA",
    "PROTEUS", "PROVIDENCIA", "PSEUDOMONAS", "SERRATIA", "STAPHYLOCOCCUS", "STENOTROPHOMONAS", "STREPTOCOCCUS",
)


def _quoted(values) -> str:
    return ", ".join("'" + str(v).replace("'", "''") + "'" for v in values)


def _linked(loader, table: str, columns: dict[str, str], where: str = "", distinct: bool = False) -> str:
    """Rows of a staged table (alias f) with culture_id attached through the loader's link."""
    select, on, _ = loader.feature_link(table)
    inner = ", ".join(f"{expr} AS {alias}" for alias, expr in columns.items())
    outer = ", ".join(f"f.{alias}" for alias in columns)
    return (
        f"SELECT k.culture_id, {outer} "
        f"FROM (SELECT {'DISTINCT ' if distinct else ''}{select}, {inner} FROM {loader.raw(table)} f {where}) f "
        f"JOIN {loader.table('culture_keys')} k ON {on}"
    )


def _days(loader, table: str) -> str:
    column = loader.resolve_column(table, loader.feature_time_columns.get(table))
    if column is None:
        raise ValueError(f"{loader.site}.{table} has no days-to-culture column")
    return f"TRY_CAST(f.{column} AS DOUBLE)"


def _antibiotic_part(loader) -> str | None:
    con, site = loader.con, loader.site
    if site == "mgb":
        table, name_column, days, mapper = "prior_abx", "drug_code", "TRY_CAST(f.last_dose_to_culture AS DOUBLE)", vocab.prior_abx_code
    else:
        table, name_column, mapper = "antibiotic_subtype_exposure", "medication_name", vocab.medication_code
        if table not in loader.staged_tables():
            return None
        days = _days(loader, table)
    if table not in loader.staged_tables():
        return None
    lookup = register_value_lookup(con, f"{site}_lookup_prior_antibiotic", loader.distinct(table, name_column), mapper)
    src = _linked(loader, table, {"drug": f"f.{name_column}"},
                  where=f"WHERE {days} BETWEEN {PRIOR_MIN_DAYS} AND {ANTIBIOTIC_WINDOW_DAYS}", distinct=True)
    flags = ", ".join(f"MAX(CAST(l.value IN ({_quoted(codes)}) AS INTEGER)) AS {name}" for name, codes in ANTIBIOTIC_CLASSES.items())
    return f"SELECT s.culture_id, {flags} FROM ({src}) s JOIN {lookup} l ON l.raw = s.drug GROUP BY s.culture_id"


def _resistance_part(loader, interpretation: str) -> str | None:
    con, site = loader.con, loader.site
    if site == "mgb":
        table = "prior_micro"
        if table not in loader.staged_tables():
            return None
        pheno = "AST_pheno" if interpretation == "lab" else "CLSI_2022_pheno"
        days = "TRY_CAST(f.prior_AST_time_to_culture AS DOUBLE)"
        columns = {"organism": "f.organism", "code": "f.drug_code"}
        where = f"WHERE f.{pheno} = 'Resistant' AND {days} >= {PRIOR_MIN_DAYS}"
        code_expr, abx_join = "s.code", ""
    else:
        table = "microbial_resistance"
        if table not in loader.staged_tables():
            return None
        # Every row in this table is a resistant result.
        abx = register_value_lookup(con, f"{site}_lookup_resistance_antibiotic", loader.distinct(table, "antibiotic"), vocab.antibiotic_code)
        columns = {"organism": "f.organism", "code": "f.antibiotic"}
        where = f"WHERE {_days(loader, table)} >= {PRIOR_MIN_DAYS}"
        code_expr, abx_join = "a.value", f"LEFT JOIN {abx} a ON a.raw = s.code"
    org = organisms.register_lookup(con, loader.distinct(table, "organism"), table=f"{site}_lookup_resistance_organism")
    src = _linked(loader, table, columns, where=where, distinct=True)
    vre_species = " OR ".join(f"starts_with(o.organism_std, '{s}')" for s in VRE_SPECIES)
    return f"""
        SELECT s.culture_id,
               MAX(CAST(COALESCE(o.is_enterobacterales AND {code_expr} IN ({_quoted(ESBL_AGENTS)}), FALSE) AS INTEGER)) AS prior_ESBL,
               MAX(CAST(COALESCE((starts_with(o.organism_std, 'STAPHYLOCOCCUS AUREUS') AND {code_expr} IN ({_quoted(MRSA_AGENTS)}))
                                 OR regexp_matches(s.organism, '{MRSA_NAME_PATTERN}'), FALSE) AS INTEGER)) AS prior_MRSA,
               MAX(CAST(COALESCE(o.is_enterobacterales AND {code_expr} IN ({_quoted(CARBAPENEMS)}), FALSE) AS INTEGER)) AS prior_CRE,
               MAX(CAST(COALESCE(({vre_species}) AND {code_expr} = 'VAN', FALSE) AS INTEGER)) AS prior_VRE
        FROM ({src}) s
        LEFT JOIN {org} o ON o.organism_raw = s.organism
        {abx_join}
        GROUP BY s.culture_id
    """


def _organism_part(loader) -> str | None:
    table, name_column = ("prior_org", "prior_org") if loader.site == "mgb" else ("prior_infecting_organism", "prior_organism")
    if table not in loader.staged_tables():
        return None
    days = _days(loader, table)
    src = _linked(loader, table, {"genus": f"upper(split_part(trim(f.{name_column}), ' ', 1))", "days": days},
                  where=f"WHERE {days} >= {PRIOR_MIN_DAYS}", distinct=True)
    return (
        f"SELECT s.culture_id, COUNT(DISTINCT s.genus) AS num_prior_orgs, MIN(s.days) AS days_since_prior_org "
        f"FROM ({src}) s WHERE s.genus IN ({_quoted(COMMON_PRIOR_GENERA)}) GROUP BY s.culture_id"
    )


def _procedure_part(loader) -> str | None:
    table = "prior_procedures"
    if table not in loader.staged_tables():
        return None
    days = _days(loader, table)
    src = _linked(loader, table, {"procedure": "f.procedure_description", "days": days},
                  where=f"WHERE {days} >= {PRIOR_MIN_DAYS}", distinct=True)
    flags = ", ".join(f"MAX(CAST(s.procedure = '{value}' AS INTEGER)) AS {name}" for name, value in PROCEDURES.items())
    recent = f"MAX(CAST(s.days <= {PROCEDURE_WINDOW_DAYS} AND s.procedure IN ({_quoted(PROCEDURES.values())}) AS INTEGER)) AS any_proc_30d"
    return f"SELECT s.culture_id, {flags}, {recent} FROM ({src}) s GROUP BY s.culture_id"


def _nursing_home_part(loader) -> str | None:
    table = "nursing_home_visits"
    if table not in loader.staged_tables():
        return None
    days = _days(loader, table)
    src = _linked(loader, table, {"days": days}, where=f"WHERE {days} BETWEEN 0 AND {NURSING_HOME_WINDOW_DAYS}", distinct=True)
    return f"SELECT s.culture_id, 1 AS nursing_home_90d FROM ({src}) s GROUP BY s.culture_id"


def _comorbidity_part(loader) -> str | None:
    table = "comorbidity"
    if table not in loader.staged_tables():
        return None
    con, site = loader.con, loader.site
    column = COMORBIDITY_COLUMN.get(site, "comorbidity_component")
    lookup = f"{site}_lookup_comorbidity"
    con.execute(f"CREATE OR REPLACE TEMP TABLE {lookup} (raw VARCHAR, {', '.join(f'{k} BOOLEAN' for k in COMORBIDITY_RULES)})")
    rows = [(n, *[bool(re.search(p, n.lower())) for p in COMORBIDITY_RULES.values()]) for n in loader.distinct(table, column) if n is not None]
    if rows:
        con.executemany(f"INSERT INTO {lookup} VALUES ({', '.join('?' for _ in range(len(COMORBIDITY_RULES) + 1))})", rows)
    src = _linked(loader, table, {"component": f"f.{column}"}, where=f"WHERE f.{column} IS NOT NULL", distinct=True)
    flags = ", ".join(f"MAX(CAST(l.{name} AS INTEGER)) AS {name}" for name in COMORBIDITY_RULES)
    return (
        f"SELECT s.culture_id, COUNT(DISTINCT s.component) AS elixhauser_count, {flags} "
        f"FROM ({src}) s JOIN {lookup} l ON l.raw = s.component GROUP BY s.culture_id"
    )


def build_features(loader, built: dict, interpretation: str = "lab") -> str:
    """Build <site>_features_<interpretation> with culture_id and FEATURES_48, one row per culture."""
    if interpretation not in ("lab", "clsi2022"):
        raise ValueError("interpretation must be 'lab' or 'clsi2022'")
    con, site = loader.con, loader.site
    parts = {
        "abx": _antibiotic_part(loader),
        "res": _resistance_part(loader, interpretation),
        "org": _organism_part(loader),
        "proc": _procedure_part(loader),
        "nh": _nursing_home_part(loader),
        "com": _comorbidity_part(loader),
    }
    owners = {
        "abx": tuple(ANTIBIOTIC_CLASSES),
        "res": ("prior_ESBL", "prior_MRSA", "prior_CRE", "prior_VRE"),
        "org": ("num_prior_orgs", "days_since_prior_org"),
        "proc": (*PROCEDURES, "any_proc_30d"),
        "nh": ("nursing_home_90d",),
        "com": ("elixhauser_count", *COMORBIDITY_RULES),
    }
    joined = "\n".join(f"LEFT JOIN ({sql}) {alias} USING (culture_id)" for alias, sql in parts.items() if sql)
    history = ",\n".join(
        f"COALESCE({alias}.{name}, 0) AS {name}" if parts[alias] else f"0 AS {name}"
        for alias, names in owners.items() for name in names
    )
    age_case = " ".join(f"WHEN '{band}' THEN {i}" for i, band in enumerate(vocab.AGE_BANDS))
    base = f"""
        SELECT e.culture_id,
               CASE e.age_band {age_case} ELSE 4 END AS age_encoded,
               COALESCE(CAST(e.sex = 'm' AS INTEGER), 0) AS gender_male,
               COALESCE(e.adi_score, 0) AS adi_score_clean,
               COALESCE(CAST(e.adi_state_rank >= 8 AS INTEGER), 0) AS adi_high,
               CAST(e.adi_missing AS INTEGER) AS adi_missing,
               CAST(COALESCE(e.setting_raw_inpatient, FALSE) AS INTEGER) AS hosp_ward_IP,
               CAST(COALESCE(e.setting_raw_outpatient, FALSE) AS INTEGER) AS hosp_ward_OP,
               CAST(COALESCE(e.setting_raw_ed, FALSE) AS INTEGER) AS hosp_ward_ER,
               CAST(COALESCE(e.setting_raw_urgent_care, FALSE) AS INTEGER) AS hosp_ward_UC,
               CAST(e.specimen_type = 'urine' AS INTEGER) AS specimen_urine,
               CAST(e.specimen_type = 'blood' AS INTEGER) AS specimen_blood,
               CAST(e.specimen_type = 'respiratory' AS INTEGER) AS specimen_respiratory,
               {history}
        FROM {built['culture_episodes']} e
        {joined}
    """
    derived = {
        "age_65plus": "CAST(age_encoded >= 5 AS INTEGER)",
        "fq_and_ceph3": "CAST(T_fq_90d > 0 AND T_ceph3_90d > 0 AS INTEGER)",
        "multi_prior_resistant": "prior_ESBL + prior_CRE + prior_MRSA + prior_VRE",
        "any_immunocomp": "CAST(comorbid_lymphoma = 1 OR comorbid_metastatic_cancer = 1 AS INTEGER)",
        "any_invasive_proc": "CAST(proc_dialysis = 1 OR proc_mechvent = 1 OR proc_urethral = 1 OR proc_cvc = 1 AS INTEGER)",
        "prior_esbl_and_ceph3": "CAST(prior_ESBL > 0 AND T_ceph3_90d > 0 AS INTEGER)",
        "hospital_multiabx": "CAST(hosp_ward_IP > 0 AND (T_fq_90d + T_ceph3_90d + T_carb_90d + T_glyco_90d) >= 2 AS INTEGER)",
        "renal_urine": "CAST(comorbid_renal_failure > 0 AND specimen_urine = 1 AS INTEGER)",
        "nursing_home_flag": "CAST(nursing_home_90d > 0 AS INTEGER)",
    }
    select = ", ".join(
        f"CAST({derived.get(name, name)} AS {'DOUBLE' if name in DOUBLE_FEATURES else 'INTEGER'}) AS {name}" for name in FEATURES_48
    )
    table = f"{site}_features_{interpretation}"
    con.execute(f"CREATE OR REPLACE TEMP TABLE {table} AS SELECT culture_id, {select} FROM ({base}) b")
    nulls = con.execute(f"SELECT {' + '.join(f'COUNT(*) - COUNT({n})' for n in FEATURES_48)} FROM {table}").fetchone()[0]
    if nulls:
        raise ValueError(f"{table} has NULL feature values")
    return table


def summarize(con, table: str) -> dict:
    n = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    binary = [f for f in FEATURES_48 if f not in NON_BINARY_FEATURES]
    sums = con.execute(f"SELECT {', '.join(f'SUM({f})' for f in binary)} FROM {table}").fetchone()
    out = {"cultures": cell(n), "binary_prevalence": {f: rate(int(s), n) for f, s in zip(binary, sums)}, "numeric": {}}
    for f in sorted(NON_BINARY_FEATURES):
        mean, p50, nonzero = con.execute(f"SELECT AVG({f}), MEDIAN({f}), COUNT(*) FILTER (WHERE {f} <> 0) FROM {table}").fetchone()
        out["numeric"][f] = {"mean": round(float(mean), 4), "median": float(p50), "nonzero_share": rate(nonzero, n)}
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the harmonized 48-feature table for one site.")
    parser.add_argument("--site", required=True)
    parser.add_argument("--interpretation", choices=("lab", "clsi2022"), default="lab")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args(argv)
    if args.interpretation == "clsi2022" and args.site != "mgb":
        parser.error("clsi2022 is only available at mgb")

    cfg = get_site(args.site)
    run = RunManifest(f"features_{args.site}", seed=args.seed, args=vars(args))
    for table, record in sorted(json.loads((cfg.parquet_dir / "_conversion.json").read_text())["tables"].items()):
        run.add_input(f"{args.site}.{table}", cfg.csv_path(table), record["sha256"])
    loader = get_loader(args.site, connect(cfg.tmp_dir))
    built = loader.build()
    table = build_features(loader, built, args.interpretation)

    out_dir = cfg.parquet_dir.parent / "features" / f"harmonized48_{args.interpretation}"
    tmp_dir = out_dir.with_name(out_dir.name + ".tmp")
    shutil.rmtree(tmp_dir, ignore_errors=True)
    tmp_dir.mkdir(parents=True)
    loader.con.execute(f"COPY (SELECT * FROM {table} ORDER BY culture_id) TO '{tmp_dir / 'features.parquet'}' (FORMAT parquet, COMPRESSION zstd)")
    shutil.rmtree(out_dir, ignore_errors=True)
    tmp_dir.rename(out_dir)

    summary = summarize(loader.con, table)
    run.count("cultures", loader.con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
    run.set("output", str(out_dir))
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"features: wrote {out_dir} and {run.dir}")
    loader.con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
