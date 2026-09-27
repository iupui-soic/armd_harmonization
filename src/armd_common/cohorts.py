"""Estimand cohorts for Paper A (plan sections 0, 5.1 and 5.2).

E1  Status quo. Cultures that grew Enterobacterales and had at least one sentinel agent reported as
    S, I or R. Label: labels.esbl_phenotype_composite.
E2  Decision time. Every ordered culture. Negative cultures and cultures without Enterobacterales
    are Y = 0. Label: labels.esbl_decision_time. Cultures with Enterobacterales but no sentinel
    result, and cultures whose isolates are unidentified, are not evaluable and leave the cohort.
E3  Population. Not a cohort. ARMD holds cultured patients only, so E3 is analysed as bounds and a
    tipping point on E2 results.

A cohort table has one row per evaluable unit of analysis with identifiers, partition and fold, the
label and reason, subgroup columns and the 48 harmonized features. The repeat rule and the unit are
set by the spec (episodes.py). The defaults, release and culture, give one row per culture as
released. Each build also returns the counts at every filter step for the TRIPOD+AI flow diagram. The test partition is sealed, so asking for it raises
unless the caller is analysis/03_evaluate_ladder.py.
"""

from __future__ import annotations

import uuid
from dataclasses import asdict, dataclass
from pathlib import Path

from . import episodes as units_mod
from . import features, implied, labels, splits
from .provenance import sha256_file

RECENT_ANTIBIOTIC_DAYS = 3
# Antifungals do not suppress bacterial growth, so they do not count as recent antibiotics.
ANTIFUNGAL_CODES = frozenset({"5FC", "AFG", "AMB", "CAS", "FLC", "ISA", "ITC", "MCZ", "MFG", "POS", "TRB", "VRC"})
SUBGROUP_COLUMNS = ("setting_primary", "specimen_type", "age_band", "sex", "prior_culture_history", "adi_quintile")
# Covariates described for each cohort, with standardized mean differences between E2 and E1.
TABLE1_COVARIATES = (
    "age_65plus", "gender_male", "hosp_ward_IP", "hosp_ward_OP", "hosp_ward_ER", "specimen_urine", "specimen_blood",
    "specimen_respiratory", "prior_culture_history", "prior_ESBL", "T_ceph3_90d", "T_fq_90d", "nursing_home_90d",
    "adi_missing", "elixhauser_count", "num_prior_orgs",
)


@dataclass(frozen=True)
class CohortSpec:
    estimand: str
    label: str = "composite"
    nonsusceptible: bool = False
    interpretation: str = "lab"
    implied_susceptibility: str | None = None
    exclude_recent_antibiotics: bool = False
    repeat_rule: str = "release"
    unit: str = "culture"

    @property
    def name(self) -> str:
        parts = [self.estimand.lower(), self.label]
        if self.nonsusceptible:
            parts.append("nonsusceptible")
        if self.interpretation != "lab":
            parts.append(self.interpretation)
        if self.implied_susceptibility:
            parts.append(f"implied_{self.implied_susceptibility}")
        if self.exclude_recent_antibiotics:
            parts.append(f"no_antibiotics_{RECENT_ANTIBIOTIC_DAYS}d")
        if self.repeat_rule != "release":
            parts.append(self.repeat_rule)
        if self.unit != "culture":
            parts.append(self.unit)
        return "_".join(parts)

    def validate(self, site: str) -> None:
        problems = []
        if self.estimand not in ("E1", "E2"):
            problems.append("estimand must be E1 or E2")
        if self.label not in ("composite", "idsa", "enzyme"):
            problems.append("label must be composite, idsa or enzyme")
        if self.estimand == "E2" and self.label != "composite":
            problems.append("E2 is defined for the composite label only")
        if self.label == "enzyme" and site != "mgb":
            problems.append("the enzyme label exists at MGB only")
        if self.label != "composite" and self.nonsusceptible:
            problems.append("nonsusceptible applies to the composite label; the IDSA label is non-susceptibility by definition")
        if self.interpretation not in ("lab", "clsi2022"):
            problems.append("interpretation must be lab or clsi2022")
        if self.interpretation == "clsi2022" and site != "mgb":
            problems.append("the CLSI 2022 reinterpretation exists at MGB only")
        if self.implied_susceptibility not in (None, "strict", "full"):
            problems.append("implied_susceptibility must be None, strict or full")
        if self.implied_susceptibility and (self.interpretation != "lab" or self.label == "enzyme"):
            problems.append("implied susceptibility applies to lab S, I and R results only")
        if self.repeat_rule not in units_mod.REPEAT_RULES:
            problems.append(f"repeat_rule must be one of {', '.join(units_mod.REPEAT_RULES)}")
        if self.unit not in units_mod.UNITS:
            problems.append(f"unit must be one of {', '.join(units_mod.UNITS)}")
        if problems:
            raise ValueError("; ".join(problems))


PRIMARY_SPECS = (CohortSpec("E1"), CohortSpec("E2"))


@dataclass(frozen=True)
class CohortResult:
    table: str
    spec: CohortSpec
    flow: list[dict]
    provenance: dict


def recent_antibiotic_cultures(loader) -> str:
    """TEMP table of cultures with an antibacterial record 1 to 3 days before the order.

    Exposure times are whole days, so a day 0 record cannot be placed before the order. Stanford's
    prior_med and MGB's last doses on or after the order day include treatment started after the
    culture was drawn, and UTSW and ECUH hold no day 0 records. Every site therefore counts days 1 to 3
    only, from the table the antibiotic features use.
    """
    con, site = loader.con, loader.site
    name = f"{site}_recent_antibiotic_cultures"
    source = features.antibiotic_source(loader)
    if source is None:
        con.execute(f"CREATE OR REPLACE TEMP TABLE {name} (culture_id VARCHAR)")
        return name
    table, name_column, days, lookup = source
    where = f"WHERE {days} BETWEEN {features.PRIOR_MIN_DAYS} AND {RECENT_ANTIBIOTIC_DAYS}"
    src = features.linked_rows(loader, table, {"drug": f"f.{name_column}"}, where=where, distinct=True)
    antifungals = ", ".join(f"'{c}'" for c in sorted(ANTIFUNGAL_CODES))
    con.execute(
        f"CREATE OR REPLACE TEMP TABLE {name} AS SELECT DISTINCT s.culture_id FROM ({src}) s "
        f"JOIN {lookup} l ON l.raw = s.drug WHERE l.value IS NOT NULL AND l.value NOT IN ({antifungals})"
    )
    return name


def _label(con, built: dict, spec: CohortSpec) -> labels.LabelResult:
    episodes, ast = built["culture_episodes"], built["isolate_ast"]
    if spec.implied_susceptibility:
        ast = implied.apply_implied_susceptibility(con, ast, implied.stanford_rules(con), strict=spec.implied_susceptibility == "strict")
    sir_column = "sir" if spec.interpretation == "lab" else "sir_clsi2022"
    if spec.label == "enzyme":
        return labels.esbl_enzyme_class(con, episodes, built["enzyme_tests"])
    if spec.label == "idsa":
        return labels.esbl_phenotype_idsa(con, episodes, ast, sir_column=sir_column)
    if spec.estimand == "E1":
        return labels.esbl_phenotype_composite(con, episodes, ast, nonsusceptible=spec.nonsusceptible, sir_column=sir_column)
    return labels.esbl_decision_time(con, episodes, ast, nonsusceptible=spec.nonsusceptible, sir_column=sir_column)


SETTING_PRECEDENCE = ("inpatient", "ed", "urgent_care", "outpatient", "unknown")


def unit_feature_sql(name: str, alias: str = "m") -> str:
    """How a unit combines one feature across its cultures.

    Flags and counts take the maximum, so a unit has a history if any culture records it.
    days_since_prior_org is 0 when there is no prior organism, so the unit takes the most recent
    non-zero value instead.
    """
    if name == "days_since_prior_org":
        return f"COALESCE(MIN(NULLIF({alias}.{name}, 0)), 0)"
    return f"MAX({alias}.{name})"


def build_cohort(
    loader,
    built: dict,
    spec: CohortSpec,
    partitions: tuple[str, ...] = ("derivation",),
    manifest: Path | None = None,
    features_table: str | None = None,
) -> CohortResult:
    """Materialize <site>_cohort_<spec name> with one row per evaluable unit in the partitions.

    Features are measured from the order date, so they mostly agree within a unit. unit_feature_sql
    combines them, which makes specimen and setting flags mean "any culture in the unit".
    specimen_type is 'multiple' for an order event with more than one specimen type, and
    setting_primary takes the first setting present in SETTING_PRECEDENCE.
    """
    spec.validate(loader.site)
    con, site = loader.con, loader.site
    episodes = built["culture_episodes"]
    manifest = manifest or splits.manifest_path(site)
    tag = uuid.uuid4().hex[:12]
    split = splits.register_manifest(con, site, table=f"_split_{tag}", include_test="test" in partitions, path=manifest)
    label = _label(con, built, spec)
    kept = units_mod.kept_cultures(con, episodes, spec.repeat_rule, f"_kept_{tag}")
    units = units_mod.unit_labels(con, label.frame.alias, episodes, kept, spec.unit, spec.estimand, f"_units_{tag}")
    feats = features_table or features.build_features(loader, built, spec.interpretation)
    recent = recent_antibiotic_cultures(loader) if spec.exclude_recent_antibiotics else None
    in_partitions = ", ".join(f"'{p}'" for p in partitions)
    unit_key = units_mod.unit_key_sql(spec.unit)

    split_join = f"JOIN {split} s ON s.patient_id_hash = sha256(u.patient_id)"
    flow: list[dict] = []

    def record(step: str, where: str, with_split: bool) -> None:
        # Steps before the partition filter include sealed test patients, so they count no labels.
        join = split_join if with_split else ""
        label_count = "COUNT(*) FILTER (WHERE u.label = 1)" if with_split else "NULL"
        n_units, n_cultures, patients, positive = con.execute(
            f"SELECT COUNT(*), SUM(u.n_cultures), COUNT(DISTINCT u.patient_id), {label_count} FROM {units} u {join} WHERE {where}"
        ).fetchone()
        flow.append({"step": step, "units": n_units, "cultures": n_cultures or 0, "patients": patients, "label_positive": positive})

    n, patients = con.execute(f"SELECT COUNT(*), COUNT(DISTINCT patient_id) FROM {episodes}").fetchone()
    flow.append({"step": "ordered cultures", "units": None, "cultures": n, "patients": patients, "label_positive": None})
    if spec.repeat_rule != "release":
        record(f"index orders under {spec.repeat_rule}", "TRUE", with_split=False)
    conditions = ["u.evaluable", f"s.partition IN ({in_partitions})"]
    record(f"evaluable for {spec.estimand}", "u.evaluable", with_split=False)
    record(f"in partition {', '.join(partitions)}", " AND ".join(conditions), with_split=True)
    if recent:
        conditions.append(f"u.unit_id NOT IN (SELECT {unit_key} FROM {episodes} e JOIN {recent} r USING (culture_id))")
        record(f"no antibacterial record {features.PRIOR_MIN_DAYS} to {RECENT_ANTIBIOTIC_DAYS} days before the order", " AND ".join(conditions), with_split=True)

    cuts = con.execute(
        f"SELECT quantile_cont(e.adi_score, [0.2, 0.4, 0.6, 0.8]) FROM {episodes} e "
        f"JOIN {split} s ON s.patient_id_hash = sha256(e.patient_id) WHERE s.partition = 'derivation' AND e.adi_score IS NOT NULL"
    ).fetchone()[0]
    adi = "MAX(m.adi_score)"
    if cuts:
        bands = " ".join(f"WHEN {adi} <= {cut} THEN 'q{i}'" for i, cut in enumerate(cuts, start=1))
        adi_quintile = f"CASE WHEN {adi} IS NULL THEN 'missing' {bands} ELSE 'q5' END"
    else:
        adi_quintile = "'missing'"
    setting_rank = "CASE m.setting_primary " + " ".join(f"WHEN '{v}' THEN {i}" for i, v in enumerate(SETTING_PRECEDENCE)) + " END"

    table = f"{site}_cohort_{spec.name}"
    feature_columns = ", ".join(f"f.{name}" for name in features.FEATURES_48)
    feature_aggregates = ", ".join(f"{unit_feature_sql(name)} AS {name}" for name in features.FEATURES_48)
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE {table} AS
        WITH m AS (
            SELECT {unit_key} AS unit_id, e.culture_id, e.setting_primary, e.specimen_type, e.age_band, e.sex, e.adi_score,
                   {feature_columns}
            FROM {episodes} e
            JOIN {kept} k USING (culture_id)
            JOIN {feats} f USING (culture_id)
        )
        SELECT '{site}' AS site, u.patient_id, u.unit_id, MIN(m.culture_id) AS culture_id, u.n_cultures,
               s.partition, s.fold, u.label, u.reason,
               arg_min(m.setting_primary, {setting_rank}) AS setting_primary,
               CASE WHEN COUNT(DISTINCT m.specimen_type) = 1 THEN MIN(m.specimen_type) ELSE 'multiple' END AS specimen_type,
               MIN(m.age_band) AS age_band, MIN(m.sex) AS sex,
               CAST(MAX(m.num_prior_orgs) > 0 OR MAX(m.multi_prior_resistant) > 0 AS INTEGER) AS prior_culture_history,
               {adi_quintile} AS adi_quintile,
               {feature_aggregates}
        FROM {units} u
        {split_join}
        JOIN m USING (unit_id)
        WHERE {' AND '.join(conditions)}
        GROUP BY u.patient_id, u.unit_id, u.n_cultures, s.partition, s.fold, u.label, u.reason
    """)
    provenance = {
        "spec": asdict(spec),
        "name": spec.name,
        "partitions": list(partitions),
        "label": label.provenance,
        "repeat_rule": spec.repeat_rule,
        "unit": spec.unit,
        "split_manifest": {"path": str(manifest), "sha256": sha256_file(manifest)},
        "features_table": feats,
        "adi_quintile_cuts": list(cuts) if cuts else None,
    }
    return CohortResult(table, spec, flow, provenance)
