"""Loaders for releases that follow the ARMD-Stanford layout: Stanford, UTSW and ECUH.

UTSW and ECUH reuse the Stanford file structure. A column that a release lacks, such as the ICU
ward flag or ordering_mode, becomes NULL instead of failing the build.

Stanford specifics that matter for harmonization:
  * One row per culture in demographics, ward_info and adi_scores, and one encounter per culture.
  * Positivity is was_positive. Negative cultures have a single row with no organism or result.
  * Sex is released as a masked 0 or 1 code, so it goes to sex_code_masked and sex stays NULL.
  * adi_score is an integer national percentile.
"""

from __future__ import annotations

from .. import vocab
from .base import SiteLoader, register_value_lookup

WARD_FLAGS = {
    "inpatient": "hosp_ward_IP",
    "outpatient": "hosp_ward_OP",
    "ed": "hosp_ward_ER",
    "icu": "hosp_ward_ICU",
    "urgent_care": "hosp_ward_UC",
    "day_surgery": "hosp_ward_day_surg",
}


class DryadStandardLoader(SiteLoader):
    order_time_column = "order_time_jittered_utc"
    order_time_resolution = "second"
    sex_masked = True
    positive_row_sql = "TRY_CAST(was_positive AS INTEGER) = 1"
    feature_time_columns = {
        "demographics": None,
        "ward_info": None,
        "adi_scores": None,
        "vitals": None,
        "labs": "Period_Day",
        "antibiotic_class_exposure": "time_to_culturetime",
        "antibiotic_subtype_exposure": "medication_time_to_cultureTime",
        "comorbidity": "comorbidity_component_start_days_culture",
        "microbial_resistance": "resistant_time_to_culturetime",
        "nursing_home_visits": "nursing_home_visit_culture",
        "prior_infecting_organism": "prior_infecting_organism_days_to_culutre",
        "prior_med": "medication_time_to_culturetime",
        "prior_procedures": "procedure_time_to_culturetime",
    }

    def culture_stage_sql(self) -> str:
        cohort_columns = set(self.columns("cohort"))
        ward_columns = set(self.columns("ward_info"))
        mode = "ANY_VALUE(ordering_mode)" if "ordering_mode" in cohort_columns else "CAST(NULL AS VARCHAR)"
        flags = ",\n                       ".join(
            f"MAX(TRY_CAST({col} AS INTEGER)) = 1 AS {key}" if col in ward_columns else f"CAST(NULL AS BOOLEAN) AS {key}"
            for key, col in WARD_FLAGS.items()
        )
        return f"""
            WITH c AS (
                SELECT anon_id, order_proc_id_coded,
                       COUNT(DISTINCT pat_enc_csn_id_coded) AS n_encounters,
                       CAST(MIN({self.order_time_column}) AS TIMESTAMP) AS order_time,
                       ANY_VALUE(culture_description) AS specimen_raw,
                       BOOL_OR({self.positive_row_sql}) AS culture_positive,
                       BOOL_OR({self.positive_row_sql} AND organism IS NULL) AS any_unidentified,
                       {mode} AS ordering_mode_raw
                FROM {self.raw('cohort')}
                GROUP BY anon_id, order_proc_id_coded
            ),
            w AS (
                SELECT anon_id, order_proc_id_coded,
                       {flags}
                FROM {self.raw('ward_info')}
                GROUP BY anon_id, order_proc_id_coded
            ),
            d AS (
                SELECT anon_id, order_proc_id_coded, ANY_VALUE(age) AS age_raw, CAST(ANY_VALUE(gender) AS VARCHAR) AS gender_raw
                FROM {self.raw('demographics')}
                GROUP BY anon_id, order_proc_id_coded
            ),
            a AS (
                SELECT anon_id, order_proc_id_coded,
                       TRY_CAST(ANY_VALUE(adi_score) AS DOUBLE) AS adi_score,
                       TRY_CAST(ANY_VALUE(adi_state_rank) AS DOUBLE) AS adi_state_rank
                FROM {self.raw('adi_scores')}
                GROUP BY anon_id, order_proc_id_coded
            )
            SELECT c.*, w.inpatient, w.outpatient, w.ed, w.icu, w.urgent_care, w.day_surgery,
                   d.age_raw, d.gender_raw, a.adi_score, a.adi_state_rank
            FROM c
            LEFT JOIN w USING (anon_id, order_proc_id_coded)
            LEFT JOIN d USING (anon_id, order_proc_id_coded)
            LEFT JOIN a USING (anon_id, order_proc_id_coded)
        """

    def isolate_ast_sql(self, org: str) -> str:
        abx = self.lookup("antibiotic", "cohort", "antibiotic", vocab.antibiotic_code)
        sir = register_value_lookup(self.con, f"{self.site}_lookup_sir", self.distinct("cohort", "susceptibility"), vocab.sir)
        return f"""
            SELECT '{self.site}' AS site,
                   {self.patient_id_sql('c')} AS patient_id,
                   {self.culture_id_sql('c')} AS culture_id,
                   c.organism AS organism_raw,
                   l.organism_std, l.organism_group, l.is_enterobacterales, l.idsa_esbl_species, l.organism_ambiguous,
                   c.antibiotic AS antibiotic_raw,
                   ab.value AS antibiotic_code,
                   c.susceptibility AS sir_raw,
                   si.value AS sir,
                   NULL AS sir_clsi2022_raw,
                   NULL AS sir_clsi2022,
                   NULL AS mic_inequality,
                   NULL AS mic_value,
                   NULL AS mic_value2,
                   NULL AS ast_method,
                   FALSE AS preliminary,
                   COUNT(*) AS n_source_rows
            FROM {self.raw('cohort')} c
            LEFT JOIN {org} l ON l.organism_raw = c.organism
            LEFT JOIN {abx} ab ON ab.raw = c.antibiotic
            LEFT JOIN {sir} si ON si.raw = c.susceptibility
            WHERE {self.positive_row_sql} AND c.antibiotic IS NOT NULL
            GROUP BY ALL
        """


class StanfordLoader(DryadStandardLoader):
    site = "stanford"


class UTSWLoader(DryadStandardLoader):
    site = "utsw"


class ECUHLoader(DryadStandardLoader):
    site = "ecuh"
