"""Loader for ARMD-MGB (PhysioNet 1.0.0).

MGB differs from the Stanford layout in ways that change the SQL:
  * Cohort and feature rows repeat once per linked encounter, so results are collapsed before use.
  * order_proc_id_coded is reused across patients, so the culture key includes anon_id.
  * The time column is a date only. The maintainers stated on 2026-10-07 that it is the specimen
    collection date, not the order date, despite its name.
  * Each AST row carries the lab interpretation (AST_pheno) and a CLSI 2022 reinterpretation.
  * A preliminary result is dropped when a final result exists for the same organism and agent.
  * Enzyme tests (ESBL, carbapenemase, beta-lactamase, PBP2a) are their own rows.
  * Sex is disclosed, and there is no ICU flag.
  * comorbidity, nursing_home_visits and prior_procedures carry encounter and culture IDs from a
    different mapping than the cohort. No ordering of the three ID columns links them, yet every
    (anon_id, order date) pair in those files exists in the cohort. They link by patient and order
    date instead. Their contents are relative to the order date, so same-day cultures share values.
"""

from __future__ import annotations

from .. import vocab
from .base import SiteLoader, register_value_lookup


class MGBLoader(SiteLoader):
    site = "mgb"
    order_time_columns = ("order_time_jittered_utc_shifted",)
    order_time_resolution = "day"
    sex_masked = False
    positive_row_sql = "neg_cx IS NULL"
    feature_time_columns = {
        "demographics": None,
        "ward_info": None,
        "adi_scores": None,
        "comorbidity": None,
        "nursing_home_visits": ("nursing_home_visit_culture",),
        "prior_abx": ("last_dose_to_culture",),
        "prior_micro": ("prior_AST_time_to_culture",),
        "prior_org": ("prior_org_days_to_culture",),
        "prior_procedures": ("procedure_days_culture",),
    }
    date_linked_tables = ("comorbidity", "nursing_home_visits", "prior_procedures")

    def culture_stage_sql(self) -> str:
        return f"""
            WITH c AS (
                SELECT anon_id, order_proc_id_coded,
                       COUNT(DISTINCT pat_enc_csn_id_coded) AS n_encounters,
                       CAST(MIN(order_time_jittered_utc_shifted) AS TIMESTAMP) AS order_time,
                       ANY_VALUE(culture_description) AS specimen_raw,
                       BOOL_OR(neg_cx IS NULL) AS culture_positive,
                       BOOL_OR(neg_cx IS NULL AND organism IS NULL) AS any_unidentified,
                       CAST(NULL AS VARCHAR) AS ordering_mode_raw
                FROM {self.raw('cohort')}
                GROUP BY anon_id, order_proc_id_coded
            ),
            w AS (
                SELECT anon_id, order_proc_id_coded,
                       MAX(hosp_ward_IP) = 1 AS inpatient,
                       MAX(hosp_ward_OP) = 1 AS outpatient,
                       MAX(hosp_ward_ER) = 1 AS ed,
                       CAST(NULL AS BOOLEAN) AS icu,
                       MAX(hosp_ward_UC) = 1 AS urgent_care,
                       MAX(hosp_ward_day_surg) = 1 AS day_surgery
                FROM {self.raw('ward_info')}
                GROUP BY anon_id, order_proc_id_coded
            ),
            d AS (
                SELECT anon_id, order_proc_id_coded, ANY_VALUE(age) AS age_raw, ANY_VALUE(gender) AS gender_raw
                FROM {self.raw('demographics')}
                GROUP BY anon_id, order_proc_id_coded
            ),
            a AS (
                SELECT anon_id, order_proc_id_coded,
                       ANY_VALUE(adi_score) AS adi_score,
                       ANY_VALUE(adi_state_rank) AS adi_state_rank
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
        sir = register_value_lookup(
            self.con, "mgb_lookup_sir", self.distinct("cohort", "AST_pheno") + self.distinct("cohort", "CLSI_2022_pheno"), vocab.sir
        )
        return f"""
            WITH src AS (
                SELECT anon_id, order_proc_id_coded, organism, AST_code, antibiotic, AST_pheno, CLSI_2022_pheno,
                       AST_inequality, AST_val1, AST_val2, AST_panel,
                       COALESCE(prelim_AST = 'X', FALSE) AS preliminary,
                       COUNT(*) AS n_source_rows
                FROM {self.raw('cohort')}
                -- Enzyme tests also carry an AST_code (ESBL, CARBA, BLAC, MEC) and belong to enzyme_tests.
                WHERE neg_cx IS NULL AND AST_code IS NOT NULL AND enzyme_class IS NULL
                GROUP BY ALL
            ),
            flagged AS (
                SELECT *, BOOL_OR(NOT preliminary) OVER (PARTITION BY anon_id, order_proc_id_coded, organism, AST_code) AS has_final
                FROM src
            )
            SELECT 'mgb' AS site,
                   {self.patient_id_sql('f')} AS patient_id,
                   {self.culture_id_sql('f')} AS culture_id,
                   f.organism AS organism_raw,
                   l.organism_std, l.organism_group, l.is_enterobacterales, l.idsa_esbl_species, l.organism_ambiguous,
                   f.antibiotic AS antibiotic_raw,
                   f.AST_code AS antibiotic_code,
                   f.AST_pheno AS sir_raw,
                   s1.value AS sir,
                   f.CLSI_2022_pheno AS sir_clsi2022_raw,
                   s2.value AS sir_clsi2022,
                   f.AST_inequality AS mic_inequality,
                   f.AST_val1 AS mic_value,
                   f.AST_val2 AS mic_value2,
                   f.AST_panel AS ast_method,
                   f.preliminary,
                   f.n_source_rows
            FROM flagged f
            LEFT JOIN {org} l ON l.organism_raw = f.organism
            LEFT JOIN {sir} s1 ON s1.raw = f.AST_pheno
            LEFT JOIN {sir} s2 ON s2.raw = f.CLSI_2022_pheno
            WHERE NOT f.preliminary OR NOT f.has_final
        """

    def enzyme_tests_sql(self, org: str) -> str:
        return f"""
            WITH src AS (
                SELECT anon_id, order_proc_id_coded, organism, enzyme_class, enzyme, AST_pheno, AST_panel,
                       COALESCE(prelim_AST = 'X', FALSE) AS preliminary,
                       COUNT(*) AS n_source_rows
                FROM {self.raw('cohort')}
                WHERE neg_cx IS NULL AND enzyme_class IS NOT NULL
                GROUP BY ALL
            )
            SELECT 'mgb' AS site,
                   {self.patient_id_sql('f')} AS patient_id,
                   {self.culture_id_sql('f')} AS culture_id,
                   f.organism AS organism_raw,
                   l.organism_std, l.organism_group, l.is_enterobacterales, l.idsa_esbl_species, l.organism_ambiguous,
                   f.enzyme_class,
                   f.enzyme,
                   f.AST_pheno AS enzyme_result_raw,
                   CASE lower(f.AST_pheno) WHEN 'positive' THEN 'positive' WHEN 'negative' THEN 'negative' ELSE 'other' END AS enzyme_result,
                   f.AST_panel AS ast_method,
                   f.preliminary,
                   f.n_source_rows
            FROM src f
            LEFT JOIN {org} l ON l.organism_raw = f.organism
        """
