"""Loader for ARMD-MGB 2.0.0 (PhysioNet doi:10.13026/1smb-ts67).

The 2.0.0 release follows the ARMD-Stanford layout, so it reuses the Dryad loader. What differs:
  * The release was de-identified again. anon_id and every date offset differ from 1.0.0, so the two
    releases cannot be linked and 2.0.0 is its own site with its own split.
  * The 14-day repeat rule is applied at source. A culture is kept if it is the patient's first or
    was collected more than 14 whole days after the patient's previous blood, urine or respiratory
    culture, kept or not. Same-day cultures are therefore dropped at source.
  * The time column is order_time_jittered_utc_shifted, the collection time with a per-patient shift.
  * Susceptibility is one three-valued column. The MIC, method, CLSI 2022 and enzyme fields of 1.0.0
    are gone, so there is no CLSI 2022 interpretation and no enzyme_tests table.
  * Sex is disclosed as Female, Male or Others, and ages are released in bands.
  * comorbidity, nursing_home_visits and prior_procedures share the cohort's IDs, so every feature
    table links by culture. A culture that spans two encounters appears twice and is collapsed.
  * ward_type has an all-NULL ICU flag and no urgent care or day surgery flags. The release says
    those settings are identified through ordering_mode, so they are derived from that column.
"""

from __future__ import annotations

from .dryad import DryadStandardLoader


def _mode_in(values) -> str:
    return "lower(trim(ordering_mode)) IN (" + ", ".join("'" + v + "'" for v in values) + ")"


class MGB2Loader(DryadStandardLoader):
    site = "mgb2"
    order_time_columns = ("order_time_jittered_utc_shifted",)
    order_time_resolution = "second"
    sex_masked = False
    feature_time_columns = {
        **DryadStandardLoader.feature_time_columns,
        "comorbidity": None,
        "prior_procedures": ("procedure_days_culture",),
    }
    # ordering_mode holds Day Surgery, ED Observation, Emergency, Inpatient, Observation, Outpatient,
    # Psych Inpatient, Surgery Admit and Urgent Care (staged release, 2026-11-03). Surgery Admit is an
    # admission for surgery, so only Day Surgery marks day surgery.
    cohort_setting_sql = {
        "urgent_care": _mode_in(("urgent care",)),
        "day_surgery": _mode_in(("day surgery",)),
    }
