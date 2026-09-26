"""Repeat cultures and the unit of analysis.

docs/HARMONIZATION.md, "Repeat cultures and the unit of analysis", gives the evidence for every
option here. Both choices apply before any label is read.

Repeat rule, applied to order dates:
  release            every culture as released
  chained_14d        an order date is an index if the patient has no culture 1 to 13 days earlier,
                     whether or not that earlier culture was kept. This reproduces the Stanford, UTSW
                     and ECUH releases exactly, so it changes only MGB.
  first_per_patient  each patient's first order date only

Unit of analysis:
  culture        one row per culture, as released
  specimen_day   one row per patient, specimen type and order date, so blood culture sets count once
  order_event    one row per patient and order date, pooling specimen types

The order date is the finest grain every site shares. MGB records dates only, and at the Dryad sites
every same-day culture of a patient shares one order time.

Labels for a unit with several cultures: the unit is 1 if any culture is 1, and 0 only if every
culture in it is evaluable. Otherwise it is not evaluable. Under an E1 label, cultures that are
negative or grew no target organism are not E1 cultures, so they leave the unit instead of making it
non-evaluable.
"""

from __future__ import annotations

from .labels import REASONS

REPEAT_RULES = ("release", "chained_14d", "first_per_patient")
UNITS = ("culture", "specimen_day", "order_event")
REPEAT_WINDOW_DAYS = 14

# E1 reasons that mean a culture is not an E1 culture at all.
OUTSIDE_E1 = ("culture_negative", "organism_absent")
# When several cultures could name a unit's reason, the first in this order wins.
REASON_PRIORITY = ("resistant", "not_resistant", "not_tested", "organism_unidentified", "organism_absent", "culture_negative")
if set(REASON_PRIORITY) != set(REASONS):
    raise ImportError("REASON_PRIORITY must list every label reason")


def order_date_sql(alias: str = "e") -> str:
    return f"CAST({alias}.order_time AS DATE)"


def unit_key_sql(unit: str, alias: str = "e") -> str:
    if unit not in UNITS:
        raise ValueError(f"unit must be one of {UNITS}")
    day = f"strftime({order_date_sql(alias)}, '%Y-%m-%d')"
    return {
        "culture": f"{alias}.culture_id",
        "specimen_day": f"{alias}.patient_id || ':' || {day} || ':' || {alias}.specimen_type",
        "order_event": f"{alias}.patient_id || ':' || {day}",
    }[unit]


def kept_cultures(con, episodes: str, rule: str, name: str) -> str:
    """TEMP table `name` with the culture_id of every culture the repeat rule keeps."""
    if rule not in REPEAT_RULES:
        raise ValueError(f"repeat rule must be one of {REPEAT_RULES}")
    if rule == "release":
        con.execute(f"CREATE OR REPLACE TEMP TABLE {name} AS SELECT culture_id FROM {episodes}")
        return name
    missing = con.execute(f"SELECT COUNT(*) - COUNT(order_time) FROM {episodes}").fetchone()[0]
    if missing:
        raise ValueError("every culture needs an order time before a repeat rule can apply")
    days = f"SELECT DISTINCT patient_id, {order_date_sql()} AS order_date FROM {episodes} e"
    if rule == "chained_14d":
        # The gap runs to the previous order date of any culture, kept or not, so the window chains.
        index = f"""
            SELECT patient_id, order_date FROM (
                SELECT patient_id, order_date,
                       order_date - LAG(order_date) OVER (PARTITION BY patient_id ORDER BY order_date) AS gap
                FROM ({days})
            ) WHERE gap IS NULL OR gap >= {REPEAT_WINDOW_DAYS}"""
    else:
        index = f"SELECT patient_id, MIN(order_date) AS order_date FROM ({days}) GROUP BY patient_id"
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE {name} AS
        SELECT e.culture_id FROM {episodes} e
        JOIN ({index}) i ON i.patient_id = e.patient_id AND i.order_date = {order_date_sql()}
    """)
    return name


def unit_labels(con, labels_table: str, episodes: str, kept: str, unit: str, estimand: str, name: str) -> str:
    """TEMP table `name` with one row per unit of kept cultures.

    Columns: unit_id, patient_id, n_cultures, label, evaluable, reason. With unit='culture' every
    unit is one culture and label, evaluable and reason equal the culture's own.
    """
    if estimand not in ("E1", "E2"):
        raise ValueError("estimand must be E1 or E2")
    outside = f"l.reason IN ({', '.join(repr(r) for r in OUTSIDE_E1)})" if estimand == "E1" else "FALSE"
    rank = "CASE l.reason " + " ".join(f"WHEN '{r}' THEN {i}" for i, r in enumerate(REASON_PRIORITY)) + " END"
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE {name} AS
        WITH members AS (
            SELECT {unit_key_sql(unit)} AS unit_id, e.patient_id, l.label, l.evaluable, l.reason,
                   {outside} AS outside, {rank} AS rank
            FROM {episodes} e
            JOIN {kept} k USING (culture_id)
            JOIN {labels_table} l USING (culture_id)
        ),
        units AS (
            SELECT unit_id, MIN(patient_id) AS patient_id, COUNT(*) AS n_cultures,
                   COALESCE(BOOL_OR(label = 1), FALSE) AS hit,
                   COUNT(*) FILTER (WHERE NOT outside) AS n_inside,
                   COALESCE(BOOL_AND(evaluable) FILTER (WHERE NOT outside), FALSE) AS inside_evaluable,
                   arg_min(reason, rank) FILTER (WHERE NOT outside) AS inside_reason,
                   arg_min(reason, rank) FILTER (WHERE NOT outside AND NOT evaluable) AS blocking_reason,
                   arg_min(reason, rank) FILTER (WHERE outside) AS outside_reason
            FROM members
            GROUP BY unit_id
        )
        SELECT unit_id, patient_id, n_cultures,
               CAST(CASE WHEN hit THEN 1 WHEN n_inside > 0 AND inside_evaluable THEN 0 END AS TINYINT) AS label,
               hit OR (n_inside > 0 AND inside_evaluable) AS evaluable,
               CASE WHEN hit THEN 'resistant'
                    WHEN n_inside = 0 THEN outside_reason
                    WHEN inside_evaluable THEN inside_reason
                    ELSE blocking_reason END AS reason
        FROM units
    """)
    return name
