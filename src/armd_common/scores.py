"""Published and guideline risk rules for ESBL, approximated from the harmonized features.

Each function takes `f`, a mapping from feature name (features.FEATURES_48) to an array, and returns an
array with one value per unit. `columns` builds that mapping from a feature matrix.

tumbarello_points  The Tumbarello 2011 score (Antimicrob Agents Chemother 55:3485), approximated:
  previous hospitalization within 12 months, 3 points   proxy: prior culture history
  admission from another healthcare facility, 3 points  proxy: nursing home visit in the prior 90 days
  beta-lactam or fluoroquinolone within 3 months, 2     third-generation cephalosporin, beta-lactam/
                                                        beta-lactamase inhibitor, carbapenem or
                                                        fluoroquinolone in the prior 90 days
  Charlson score of 4 or more, 2 points                 approximate Charlson from comorbidity flags
  age 70 or older, 2 points                             proxy: age 65 or older, the nearest band
  recent urinary catheterization, 2 points              a prior urethral catheter procedure
  High risk is TUMBARELLO_CUTOFF points or more, as published.

prior_esbl                   R1: a prior third-generation-cephalosporin-resistant Enterobacterales
prior_esbl_or_recent_blfq    R2: R1, or beta-lactam or fluoroquinolone treatment in the prior 90 days
minimal_design               R3's predictors: specimen type, prior ESBL and that recent treatment
"""

from __future__ import annotations

import numpy as np

from .features import FEATURES_48

RECENT_BLFQ = ("T_ceph3_90d", "T_esp_90d", "T_carb_90d", "T_fq_90d")
CHARLSON_WEIGHTS = {
    "comorbid_heart_failure": 1,
    "comorbid_liver_disease": 1,
    "comorbid_renal_failure": 2,
    "comorbid_lymphoma": 2,
    "comorbid_metastatic_cancer": 6,
}
TUMBARELLO_CUTOFF = 3
TUMBARELLO_MAX = 14
MINIMAL_PREDICTORS = ("specimen_urine", "specimen_blood", "specimen_respiratory", "prior_ESBL", "recent_blfq")


def columns(X: np.ndarray, names=FEATURES_48) -> dict:
    return {name: X[:, i] for i, name in enumerate(names)}


def _flag(v) -> np.ndarray:
    return np.asarray(v, dtype=float) > 0


def recent_blfq(f: dict) -> np.ndarray:
    return np.any([_flag(f[c]) for c in RECENT_BLFQ], axis=0)


def prior_culture_history(f: dict) -> np.ndarray:
    return _flag(f["num_prior_orgs"]) | _flag(f["multi_prior_resistant"])


def charlson_approx(f: dict) -> np.ndarray:
    return sum(w * _flag(f[c]).astype(int) for c, w in CHARLSON_WEIGHTS.items())


def tumbarello_points(f: dict) -> np.ndarray:
    return (3 * prior_culture_history(f)
            + 3 * _flag(f["nursing_home_90d"])
            + 2 * recent_blfq(f)
            + 2 * (charlson_approx(f) >= 4)
            + 2 * _flag(f["age_65plus"])
            + 2 * _flag(f["proc_urethral"])).astype(int)


def prior_esbl(f: dict) -> np.ndarray:
    return _flag(f["prior_ESBL"])


def prior_esbl_or_recent_blfq(f: dict) -> np.ndarray:
    return prior_esbl(f) | recent_blfq(f)


def minimal_design(f: dict) -> np.ndarray:
    """R3's predictor matrix, in MINIMAL_PREDICTORS order."""
    return np.column_stack([
        _flag(f["specimen_urine"]), _flag(f["specimen_blood"]), _flag(f["specimen_respiratory"]),
        prior_esbl(f), recent_blfq(f),
    ]).astype(float)
