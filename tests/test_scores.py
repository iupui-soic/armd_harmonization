"""Tumbarello approximation and the simple rules, on hand-built feature rows."""

import numpy as np

from armd_common import scores
from armd_common.features import FEATURES_48


def rows(*overrides):
    """One feature row per dict of overrides, every other feature 0."""
    X = np.zeros((len(overrides), len(FEATURES_48)), dtype=np.float32)
    for i, o in enumerate(overrides):
        for name, v in o.items():
            X[i, FEATURES_48.index(name)] = v
    return scores.columns(X)


def test_tumbarello_items_and_proxies():
    f = rows(
        {},
        {"num_prior_orgs": 2},                                   # prior culture history: 3
        {"nursing_home_90d": 1},                                 # facility: 3
        {"T_fq_90d": 1},                                         # recent beta-lactam or fluoroquinolone: 2
        {"T_glyco_90d": 1},                                      # vancomycin is not a beta-lactam: 0
        {"comorbid_metastatic_cancer": 1},                       # Charlson 6: 2
        {"comorbid_heart_failure": 1, "comorbid_renal_failure": 1},  # Charlson 3: 0
        {"age_65plus": 1},                                       # age proxy: 2
        {"proc_urethral": 1},                                    # catheter: 2
        {"multi_prior_resistant": 1, "nursing_home_90d": 1, "T_esp_90d": 1, "comorbid_lymphoma": 1,
         "comorbid_renal_failure": 1, "age_65plus": 1, "proc_urethral": 1},  # every item
    )
    assert scores.tumbarello_points(f).tolist() == [0, 3, 3, 2, 0, 2, 0, 2, 2, 14]
    assert scores.TUMBARELLO_MAX == 14


def test_rules():
    f = rows({}, {"prior_ESBL": 1}, {"T_ceph3_90d": 1}, {"T_carb_90d": 1, "prior_ESBL": 1})
    assert scores.prior_esbl(f).tolist() == [False, True, False, True]
    assert scores.prior_esbl_or_recent_blfq(f).tolist() == [False, True, True, True]


def test_minimal_design():
    f = rows({"specimen_urine": 1, "prior_ESBL": 1}, {"specimen_blood": 1, "T_fq_90d": 1})
    assert scores.minimal_design(f).tolist() == [[1, 0, 0, 1, 0], [0, 1, 0, 0, 1]]
    assert len(scores.MINIMAL_PREDICTORS) == 5
