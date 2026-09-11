from armd_common import sex_inference as si

REFERENCE = {
    ("urine", "18-24"): {"f": 900, "m": 100},
    ("urine", "85-89"): {"f": 700, "m": 300},
    ("blood", "75-84"): {"f": 450, "m": 550},
    ("respiratory", "65-74"): {"f": 400, "m": 600},
}


def test_compare_patterns_prefers_the_matching_mapping():
    code0_female = {k: {"0": v["f"], "1": v["m"]} for k, v in REFERENCE.items()}
    out = si.compare_patterns(code0_female, REFERENCE)
    assert out["strata_compared"] == 4
    assert out["distance_if_code0_female"] < out["distance_if_code0_male"]
    assert out["correlation_code0_share_with_mgb_female_share"] > 0.99

    code0_male = {k: {"0": v["m"], "1": v["f"]} for k, v in REFERENCE.items()}
    out = si.compare_patterns(code0_male, REFERENCE)
    assert out["distance_if_code0_female"] > out["distance_if_code0_male"]


def test_compare_patterns_skips_small_strata():
    reference = {("urine", "18-24"): {"f": 50, "m": 10}}
    assert si.compare_patterns({("urine", "18-24"): {"0": 60, "1": 5}}, reference) == {"strata_compared": 0}


def test_decide_requires_both_lines_to_agree():
    evidence = {
        "0": {"patients": 1000, "male_dx": 1, "female_dx": 300, "conflicting_codes": 0},
        "1": {"patients": 1000, "male_dx": 200, "female_dx": 2, "conflicting_codes": 0},
    }
    agree = si.decide(evidence, {"distance_if_code0_female": 0.1, "distance_if_code0_male": 0.3})
    assert agree["evidence_agrees"] and agree["proposed_mapping"] == {"0": "f", "1": "m"}

    disagree = si.decide(evidence, {"distance_if_code0_female": 0.3, "distance_if_code0_male": 0.1})
    assert not disagree["evidence_agrees"] and disagree["proposed_mapping"] is None

    mixed = {
        "0": {"patients": 1000, "male_dx": 100, "female_dx": 300, "conflicting_codes": 0},
        "1": {"patients": 1000, "male_dx": 50, "female_dx": 200, "conflicting_codes": 0},
    }
    assert si.decide(mixed, {"distance_if_code0_female": 0.1, "distance_if_code0_male": 0.3})["female_code_by_diagnoses"] is None
