"""Checks for analysis/mimic/code_lists.yaml: shape, code format and known memberships.

Matching is by dotless prefix, as in analysis/mimic/build.py: a code matches a list when it starts
with one of the list's entries.
"""
import re
from pathlib import Path

import pytest
import yaml

CODE_LISTS = Path(__file__).resolve().parents[1] / "analysis" / "mimic" / "code_lists.yaml"
CODE = re.compile(r"^[A-Z0-9]+$")

CHARLSON_1987_WEIGHTS = {
    "myocardial_infarction": 1,
    "congestive_heart_failure": 1,
    "peripheral_vascular_disease": 1,
    "cerebrovascular_disease": 1,
    "dementia": 1,
    "chronic_pulmonary_disease": 1,
    "rheumatic_disease": 1,
    "peptic_ulcer_disease": 1,
    "mild_liver_disease": 1,
    "diabetes_without_chronic_complication": 1,
    "diabetes_with_chronic_complication": 2,
    "hemiplegia_or_paraplegia": 2,
    "renal_disease": 2,
    "any_malignancy": 2,
    "moderate_or_severe_liver_disease": 3,
    "metastatic_solid_tumor": 6,
    "aids_hiv": 6,
}


@pytest.fixture(scope="module")
def lists():
    return yaml.safe_load(CODE_LISTS.read_text())


def matches(code, prefixes):
    return any(code.startswith(p) for p in prefixes)


def entries(lists, path, key):
    node = lists
    for part in path.split("."):
        node = node[part]
    return node[key]


def every_list(lists):
    """Yield (label, kind, entries) for each code list. kind is dx9, dx10, pr9 or pcs."""
    for section in ("elixhauser_quan2005", "charlson_quan2005"):
        for name, cat in lists[section].items():
            yield f"{section}.{name}.icd9", "dx9", cat["icd9"]
            yield f"{section}.{name}.icd10", "dx10", cat["icd10"]
    yield "heart_failure_ccsr.icd9", "dx9", lists["heart_failure_ccsr"]["icd9"]
    yield "heart_failure_ccsr.icd10", "dx10", lists["heart_failure_ccsr"]["icd10"]
    for name, proc in lists["procedures"].items():
        yield f"procedures.{name}.icd9", "pr9", proc["icd9"]
        yield f"procedures.{name}.icd10pcs", "pcs", proc["icd10pcs"]


def test_top_level_sections(lists):
    assert set(lists) == {"elixhauser_quan2005", "charlson_quan2005", "heart_failure_ccsr", "procedures"}


def test_elixhauser_has_31_categories(lists):
    elix = lists["elixhauser_quan2005"]
    assert len(elix) == 31
    for name, cat in elix.items():
        assert set(cat) == {"icd9", "icd10"}, name


def test_charlson_has_17_categories_with_original_weights(lists):
    charl = lists["charlson_quan2005"]
    assert len(charl) == 17
    for name, cat in charl.items():
        assert set(cat) == {"icd9", "icd10", "weight"}, name
    assert {name: cat["weight"] for name, cat in charl.items()} == CHARLSON_1987_WEIGHTS


def test_heart_failure_and_procedure_shapes(lists):
    assert set(lists["heart_failure_ccsr"]) == {"icd9", "icd10"}
    procs = lists["procedures"]
    assert set(procs) == {"cvc", "mechvent", "dialysis", "urethral_catheter", "surgical"}
    for name, proc in procs.items():
        assert set(proc) == {"icd9", "icd10pcs"}, name


def test_category_names_used_by_the_build(lists):
    elix = lists["elixhauser_quan2005"]
    for name in ("liver_disease", "lymphoma", "metastatic_cancer", "obesity", "renal_failure"):
        assert name in elix
    charl = lists["charlson_quan2005"]
    for name in ("mild_liver_disease", "moderate_or_severe_liver_disease", "diabetes_without_chronic_complication",
                 "diabetes_with_chronic_complication", "any_malignancy", "metastatic_solid_tumor"):
        assert name in charl


def test_every_entry_is_a_dotless_uppercase_alphanumeric_string(lists):
    for label, _, codes in every_list(lists):
        assert isinstance(codes, list) and codes, label
        for code in codes:
            assert isinstance(code, str), (label, code)
            assert CODE.match(code), (label, code)
        assert len(codes) == len(set(codes)), f"duplicate entries in {label}"


def test_entry_lengths_fit_the_code_system(lists):
    # ICD-9-CM diagnosis codes have 3 to 5 characters, ICD-9-CM procedure codes 3 to 4 (2-digit chapter
    # prefixes are allowed), ICD-10-CM codes 3 to 7 and ICD-10-PCS codes exactly 7 (shorter = prefix).
    limits = {"dx9": (3, 5), "pr9": (2, 4), "dx10": (3, 7), "pcs": (1, 7)}
    for label, kind, codes in every_list(lists):
        low, high = limits[kind]
        for code in codes:
            assert low <= len(code) <= high, (label, code)


MEMBERS = [
    ("elixhauser_quan2005.congestive_heart_failure", "icd10", "I509"),
    ("elixhauser_quan2005.congestive_heart_failure", "icd9", "4280"),
    ("charlson_quan2005.congestive_heart_failure", "icd9", "4280"),
    ("elixhauser_quan2005.metastatic_cancer", "icd10", "C787"),
    ("charlson_quan2005.metastatic_solid_tumor", "icd10", "C787"),
    ("elixhauser_quan2005.renal_failure", "icd10", "N186"),
    ("charlson_quan2005.renal_disease", "icd10", "N186"),
    ("elixhauser_quan2005.obesity", "icd10", "E6601"),
    ("elixhauser_quan2005.lymphoma", "icd10", "C8190"),
    ("elixhauser_quan2005.liver_disease", "icd10", "K7031"),
    ("charlson_quan2005.mild_liver_disease", "icd10", "K7031"),
    ("heart_failure_ccsr", "icd10", "I509"),
    ("heart_failure_ccsr", "icd10", "I110"),
    ("heart_failure_ccsr", "icd10", "I130"),
    ("heart_failure_ccsr", "icd10", "I132"),
    ("heart_failure_ccsr", "icd9", "4280"),
    ("heart_failure_ccsr", "icd9", "40291"),
    ("procedures.mechvent", "icd10pcs", "5A1955Z"),
    ("procedures.mechvent", "icd9", "9672"),
    ("procedures.dialysis", "icd10pcs", "5A1D70Z"),
    ("procedures.dialysis", "icd10pcs", "5A1D90Z"),
    ("procedures.dialysis", "icd9", "3995"),
    ("procedures.urethral_catheter", "icd10pcs", "0T9B70Z"),
    ("procedures.urethral_catheter", "icd9", "5794"),
    ("procedures.cvc", "icd10pcs", "02HV33Z"),
    ("procedures.cvc", "icd10pcs", "0JH63XZ"),
    ("procedures.cvc", "icd9", "3893"),
    ("procedures.surgical", "icd10pcs", "0DTJ4ZZ"),
    ("procedures.surgical", "icd9", "4701"),
]

NON_MEMBERS = [
    ("heart_failure_ccsr", "icd10", "I119"),  # hypertensive heart disease without heart failure
    ("heart_failure_ccsr", "icd10", "I10"),
    ("heart_failure_ccsr", "icd9", "40200"),  # hypertensive heart disease without heart failure
    ("elixhauser_quan2005.congestive_heart_failure", "icd9", "4019"),
    ("elixhauser_quan2005.metastatic_cancer", "icd10", "C189"),  # primary colon cancer
    ("elixhauser_quan2005.solid_tumor_without_metastasis", "icd10", "C787"),
    ("procedures.mechvent", "icd10pcs", "5A09357"),  # CPAP, non-invasive
    ("procedures.mechvent", "icd9", "9390"),  # CPAP
    ("procedures.urethral_catheter", "icd10pcs", "0T9B30Z"),  # suprapubic
    ("procedures.cvc", "icd10pcs", "03HY32Z"),  # arterial monitoring device
    ("procedures.cvc", "icd10pcs", "05HB33Z"),  # basilic vein, a midline catheter
    ("procedures.cvc", "icd9", "3891"),  # arterial catheterization
    ("procedures.dialysis", "icd9", "3927"),  # dialysis access creation
    ("procedures.surgical", "icd10pcs", "5A1955Z"),
    ("procedures.surgical", "icd9", "9672"),
    ("procedures.surgical", "icd9", "0066"),  # chapter 00
]


@pytest.mark.parametrize("path, key, code", MEMBERS, ids=[f"{p}.{k}:{c}" for p, k, c in MEMBERS])
def test_known_member(lists, path, key, code):
    assert matches(code, entries(lists, path, key))


@pytest.mark.parametrize("path, key, code", NON_MEMBERS, ids=[f"{p}.{k}:{c}" for p, k, c in NON_MEMBERS])
def test_known_non_member(lists, path, key, code):
    assert not matches(code, entries(lists, path, key))
