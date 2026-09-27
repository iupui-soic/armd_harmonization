"""analysis/05_e3_bounds.py on a synthetic 03 summary."""

import importlib.util

import pytest

from armd_common import bounds
from armd_common.config import REPO_ROOT
from armd_common.qc import WITHHELD

spec = importlib.util.spec_from_file_location("e3", REPO_ROOT / "analysis" / "05_e3_bounds.py")
e3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e3)

COUNTS = {"tp": 20, "fp": 60, "tn": 900, "fn": 20}
EVALUATION = {"sites": {"site": {"grid": {
    "xgboost:E2->E2": {"threshold_90_sensitivity": 0.02, "confusion": {"0.20": COUNTS, "90_sensitivity": COUNTS}},
    "xgboost:E1->E2": {"threshold_90_sensitivity": 0.3,
                       "confusion": {"0.20": {k: WITHHELD for k in COUNTS}, "90_sensitivity": COUNTS}},
    "xgboost:E1->E1": {"threshold_90_sensitivity": 0.3, "confusion": {"0.20": COUNTS}},
}}}}


def test_bounds_for_e2_cells_only():
    out = e3.e3_bounds(EVALUATION)
    assert set(out["site"]) == {"xgboost:E2->E2", "xgboost:E1->E2"}
    b = out["site"]["xgboost:E2->E2"]["0.20"]
    assert b["computable"] and set(b["by_uncultured_share"]) == {"0.10", "0.25", "0.50", "0.75"}
    assert b["by_uncultured_share"]["0.50"]["manski"] == bounds.manski_bounds(20, 60, 900, 20, 0.5)
    assert b["by_uncultured_share"]["0.50"]["tipping"]["relative_risk"] == pytest.approx(0.6)
    # The operating point uses the cell's own threshold.
    assert out["site"]["xgboost:E2->E2"]["90_sensitivity"]["threshold"] == 0.02


def test_withheld_counts_leave_bounds_uncomputed():
    out = e3.e3_bounds(EVALUATION)
    assert out["site"]["xgboost:E1->E2"]["0.20"] == {"threshold": 0.2, "computable": False}
    assert "withheld" in e3.render({"sites": out})
