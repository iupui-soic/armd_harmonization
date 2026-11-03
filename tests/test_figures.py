"""analysis/06_figures.py on synthetic summaries shaped like the real ones."""

import importlib.util
import json

import pytest

from armd_common.config import REPO_ROOT

pytest.importorskip("matplotlib")
pytest.importorskip("xgboost")


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "analysis" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


figures = load("figures", "06_figures.py")
bounds_script = load("e3", "05_e3_bounds.py")
from test_evaluate_ladder import evaluate, synthetic_test, trained  # noqa: E402,F401  (trained is a fixture)

FLOW = lambda e, n_eval: [  # noqa: E731
    {"step": "ordered cultures", "units": None, "cultures": 750000, "patients": 280000, "label_positive": None, "prevalence": None},
    {"step": f"evaluable for {e}", "units": n_eval, "cultures": n_eval, "patients": 200000, "label_positive": None, "prevalence": None},
    {"step": "in partition derivation", "units": n_eval // 2, "cultures": n_eval // 2, "patients": 1000, "label_positive": 500, "prevalence": 0.1},
]
COHORTS = {"sites": {s: {"cohorts": {"e1_composite": {"flow": FLOW("E1", 78000)}, "e2_composite": {"flow": FLOW("E2", 744000)}}}
                     for s in ("stanford", "mgb")}}
CASCADE = {"sites": {s: {"reporting": {"sir": {**{a: {"reported_sir": 0.5 + 0.1 * i} for i, a in enumerate(("CRO", "CAZ", "FEP", "TZP"))},
                                               "any_sentinel_reported_sir": 0.95}}} for s in ("stanford", "mgb")}}


@pytest.fixture(scope="module")
def evaluation(trained):
    result = evaluate.evaluate_site(synthetic_test(), trained, n_boot=10, seed=1)
    # Round-trip through JSON, as the figures script reads the saved summary.
    return json.loads(json.dumps({"sites": {"stanford": result}}, default=str))


def names(paths):
    return sorted(p.name for p in paths)


def test_every_figure_writes_pdf_png_and_csv(tmp_path, evaluation):
    figures.style()
    written = figures.figure_flow(COHORTS, {"E1": "e1_composite", "E2": "e2_composite"}, tmp_path)
    written += figures.figure_calibration(evaluation, tmp_path)
    written += figures.figure_decision(evaluation, tmp_path)
    written += figures.table_operating_points(evaluation, tmp_path)
    written += figures.figure_cascade(CASCADE, tmp_path)
    written += figures.figure_tipping({"sites": bounds_script.e3_bounds(evaluation)}, tmp_path)
    assert names(written) == sorted(
        [f"{stem}.{ext}" for stem in ("figure1_flow", "figure2_calibration_stanford", "figure3_decision_stanford",
                                      "figure5_cascade", "figure6_tipping") for ext in ("pdf", "png", "csv")]
        + ["table4_operating_points.csv", "table4_operating_points.md"])
    assert all(p.stat().st_size > 0 for p in written)


def test_table_view_holds_every_plotted_value(tmp_path, evaluation):
    figures.style()
    csv_path = [p for p in figures.figure_decision(evaluation, tmp_path) if p.suffix == ".csv"][0]
    lines = csv_path.read_text().splitlines()
    # Four E2 variants and two E1 variants over 40 thresholds, plus treat-all on both cohorts.
    assert len(lines) == 1 + 6 * 40 + 2 * 40


def test_withheld_flow_counts_are_marked_not_plotted(tmp_path):
    figures.style()
    cohorts = json.loads(json.dumps(COHORTS))
    cohorts["sites"]["mgb"]["cohorts"]["e1_composite"]["flow"][1]["units"] = figures.WITHHELD
    csv_path = [p for p in figures.figure_flow(cohorts, {"E1": "e1_composite", "E2": "e2_composite"}, tmp_path) if p.suffix == ".csv"][0]
    assert figures.WITHHELD in csv_path.read_text()


def test_merge_summaries_joins_sites_from_several_runs(tmp_path):
    for name, sites in (("a", ("stanford", "utsw")), ("b", ("mgb2",))):
        (tmp_path / name).mkdir()
        (tmp_path / name / "summary.json").write_text(json.dumps({"partition": "derivation", "sites": {s: {"k": name} for s in sites}}))
    merged = figures.merge_summaries([tmp_path / "a", tmp_path / "b"])
    assert set(merged["sites"]) == {"stanford", "utsw", "mgb2"} and merged["partition"] == "derivation"
    assert figures.run_paths(f"{tmp_path / 'a'},{tmp_path / 'b'}", "01_build_cohorts") == [tmp_path / "a", tmp_path / "b"]


def test_merge_summaries_selects_and_orders_sites(tmp_path):
    for name, sites in (("a", ("stanford", "mgb")), ("b", ("mgb2",))):
        (tmp_path / name).mkdir()
        (tmp_path / name / "summary.json").write_text(json.dumps({"sites": {s: {"k": name} for s in sites}}))
    paths = [tmp_path / "a", tmp_path / "b"]
    merged = figures.merge_summaries(paths, ["stanford", "mgb2"])
    assert list(merged["sites"]) == ["stanford", "mgb2"]
    with pytest.raises(SystemExit):
        figures.merge_summaries(paths, ["utsw"])
