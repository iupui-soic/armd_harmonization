"""Figures for OrderGap from saved results (plan section 5.8).

Reads only aggregate summaries, never row-level data:
  1  cohort flow, ordered cultures to E2 to E1        results/01_build_cohorts
  2  train and evaluate grid, calibration curves      results/03_evaluate_ladder
  3  decision curves on the E1 and E2 test cohorts    results/03_evaluate_ladder
  4  operating points at 90% sensitivity (a table)    results/03_evaluate_ladder
  5  sentinel agent reporting by site                 results/04_cascade_sensitivity
  6  E3 tipping point against the uncultured share    results/05_e3_bounds

A figure whose source run does not exist yet is skipped. Each figure is written as PDF and PNG with
a CSV of the plotted values beside it, which is the table view. Figures 2, 3 and 6 show the primary
model (XGBoost). The CSVs hold every model class. Withheld values are left out of the plots and
marked in the CSVs.

Colors follow the entity in every figure. Model variants and sites take fixed categorical slots,
validated for color vision deficiency with the dataviz palette checks. Reference lines are neutral.

Output: results/06_figures/<timestamp>/.

Usage: python analysis/06_figures.py [--cohorts-run ...] [--evaluate-run ...] [--cascade-run ...]
       [--bounds-run ...] [--repeat-rule release] [--unit culture] [--seed 20260911]
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402

from armd_common.cohorts import CohortSpec  # noqa: E402
from armd_common.config import DEFAULT_SEED, RESULTS_DIR  # noqa: E402
from armd_common.episodes import REPEAT_RULES, UNITS  # noqa: E402
from armd_common.provenance import RunManifest, sha256_file  # noqa: E402
from armd_common.qc import WITHHELD  # noqa: E402

SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
# Categorical slots 1 to 4 of the reference palette, which pass the adjacent-pair checks in this order.
VARIANT_COLORS = {"E2": "#2a78d6", "E1": "#eb6834", "E1_platt": "#1baf7a", "E1_isotonic": "#eda100"}
VARIANT_LABELS = {"E2": "E2-trained", "E1": "E1-trained", "E1_platt": "E1-trained, Platt", "E1_isotonic": "E1-trained, isotonic"}
SITE_COLORS = {"stanford": "#2a78d6", "mgb": "#eb6834", "mgb2": "#eb6834", "utsw": "#1baf7a", "ecuh": "#eda100"}
SITE_LABELS = {"stanford": "Stanford", "mgb": "MGB", "mgb2": "MGB", "utsw": "UTSW", "ecuh": "ECUH"}
PRIMARY_MODEL = "xgboost"
AGENTS = ("CRO", "CAZ", "FEP", "TZP")
TIPPING_THRESHOLD = "0.10"


def style() -> None:
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif", "font.size": 9, "text.color": INK_2,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK_2, "axes.titlecolor": INK,
        "axes.titlesize": 10, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
        "lines.linewidth": 1.5, "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
        "legend.frameon": False, "legend.fontsize": 8,
    })


def number(v):
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def save(fig, out: Path, name: str, rows: list[dict]) -> list[Path]:
    out.mkdir(parents=True, exist_ok=True)
    paths = [out / f"{name}.pdf", out / f"{name}.png"]
    fig.savefig(paths[0], bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight", dpi=200)
    plt.close(fig)
    paths.append(write_csv(out / f"{name}.csv", rows))
    return paths


def write_csv(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


# Figure 1 ----------------------------------------------------------------------------------------

def figure_flow(cohorts: dict, names: dict, out: Path) -> list[Path]:
    rows = []
    for site, s in cohorts["sites"].items():
        e1, e2 = s["cohorts"][names["E1"]]["flow"], s["cohorts"][names["E2"]]["flow"]
        steps = [("Ordered cultures", e2[0], "cultures"),
                 ("Evaluable for E2", next(x for x in e2 if x["step"] == "evaluable for E2"), "units"),
                 ("Evaluable for E1", next(x for x in e1 if x["step"] == "evaluable for E1"), "units")]
        for label, step, key in steps:
            value = step.get(key) if step.get(key) is not None else step.get("cultures")
            rows.append({"site": site, "step": label, "flow_step": step["step"], "count": value})
    sites = list(cohorts["sites"])
    fig, axes = plt.subplots(1, len(sites), figsize=(2.6 * len(sites), 2.2), sharey=True, squeeze=False)
    labels = ["Ordered cultures", "Evaluable for E2", "Evaluable for E1"]
    for ax, site in zip(axes[0], sites):
        values = [next(r["count"] for r in rows if r["site"] == site and r["step"] == label) for label in labels]
        top = max(v for v in values if number(v) is not None)
        for i, v in enumerate(values):
            if number(v) is None:
                ax.text(0, i, "withheld", va="center", fontsize=7.5, color=MUTED)
                continue
            ax.barh(i, v, height=0.45, color=SITE_COLORS.get(site, VARIANT_COLORS["E2"]))
            ax.text(v + top * 0.02, i, f"{v:,}", va="center", fontsize=7.5, color=INK_2)
        ax.set_title(SITE_LABELS.get(site, site), loc="left")
        ax.set_xlim(0, top * 1.35)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: "0" if x == 0 else f"{x / 1000:,.0f}k"))
        ax.grid(axis="y", visible=False)
    axes[0][0].set_yticks(range(len(labels)), labels)
    axes[0][0].invert_yaxis()
    fig.suptitle("Figure 1. From ordered cultures to the E2 and E1 cohorts", x=0.02, ha="left", color=INK, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    return save(fig, out, "figure1_flow", rows)


# Figures 2 to 4 ----------------------------------------------------------------------------------

def cell_name(model: str, variant: str, test: str) -> str:
    return f"{model}:{variant}->{test}"


def figure_calibration(evaluation: dict, out: Path, model: str = PRIMARY_MODEL) -> list[Path]:
    rows, paths = [], []
    for site, s in evaluation["sites"].items():
        fig, axes = plt.subplots(2, 2, figsize=(5.6, 5.2), squeeze=False)
        for i, train in enumerate(("E1", "E2")):
            for j, test in enumerate(("E1", "E2")):
                ax = axes[i][j]
                variants = [train] + (["E1_platt", "E1_isotonic"] if (train, test) == ("E1", "E2") else [])
                top = 0.0
                for v in variants:
                    c = s["grid"].get(cell_name(model, v, test))
                    if not c:
                        continue
                    curve = c["curves"]["calibration_curve"]
                    x, y = curve["predicted"], curve["observed"]
                    top = max(top, max(x), max(y))
                    ax.plot(x, y, color=VARIANT_COLORS[v], label=VARIANT_LABELS[v])
                    for xv, yv in zip(x, y):
                        rows.append({"site": site, "model": model, "trained": train, "tested": test, "variant": v,
                                     "predicted": xv, "observed": yv})
                top = top * 1.05 or 1.0
                ax.plot([0, top], [0, top], color=MUTED, linewidth=0.8)
                ax.set_xlim(0, top)
                ax.set_ylim(0, top)
                ax.set_title(f"Trained on {train}, tested on {test}", loc="left", fontsize=9)
                if i == 1:
                    ax.set_xlabel("Predicted risk")
                if j == 0:
                    ax.set_ylabel("Observed risk")
        handles = [plt.Line2D([], [], color=VARIANT_COLORS[v], label=VARIANT_LABELS[v]) for v in VARIANT_COLORS]
        fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.03))
        fig.suptitle(f"Figure 2. Calibration across the train and evaluate grid, {SITE_LABELS.get(site, site)}",
                     x=0.02, ha="left", color=INK, fontsize=10)
        fig.tight_layout(rect=(0, 0.04, 1, 0.97))
        paths += save(fig, out, f"figure2_calibration_{site}", [r for r in rows if r["site"] == site])
    return paths


def figure_decision(evaluation: dict, out: Path, model: str = PRIMARY_MODEL) -> list[Path]:
    rows, paths = [], []
    for site, s in evaluation["sites"].items():
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), squeeze=False)
        for ax, test, variants in ((axes[0][0], "E1", ("E1", "E2")), (axes[0][1], "E2", ("E1", "E1_platt", "E1_isotonic", "E2"))):
            treat_all = None
            peak = 0.0
            for v in variants:
                c = s["grid"].get(cell_name(model, v, test))
                if not c:
                    continue
                dc = c["curves"]["decision_curve"]
                treat_all = c["curves"]["treat_all"]
                pts = [(float(t), nb) for t, nb in dc.items() if number(nb) is not None]
                for t, nb in dc.items():
                    rows.append({"site": site, "model": model, "tested": test, "variant": v, "threshold": t,
                                 "net_benefit": nb if number(nb) is not None else WITHHELD})
                if pts:
                    xs, ys = zip(*pts)
                    # Curves converge near zero, so the legend and the CSV carry identity, not end labels.
                    ax.plot(xs, ys, color=VARIANT_COLORS[v], label=VARIANT_LABELS[v])
                    peak = max(peak, max(ys))
            if treat_all:
                xs = [float(t) for t in treat_all]
                ys = list(treat_all.values())
                ax.plot(xs, ys, color=MUTED, linewidth=1.0, label="Treat all")
                for t, nb in treat_all.items():
                    rows.append({"site": site, "model": "treat_all", "tested": test, "variant": "", "threshold": t, "net_benefit": nb})
                peak = max(peak, max(ys))
            ax.axhline(0, color=AXIS, linewidth=0.8)
            ax.set_ylim(-peak if peak > 0 else -0.01, peak * 1.15 if peak > 0 else 0.01)
            ax.set_xlim(0, 0.4)
            ax.set_title(f"{test} test cohort", loc="left", fontsize=9)
            ax.set_xlabel("Threshold probability")
        axes[0][0].set_ylabel("Net benefit")
        handles = [plt.Line2D([], [], color=VARIANT_COLORS[v], label=VARIANT_LABELS[v]) for v in VARIANT_COLORS]
        handles.append(plt.Line2D([], [], color=MUTED, linewidth=1.0, label="Treat all"))
        fig.legend(handles=handles, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.08))
        fig.suptitle(f"Figure 3. Decision curves, {SITE_LABELS.get(site, site)}", x=0.02, ha="left", color=INK, fontsize=10)
        fig.tight_layout(rect=(0, 0.05, 1, 0.95))
        paths += save(fig, out, f"figure3_decision_{site}", [r for r in rows if r["site"] == site])
    return paths


def table_operating_points(evaluation: dict, out: Path) -> list[Path]:
    rows = []
    for site, s in evaluation["sites"].items():
        for name, c in s["grid"].items():
            row = {"site": site, "cell": name, "threshold": c["threshold_90_sensitivity"], "units": c["units"], "events": c["events"]}
            for k in ("sensitivity", "specificity", "ppv", "npv", "net_benefit"):
                key = f"{k}_at_90_sensitivity"
                point, interval = c["point"].get(key), c["interval"].get(key)
                row[k] = point
                row[f"{k}_95ci"] = interval if number(point) is not None else WITHHELD
            rows.append(row)
    lines = ["| site | cell | threshold | sensitivity | specificity | PPV | NPV | net benefit |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        cells = [f"{r[k]:.3f}" if number(r[k]) is not None else "withheld" for k in ("sensitivity", "specificity", "ppv", "npv", "net_benefit")]
        lines.append(f"| {r['site']} | {r['cell']} | {r['threshold']:.4f} | " + " | ".join(cells) + " |")
    out.mkdir(parents=True, exist_ok=True)
    (out / "table4_operating_points.md").write_text("\n".join(lines) + "\n")
    return [write_csv(out / "table4_operating_points.csv", rows), out / "table4_operating_points.md"]


# Figure 5 ----------------------------------------------------------------------------------------

def figure_cascade(cascade: dict, out: Path) -> list[Path]:
    rows = []
    sites = list(cascade["sites"])
    items = list(AGENTS) + ["any_sentinel_reported_sir"]
    for site in sites:
        r = cascade["sites"][site]["reporting"]["sir"]
        for agent in items:
            value = r[agent] if agent == "any_sentinel_reported_sir" else r[agent]["reported_sir"]
            rows.append({"site": site, "agent": "Any sentinel" if agent == items[-1] else agent,
                         "share_reported": value if number(value) is not None else WITHHELD})
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    width = 0.8 / len(sites)
    labels = list(AGENTS) + ["Any sentinel"]
    for k, site in enumerate(sites):
        for i, label in enumerate(labels):
            v = next(r["share_reported"] for r in rows if r["site"] == site and r["agent"] == label)
            if number(v) is not None:
                ax.bar(i + (k - (len(sites) - 1) / 2) * width, v, width=width * 0.85, color=SITE_COLORS.get(site))
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.set_ylabel("Enterobacterales cultures with a result")
    ax.grid(axis="x", visible=False)
    handles = [matplotlib.patches.Patch(color=SITE_COLORS.get(s), label=SITE_LABELS.get(s, s)) for s in sites]
    ax.legend(handles=handles, loc="upper left", ncol=len(sites), bbox_to_anchor=(0, 1.12))
    fig.suptitle("Figure 5. Sentinel agents reported as S, I or R, derivation partition", x=0.02, y=1.06, ha="left", color=INK, fontsize=10)
    return save(fig, out, "figure5_cascade", rows)


# Figure 6 ----------------------------------------------------------------------------------------

def figure_tipping(bounds_summary: dict, out: Path, model: str = PRIMARY_MODEL, threshold: str = TIPPING_THRESHOLD) -> list[Path]:
    rows = []
    sites = list(bounds_summary["sites"])
    fig, axes = plt.subplots(1, len(sites), figsize=(max(6.4, 2.7 * len(sites)), 2.8), sharey=True, squeeze=False)
    for ax, site in zip(axes[0], sites):
        any_line = False
        for v in VARIANT_COLORS:
            b = bounds_summary["sites"][site].get(cell_name(model, v, "E2"), {}).get(threshold)
            if not b or not b["computable"]:
                continue
            pts = []
            for u, entry in b["by_uncultured_share"].items():
                rr = entry["tipping"]["relative_risk"]
                rows.append({"site": site, "model": model, "variant": v, "threshold": threshold, "uncultured_share": u,
                             "tipping_relative_risk": rr, "net_benefit_cultured": entry["tipping"]["net_benefit_cultured"]})
                if rr is not None:
                    pts.append((float(u), rr))
            if pts:
                any_line = True
                xs, ys = zip(*pts)
                ax.plot(xs, ys, color=VARIANT_COLORS[v], marker="o", markersize=4, markeredgecolor=SURFACE,
                        markeredgewidth=1.0, label=VARIANT_LABELS[v])
        ax.axhline(1, color=MUTED, linewidth=0.8)
        ax.annotate("equal risk", (1, 1), xycoords=("axes fraction", "data"), xytext=(-2, 3), textcoords="offset points",
                    ha="right", va="bottom", fontsize=7, color=MUTED)
        ax.set_ylim(bottom=0)
        if not any_line:
            ax.text(0.5, 0.5, "No sign change\nin the feasible range", transform=ax.transAxes, ha="center", va="center",
                    fontsize=7.5, color=MUTED)
        ax.set_title(SITE_LABELS.get(site, site), loc="left")
        ax.set_xlabel("Share of events uncultured")
    axes[0][0].set_ylabel("Tipping relative risk")
    handles = [plt.Line2D([], [], color=VARIANT_COLORS[v], marker="o", markersize=4, label=VARIANT_LABELS[v]) for v in VARIANT_COLORS]
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.1))
    fig.suptitle(f"Figure 6. Tipping relative risk for net benefit at threshold {threshold}", x=0.02, ha="left", color=INK, fontsize=10)
    fig.tight_layout(rect=(0, 0.06, 1, 0.92))
    return save(fig, out, "figure6_tipping", rows)


# Main --------------------------------------------------------------------------------------------

def latest_run(script: str) -> Path | None:
    folder = RESULTS_DIR / script
    runs = sorted(p for p in folder.iterdir() if (p / "summary.json").exists()) if folder.is_dir() else []
    return runs[-1] if runs else None


def run_paths(spec: str | None, script: str) -> list[Path]:
    """Run folders named in a comma-separated list, or the latest run. Several runs are merged by site,
    later runs winning, so a site built later joins the frozen runs of the others."""
    if spec:
        return [Path(p) for p in spec.split(",")]
    latest = latest_run(script)
    return [latest] if latest else []


def merge_summaries(paths: list[Path], sites: list[str] | None = None) -> dict:
    merged: dict = {}
    for path in paths:
        summary = json.loads((path / "summary.json").read_text())
        merged = {**merged, **summary, "sites": {**merged.get("sites", {}), **summary.get("sites", {})}}
    if sites is not None:
        missing = [s for s in sites if s not in merged["sites"]]
        if missing:
            raise SystemExit(f"no results for {missing} in {[str(p) for p in paths]}")
        merged["sites"] = {s: merged["sites"][s] for s in sites}
    return merged


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for script in ("cohorts", "evaluate", "cascade", "bounds"):
        parser.add_argument(f"--{script}-run", help="a results run folder, or several joined by commas and merged by site (default: the latest)")
    parser.add_argument("--sites", help="comma-separated sites to plot, in order (default: every site in the runs)")
    parser.add_argument("--repeat-rule", default="release", choices=REPEAT_RULES)
    parser.add_argument("--unit", default="culture", choices=UNITS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    style()
    run = RunManifest("06_figures", seed=args.seed, args=vars(args))
    out = run.dir
    sources = {
        "cohorts": run_paths(args.cohorts_run, "01_build_cohorts"),
        "evaluate": run_paths(args.evaluate_run, "03_evaluate_ladder"),
        "cascade": run_paths(args.cascade_run, "04_cascade_sensitivity"),
        "bounds": run_paths(args.bounds_run, "05_e3_bounds"),
    }
    loaded = {}
    for key, paths in sources.items():
        if not paths:
            print(f"figures: no {key} run yet, skipping its figures")
            continue
        for path in paths:
            run.add_input(f"{key}.summary:{path.name}", path / "summary.json", sha256_file(path / "summary.json"))
        loaded[key] = merge_summaries(paths, args.sites.split(",") if args.sites else None)
    written = []
    if "cohorts" in loaded:
        names = {e: CohortSpec(e, repeat_rule=args.repeat_rule, unit=args.unit).name for e in ("E1", "E2")}
        written += figure_flow(loaded["cohorts"], names, out)
    if "evaluate" in loaded:
        written += figure_calibration(loaded["evaluate"], out)
        written += figure_decision(loaded["evaluate"], out)
        written += table_operating_points(loaded["evaluate"], out)
    if "cascade" in loaded:
        written += figure_cascade(loaded["cascade"], out)
    if "bounds" in loaded:
        written += figure_tipping(loaded["bounds"], out)
    run.set("outputs", [p.name for p in written])
    run.write()
    print(f"figures: wrote {len(written)} files to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
