"""Bounds for the population estimand E3 (plan section 5.6).

ARMD holds cultured events only, so E3 is not identified. This script takes the confusion counts that
analysis/03_evaluate_ladder.py publishes for each model on the E2 test cohort and reports, for each
share u of empiric prescribing events that has no culture:
- Manski worst-case bounds on population prevalence, PPV and NPV (bounds.manski_bounds);
- the relative risk of resistance among uncultured events at which net benefit changes sign, under
  the assumptions stated in bounds.tipping_relative_risk.

Both are computed at every net-benefit threshold of the plan and at the 90%-sensitivity operating
point. Counts that 03 withheld leave the bounds for that threshold uncomputed. The uncultured shares
need justification from the literature before the paper uses them.

It reads only the aggregate results/03_evaluate_ladder/<timestamp>/summary.json, never row-level data.
results/05_e3_bounds/<timestamp>/ receives manifest.json, summary.json and report.md.

Usage: python analysis/05_e3_bounds.py [--evaluate-run results/03_evaluate_ladder/<timestamp>] [--seed 20260911]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from armd_common import bounds
from armd_common.config import DEFAULT_SEED, RESULTS_DIR
from armd_common.provenance import RunManifest, sha256_file

UNCULTURED_SHARES = (0.10, 0.25, 0.50)


def thresholds_of(cell: dict) -> dict:
    """Threshold label -> threshold value, from a 03 grid cell."""
    out = {}
    for name in cell["confusion"]:
        out[name] = cell["threshold_90_sensitivity"] if name == "90_sensitivity" else float(name)
    return out


def cell_bounds(cell: dict, shares=UNCULTURED_SHARES) -> dict:
    out = {}
    for name, t in thresholds_of(cell).items():
        c = cell["confusion"][name]
        if not all(isinstance(c[k], int) for k in ("tp", "fp", "tn", "fn")):
            out[name] = {"threshold": t, "computable": False}
            continue
        out[name] = {
            "threshold": t,
            "computable": True,
            "by_uncultured_share": {
                f"{u:.2f}": {"manski": bounds.manski_bounds(c["tp"], c["fp"], c["tn"], c["fn"], u),
                             "tipping": bounds.tipping_relative_risk(c["tp"], c["fp"], c["tn"], c["fn"], u, t)}
                for u in shares
            },
        }
    return out


def e3_bounds(evaluation: dict, shares=UNCULTURED_SHARES) -> dict:
    """Bounds for every cell of every site evaluated on the E2 test cohort."""
    out = {}
    for site, s in evaluation["sites"].items():
        out[site] = {name: cell_bounds(c, shares) for name, c in s["grid"].items() if name.endswith("->E2")}
    return out


def render(summary: dict) -> str:
    lines = ["# 05_e3_bounds", "", "Bounds on the E2 test cohort for uncultured shares u. PPV bounds at each threshold, "
             "and the tipping relative risk.", ""]
    for site, cells in summary["sites"].items():
        lines += [f"## {site}", "", "| cell | threshold | u | PPV lower | PPV upper | tipping RR |", "|---|---|---|---|---|---|"]
        for name, by_threshold in cells.items():
            for label, b in by_threshold.items():
                if not b["computable"]:
                    lines.append(f"| {name} | {label} | all | withheld | withheld | withheld |")
                    continue
                for u, v in b["by_uncultured_share"].items():
                    lo, hi = v["manski"]["ppv"]
                    rr = v["tipping"]["relative_risk"]
                    lines.append(f"| {name} | {label} | {u} | {lo:.3f} | {hi:.3f} | {'none' if rr is None else f'{rr:.2f}'} |")
        lines.append("")
    return "\n".join(lines) + "\n"


def latest_run(script: str) -> Path:
    runs = sorted(p for p in (RESULTS_DIR / script).iterdir() if (p / "summary.json").exists()) if (RESULTS_DIR / script).is_dir() else []
    if not runs:
        raise SystemExit(f"no {script} run with a summary.json")
    return runs[-1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evaluate-run", help="results/03_evaluate_ladder/<timestamp> (default: the latest)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    source = Path(args.evaluate_run) if args.evaluate_run else latest_run("03_evaluate_ladder")
    evaluation = json.loads((source / "summary.json").read_text())
    run = RunManifest("05_e3_bounds", seed=args.seed, args=vars(args))
    run.add_input("03_evaluate_ladder.summary", source / "summary.json", sha256_file(source / "summary.json"))
    summary = {"evaluate_run": str(source), "uncultured_shares": list(UNCULTURED_SHARES), "sites": e3_bounds(evaluation)}
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (run.dir / "report.md").write_text(render(summary))
    print(f"bounds: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
