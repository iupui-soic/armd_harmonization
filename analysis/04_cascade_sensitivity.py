"""Cascade and selective reporting (layer T), on the derivation partition (plan section 5.5).

For each staged site, among cultures that grew Enterobacterales in the derivation partition:
- reporting: the share with each sentinel agent reported as S, I or R;
- non-evaluable cultures: those without any sentinel result, compared with evaluable cultures on
  setting, specimen, organism group and prior history, with standardized mean differences;
- MGB cross-check: ESBL enzyme test positivity and the CLSI 2022 composite phenotype, among cultures
  with and without a sentinel result under the laboratory interpretation;
- implied susceptibility: E1 evaluable cultures and prevalence with no rules, strict rules and all
  rules (implied.py), to judge whether imputation changes the estimand.

Everything is per culture as released, because reporting practice applies to each culture. The test
partition is not read and no model is trained. results/04_cascade_sensitivity/<timestamp>/ receives
manifest.json, summary.json and report.md. Counts from 1 to 10 are suppressed, and so are counts that
would expose one by subtraction.

Usage: python analysis/04_cascade_sensitivity.py [--sites stanford,mgb] [--seed 20260911]
"""

from __future__ import annotations

import argparse
import json
import math

from armd_common import features, implied, labels, splits
from armd_common.config import DEFAULT_SEED, MIN_CELL, SITES, get_site
from armd_common.db import connect
from armd_common.provenance import RunManifest, sha256_file
from armd_common.qc import cell, publish_flows, rate, sentinel_reporting
from armd_common.sites import get_loader

ENTEROBACTERALES_GROUPS = ("ecoli", "klebsiella", "proteus", "enterobacter", "other_enterobacterales")
COMPARED = {
    "hosp_ward_IP": "f.hosp_ward_IP",
    "hosp_ward_OP": "f.hosp_ward_OP",
    "hosp_ward_ER": "f.hosp_ward_ER",
    "specimen_urine": "f.specimen_urine",
    "specimen_blood": "f.specimen_blood",
    "specimen_respiratory": "f.specimen_respiratory",
    "prior_culture_history": "CAST(f.num_prior_orgs > 0 OR f.multi_prior_resistant > 0 AS INTEGER)",
    "prior_ESBL": "f.prior_ESBL",
    "num_prior_orgs": "f.num_prior_orgs",
    **{f"organism_{g}": f"CAST(e.organism_group = '{g}' AS INTEGER)" for g in ENTEROBACTERALES_GROUPS},
}
COUNT_COVARIATES = ("num_prior_orgs",)
IMPUTATION = (("none", None), ("strict", True), ("all rules", False))


def derivation_episodes(con, site: str, episodes: str, manifest=None) -> str:
    split = splits.register_manifest(con, site, table=f"{site}_split_cascade", path=manifest)
    name = f"{site}_derivation_episodes"
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE {name} AS
        SELECT e.* FROM {episodes} e JOIN {split} s ON s.patient_id_hash = sha256(e.patient_id)
        WHERE s.partition = 'derivation'
    """)
    return name


def smd(m1, v1, m2, v2):
    if None in (m1, m2, v1, v2):
        return None
    pooled = math.sqrt((v1 + v2) / 2)
    return round((m2 - m1) / pooled, 3) if pooled > 0 else None


def evaluable_comparison(con, episodes: str, ast: str, feats: str) -> dict:
    """Enterobacterales cultures without a sentinel result against those with one."""
    label = labels.esbl_phenotype_composite(con, episodes, ast)
    exprs = ", ".join(f"AVG({sql}), VAR_SAMP({sql})" for sql in COMPARED.values())
    rows = dict((r[0], r[1:]) for r in con.execute(f"""
        SELECT l.reason = 'not_tested' AS not_tested, COUNT(*), {exprs}
        FROM {episodes} e JOIN {label.frame.alias} l USING (culture_id) JOIN {feats} f USING (culture_id)
        WHERE l.reason IN ('resistant', 'not_resistant', 'not_tested')
        GROUP BY 1
    """).fetchall())
    evaluable, not_tested = rows.get(False), rows.get(True)
    out = {"evaluable_cultures": cell(evaluable[0] if evaluable else 0),
           "not_tested_cultures": cell(not_tested[0] if not_tested else 0),
           "not_tested_share": rate(not_tested[0] if not_tested else 0,
                                    (evaluable[0] if evaluable else 0) + (not_tested[0] if not_tested else 0)),
           "covariates": {}}
    small_group = any(0 < (r[0] if r else 0) < MIN_CELL for r in (evaluable, not_tested))
    for i, name in enumerate(COMPARED):
        if not evaluable or not not_tested or small_group:
            out["covariates"][name] = None
            continue
        m1, v1 = evaluable[1 + 2 * i], evaluable[2 + 2 * i]
        m2, v2 = not_tested[1 + 2 * i], not_tested[2 + 2 * i]
        # The mean of a flag times the group size is a count, so it obeys the small-cell rule.
        if name not in COUNT_COVARIATES and any(
            m is not None and (0 < round(m * n) < MIN_CELL or 0 < n - round(m * n) < MIN_CELL)
            for m, n in ((m1, evaluable[0]), (m2, not_tested[0]))
        ):
            out["covariates"][name] = None
            continue
        out["covariates"][name] = {"evaluable_mean": None if m1 is None else round(m1, 4),
                                   "not_tested_mean": None if m2 is None else round(m2, 4),
                                   "smd_not_tested_minus_evaluable": smd(m1, v1, m2, v2)}
    return out


def mgb_cross_check(con, episodes: str, ast: str, enzymes: str) -> dict:
    """ESBL enzyme tests and the CLSI 2022 phenotype, by whether a sentinel agent had a lab result."""
    clsi = labels.esbl_phenotype_composite(con, episodes, ast, sir_column="sir_clsi2022")
    lab = labels.esbl_phenotype_composite(con, episodes, ast)
    rows = con.execute(f"""
        WITH enz AS (
            SELECT culture_id, BOOL_OR(enzyme_result = 'positive') AS positive
            FROM {enzymes}
            WHERE is_enterobacterales AND enzyme_class = 'ESBL' AND enzyme_result IN ('positive', 'negative')
            GROUP BY culture_id
        )
        SELECT l.reason <> 'not_tested' AS lab_sentinel_reported,
               COUNT(*),
               COUNT(z.culture_id), COUNT(*) FILTER (WHERE z.positive),
               COUNT(*) FILTER (WHERE c.evaluable), COUNT(*) FILTER (WHERE c.label = 1)
        FROM {episodes} e
        JOIN {lab.frame.alias} l USING (culture_id)
        JOIN {clsi.frame.alias} c USING (culture_id)
        LEFT JOIN enz z USING (culture_id)
        WHERE l.reason IN ('resistant', 'not_resistant', 'not_tested')
        GROUP BY 1 ORDER BY 1
    """).fetchall()
    out = {}
    for reported, n, tested, positive, clsi_evaluable, clsi_resistant in rows:
        out["sentinel_reported" if reported else "no_sentinel_result"] = {
            "cultures": cell(n),
            "esbl_tested": cell(tested),
            "esbl_positive_share": rate(positive, tested),
            "clsi2022_evaluable": cell(clsi_evaluable),
            "clsi2022_resistant_share": rate(clsi_resistant, clsi_evaluable),
        }
    return out


def imputation(con, episodes: str, ast: str) -> dict:
    """E1 evaluable cultures and resistant cultures with no, strict and all implied-susceptibility rules."""
    rules = implied.stanford_rules(con)
    flows = {}
    for name, strict in IMPUTATION:
        source = ast if strict is None else implied.apply_implied_susceptibility(con, ast, rules, strict=strict)
        label = labels.esbl_phenotype_composite(con, episodes, source)
        n, resistant = con.execute(
            f"SELECT COUNT(*) FILTER (WHERE evaluable), COUNT(*) FILTER (WHERE label = 1) FROM {label.frame.alias}"
        ).fetchone()
        flows[name] = [{"step": "evaluable for E1", "units": None, "cultures": n, "patients": None, "label_positive": resistant}]
    return {name: flow[0] for name, flow in publish_flows(flows).items()}


def run_site(site: str, run: RunManifest | None = None, loader=None, manifest=None) -> dict:
    if loader is None:
        cfg = get_site(site)
        conversion = json.loads((cfg.parquet_dir / "_conversion.json").read_text())
        if run is not None:
            for table, record in sorted(conversion["tables"].items()):
                run.add_input(f"{site}.{table}", cfg.csv_path(table), record["sha256"])
            run.add_input(f"{site}.split_manifest", splits.manifest_path(site), sha256_file(splits.manifest_path(site)))
        loader = get_loader(site, connect(cfg.tmp_dir))
    built = loader.build()
    con = loader.con
    episodes = derivation_episodes(con, site, built["culture_episodes"], manifest)
    feats = features.build_features(loader, built, "lab")
    out = {"reporting": {"sir": sentinel_reporting(con, {"culture_episodes": episodes, "isolate_ast": built["isolate_ast"]})}}
    if site == "mgb":
        out["reporting"]["sir_clsi2022"] = sentinel_reporting(
            con, {"culture_episodes": episodes, "isolate_ast": built["isolate_ast"]}, "sir_clsi2022")
    out["not_tested_vs_evaluable"] = evaluable_comparison(con, episodes, built["isolate_ast"], feats)
    enzymes = built.get("enzyme_tests")
    if enzymes and con.execute(f"SELECT COUNT(*) FROM {enzymes}").fetchone()[0]:
        out["mgb_cross_check"] = mgb_cross_check(con, episodes, built["isolate_ast"], built["enzyme_tests"])
    out["implied_susceptibility"] = imputation(con, episodes, built["isolate_ast"])
    if run is not None:
        run.count(f"{site}.derivation cultures", con.execute(f"SELECT COUNT(*) FROM {episodes}").fetchone()[0])
    return out


def render(summary: dict) -> str:
    lines = ["# 04_cascade_sensitivity", "", "Derivation partition, per culture as released.", ""]
    for site, s in summary["sites"].items():
        lines += [f"## {site}", "", f"Sentinel reporting among Enterobacterales cultures: {s['reporting']['sir']}", ""]
        c = s["not_tested_vs_evaluable"]
        lines += [f"Not tested: {c['not_tested_cultures']} of {c['evaluable_cultures']} evaluable plus not tested "
                  f"(share {c['not_tested_share']}).", "",
                  "| covariate | evaluable mean | not tested mean | SMD (not tested minus evaluable) |", "|---|---|---|---|"]
        for name, v in c["covariates"].items():
            lines.append(f"| {name} | " + (" | ".join(str(v[k]) for k in ("evaluable_mean", "not_tested_mean", "smd_not_tested_minus_evaluable"))
                                          if v else "withheld | withheld | withheld") + " |")
        if "mgb_cross_check" in s:
            lines += ["", "MGB cross-check:", ""] + [f"- {k}: {v}" for k, v in s["mgb_cross_check"].items()]
        lines += ["", "Implied susceptibility, E1:", ""] + [f"- {k}: {v}" for k, v in s["implied_susceptibility"].items()]
        lines.append("")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", help="comma-separated sites (default: every staged site)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    sites = args.sites.split(",") if args.sites else [s for s in SITES if (get_site(s).parquet_dir / "cohort").is_dir()]
    run = RunManifest("04_cascade_sensitivity", seed=args.seed, args=vars(args))
    summary = {"partition": "derivation", "sites": {}}
    for site in sites:
        print(f"cascade: {site}", flush=True)
        summary["sites"][site] = run_site(site, run)
    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (run.dir / "report.md").write_text(render(summary))
    print(f"cascade: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
