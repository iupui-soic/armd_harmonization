"""Infer which masked sex code means female at Stanford, UTSW and ECUH.

Those releases publish sex as 0 or 1 without saying which is which. MGB discloses sex and serves as
the reference. Two independent lines of evidence are computed for each masked site.

1. Sex-specific diagnoses. For patients with each code, the share with any male-only diagnosis
   (prostate hyperplasia, prostate, testis or penis cancer, erectile dysfunction) and the share
   with any female-only diagnosis (female reproductive cancers, endometriosis, ovarian cyst,
   menstrual disorders, pregnancy). The same shares by disclosed sex at MGB show how clean the
   signal is. All four releases use the same CCSR category names.
2. Specimen and age pattern. Within each specimen type and age band, the share of code 0 is
   compared with MGB's female share in the same stratum, under both possible mappings.

A mapping is proposed only when both lines agree. It is then recorded by hand in
mappings/sex_codes.yaml with this run's evidence, and it stays marked as inferred.

Usage: python -m armd_common.sex_inference [--sites stanford,utsw,ecuh]
Writes results/sex_inference/<timestamp>/ with manifest.json, summary.json and report.md.
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from .config import DEFAULT_SEED, get_site
from .db import connect
from .provenance import RunManifest
from .qc import cell, rate
from .sites import get_loader

MASKED_SITES = ("stanford", "utsw", "ecuh")

MALE_ONLY = (
    "Hyperplasia of prostate",
    "Male reproductive system cancers - prostate",
    "Male reproductive system cancers - testis",
    "Male reproductive system cancers - penis",
    "Erectile dysfunction",
)

FEMALE_ONLY = (
    "Female reproductive system cancers - cervix",
    "Female reproductive system cancers - endometrium",
    "Female reproductive system cancers - ovary",
    "Female reproductive system cancers - uterus",
    "Female reproductive system cancers - vagina",
    "Female reproductive system cancers - vulva",
    "Endometriosis",
    "Benign ovarian cyst",
    "Menstrual disorders",
    "Uncomplicated pregnancy, delivery or puerperium",
    "Supervision of high-risk pregnancy",
    "Ectopic pregnancy and complications of ectopic pregnancy",
    "Spontaneous abortion and complications of spontaneous abortion",
    "Induced abortion and complications of termination of pregnancy",
    "Multiple gestation",
    "Prolonged pregnancy",
    "Obstetric history affecting care in pregnancy",
    "Other specified complications in pregnancy",
    "Maternal care for abnormality of pelvic organs",
    "Maternal care related to disorders of the placenta and placental implantation",
    "Maternal care related to fetal conditions",
    "Maternal intrauterine infection",
    "Complications specified during the puerperium",
    "OB-related trauma to perineum and vulva",
)

COMPONENT_COLUMN = {"stanford": "comorbidity_component", "utsw": "comorbidity_component", "ecuh": "comorbidity_component", "mgb": "category", "mgb2": "category"}

# Strata smaller than this at either site are left out of the pattern comparison.
MIN_STRATUM = 100


def _quoted(values) -> str:
    return ", ".join("'" + v.replace("'", "''") + "'" for v in values)


def diagnosis_evidence(con, loader, built: dict, code_expr: str) -> dict:
    """Per code: patients, and how many have a male-only or female-only diagnosis. Raw counts."""
    column = COMPONENT_COLUMN[loader.site]
    select, on, _ = loader.feature_link("comorbidity")
    male, female = _quoted(MALE_ONLY), _quoted(FEMALE_ONLY)
    rows = con.execute(f"""
        WITH dx AS (
            SELECT k.culture_id,
                   BOOL_OR(f.component IN ({male})) AS male_dx,
                   BOOL_OR(f.component IN ({female})) AS female_dx
            FROM (
                SELECT {select}, f.{column} AS component
                FROM {loader.raw('comorbidity')} f
                WHERE f.{column} IN ({male}, {female})
            ) f
            JOIN {built['culture_keys']} k ON {on}
            GROUP BY k.culture_id
        ),
        patients AS (
            SELECT e.patient_id,
                   MIN({code_expr}) AS code,
                   COUNT(DISTINCT {code_expr}) AS n_codes,
                   COALESCE(BOOL_OR(dx.male_dx), FALSE) AS male_dx,
                   COALESCE(BOOL_OR(dx.female_dx), FALSE) AS female_dx
            FROM {built['culture_episodes']} e
            LEFT JOIN dx USING (culture_id)
            GROUP BY e.patient_id
        )
        SELECT code, COUNT(*), COUNT(*) FILTER (WHERE male_dx), COUNT(*) FILTER (WHERE female_dx),
               COUNT(*) FILTER (WHERE n_codes > 1)
        FROM patients
        GROUP BY code
    """).fetchall()
    return {code: {"patients": n, "male_dx": m, "female_dx": f, "conflicting_codes": c} for code, n, m, f, c in rows}


def stratum_counts(con, built: dict, code_expr: str) -> dict:
    rows = con.execute(f"""
        SELECT specimen_type, age_band, {code_expr} AS code, COUNT(*)
        FROM {built['culture_episodes']} e
        WHERE {code_expr} IS NOT NULL AND age_band IS NOT NULL
        GROUP BY ALL
    """).fetchall()
    out: dict = {}
    for specimen, age, code, n in rows:
        out.setdefault((specimen, age), {})[code] = n
    return out


def compare_patterns(site_strata: dict, reference_strata: dict) -> dict:
    """Weighted distance between code 0 share and MGB female share under both mappings."""
    shares, female, weights = [], [], []
    for key, counts in sorted(site_strata.items()):
        ref = reference_strata.get(key, {})
        n_site = counts.get("0", 0) + counts.get("1", 0)
        n_ref = ref.get("f", 0) + ref.get("m", 0)
        if n_site < MIN_STRATUM or n_ref < MIN_STRATUM:
            continue
        shares.append(counts.get("0", 0) / n_site)
        female.append(ref.get("f", 0) / n_ref)
        weights.append(n_site)
    if len(shares) < 3:
        return {"strata_compared": len(shares)}
    s, f, w = np.array(shares), np.array(female), np.array(weights, dtype=float)
    return {
        "strata_compared": len(shares),
        "distance_if_code0_female": round(float(np.average(np.abs(s - f), weights=w)), 4),
        "distance_if_code0_male": round(float(np.average(np.abs(s - (1 - f)), weights=w)), 4),
        "correlation_code0_share_with_mgb_female_share": round(float(np.corrcoef(s, f)[0, 1]), 4),
    }


def dx_rates(evidence: dict, codes: tuple[str, str]) -> dict:
    return {
        code: {
            "patients": cell(evidence[code]["patients"]),
            "male_dx_share": rate(evidence[code]["male_dx"], evidence[code]["patients"]),
            "female_dx_share": rate(evidence[code]["female_dx"], evidence[code]["patients"]),
            "patients_with_conflicting_codes": cell(evidence[code]["conflicting_codes"]),
        }
        for code in codes
        if code in evidence
    }


def decide(evidence: dict, pattern: dict) -> dict:
    def share(code, key):
        e = evidence.get(code)
        return e[key] / e["patients"] if e and e["patients"] else None

    m0, m1, f0, f1 = share("0", "male_dx"), share("1", "male_dx"), share("0", "female_dx"), share("1", "female_dx")
    by_dx = None
    if None not in (m0, m1, f0, f1):
        if m0 < m1 and f0 > f1:
            by_dx = "0"
        elif m1 < m0 and f1 > f0:
            by_dx = "1"
    by_pattern = None
    if "distance_if_code0_female" in pattern:
        by_pattern = "0" if pattern["distance_if_code0_female"] < pattern["distance_if_code0_male"] else "1"
    agree = by_dx is not None and by_dx == by_pattern
    return {
        "female_code_by_diagnoses": by_dx,
        "female_code_by_pattern": by_pattern,
        "evidence_agrees": agree,
        "proposed_mapping": ({by_dx: "f", "1" if by_dx == "0" else "0": "m"} if agree else None),
    }


def render(summary: dict) -> str:
    ref = summary["reference_mgb"]["diagnoses"]
    lines = [
        "# Masked sex code inference",
        "",
        "MGB reference, by disclosed sex:",
        "",
        "| sex | patients | male-only diagnosis | female-only diagnosis |",
        "|---|---|---|---|",
    ]
    lines += [f"| {k} | {v['patients']} | {v['male_dx_share']} | {v['female_dx_share']} |" for k, v in ref.items()]
    for site, s in summary["sites"].items():
        d, p = s["decision"], s["pattern"]
        lines += [
            "",
            f"## {site}",
            "",
            "| code | patients | male-only diagnosis | female-only diagnosis | patients with conflicting codes |",
            "|---|---|---|---|---|",
        ]
        lines += [
            f"| {k} | {v['patients']} | {v['male_dx_share']} | {v['female_dx_share']} | {v['patients_with_conflicting_codes']} |"
            for k, v in s["diagnoses"].items()
        ]
        lines += [
            "",
            f"Pattern over {p.get('strata_compared')} specimen and age strata: distance {p.get('distance_if_code0_female')} if code 0 is "
            f"female, {p.get('distance_if_code0_male')} if code 0 is male, correlation {p.get('correlation_code0_share_with_mgb_female_share')}.",
            "",
            f"Female code by diagnoses: {d['female_code_by_diagnoses']}. By pattern: {d['female_code_by_pattern']}. "
            f"Agree: {d['evidence_agrees']}. Proposed mapping: {d['proposed_mapping']}.",
        ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Infer the masked sex code mapping from aggregate evidence.")
    parser.add_argument("--sites", default=",".join(MASKED_SITES))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args(argv)
    run = RunManifest("sex_inference", seed=args.seed, args=vars(args))

    def open_site(site):
        cfg = get_site(site)
        for table, record in sorted(json.loads((cfg.parquet_dir / "_conversion.json").read_text())["tables"].items()):
            run.add_input(f"{site}.{table}", cfg.csv_path(table), record["sha256"])
        loader = get_loader(site, connect(cfg.tmp_dir, threads=4, memory_limit="24GB"))
        return loader, loader.build()

    loader, built = open_site("mgb")
    reference_dx = diagnosis_evidence(loader.con, loader, built, "e.sex")
    reference_strata = stratum_counts(loader.con, built, "e.sex")
    loader.con.close()
    summary = {"reference_mgb": {"diagnoses": dx_rates(reference_dx, ("f", "m"))}, "sites": {}}

    for site in args.sites.split(","):
        loader, built = open_site(site)
        evidence = diagnosis_evidence(loader.con, loader, built, "e.sex_code_masked")
        pattern = compare_patterns(stratum_counts(loader.con, built, "e.sex_code_masked"), reference_strata)
        decision = decide(evidence, pattern)
        summary["sites"][site] = {"diagnoses": dx_rates(evidence, ("0", "1")), "pattern": pattern, "decision": decision}
        run.set(f"{site}_decision", decision)
        loader.con.close()
        print(f"{site}: {decision}", flush=True)

    run.write()
    (run.dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (run.dir / "report.md").write_text(render(summary))
    print(f"sex_inference: wrote {run.dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
