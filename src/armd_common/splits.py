"""Patient-level derivation and test split with inner folds.

Every culture follows its patient, so assigning folds to patients is leakage-safe by construction.
With one row per patient, StratifiedGroupKFold reduces to StratifiedKFold, which is what runs here.

One of five outer folds becomes the 20% test partition. The remaining patients get five inner folds.
Splits are generated once per site with a fixed seed and written as a gzip CSV manifest of
(site, patient_id_hash, partition, fold). patient_id_hash is the SHA-256 of the site-prefixed
patient_id, so anyone holding the same release can rebuild the mapping.

The test partition is sealed (CLAUDE.md rule 3). read_manifest drops test rows unless the caller
is analysis/03_evaluate_ladder.py.

Usage: python -m armd_common.splits --site stanford [--seed 20260911]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import sys
from pathlib import Path

import numpy as np
from sklearn.model_selection import StratifiedKFold

from .config import DEFAULT_SEED, REPO_ROOT, SITES

OUTER_FOLDS = 5
INNER_FOLDS = 5
MANIFEST_DIR = REPO_ROOT / "benchmarks" / "manifests"
SEALED_READER = "03_evaluate_ladder.py"
COLUMNS = ("site", "patient_id_hash", "partition", "fold")

# Patient strata: 0 has no evaluable E1 culture, 1 has one but never resistant, 2 was ever resistant.
STRATA_SQL = """
    SELECT e.patient_id,
           CASE WHEN BOOL_OR(l.label = 1) THEN 2 WHEN BOOL_OR(l.evaluable) THEN 1 ELSE 0 END AS stratum
    FROM {episodes} e JOIN {labels} l USING (culture_id)
    GROUP BY e.patient_id
"""


class SealedPartitionError(PermissionError):
    pass


def patient_hash(patient_id: str) -> str:
    return hashlib.sha256(patient_id.encode("utf-8")).hexdigest()


def assign(patient_ids, strata, seed: int = DEFAULT_SEED) -> list[tuple[str, str, int | None]]:
    """Return (patient_id, partition, fold) per patient. Input order does not affect the result."""
    ids = np.asarray(patient_ids, dtype=object)
    y = np.asarray(strata)
    if len(set(ids.tolist())) != len(ids):
        raise ValueError("patient_ids must be unique")
    order = np.argsort(ids.astype(str), kind="stable")
    ids, y = ids[order], y[order]

    outer = StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=seed)
    derivation_idx, test_idx = next(outer.split(ids, y))
    inner = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True, random_state=seed + 1)
    fold = np.zeros(len(ids), dtype=int)
    for k, (_, held_out) in enumerate(inner.split(ids[derivation_idx], y[derivation_idx]), start=1):
        fold[derivation_idx[held_out]] = k

    test = np.zeros(len(ids), dtype=bool)
    test[test_idx] = True
    return [(pid, "test" if is_test else "derivation", None if is_test else int(f)) for pid, is_test, f in zip(ids.tolist(), test, fold)]


def manifest_path(site: str) -> Path:
    return MANIFEST_DIR / f"{site}.csv.gz"


def write_manifest(site: str, assignments, path: Path | None = None, force: bool = False) -> Path:
    path = path or manifest_path(site)
    if path.exists() and not force:
        raise FileExistsError(f"{path} is frozen; pass force=True only if the split must be regenerated")
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted((site, patient_hash(pid), partition, "" if fold is None else fold) for pid, partition, fold in assignments)
    with gzip.open(path, "wt", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    return path


def _caller_may_read_test() -> bool:
    return Path(sys.argv[0]).name == SEALED_READER


def read_manifest(site: str, include_test: bool = False, path: Path | None = None) -> list[dict]:
    if include_test and not _caller_may_read_test():
        raise SealedPartitionError(f"only analysis/{SEALED_READER} may read the test partition")
    with gzip.open(path or manifest_path(site), "rt", newline="") as fh:
        rows = list(csv.DictReader(fh))
    return [r for r in rows if include_test or r["partition"] != "test"]


def register_manifest(con, site: str, table: str | None = None, include_test: bool = False, path: Path | None = None) -> str:
    """Load a manifest into DuckDB so it joins to culture_episodes on sha256(patient_id)."""
    table = table or f"{site}_split"
    rows = read_manifest(site, include_test=include_test, path=path)
    con.execute(f"CREATE OR REPLACE TEMP TABLE {table} (site VARCHAR, patient_id_hash VARCHAR, partition VARCHAR, fold INTEGER)")
    con.executemany(
        f"INSERT INTO {table} VALUES (?, ?, ?, ?)",
        [(r["site"], r["patient_id_hash"], r["partition"], int(r["fold"]) if r["fold"] else None) for r in rows],
    )
    return table


def main(argv: list[str] | None = None) -> int:
    from . import labels
    from .config import get_site
    from .db import connect
    from .sites import get_loader

    parser = argparse.ArgumentParser(description="Generate and freeze a site's split manifest.")
    parser.add_argument("--site", required=True, choices=SITES)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--force", action="store_true", help="overwrite a frozen manifest")
    args = parser.parse_args(argv)

    loader = get_loader(args.site, connect(get_site(args.site).tmp_dir))
    built = loader.build()
    e1 = labels.esbl_phenotype_composite(loader.con, built["culture_episodes"], built["isolate_ast"])
    strata = loader.con.execute(STRATA_SQL.format(episodes=built["culture_episodes"], labels=e1.frame.alias)).fetchall()
    path = write_manifest(args.site, assign([p for p, _ in strata], [s for _, s in strata], seed=args.seed), force=args.force)
    print(f"wrote {path} for {len(strata)} patients")
    return 0


if __name__ == "__main__":
    sys.exit(main())
