"""File hashing and run manifests.

Every script records the git SHA, seed, input hashes, package versions and row counts at each
filter step. Counts are aggregates, so manifests never hold row-level data.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from .config import DEFAULT_SEED, REPO_ROOT, RESULTS_DIR

_PACKAGES = ("duckdb", "polars", "pyarrow", "numpy", "scikit-learn", "xgboost", "statsmodels", "pyyaml")


def sha256_file(path: Path, cache_path: Path | None = None, chunk: int = 1 << 24) -> str:
    """Hash a file. A cache keyed by size and mtime avoids rehashing multi-GB releases."""
    path = Path(path)
    stat = path.stat()
    key = f"{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}"
    cache = {}
    if cache_path is not None and Path(cache_path).exists():
        cache = json.loads(Path(cache_path).read_text())
        if key in cache:
            return cache[key]
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    digest = h.hexdigest()
    if cache_path is not None:
        cache[key] = digest
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        Path(cache_path).write_text(json.dumps(cache, indent=1))
    return digest


def git_state() -> dict:
    def run(*args):
        return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True).stdout.strip()

    sha = run("rev-parse", "HEAD") or None
    dirty = bool(run("status", "--porcelain", "--untracked-files=no"))
    return {"sha": sha, "dirty": dirty}


def require_git_tag(tag: str, pinned: Path | None = None) -> dict:
    """Stop the run unless the repository has `tag`.

    `pinned` is a local file the run depends on but the repository does not track, such as an
    analysis plan. Its SHA-256 goes into the returned record, so the manifest pins its contents.
    """
    tags = subprocess.run(["git", "tag", "-l", tag], cwd=REPO_ROOT, capture_output=True, text=True).stdout.split()
    if tag not in tags:
        raise SystemExit(f"refusing to run: the git tag {tag} does not exist")
    return {"tag": tag, "pinned_file": None if pinned is None else Path(pinned).name,
            "pinned_sha256": sha256_file(pinned) if pinned is not None and Path(pinned).exists() else None}


def package_versions() -> dict:
    versions = {"python": sys.version.split()[0]}
    for name in _PACKAGES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


class RunManifest:
    """Collects provenance for one script run and writes results/<script>/<timestamp>/manifest.json."""

    def __init__(self, script: str, seed: int = DEFAULT_SEED, args: dict | None = None):
        self.started = datetime.now(timezone.utc)
        self.script = script
        self.dir = RESULTS_DIR / script / self.started.strftime("%Y%m%dT%H%M%SZ")
        self.data = {
            "script": script,
            "started_utc": self.started.isoformat(),
            "seed": seed,
            "args": args or {},
            "git": git_state(),
            "packages": package_versions(),
            "inputs": {},
            "row_counts": [],
        }

    def add_input(self, label: str, path: Path, digest: str) -> None:
        self.data["inputs"][label] = {"path": str(path), "sha256": digest}

    def count(self, step: str, n: int, **context) -> None:
        self.data["row_counts"].append({"step": step, "n": int(n), **context})

    def set(self, key: str, value) -> None:
        self.data[key] = value

    def write(self) -> Path:
        self.dir.mkdir(parents=True, exist_ok=True)
        self.data["finished_utc"] = datetime.now(timezone.utc).isoformat()
        out = self.dir / "manifest.json"
        out.write_text(json.dumps(self.data, indent=2, default=str))
        return out
