"""Download Dryad-hosted ARMD releases into the data root.

Dryad requires an API token for file downloads. Credentials come from .env as DRYAD_CLIENT_ID and
DRYAD_CLIENT_SECRET. Files stream to <name>.part, resume after interruption, and are renamed into
place only once size and SHA-256 match data_sources.yaml. Tokens and signed URLs are never printed.

Usage: python -m armd_common.download --site utsw [--workers 4] [--files README.md,...]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from .config import SITES, get_site, load_dotenv, load_sources

API = "https://datadryad.org"
USER_AGENT = "armd-common-downloader"
CHUNK = 1 << 22
MAX_ATTEMPTS = 6


class DownloadError(RuntimeError):
    pass


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


_no_redirect = urllib.request.build_opener(_NoRedirect)


class DryadClient:
    def __init__(self, client_id: str, client_secret: str):
        self._client_id = client_id
        self._client_secret = client_secret
        self._token: str | None = None
        self._expires = 0.0
        self._lock = threading.Lock()

    def token(self) -> str:
        with self._lock:
            if self._token is None or time.time() > self._expires - 600:
                body = urllib.parse.urlencode(
                    {"client_id": self._client_id, "client_secret": self._client_secret, "grant_type": "client_credentials"}
                ).encode()
                req = urllib.request.Request(
                    f"{API}/oauth/token",
                    data=body,
                    headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": USER_AGENT},
                )
                try:
                    with urllib.request.urlopen(req, timeout=60) as resp:
                        payload = json.load(resp)
                except urllib.error.HTTPError as err:
                    raise DownloadError(f"Dryad token request failed with HTTP {err.code}; check the .env credentials") from None
                self._token = payload["access_token"]
                self._expires = time.time() + float(payload.get("expires_in", 36000))
            return self._token

    def get_json(self, path: str) -> dict:
        req = urllib.request.Request(API + path, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.load(resp)

    def list_files(self, doi: str, version: str) -> dict[str, str]:
        """Map file name to download href for one pinned dataset version."""
        encoded = urllib.parse.quote(f"doi:{doi}", safe="")
        versions = self.get_json(f"/api/v2/datasets/{encoded}/versions")["_embedded"]["stash:versions"]
        pinned = [v for v in versions if str(v.get("versionNumber")) == str(version)]
        if not pinned:
            raise DownloadError(f"version {version} of {doi} is not public on Dryad; data_sources.yaml may be stale")
        files: dict[str, str] = {}
        path = pinned[0]["_links"]["stash:files"]["href"]
        while path:
            page = self.get_json(path)
            for item in page.get("_embedded", {}).get("stash:files", []):
                files[item["path"]] = item["_links"]["stash:download"]["href"]
            path = page.get("_links", {}).get("next", {}).get("href")
        return files

    def resolve(self, href: str) -> tuple[str, bool]:
        """Return the URL that serves the bytes and whether it needs the bearer token."""
        req = urllib.request.Request(API + href, headers={"Authorization": f"Bearer {self.token()}", "User-Agent": USER_AGENT})
        try:
            resp = _no_redirect.open(req, timeout=120)
        except urllib.error.HTTPError as err:
            if err.code in (301, 302, 303, 307, 308) and err.headers.get("Location"):
                return err.headers["Location"], False
            raise DownloadError(f"Dryad refused {href} with HTTP {err.code}") from None
        resp.close()
        return API + href, True


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(CHUNK):
            h.update(block)
    return h.hexdigest()


def download_file(client: DryadClient, href: str, dest: Path, expected: dict) -> dict:
    if dest.exists():
        if dest.stat().st_size == expected["bytes"] and _sha256(dest) == expected["sha256"]:
            return {"file": dest.name, "status": "present"}
        raise DownloadError(f"{dest.name} exists but does not match data_sources.yaml; move it aside first")
    part = dest.with_name(dest.name + ".part")
    started = time.time()
    fetched = 0
    for attempt in range(MAX_ATTEMPTS):
        offset = part.stat().st_size if part.exists() else 0
        if offset > expected["bytes"]:
            part.unlink()
            offset = 0
        if offset == expected["bytes"]:
            break
        url, needs_auth = client.resolve(href)
        headers = {"User-Agent": USER_AGENT}
        if needs_auth:
            headers["Authorization"] = f"Bearer {client.token()}"
        if offset:
            headers["Range"] = f"bytes={offset}-"
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=300) as resp:
                mode = "ab" if offset and resp.status == 206 else "wb"
                with part.open(mode) as fh:
                    while block := resp.read(CHUNK):
                        fh.write(block)
                        fetched += len(block)
        except urllib.error.HTTPError as err:
            print(f"  {dest.name}: HTTP {err.code} on attempt {attempt + 1}, retrying", flush=True)
            time.sleep(5 * 2**attempt)
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as err:
            print(f"  {dest.name}: {type(err).__name__} on attempt {attempt + 1}, retrying", flush=True)
            time.sleep(5 * 2**attempt)
    if not part.exists() or part.stat().st_size != expected["bytes"]:
        raise DownloadError(f"{dest.name}: incomplete after {MAX_ATTEMPTS} attempts, rerun to resume")
    digest = _sha256(part)
    if digest != expected["sha256"]:
        part.unlink()
        raise DownloadError(f"{dest.name}: SHA-256 mismatch, partial file removed")
    part.rename(dest)
    seconds = time.time() - started
    return {"file": dest.name, "status": "downloaded", "bytes": expected["bytes"], "seconds": round(seconds, 1),
            "mb_per_s": round(fetched / 1e6 / max(seconds, 1e-6), 1)}


def main(argv: list[str] | None = None) -> int:
    dryad_sites = sorted(s for s, v in load_sources().items() if v["repository"] == "Dryad")
    parser = argparse.ArgumentParser(description="Download a Dryad-hosted ARMD release.")
    parser.add_argument("--site", required=True, choices=dryad_sites)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--files", help="comma-separated file names (default: every file in the release)")
    args = parser.parse_args(argv)

    load_dotenv()
    client_id, client_secret = os.environ.get("DRYAD_CLIENT_ID"), os.environ.get("DRYAD_CLIENT_SECRET")
    if not client_id or not client_secret:
        parser.error("set DRYAD_CLIENT_ID and DRYAD_CLIENT_SECRET in .env (see .env.example)")

    cfg = get_site(args.site)
    if not cfg.raw_dir.is_dir():
        parser.error(f"{cfg.raw_dir} does not exist; create it under the data root first")
    source = cfg.source
    expected = load_sources()[args.site]["files"]
    wanted = args.files.split(",") if args.files else sorted(expected)
    if unknown := sorted(set(wanted) - set(expected)):
        parser.error(f"not in data_sources.yaml for {args.site}: {unknown}")

    client = DryadClient(client_id, client_secret)
    remote = client.list_files(source["doi"], source["version"])
    if absent := sorted(set(wanted) - set(remote)):
        raise DownloadError(f"missing from Dryad version {source['version']}: {absent}")

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(download_file, client, remote[name], cfg.raw_dir / name, expected[name]): name for name in wanted}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(json.dumps(result), flush=True)

    log_path = cfg.raw_dir / "_download.json"
    log = json.loads(log_path.read_text()) if log_path.exists() else {"site": args.site, "events": []}
    log["source"] = source
    log["events"].append({"finished_utc": datetime.now(timezone.utc).isoformat(), "results": sorted(results, key=lambda r: r["file"])})
    log_path.write_text(json.dumps(log, indent=2))
    print(f"logged to {log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
