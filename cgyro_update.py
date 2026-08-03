"""
Small, dependency-free update checker for the CGYRO comparison tool.

The checker prefers the latest GitHub release and falls back to the VERSION
file on the default branch.  It is deliberately synchronous so callers can
choose the appropriate execution model (the Tk UI runs it in a worker thread).
"""

import json
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


REPOSITORY_URL = "https://github.com/tdxy-liu/cgyro-plot-tools"
RELEASE_API_URL = f"https://api.github.com/repos/tdxy-liu/cgyro-plot-tools/releases/latest"
VERSION_URL = "https://raw.githubusercontent.com/tdxy-liu/cgyro-plot-tools/main/VERSION"
DEFAULT_VERSION = "0.2.0"
DEFAULT_TIMEOUT = 5.0

_VERSION_PATTERN = re.compile(
    r"(?<!\d)v?(\d+)\.(\d+)(?:\.(\d+))?(?:[-+][0-9A-Za-z.-]+)?"
)


class UpdateCheckError(RuntimeError):
    """Raised when the remote version cannot be retrieved or interpreted."""


@dataclass(frozen=True)
class UpdateInfo:
    """Result of an update check."""

    current_version: str
    latest_version: str
    update_available: bool
    release_url: str = REPOSITORY_URL
    release_name: str = ""
    published_at: str = ""
    source: str = ""


def normalize_version(value):
    """Return a comparable ``major.minor.patch`` version or ``None``."""
    match = _VERSION_PATTERN.search(str(value or ""))
    if not match:
        return None
    major, minor, patch = match.groups()
    return f"{int(major)}.{int(minor)}.{int(patch or 0)}"


def _version_key(value):
    normalized = normalize_version(value)
    if normalized is None:
        return None
    return tuple(int(part) for part in normalized.split("."))


def _read_local_version():
    version_path = Path(__file__).with_name("VERSION")
    try:
        version = normalize_version(version_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError):
        version = None
    return version or DEFAULT_VERSION


APP_VERSION = _read_local_version()


def _fetch_text(url, timeout=DEFAULT_TIMEOUT, allow_not_found=False):
    request = Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": "CGYRO-Comparison-Tool-Update-Checker",
        },
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        if allow_not_found and exc.code == 404:
            return None
        raise UpdateCheckError(f"Update service returned HTTP {exc.code}.") from exc
    except (URLError, TimeoutError, OSError) as exc:
        raise UpdateCheckError("Could not connect to the update service.") from exc


def _release_info(payload, current_version):
    try:
        release = json.loads(payload)
    except (TypeError, ValueError) as exc:
        raise UpdateCheckError("The update service returned invalid data.") from exc

    if not isinstance(release, dict):
        raise UpdateCheckError("The update service returned an unexpected response.")

    latest_version = normalize_version(release.get("tag_name") or release.get("name"))
    if latest_version is None:
        return None

    current_key = _version_key(current_version)
    latest_key = _version_key(latest_version)
    if current_key is None or latest_key is None:
        raise UpdateCheckError("Could not compare the local and remote versions.")

    return UpdateInfo(
        current_version=normalize_version(current_version) or str(current_version),
        latest_version=latest_version,
        update_available=latest_key > current_key,
        release_url=release.get("html_url") or REPOSITORY_URL,
        release_name=str(release.get("name") or ""),
        published_at=str(release.get("published_at") or ""),
        source="GitHub release",
    )


def _manifest_info(payload, current_version):
    latest_version = normalize_version(payload)
    if latest_version is None:
        raise UpdateCheckError("The remote VERSION file contains an invalid version.")

    current_key = _version_key(current_version)
    latest_key = _version_key(latest_version)
    if current_key is None or latest_key is None:
        raise UpdateCheckError("Could not compare the local and remote versions.")

    return UpdateInfo(
        current_version=normalize_version(current_version) or str(current_version),
        latest_version=latest_version,
        update_available=latest_key > current_key,
        release_url=REPOSITORY_URL,
        source="GitHub VERSION file",
    )


def check_for_updates(current_version=APP_VERSION, timeout=DEFAULT_TIMEOUT):
    """Check GitHub for a newer version.

    ``UpdateCheckError`` is raised only after both supported endpoints fail.
    The function performs no UI work and can safely be called from a worker.
    """
    current_version = normalize_version(current_version)
    if current_version is None:
        raise UpdateCheckError("The local application version is invalid.")

    errors = []

    try:
        release_payload = _fetch_text(RELEASE_API_URL, timeout=timeout, allow_not_found=True)
        if release_payload:
            release_result = _release_info(release_payload, current_version)
            if release_result is not None:
                return release_result
    except UpdateCheckError as exc:
        errors.append(exc)

    try:
        manifest_payload = _fetch_text(VERSION_URL, timeout=timeout)
        return _manifest_info(manifest_payload, current_version)
    except UpdateCheckError as exc:
        errors.append(exc)

    if errors:
        raise UpdateCheckError(
            "Unable to check for updates. Please verify your internet connection."
        ) from errors[-1]
    raise UpdateCheckError("No valid remote version was found.")
