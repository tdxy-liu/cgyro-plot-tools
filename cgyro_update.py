"""
Small, dependency-free update checker for the CGYRO comparison tool.

The checker prefers the latest GitHub release and falls back to the VERSION
file on the default branch.  It is deliberately synchronous so callers can
choose the appropriate execution model (the Tk UI runs it in a worker thread).
"""

import json
import re
import subprocess
import sys
import time
import argparse
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


REPOSITORY_URL = "https://github.com/tdxy-liu/cgyro-plot-tools"
RELEASE_API_URL = f"https://api.github.com/repos/tdxy-liu/cgyro-plot-tools/releases/latest"
VERSION_URL = "https://raw.githubusercontent.com/tdxy-liu/cgyro-plot-tools/main/VERSION"
DEFAULT_VERSION = "0.2.9"
DEFAULT_TIMEOUT = 5.0

_VERSION_PATTERN = re.compile(
    r"(?<!\d)v?(\d+)\.(\d+)(?:\.(\d+))?(?:[-+][0-9A-Za-z.-]+)?"
)


class UpdateCheckError(RuntimeError):
    """Raised when the remote version cannot be retrieved or interpreted."""


class GitUpdateError(RuntimeError):
    """Raised when a safe fast-forward Git update cannot be completed."""


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


@dataclass(frozen=True)
class GitUpdateResult:
    """Summary of a command-line Git update."""

    repo_dir: str
    remote: str
    branch: str
    updated: bool
    previous_commit: str
    current_commit: str
    output: str = ""


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


def _cache_busted_url(url):
    """Avoid stale CDN responses when reading the GitHub version manifest."""
    separator = "&" if "?" in str(url) else "?"
    return f"{url}{separator}cb={time.time_ns()}"


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
        manifest_payload = _fetch_text(_cache_busted_url(VERSION_URL), timeout=timeout)
        return _manifest_info(manifest_payload, current_version)
    except UpdateCheckError as exc:
        errors.append(exc)

    if errors:
        raise UpdateCheckError(
            "Unable to check for updates. Please verify your internet connection."
        ) from errors[-1]
    raise UpdateCheckError("No valid remote version was found.")


def _run_git(repo_dir, args, timeout):
    """Run one Git command and return its completed process safely."""
    # Use subprocess' working-directory support instead of ``git -C``.
    # Some OMFIT installations ship an older Git that does not understand
    # the global ``-C`` option.
    command = ["git", *[str(arg) for arg in args]]
    try:
        return subprocess.run(
            command,
            cwd=str(repo_dir),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise GitUpdateError("Git was not found on PATH.") from exc
    except (OSError, subprocess.SubprocessError) as exc:
        raise GitUpdateError(f"Could not run Git: {exc}") from exc


def _git_output(result):
    """Return trimmed command output for errors and CLI display."""
    return str(result.stdout or "").strip()


def update_from_git(
    repo_dir=None,
    remote="origin",
    branch="main",
    timeout=60.0,
):
    """Fast-forward a clean checkout from GitHub.

    The default repository is the directory containing this module.  The
    update is intentionally conservative: it requires the requested branch,
    refuses a dirty worktree, and runs ``git pull --ff-only`` so local commits
    or edits are never overwritten silently.
    """
    if repo_dir is None:
        repo_dir = Path(__file__).resolve().parent
    else:
        repo_dir = Path(repo_dir).expanduser().resolve()

    remote = str(remote or "origin").strip()
    branch = str(branch or "main").strip()
    if not remote:
        raise GitUpdateError("The Git remote name cannot be empty.")
    if not branch:
        raise GitUpdateError("The Git branch name cannot be empty.")

    repo_check = _run_git(repo_dir, ["rev-parse", "--show-toplevel"], timeout)
    if repo_check.returncode != 0:
        detail = _git_output(repo_check)
        raise GitUpdateError(
            "The application directory is not a Git checkout."
            + (f"\n{detail}" if detail else "")
        )
    repo_root = Path(_git_output(repo_check) or str(repo_dir)).resolve()

    branch_check = _run_git(repo_root, ["symbolic-ref", "HEAD"], timeout)
    branch_ref = _git_output(branch_check)
    branch_prefix = "refs/heads/"
    current_branch = (
        branch_ref[len(branch_prefix):]
        if branch_ref.startswith(branch_prefix)
        else ""
    )
    if branch_check.returncode != 0 or current_branch != branch:
        current = current_branch or "detached HEAD"
        raise GitUpdateError(
            f"The current Git branch is '{current}', but the requested update branch is '{branch}'."
        )

    status_check = _run_git(
        repo_root,
        ["status", "--porcelain", "--untracked-files=normal"],
        timeout,
    )
    if status_check.returncode != 0:
        detail = _git_output(status_check)
        raise GitUpdateError(
            "Could not inspect the Git worktree."
            + (f"\n{detail}" if detail else "")
        )
    if _git_output(status_check):
        raise GitUpdateError(
            "Local Git changes were found; commit or stash them before updating."
        )

    before_check = _run_git(repo_root, ["rev-parse", "HEAD"], timeout)
    if before_check.returncode != 0:
        detail = _git_output(before_check)
        raise GitUpdateError(
            "Could not determine the current Git commit."
            + (f"\n{detail}" if detail else "")
        )
    previous_commit = _git_output(before_check)

    pull_result = _run_git(
        repo_root,
        ["pull", "--ff-only", remote, branch],
        timeout,
    )
    output = _git_output(pull_result)
    if pull_result.returncode != 0:
        raise GitUpdateError(output or "The Git update failed.")

    after_check = _run_git(repo_root, ["rev-parse", "HEAD"], timeout)
    if after_check.returncode != 0:
        detail = _git_output(after_check)
        raise GitUpdateError(
            "The update completed, but the new Git commit could not be verified."
            + (f"\n{detail}" if detail else "")
        )
    current_commit = _git_output(after_check)

    return GitUpdateResult(
        repo_dir=str(repo_root),
        remote=remote,
        branch=branch,
        updated=previous_commit != current_commit,
        previous_commit=previous_commit,
        current_commit=current_commit,
        output=output,
    )


def _build_cli_parser():
    """Build the standalone command-line update parser."""
    parser = argparse.ArgumentParser(
        description="Check for or safely fast-forward the CGYRO Comparison Tool from GitHub."
    )
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument(
        "--check-update",
        "--check",
        dest="check_update",
        action="store_true",
        help="check the latest GitHub version without changing files",
    )
    actions.add_argument(
        "--update",
        action="store_true",
        help="pull the selected branch with git pull --ff-only",
    )
    parser.add_argument(
        "--repo",
        default=None,
        help="Git checkout to update (default: the directory containing this script)",
    )
    parser.add_argument("--remote", default="origin", help="Git remote name (default: origin)")
    parser.add_argument("--branch", default="main", help="Git branch (default: main)")
    parser.add_argument(
        "--timeout",
        type=float,
        default=60.0,
        help="timeout in seconds for each Git/network operation (default: 60)",
    )
    return parser


def _run_cli(argv=None):
    """Run the standalone update/check command and return an exit code."""
    parser = _build_cli_parser()
    args = parser.parse_args(argv)

    if args.check_update:
        try:
            result = check_for_updates(APP_VERSION, timeout=args.timeout)
        except UpdateCheckError as exc:
            print(f"Update check failed: {exc}", file=sys.stderr)
            return 1
        print(f"Current version: {result.current_version}")
        print(f"Latest version:  {result.latest_version}")
        print(f"Update available: {'yes' if result.update_available else 'no'}")
        print(f"Source: {result.source}")
        if result.release_url:
            print(f"Release: {result.release_url}")
        return 0

    try:
        result = update_from_git(
            repo_dir=args.repo,
            remote=args.remote,
            branch=args.branch,
            timeout=args.timeout,
        )
    except GitUpdateError as exc:
        print(f"Update failed: {exc}", file=sys.stderr)
        return 1

    if result.updated:
        print(f"Updated {result.branch} from {result.previous_commit[:12]} to {result.current_commit[:12]}.")
    else:
        print(f"Already up to date at {result.current_commit[:12]}.")
    if result.output:
        print(result.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(_run_cli())
