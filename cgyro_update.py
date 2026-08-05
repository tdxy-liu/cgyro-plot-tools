"""
Small, dependency-free update checker for the CGYRO comparison tool.

The checker prefers the latest GitHub release and falls back to the VERSION
file on the default branch.  It is deliberately synchronous so callers can
choose the appropriate execution model (the Tk UI runs it in a worker thread).
"""

import json
import http.client
import os
import re
import shlex
import socket
import ssl
import subprocess
import sys
import time
import argparse
import threading
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from functools import partial
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import HTTPSHandler, ProxyHandler, Request, build_opener, urlopen


REPOSITORY_URL = "https://github.com/tdxy-liu/cgyro-plot-tools"
RELEASE_API_URL = f"https://api.github.com/repos/tdxy-liu/cgyro-plot-tools/releases/latest"
VERSION_URL = "https://raw.githubusercontent.com/tdxy-liu/cgyro-plot-tools/main/VERSION"
DEFAULT_VERSION = "0.2.13"
DEFAULT_TIMEOUT = 5.0
DEFAULT_SSH_PORT = 22
DEFAULT_LOCAL_SOCKS_PORT = 0
UPDATE_CONFIG_ENV = "CGYRO_UPDATE_CONFIG"
UPDATE_CONFIG_NAME = "update_connection.json"

_VERSION_PATTERN = re.compile(
    r"(?<!\d)v?(\d+)\.(\d+)(?:\.(\d+))?(?:[-+][0-9A-Za-z.-]+)?"
)


class UpdateCheckError(RuntimeError):
    """Raised when the remote version cannot be retrieved or interpreted."""


class GitUpdateError(RuntimeError):
    """Raised when a safe fast-forward Git update cannot be completed."""


class UpdateConnectionError(RuntimeError):
    """Raised when the configured direct/SOCKS/SSH update connection fails."""


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


@dataclass(frozen=True)
class UpdateProxyConfig:
    """Connection settings used only while checking or downloading updates.

    ``mode`` is one of ``direct``, ``socks5`` (an already running SOCKS5
    proxy), or ``ssh-socks`` (the application starts ``ssh -D`` itself).
    SSH passwords are intentionally not represented or stored; use an SSH
    key or an ssh-agent.
    """

    mode: str = "direct"
    socks_proxy: str = ""
    ssh_host: str = ""
    ssh_user: str = ""
    ssh_port: int = DEFAULT_SSH_PORT
    local_socks_port: int = DEFAULT_LOCAL_SOCKS_PORT
    identity_file: str = ""


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


def update_connection_config_path():
    """Return the per-user update connection settings path."""
    configured = os.environ.get(UPDATE_CONFIG_ENV, "").strip()
    if configured:
        return Path(configured).expanduser()
    return Path.home() / ".cgyro_comparison_tool" / UPDATE_CONFIG_NAME


def _safe_int(value, default, minimum=None, maximum=None):
    try:
        number = int(value)
    except (TypeError, ValueError):
        number = int(default)
    if minimum is not None:
        number = max(int(minimum), number)
    if maximum is not None:
        number = min(int(maximum), number)
    return number


def normalize_update_proxy_config(value=None):
    """Convert a dict or dataclass into a safe update connection config."""
    if isinstance(value, UpdateProxyConfig):
        raw = asdict(value)
    elif isinstance(value, dict):
        raw = value
    else:
        raw = {}

    mode = str(raw.get("mode", "direct") or "direct").strip().lower()
    mode = {
        "ssh": "ssh-socks",
        "ssh_socks": "ssh-socks",
        "ssh-socks5": "ssh-socks",
        "socks": "socks5",
    }.get(mode, mode)
    if mode not in ("direct", "socks5", "ssh-socks"):
        mode = "direct"

    return UpdateProxyConfig(
        mode=mode,
        socks_proxy=str(raw.get("socks_proxy", "") or "").strip(),
        ssh_host=str(raw.get("ssh_host", "") or "").strip(),
        ssh_user=str(raw.get("ssh_user", "") or "").strip(),
        ssh_port=_safe_int(raw.get("ssh_port", DEFAULT_SSH_PORT), DEFAULT_SSH_PORT, 1, 65535),
        local_socks_port=_safe_int(
            raw.get("local_socks_port", DEFAULT_LOCAL_SOCKS_PORT),
            DEFAULT_LOCAL_SOCKS_PORT,
            0,
            65535,
        ),
        identity_file=str(raw.get("identity_file", "") or "").strip(),
    )


def load_update_proxy_config():
    """Load saved update connection settings, falling back to direct access."""
    path = update_connection_config_path()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError, UnicodeError):
        return UpdateProxyConfig()
    return normalize_update_proxy_config(payload)


def save_update_proxy_config(config):
    """Persist update connection settings without storing an SSH password."""
    config = normalize_update_proxy_config(config)
    _validate_update_proxy_config(config)
    path = update_connection_config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(asdict(config), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)
    return path


def _parse_socks_proxy(proxy_url):
    """Return ``(host, port)`` for a SOCKS5 URL or host:port value."""
    value = str(proxy_url or "").strip()
    if not value:
        raise UpdateConnectionError("The SOCKS5 proxy address is empty.")
    if "://" not in value:
        value = "socks5h://" + value
    try:
        parsed = urlsplit(value)
        port = parsed.port or 1080
    except ValueError as exc:
        raise UpdateConnectionError("The SOCKS5 proxy address is invalid.") from exc
    if parsed.scheme.lower() not in ("socks5", "socks5h"):
        raise UpdateConnectionError("The update proxy must use socks5:// or socks5h://.")
    if parsed.username or parsed.password:
        raise UpdateConnectionError(
            "SOCKS5 username/password is not supported; use a local SSH tunnel instead."
        )
    if not parsed.hostname or not 1 <= int(port) <= 65535:
        raise UpdateConnectionError("The SOCKS5 proxy host or port is invalid.")
    return parsed.hostname, int(port)


def _socks_url(host, port):
    return "socks5h://{}:{}".format(host, int(port))


def _validate_update_proxy_config(config):
    config = normalize_update_proxy_config(config)
    if config.mode == "socks5":
        _parse_socks_proxy(config.socks_proxy)
    elif config.mode == "ssh-socks":
        if not config.ssh_host:
            raise UpdateConnectionError("The SSH relay host is empty.")
        if not 1 <= int(config.ssh_port) <= 65535:
            raise UpdateConnectionError("The SSH relay port is invalid.")
        if not 0 <= int(config.local_socks_port) <= 65535:
            raise UpdateConnectionError("The local SOCKS5 port is invalid.")
        if config.identity_file and not Path(config.identity_file).expanduser().is_file():
            raise UpdateConnectionError(
                "The configured SSH identity file does not exist: {}".format(
                    config.identity_file
                )
            )
    return config


def _find_free_local_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", 0))
        return int(probe.getsockname()[1])


class SshSocksTunnel:
    """Run an SSH dynamic forward for one update operation."""

    def __init__(self, config, timeout=DEFAULT_TIMEOUT):
        self.config = _validate_update_proxy_config(config)
        self.timeout = max(1.0, float(timeout))
        self.process = None
        self.local_port = None

    @property
    def proxy_url(self):
        if self.local_port is None:
            raise UpdateConnectionError("The SSH SOCKS tunnel has not started.")
        return _socks_url("127.0.0.1", self.local_port)

    def __enter__(self):
        return self.start()

    def __exit__(self, _exc_type, _exc_value, _traceback):
        self.close()

    def start(self):
        self.local_port = self.config.local_socks_port or _find_free_local_port()
        target = self.config.ssh_host
        if self.config.ssh_user:
            target = "{}@{}".format(self.config.ssh_user, self.config.ssh_host)

        command = [
            "ssh",
            "-N",
            "-D",
            "127.0.0.1:{}".format(self.local_port),
            "-p",
            str(self.config.ssh_port),
            "-o",
            "BatchMode=yes",
            "-o",
            "ExitOnForwardFailure=yes",
            "-o",
            "ServerAliveInterval=30",
            "-o",
            "ServerAliveCountMax=2",
            "-o",
            "ConnectTimeout={}".format(max(1, int(self.timeout))),
        ]
        if self.config.identity_file:
            command.extend(["-i", str(Path(self.config.identity_file).expanduser())])
        command.append(target)

        try:
            self.process = subprocess.Popen(
                command,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
            )
        except FileNotFoundError as exc:
            raise UpdateConnectionError("ssh was not found on PATH.") from exc
        except OSError as exc:
            raise UpdateConnectionError("Could not start ssh: {}".format(exc)) from exc

        deadline = time.monotonic() + self.timeout
        while time.monotonic() < deadline:
            if self.process.poll() is not None:
                detail = self._stderr_text()
                self.close()
                raise UpdateConnectionError(
                    "SSH SOCKS tunnel exited before it was ready.{}".format(
                        "\n" + detail if detail else ""
                    )
                )
            try:
                with socket.create_connection(("127.0.0.1", self.local_port), timeout=0.25):
                    return self
            except OSError:
                time.sleep(0.05)

        self.close()
        raise UpdateConnectionError(
            "Timed out waiting for the local SSH SOCKS5 port {}.".format(self.local_port)
        )

    def _stderr_text(self):
        if self.process is None or self.process.stderr is None:
            return ""
        try:
            return str(self.process.stderr.read() or "").strip()
        except (OSError, ValueError):
            return ""

    def close(self):
        process = self.process
        self.process = None
        if process is None:
            return
        try:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=2.0)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=2.0)
        except (OSError, subprocess.SubprocessError):
            pass
        finally:
            if process.stderr is not None:
                try:
                    process.stderr.close()
                except OSError:
                    pass


@contextmanager
def _update_proxy_session(proxy_config=None, timeout=DEFAULT_TIMEOUT):
    config = _validate_update_proxy_config(
        load_update_proxy_config() if proxy_config is None else proxy_config
    )
    if config.mode == "direct":
        yield None
    elif config.mode == "socks5":
        yield config.socks_proxy
    else:
        with SshSocksTunnel(config, timeout=timeout) as tunnel:
            yield tunnel.proxy_url


def _recv_exact(sock, size):
    chunks = []
    remaining = int(size)
    while remaining > 0:
        chunk = sock.recv(remaining)
        if not chunk:
            raise UpdateConnectionError("The SOCKS5 proxy closed the connection.")
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def _socks5_connect(proxy_host, proxy_port, target_host, target_port, timeout):
    """Connect to a target through a no-auth SOCKS5 proxy."""
    try:
        sock = socket.create_connection((proxy_host, int(proxy_port)), timeout=timeout)
        sock.sendall(b"\x05\x01\x00")
        greeting = _recv_exact(sock, 2)
        if greeting != b"\x05\x00":
            raise UpdateConnectionError(
                "The SOCKS5 proxy does not allow unauthenticated connections."
            )

        encoded_host = str(target_host).encode("idna")
        if len(encoded_host) > 255:
            raise UpdateConnectionError("The target hostname is too long for SOCKS5.")
        request = (
            b"\x05\x01\x00\x03"
            + bytes([len(encoded_host)])
            + encoded_host
            + int(target_port).to_bytes(2, "big")
        )
        sock.sendall(request)
        response = _recv_exact(sock, 4)
        if response[0] != 5 or response[1] != 0:
            raise UpdateConnectionError(
                "The SOCKS5 proxy could not connect to {}:{} (code {}).".format(
                    target_host, target_port, response[1]
                )
            )
        address_type = response[3]
        if address_type == 1:
            _recv_exact(sock, 4)
        elif address_type == 3:
            address_length = _recv_exact(sock, 1)[0]
            _recv_exact(sock, address_length)
        elif address_type == 4:
            _recv_exact(sock, 16)
        else:
            raise UpdateConnectionError("The SOCKS5 proxy returned an invalid address type.")
        _recv_exact(sock, 2)
        return sock
    except UpdateConnectionError:
        try:
            sock.close()
        except (UnboundLocalError, AttributeError, OSError):
            pass
        raise
    except (OSError, ValueError) as exc:
        try:
            sock.close()
        except (UnboundLocalError, AttributeError, OSError):
            pass
        raise UpdateConnectionError("Could not connect through the SOCKS5 proxy: {}".format(exc)) from exc


class _SocksHTTPSConnection(http.client.HTTPSConnection):
    """HTTPSConnection that dials the destination through SOCKS5."""

    def __init__(self, host, proxy_host, proxy_port, **kwargs):
        self._socks_proxy_host = proxy_host
        self._socks_proxy_port = proxy_port
        super().__init__(host, **kwargs)

    def connect(self):
        self.sock = _socks5_connect(
            self._socks_proxy_host,
            self._socks_proxy_port,
            self.host,
            self.port,
            self.timeout,
        )
        if self._tunnel_host:
            self._tunnel()
        self.sock = self._context.wrap_socket(self.sock, server_hostname=self.host)


class _SocksHTTPSHandler(HTTPSHandler):
    def __init__(self, proxy_host, proxy_port):
        super().__init__(context=ssl.create_default_context())
        self._socks_proxy_host = proxy_host
        self._socks_proxy_port = proxy_port

    def https_open(self, req):
        connection_factory = partial(
            _SocksHTTPSConnection,
            proxy_host=self._socks_proxy_host,
            proxy_port=self._socks_proxy_port,
        )
        return self.do_open(
            connection_factory,
            req,
            context=self._context,
            check_hostname=self._check_hostname,
        )


def _fetch_text(url, timeout=DEFAULT_TIMEOUT, allow_not_found=False, proxy_url=None):
    request = Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": "CGYRO-Comparison-Tool-Update-Checker",
        },
    )
    try:
        if proxy_url:
            proxy_host, proxy_port = _parse_socks_proxy(proxy_url)
            opener = build_opener(
                ProxyHandler({}),
                _SocksHTTPSHandler(proxy_host, proxy_port),
            )
            response_context = opener.open(request, timeout=timeout)
        else:
            response_context = urlopen(request, timeout=timeout)
        with response_context as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        if allow_not_found and exc.code == 404:
            return None
        raise UpdateCheckError(f"Update service returned HTTP {exc.code}.") from exc
    except UpdateConnectionError as exc:
        raise UpdateCheckError(str(exc)) from exc
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


def check_for_updates(
    current_version=APP_VERSION,
    timeout=DEFAULT_TIMEOUT,
    proxy_config=None,
):
    """Check GitHub for a newer version.

    ``UpdateCheckError`` is raised only after both supported endpoints fail.
    The function performs no UI work and can safely be called from a worker.
    """
    current_version = normalize_version(current_version)
    if current_version is None:
        raise UpdateCheckError("The local application version is invalid.")

    try:
        with _update_proxy_session(proxy_config, timeout=timeout) as proxy_url:
            errors = []

            try:
                release_payload = _fetch_text(
                    RELEASE_API_URL,
                    timeout=timeout,
                    allow_not_found=True,
                    proxy_url=proxy_url,
                )
                if release_payload:
                    release_result = _release_info(release_payload, current_version)
                    if release_result is not None:
                        return release_result
            except UpdateCheckError as exc:
                errors.append(exc)

            try:
                manifest_payload = _fetch_text(
                    _cache_busted_url(VERSION_URL),
                    timeout=timeout,
                    proxy_url=proxy_url,
                )
                return _manifest_info(manifest_payload, current_version)
            except UpdateCheckError as exc:
                errors.append(exc)

            if errors:
                raise UpdateCheckError(
                    "Unable to check for updates. Please verify your internet connection."
                ) from errors[-1]
            raise UpdateCheckError("No valid remote version was found.")
    except UpdateConnectionError as exc:
        raise UpdateCheckError(str(exc)) from exc


def _shell_join(parts):
    """Quote a command for the shell used by GIT_SSH_COMMAND."""
    if os.name == "nt":
        return subprocess.list2cmdline([str(part) for part in parts])
    return " ".join(shlex.quote(str(part)) for part in parts)


def _build_socks_proxy_command(proxy_host, proxy_port):
    helper = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--socks-proxy-command",
        str(proxy_host),
        str(proxy_port),
        "%h",
        "%p",
    ]
    return _shell_join(helper)


def _git_environment_for_proxy(proxy_url):
    """Return an environment that routes both Git HTTPS and SSH through SOCKS."""
    if not proxy_url:
        return None
    proxy_host, proxy_port = _parse_socks_proxy(proxy_url)
    environment = os.environ.copy()
    # Git/libcurl honors http.proxy and these standard proxy variables for
    # HTTPS remotes.  The explicit -c option below is used as an additional
    # compatibility path for older Git builds.
    for variable in ("ALL_PROXY", "all_proxy", "HTTP_PROXY", "http_proxy", "HTTPS_PROXY", "https_proxy"):
        environment[variable] = proxy_url

    proxy_command = _build_socks_proxy_command(proxy_host, proxy_port)
    ssh_command = [
        "ssh",
        "-o",
        "BatchMode=yes",
        "-o",
        "ProxyCommand={}".format(proxy_command),
    ]
    existing_ssh_command = environment.get("GIT_SSH_COMMAND", "").strip()
    if existing_ssh_command:
        ssh_command = shlex.split(existing_ssh_command) + ssh_command[1:]
    environment["GIT_SSH_COMMAND"] = _shell_join(ssh_command)
    return environment


def _pump_stdin_to_socket(sock):
    stream = getattr(sys.stdin, "buffer", sys.stdin)
    try:
        while True:
            chunk = stream.read(65536)
            if not chunk:
                break
            sock.sendall(chunk)
    except (OSError, ValueError):
        pass
    finally:
        try:
            sock.shutdown(socket.SHUT_WR)
        except OSError:
            pass


def _pump_socket_to_stdout(sock):
    stream = getattr(sys.stdout, "buffer", sys.stdout)
    try:
        while True:
            chunk = sock.recv(65536)
            if not chunk:
                break
            stream.write(chunk)
            stream.flush()
    except (OSError, ValueError):
        pass


def _run_socks_proxy_command(values):
    """Implement the SSH ProxyCommand using only the Python standard library."""
    proxy_host, proxy_port, target_host, target_port = values
    try:
        sock = _socks5_connect(
            proxy_host,
            int(proxy_port),
            target_host,
            int(target_port),
            timeout=None,
        )
    except (UpdateConnectionError, ValueError) as exc:
        print("SOCKS5 ProxyCommand failed: {}".format(exc), file=sys.stderr)
        return 1

    inbound = threading.Thread(target=_pump_stdin_to_socket, args=(sock,), daemon=True)
    outbound = threading.Thread(target=_pump_socket_to_stdout, args=(sock,), daemon=True)
    inbound.start()
    outbound.start()
    inbound.join()
    outbound.join()
    try:
        sock.close()
    except OSError:
        pass
    return 0


def _run_git(repo_dir, args, timeout, env=None):
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
            env=env,
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
    proxy_config=None,
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

    try:
        with _update_proxy_session(proxy_config, timeout=timeout) as proxy_url:
            pull_args = ["pull"]
            pull_environment = None
            if proxy_url:
                pull_args[0:0] = ["-c", "http.proxy={}".format(proxy_url)]
                pull_environment = _git_environment_for_proxy(proxy_url)
            pull_args.extend(["--ff-only", remote, branch])
            pull_result = _run_git(
                repo_root,
                pull_args,
                timeout,
                env=pull_environment,
            )
    except UpdateConnectionError as exc:
        raise GitUpdateError(str(exc)) from exc
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
    actions = parser.add_mutually_exclusive_group(required=False)
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
    connection = parser.add_argument_group("update connection")
    connection.add_argument(
        "--direct",
        action="store_true",
        help="ignore saved proxy settings and use the normal network connection",
    )
    connection.add_argument(
        "--socks5-proxy",
        default=None,
        metavar="URL",
        help="use an already running SOCKS5 proxy, e.g. socks5h://127.0.0.1:1080",
    )
    connection.add_argument(
        "--ssh-relay",
        default=None,
        metavar="HOST",
        help="start an SSH dynamic SOCKS5 tunnel through this relay host",
    )
    connection.add_argument(
        "--ssh-user",
        default="",
        metavar="USER",
        help="SSH user for --ssh-relay",
    )
    connection.add_argument(
        "--ssh-port",
        type=int,
        default=DEFAULT_SSH_PORT,
        metavar="PORT",
        help="SSH port for --ssh-relay (default: 22)",
    )
    connection.add_argument(
        "--ssh-key",
        default="",
        metavar="FILE",
        help="optional SSH private key; passwords are not stored or requested",
    )
    connection.add_argument(
        "--local-socks-port",
        type=int,
        default=DEFAULT_LOCAL_SOCKS_PORT,
        metavar="PORT",
        help="local SOCKS5 port for --ssh-relay (0 chooses a free port)",
    )
    connection.add_argument(
        "--socks-proxy-command",
        nargs=4,
        metavar=("PROXY_HOST", "PROXY_PORT", "TARGET_HOST", "TARGET_PORT"),
        help=argparse.SUPPRESS,
    )
    return parser


def _connection_config_from_cli(args):
    selected = int(bool(args.direct)) + int(args.socks5_proxy is not None) + int(args.ssh_relay is not None)
    if selected > 1:
        raise UpdateConnectionError(
            "Choose only one of --direct, --socks5-proxy, or --ssh-relay."
        )
    if args.direct:
        return UpdateProxyConfig()
    if args.socks5_proxy is not None:
        return _validate_update_proxy_config(
            UpdateProxyConfig(mode="socks5", socks_proxy=args.socks5_proxy)
        )
    if args.ssh_relay is not None:
        return _validate_update_proxy_config(
            UpdateProxyConfig(
                mode="ssh-socks",
                ssh_host=args.ssh_relay,
                ssh_user=args.ssh_user,
                ssh_port=args.ssh_port,
                local_socks_port=args.local_socks_port,
                identity_file=args.ssh_key,
            )
        )
    return load_update_proxy_config()


def _run_cli(argv=None):
    """Run the standalone update/check command and return an exit code."""
    parser = _build_cli_parser()
    args = parser.parse_args(argv)

    if args.socks_proxy_command:
        return _run_socks_proxy_command(args.socks_proxy_command)
    if not args.check_update and not args.update:
        parser.error("one of --check-update or --update is required")

    try:
        proxy_config = _connection_config_from_cli(args)
    except UpdateConnectionError as exc:
        print(f"Update connection configuration failed: {exc}", file=sys.stderr)
        return 1

    if args.check_update:
        try:
            result = check_for_updates(
                APP_VERSION,
                timeout=args.timeout,
                proxy_config=proxy_config,
            )
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
            proxy_config=proxy_config,
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
