"""Streaming FULLT extrema tracing, independent of Tk and pygacode.

Storage: (partner kx, signed partner ky, channel, diagnosed kx).
Plots use a positive diagnosed-ky convention; channel zero is native Re(T).
"""

from collections import OrderedDict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional
import importlib.util
import sys
import sysconfig

import numpy as np

try:
    from cgyro_fullt_reader import fullt_reader_for_case
except ImportError:
    from .cgyro_fullt_reader import fullt_reader_for_case


class TraceError(ValueError):
    """An invalid request or a diagnostic that cannot support the requested trace."""


class TraceCancelled(TraceError):
    pass


@dataclass(frozen=True)
class Mode:
    kx: float
    ky: float


@dataclass(frozen=True)
class TraceSettings:
    kx: float = 0.0
    ky: float = 0.2
    steps: int = 4
    direction: str = "positive"
    diagnostic: str = "asym"
    t_start: Optional[float] = None
    t_end: Optional[float] = None
    cycles: str = "stop"

    def validate(self):
        if self.direction not in ("positive", "negative"):
            raise TraceError("Direction must be positive or negative.")
        if self.diagnostic not in ("fullt", "asym"):
            raise TraceError("Diagnostic must be fullt or asym.")
        if self.cycles not in ("stop", "allow", "skip"):
            raise TraceError("Cycle policy must be stop, allow or skip.")
        if not isinstance(self.steps, int) or not 1 <= self.steps <= 20:
            raise TraceError("Choose between 1 and 20 panels/steps.")
        for name in ("kx", "ky", "t_start", "t_end"):
            value = getattr(self, name)
            if value is not None and not np.isfinite(value):
                raise TraceError(name + " must be finite.")


@dataclass
class Peak:
    kx: float
    ky: float
    value: float
    next_mode: Mode
    conjugated: bool
    ties: int = 1


@dataclass
class TracePanel:
    mode: Mode
    values: np.ndarray
    peak: Optional[Peak]
    excluded_visited: int = 0


@dataclass
class TraceResult:
    case: str
    source_file: str
    settings: TraceSettings
    kx: np.ndarray
    ky: np.ndarray
    times: np.ndarray
    time_indices: np.ndarray
    layout: dict
    source_signature: str
    panels: list = field(default_factory=list)
    stop_reason: str = ""
    stop_code: str = "completed"

    def metadata(self):
        return {
            "case": self.case, "source_file": self.source_file,
            "settings": asdict(self.settings), "layout": self.layout,
            "source_signature": self.source_signature,
            "time_indices": self.time_indices.tolist(), "times": self.times.tolist(),
            "averaging": "arithmetic mean of all selected saved Re(T) samples",
            "ky_convention": "positive diagnosed ky; signed partner branches retained",
            "negative_partner": "next (-kx_prime, -ky_prime) when ky_prime < 0",
            "stop_code": self.stop_code, "stop_reason": self.stop_reason,
            "panels": [{
                "step": i + 1, "mode": asdict(panel.mode),
                "peak": asdict(panel.peak) if panel.peak else None,
                "nonfinite_points": int(np.count_nonzero(~np.isfinite(panel.values))),
                "excluded_visited": panel.excluded_visited,
            } for i, panel in enumerate(self.panels)],
        }


def _enable_private_decoder():
    """Honor this interpreter's existing optional decoder without loading pygacode."""
    if importlib.util.find_spec("zstandard") is not None:
        return
    abi = sys.implementation.cache_tag + "-" + sysconfig.get_platform()
    location = Path(__file__).resolve().parent / ".runtime" / "python" / abi
    if (location / "zstandard" / "__init__.py").is_file():
        sys.path.insert(0, str(location))


def _cancelled(cancel):
    if cancel is not None and cancel.is_set():
        raise TraceCancelled("Trace cancelled; the previous result is retained.")


def canonical_partner(kx, ky):
    """Choose the stored half-plane representative, reversing BOTH coordinates."""
    return Mode(-float(kx), -float(ky)) if ky < 0 else Mode(float(kx), abs(float(ky)))


class FullTSource:
    """Snapshot one native reader and cache bounded time means, never full histories."""

    def __init__(self, case, settings, progress=None, cancel=None):
        settings.validate()
        _enable_private_decoder()
        suffix = ".cgyro.fullt_asym" if settings.diagnostic == "asym" else ".cgyro.fullt"
        self.reader = fullt_reader_for_case(case, suffix)
        if self.reader is None:
            raise TraceError("No native " + suffix + " diagnostic found in the selected case.")
        try:
            self.case = Path(case).resolve()
            self.settings, self.progress, self.cancel = settings, progress, cancel
            m = self.reader.layout
            if m.dx <= 0 or (m.ny > 1 and m.dy == 0):
                raise TraceError("FULLT physical grid spacing is missing or invalid.")
            self.kx = (np.arange(m.nx) - m.nx // 2) * m.dx
            self.ky = np.arange(-(m.ny - 1), m.ny) * abs(m.dy)
            self.source_x = np.array([0.]) if m.ns == 1 else self.kx.copy()
            self.source_y = np.arange(m.ny) * abs(m.dy)
            times = np.asarray(self.reader.times, dtype=float)
            if times.size == 0 or not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
                raise TraceError("Time labels must be nonempty, finite and strictly increasing.")
            end = float(times[-1]) if settings.t_end is None else float(settings.t_end)
            start = float(times[0] + .5 * (end - times[0])) if settings.t_start is None else float(settings.t_start)
            if start > end:
                raise TraceError("Time start must not exceed time end.")
            self.tids = np.flatnonzero((times >= start) & (times <= end))
            if not self.tids.size:
                raise TraceError("No saved samples fall within the requested time window.")
            self.times = times[self.tids]
            self.cache = OrderedDict()
            self.cache_bytes = 0
            self.cache_limit = 128 * 1024 ** 2
            self.initial_stat = self.reader.path.stat()
        except Exception:
            self.reader.close()
            raise

    def close(self):
        self.cache.clear()
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def indices(self, mode, *, exact=False):
        """Only the initial user request may snap to a nearby stored mode."""
        m = self.reader.layout
        if m.ns == 1 and abs(mode.kx) > 1e-10 * max(m.dx, 1.):
            raise TraceError("Only kx=0 was saved (FULL_T_KX0=1); this next kx has no diagnostic.")
        ix = int(np.argmin(np.abs(self.source_x - mode.kx)))
        iy = int(np.argmin(np.abs(self.source_y - mode.ky)))
        if exact and (not np.isclose(self.source_x[ix], mode.kx, rtol=1e-8, atol=1e-12)
                      or not np.isclose(self.source_y[iy], mode.ky, rtol=1e-8, atol=1e-12)):
            raise TraceError("The next mode is outside the saved diagnosed grid (including the Nyquist edge).")
        return ix, iy

    def mode(self, indices):
        ix, iy = indices
        return Mode(float(self.source_x[ix]), float(self.source_y[iy]))

    def mean_map(self, indices):
        ix, iy = indices
        m = self.reader.layout
        # When affordable, average all diagnosed kx at this ky in one read pass.
        # This is especially useful for radial chains through compressed blocks.
        group = m.nx * (2 * m.ny - 1) * m.ns * 8 <= self.cache_limit // 2
        key = (iy, None if group else ix)
        if key not in self.cache:
            shape = (m.nx, 2 * m.ny - 1, m.ns) if group else (m.nx, 2 * m.ny - 1)
            total = np.zeros(shape, dtype=np.float64)
            for count, tid in enumerate(self.tids, 1):
                _cancelled(self.cancel)
                block = self.reader.block(int(tid), iy)
                values = block[:, :, 0, :] if group else block[:, :, 0, ix]
                # A nonfinite sample leaves that grid point invalid, rather than
                # silently averaging a different time set at different locations.
                total += values
                del values, block
                if self.progress and (count % 20 == 0 or count == self.tids.size):
                    self.progress(f"ky={self.source_y[iy]:.5g}: averaging {count}/{self.tids.size} samples")
            total /= self.tids.size
            while self.cache and self.cache_bytes + total.nbytes > self.cache_limit:
                _, removed = self.cache.popitem(last=False)
                self.cache_bytes -= removed.nbytes
            self.cache[key] = total
            self.cache_bytes += total.nbytes
        self.cache.move_to_end(key)
        averaged = self.cache[key]
        result = averaged[:, :, ix] if group else averaged
        return np.array(result.T, copy=True)

    def check_snapshot(self):
        # FTZ can safely append beyond its committed snapshot. Native raw files
        # cannot change underneath a time mean without invalidating provenance.
        if self.reader.path.suffix != ".ftz":
            current = self.reader.path.stat()
            if (current.st_size, current.st_mtime_ns) != (self.initial_stat.st_size, self.initial_stat.st_mtime_ns):
                raise TraceError("The raw FULLT file changed during tracing; rerun after writing finishes.")


def select_peak(values, kx, ky, direction, excluded=()):
    """Pick a strict positive maximum or negative minimum from actual samples."""
    valid = np.isfinite(values) & (values > 0 if direction == "positive" else values < 0)
    excluded_count = 0
    if excluded:
        before = np.count_nonzero(valid)
        candidate_x = np.where(ky[:, None] < 0, -kx[None, :], kx[None, :])
        candidate_y = np.abs(ky[:, None])
        for old in excluded:
            visited = (np.isclose(candidate_x, old.kx, rtol=1e-8, atol=1e-12) &
                       np.isclose(candidate_y, old.ky, rtol=1e-8, atol=1e-12))
            valid &= ~visited
        excluded_count = int(before - np.count_nonzero(valid))
    if not np.any(valid):
        return None, excluded_count
    extreme = float(np.max(values[valid]) if direction == "positive" else np.min(values[valid]))
    matches = np.argwhere(valid & (values == extreme))
    # Exact ties prefer a directly stored positive partner, then small |ky|,
    # small |kx|, negative kx. The tie count is included in the exported trace.
    iy, ix = min(matches, key=lambda p: (ky[p[0]] < 0, abs(ky[p[0]]), abs(kx[p[1]]), kx[p[1]]))
    return Peak(float(kx[ix]), float(ky[iy]), extreme,
                canonical_partner(kx[ix], ky[iy]), bool(ky[iy] < 0), len(matches)), excluded_count


def trace_case(case, settings, *, progress=None, cancel=None):
    """Trace extrema of time-averaged Re(T); each returned panel is one step."""
    with FullTSource(case, settings, progress, cancel) as source:
        indices = source.indices(Mode(settings.kx, abs(settings.ky)))
        result = TraceResult(source.case.name, source.reader.path.name, settings,
                             source.kx, source.ky, source.times, source.tids,
                             asdict(source.reader.layout), repr(source.reader.signature))
        visited = set()
        modes = []
        for step in range(settings.steps):
            _cancelled(cancel)
            mode = source.mode(indices)
            if progress:
                progress(f"Step {step + 1}/{settings.steps}: kx={mode.kx:.5g}, ky={mode.ky:.5g}")
            visited.add(indices)
            modes.append(mode)
            values = source.mean_map(indices)
            peak, excluded = select_peak(values, source.kx, source.ky, settings.direction,
                                         modes if settings.cycles == "skip" else ())
            result.panels.append(TracePanel(mode, values, peak, excluded))
            if peak is None:
                result.stop_code = "no_extremum"
                result.stop_reason = "No eligible finite " + settings.direction + " T remains."
                break
            if step + 1 == settings.steps:
                result.stop_reason = f"Completed {settings.steps} requested panels."
                break
            try:
                following = source.indices(peak.next_mode, exact=True)
            except TraceError as exc:
                result.stop_code, result.stop_reason = "unavailable_mode", str(exc)
                break
            if settings.cycles == "stop" and following in visited:
                result.stop_code = "cycle"
                result.stop_reason = "The selected extremum revisits an earlier mode. Use Allow repeats or Skip visited to continue."
                break
            indices = following
        source.check_snapshot()
        return result
