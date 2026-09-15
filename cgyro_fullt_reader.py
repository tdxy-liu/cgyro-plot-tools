# SPDX-License-Identifier: Apache-2.0
# Distributed as the Comparison tool's pinned FTZ v1 reader (2026-09-15).
# Renamed from the customized GACODE reader; see third_party/gacode/README.md.
"""FTZ v1: lossless FULL_T_ASYM blocks, lazy slices and committed-record indexes.

Native block shape is (kx_prime, signed_ky_prime, channel, diagnosed_kx), Fortran
order. A block is independent for each (output counter, diagnosed_ky).
No physics, sign or normalization transformations are performed here.
"""
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import binascii
import math
import os
import struct
import uuid
import warnings
import numpy as np

HS, BS, FS = 128, 64, 40


class FTZError(ValueError):
    """Missing, incompatible or corrupt FTZ data; never use a synthetic fallback."""


def _crc(b):
    return binascii.crc32(b) & 0xffffffff


def _zstd():
    try:
        import zstandard
    except ImportError as exc:
        raise FTZError("Compressed FULL_T_ASYM requires zstandard in the Python "
                       "environment running this tool: python -m pip install -r "
                       "requirements-compression.txt. See FULLT_COMPRESSED_READING.md "
                       "for offline installation; no full-file extraction is needed.") from exc
    return zstandard


@dataclass(frozen=True)
class Layout:
    nx: int
    ny: int
    nc: int
    ns: int
    word: int
    dx: float = 0.0
    dy: float = 0.0

    def __post_init__(self):
        if (self.nx < 2 or self.nx % 2 or self.ny < 1 or self.nc not in (1, 2)
                or self.ns not in (1, self.nx) or self.word not in (4, 8)
                or not math.isfinite(self.dx) or not math.isfinite(self.dy)):
            raise FTZError("Unsupported FTZ grid/precision/layout")
        if self.nbytes >= 2**63 - BS:
            raise FTZError("FTZ block size overflow")

    @property
    def shape(self):
        return (self.nx, 2*self.ny-1, self.nc, self.ns)

    @property
    def nbytes(self):
        return math.prod(self.shape)*self.word

    @property
    def dtype(self):
        return np.dtype("<f%d" % self.word)

    @property
    def uint(self):
        return np.dtype("<u%d" % self.word)

    def header(self, identity=None, index=False):
        b = bytearray(HS)
        b[:8] = b"CGFTI01\0" if index else b"CGFTZ01\0"
        struct.pack_into("<10I", b, 8, 1, HS, self.nx, self.ny, self.nc,
                         self.ns, self.word, self.nx//2, 1, 3)
        b[48:64] = identity if identity is not None else uuid.uuid4().bytes
        struct.pack_into("<dd", b, 64, self.dx, self.dy)
        struct.pack_into("<I", b, 124, _crc(b[:124]))
        return bytes(b)

    @classmethod
    def from_header(cls, b, index=False):
        magic = b"CGFTI01\0" if index else b"CGFTZ01\0"
        if len(b) != HS or b[:8] != magic or _crc(b[:124]) != struct.unpack_from("<I", b, 124)[0]:
            raise FTZError("Invalid FTZ header magic/checksum")
        ver, hs, nx, ny, nc, ns, word, center, endian, level = struct.unpack_from("<10I", b, 8)
        if (ver, hs, center, endian, level) != (1, HS, nx//2, 1, 3):
            raise FTZError("Unsupported FTZ version/layout")
        return cls(nx, ny, nc, ns, word, *struct.unpack_from("<dd", b, 64))


@lru_cache(maxsize=4)
def _pairs(nx, ny, nc, ns, ky):
    """Conservative complementary-index mapping; no periodic boundary aliasing."""
    nq, center = 2*ny-1, nx//2
    i, j = np.meshgrid(np.arange(nx), np.arange(nq), indexing="ij")
    i, j = i.ravel(order="F"), j.ravel(order="F")
    left, right = [], []
    for k in range(ns):
        kp = 0 if ns == 1 else k-center
        bi, bj = kp-(i-center)+center, ky+2*(ny-1)-j
        valid = (bi >= 0) & (bi < nx) & (bj >= 0) & (bj < nq)
        a, b = i+nx*j, bi+nx*bj
        valid &= a < b
        for c in range(nc):
            base = nx*nq*(c+nc*k)
            left.append(a[valid]+base)
            right.append(b[valid]+base)
    return np.concatenate(left), np.concatenate(right)


def _predict(words, layout, ky):
    a, b = _pairs(layout.nx, layout.ny, layout.nc, layout.ns, ky)
    words[b] ^= words[a]


def block_info(header, layout, ky=None, step=None):
    if len(header) != BS or header[:8] != b"FTZBLK1\0":
        raise FTZError("Invalid FTZ block magic")
    s, q, mode, raw_n, payload_n, raw_crc, payload_crc, tv, hc = struct.unpack_from("<QIIQQIIdI", header, 8)
    if (hc != _crc(header[:56]) or mode not in (0, 1) or raw_n != layout.nbytes
            or payload_n > raw_n or not 0 <= q < layout.ny or s < 1
            or not math.isfinite(tv)
            or (ky is not None and q != ky) or (step is not None and s != step)):
        raise FTZError("Invalid/corrupt FTZ block metadata")
    return s, q, mode, payload_n, raw_crc, payload_crc, tv


def encode_block(array, layout, ky, step, time):
    if not 0 <= ky < layout.ny or step < 1 or not math.isfinite(time):
        raise FTZError("Invalid block identity")
    raw = array if isinstance(array, bytes) else np.asarray(array, dtype=layout.dtype).tobytes(order="F")
    if len(raw) != layout.nbytes:
        raise FTZError("Wrong FTZ raw block length")
    words = np.frombuffer(raw, dtype=layout.uint).copy()
    _predict(words, layout, ky)
    mask = words != 0
    packed = np.packbits(mask, bitorder="little").tobytes() + words[mask].tobytes()
    payload = _zstd().ZstdCompressor(level=3).compress(packed)
    mode = 1
    if len(payload) >= len(raw):
        mode, payload = 0, raw
    b = bytearray(BS)
    b[:8] = b"FTZBLK1\0"
    struct.pack_into("<QIIQQIId", b, 8, step, ky, mode, len(raw), len(payload), _crc(raw), _crc(payload), time)
    struct.pack_into("<I", b, 56, _crc(b[:56]))
    return bytes(b)+payload


def decode_block(block, layout, ky=None, step=None):
    s, q, mode, pn, raw_crc, payload_crc, tv = block_info(block[:BS], layout, ky, step)
    payload = block[BS:]
    if len(payload) != pn or _crc(payload) != payload_crc:
        raise FTZError("FTZ compressed payload checksum/length mismatch")
    if mode == 0:
        if pn != layout.nbytes:
            raise FTZError("FTZ raw block length mismatch")
        raw = payload
    else:
        nw = layout.nbytes//layout.word
        mn = (nw+7)//8
        # Keep missing-decoder errors actionable; they are not corrupt data.
        zstd = _zstd()
        try:
            frame_size = zstd.frame_content_size(payload)
            if frame_size > layout.nbytes+mn:
                raise FTZError("FTZ frame exceeds the expected block allocation")
            packed = zstd.ZstdDecompressor().decompress(
                payload, max_output_size=layout.nbytes+mn, allow_extra_data=False)
        except Exception as exc:
            raise FTZError("Invalid FTZ Zstandard frame") from exc
        if len(packed) < mn or len(packed) > layout.nbytes+mn:
            raise FTZError("Invalid FTZ packed size")
        mask = np.unpackbits(np.frombuffer(packed[:mn], dtype=np.uint8), bitorder="little")[:nw].astype(bool)
        if len(packed) != mn+int(mask.sum())*layout.word:
            raise FTZError("Invalid FTZ mask/word count")
        words = np.zeros(nw, dtype=layout.uint)
        words[mask] = np.frombuffer(packed[mn:], dtype=layout.uint)
        _predict(words, layout, q)
        raw = words.tobytes()
    if _crc(raw) != raw_crc:
        raise FTZError("Restored FULL_T_ASYM checksum mismatch")
    return np.frombuffer(raw, dtype=layout.dtype).reshape(layout.shape, order="F")


def _footer(step, ny, start):
    b = bytearray(FS)
    b[:8] = b"FTZEND1\0"
    struct.pack_into("<QI", b, 8, step, ny)
    struct.pack_into("<Q", b, 24, start)
    struct.pack_into("<I", b, 32, _crc(b[:32]))
    return bytes(b)


def _index_record(step, time, start, end, entries):
    b = bytearray(48+24*len(entries))
    b[:8] = b"FTZIDX1\0"
    struct.pack_into("<QdQQ", b, 8, step, time, start, end)
    for q, entry in enumerate(entries):
        struct.pack_into("<QQII", b, 40+q*24, *entry)
    struct.pack_into("<I", b, len(b)-8, _crc(b[:-8]))
    b[-4:] = b"CMIT"
    return bytes(b)


def _sync(f):
    f.flush()
    os.fsync(f.fileno())


class FTZReader:
    """Read only complete index commits. Payloads are decoded on demand."""

    def __init__(self, path):
        self.path = Path(path)
        self.index_path = self.path.with_suffix(".fti")
        self.refresh()

    def refresh(self):
        with self.path.open("rb") as df, self.index_path.open("rb") as ix:
            h, hi = df.read(HS), ix.read(HS)
            layout = Layout.from_header(h)
            Layout.from_header(hi, index=True)
            if h[8:124] != hi[8:124]:
                raise FTZError("FTZ data/index identity mismatch")
            rn = 48+24*layout.ny
            # Snapshot the index size; a concurrently appended record is read next refresh.
            index_size = os.fstat(ix.fileno()).st_size
            data_size = os.fstat(df.fileno()).st_size
            records, cursor = [], HS
            while ix.tell()+rn <= index_size:
                b = ix.read(rn)
                # Short trailing records are ignored by the loop bound. A full
                # record with a damaged marker is ambiguous: fail conservatively
                # rather than silently discarding possibly committed history.
                if (b[:8] != b"FTZIDX1\0" or b[-4:] != b"CMIT"
                        or _crc(b[:-8]) != struct.unpack_from("<I", b, rn-8)[0]):
                    raise FTZError("Corrupt committed FTZ index")
                step, tv, start, end = struct.unpack_from("<QdQQ", b, 8)
                if step != len(records)+1 or start != cursor or not math.isfinite(tv):
                    raise FTZError("FTZ index counter/time/offset mismatch")
                entries = []
                for q in range(layout.ny):
                    e = struct.unpack_from("<QQII", b, 40+24*q)
                    if e[0] != cursor or not BS <= e[1] <= BS+layout.nbytes:
                        raise FTZError("FTZ block offsets/lengths inconsistent")
                    cursor += e[1]
                    entries.append(e)
                if end != cursor+FS or end > data_size:
                    raise FTZError("Committed FTZ data is truncated")
                df.seek(cursor)
                if df.read(FS) != _footer(step, layout.ny, start):
                    raise FTZError("FTZ data commit footer mismatch")
                cursor = end
                records.append((step, tv, start, end, tuple(entries)))
            self.layout, self.identity, self.records = layout, h[48:64].hex(), records
            self.committed_end = cursor
            self.times = np.asarray([r[1] for r in records])
            self.signature = (str(self.path.resolve()), self.identity, len(records), cursor,
                              index_size, os.fstat(ix.fileno()).st_mtime_ns,
                              data_size, os.fstat(df.fileno()).st_mtime_ns)
        return self

    @property
    def nt(self):
        return len(self.records)

    def block(self, time_index, ky):
        if not 0 <= time_index < self.nt or not 0 <= ky < self.layout.ny:
            raise IndexError("FTZ time or diagnosed ky outside committed data")
        step, tv, _, _, entries = self.records[time_index]
        offset, length, rc, pc = entries[ky]
        with self.path.open("rb") as f:
            f.seek(offset)
            block = f.read(length)
        info = block_info(block[:BS], self.layout, ky, step)
        if (info[4], info[5], info[6]) != (rc, pc, tv):
            raise FTZError("FTZ block/index mismatch")
        return decode_block(block, self.layout, ky, step)

    def read(self, ky, kx=None, time_indices=None):
        """Return (kx_prime, ky_prime, channel, diagnosed_kx, time)."""
        tids = list(range(self.nt)) if time_indices is None else list(np.asarray(time_indices, dtype=int).ravel())
        ks = list(range(self.layout.ns)) if kx is None else [int(kx)]
        if any(not 0 <= k < self.layout.ns for k in ks):
            raise IndexError("FTZ diagnosed kx outside saved range")
        out = np.empty(self.layout.shape[:3]+(len(ks), len(tids)), dtype=self.layout.dtype)
        for j, t in enumerate(tids):
            out[..., j] = self.block(t, int(ky))[..., ks]
        return out

    def legacy(self, ky=None, kx=None, time_indices=None):
        """Explicit float64 legacy layout: [ri,K_y,K_x,kx_prime,ky_prime,channel,time]."""
        qs = range(self.layout.ny) if ky is None else [int(ky)]
        tids = list(range(self.nt)) if time_indices is None else list(time_indices)
        ns = self.layout.ns if kx is None else 1
        out = np.zeros((2, len(qs), ns, self.layout.nx, 2*self.layout.ny-1,
                        self.layout.nc, len(tids)), dtype=np.float64)
        for j, q in enumerate(qs):
            out[0, j] = self.read(q, kx, tids).transpose(3, 0, 1, 2, 4)
        return out

    def verify(self):
        for t in range(self.nt):
            for q in range(self.layout.ny):
                self.block(t, q)
        return {"records": self.nt, "blocks": self.nt*self.layout.ny,
                "identity": self.identity, "committed_bytes": self.committed_end}

    def close(self):
        # Payload files are opened only inside block()/refresh().
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class FTZWriter:
    """Serial streaming converter/test writer; CGYRO uses the MPI writer."""

    def __init__(self, path, layout, identity=None):
        self.path, self.layout = Path(path), layout
        self.index_path = self.path.with_suffix(".fti")
        self.header = layout.header(identity)
        self.df = self.path.open("xb")
        try:
            self.ix = self.index_path.open("xb")
        except Exception:
            self.df.close()
            raise
        self.df.write(self.header)
        self.ix.write(layout.header(self.header[48:64], index=True))
        _sync(self.df)
        _sync(self.ix)
        self.step = 0
        self.failed = False

    def append(self, blocks, time, verify=True):
        """Consume exactly ny arrays, keeping only the current block in memory."""
        if self.failed:
            raise FTZError("Writer failed previously; inspect the partial files")
        start, step, entries = self.df.tell(), self.step+1, []
        try:
            count = 0
            for q, array in enumerate(blocks):
                if q >= self.layout.ny:
                    raise FTZError("Too many diagnosed ky blocks")
                raw = np.asarray(array, dtype=self.layout.dtype).tobytes(order="F")
                block = encode_block(raw, self.layout, q, step, time)
                if verify and decode_block(block, self.layout, q, step).tobytes(order="F") != raw:
                    raise FTZError("Conversion is not bit-exact")
                info = block_info(block[:BS], self.layout, q, step)
                entries.append((self.df.tell(), len(block), info[4], info[5]))
                self.df.write(block)
                count += 1
            if count != self.layout.ny:
                raise FTZError("Missing diagnosed ky blocks")
            self.df.write(_footer(step, self.layout.ny, start))
            _sync(self.df)
            self.ix.write(_index_record(step, time, start, self.df.tell(), entries))
            _sync(self.ix)
            self.step = step
        except Exception:
            self.failed = True
            raise

    def close(self):
        self.df.close()
        self.ix.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def input_values(case):
    """Generated run inputs, then explicit current input overrides."""
    result = {}
    case = Path(case)
    for name in ("input.cgyro.gen", "input.cgyro"):
        path = case/name
        if not path.exists():
            continue
        for line in path.read_text(errors="replace").splitlines():
            line = line.split("#", 1)[0].split("!", 1)[0].strip()
            if not line:
                continue
            if "=" in line:
                key, value = line.split("=", 1)
            else:
                fields = line.split()
                if len(fields) != 2:
                    continue
                value, key = fields
            try:
                result[key.strip()] = float(value.strip().replace("D", "e").replace("d", "e"))
            except ValueError:
                continue
    return result


def selected_path(case):
    """Select raw vs compressed explicitly; never silently mask missing history."""
    case = Path(case)
    requested = input_values(case).get("FULL_T_ASYM_COMPRESSION", 0)
    if requested not in (0, 1):
        raise FTZError("FULL_T_ASYM_COMPRESSION must be 0 or 1")
    flag = int(requested)
    bases = [case/"bin.cgyro.fullt_asym", case/"bin"/"bin.cgyro.fullt_asym"]
    raw = [p for p in bases if p.is_file()]
    packed = [Path(str(p)+".ftz") for p in bases if Path(str(p)+".ftz").is_file()]
    candidates = packed if flag else raw
    if not candidates and (packed or raw or flag):
        raise FTZError("FULL_T_ASYM format conflicts with FULL_T_ASYM_COMPRESSION; "
                       "select the correct input flag or convert existing history first")
    if len(candidates) > 1:
        raise FTZError("Ambiguous FULL_T_ASYM files in case and case/bin")
    return candidates[0] if candidates else None


def case_grid(case):
    """Read the native grid, never infer physical axes from a plot label."""
    g = np.fromfile(Path(case)/"out.cgyro.grids", sep=" ")
    if g.size < 11 or not np.all(np.isfinite(g)):
        raise FTZError("Missing or invalid case grid metadata")
    ny, nx, ntheta, ne, nxi, box = (int(g[k]) for k in (0, 3, 4, 5, 6, 7))
    if min(ny, nx, ntheta, ne, nxi, box) < 1 or g[8] <= 0:
        raise FTZError("Invalid case grid dimensions/length")
    # CGYRO warns, but permits this grid. Native thetab is allocated as
    # (n_theta, n_radial/box_size), with integer division before multiplication.
    # Keep that layout; do not round up or change the simulation's radial grid.
    if nx % box:
        warnings.warn(
            f"CGYRO resolution warning: N_RADIAL={nx} is not a multiple of BOX_SIZE={box}; "
            "reading the native grid with integer-division thetab length.",
            RuntimeWarning, stacklevel=2)
    mark = 11+nx+ntheta+ne+nxi+ntheta*(nx//box)
    ky = g[mark:mark+ny]
    if ky.size != ny or not np.array_equal(g[11:11+nx], np.arange(nx)-nx//2):
        raise FTZError("Unsupported case radial grid layout")
    dy = float(ky[1]) if ny > 1 else float(ky[0])
    if not np.allclose(ky, np.arange(ny)*dy, rtol=2e-6, atol=1e-12):
        raise FTZError("Unsupported case toroidal grid layout")
    return nx, ny, float(2*np.pi/g[8]), dy


def case_times(case):
    """Read complete time-label lines; a live, incomplete final line is ignored."""
    lines = (Path(case)/"out.cgyro.time").read_text().splitlines(keepends=True)
    values = []
    for i, line in enumerate(lines):
        if not line.strip():
            continue
        if i == len(lines)-1 and not line.endswith(("\n", "\r")):
            # Files written by CGYRO terminate each record with a newline.
            break
        try:
            value = float(line.split()[0].replace("D", "e").replace("d", "e"))
        except (ValueError, IndexError) as exc:
            raise FTZError("Invalid case time labels") from exc
        if not math.isfinite(value):
            raise FTZError("Non-finite case time label")
        values.append(value)
    return np.asarray(values)


def _check_time_prefix(actual, expected):
    actual, expected = np.asarray(actual).ravel(), np.asarray(expected).ravel()
    n = min(len(actual), len(expected))
    # out.cgyro.time uses rounded printed labels; FTZ stores the native double.
    # Unequal lengths are normal while a run is appending or after a rewind.
    if not np.allclose(actual[:n], expected[:n], rtol=5e-5, atol=1e-10):
        raise FTZError("FULLT time prefix differs from the case; reload the correct case")


def validate_loaded_case(reader, data):
    """Reject stale GUI axes/time labels before selecting a native slice."""
    m = reader.layout
    if (int(data.n_radial), int(data.n_n)) != (m.nx, m.ny):
        raise FTZError("FULLT dimensions differ from the loaded case")
    if hasattr(data, "ky"):
        ky = np.asarray(data.ky).ravel()
        if ky.size != m.ny or not np.allclose(ky, np.arange(m.ny)*m.dy, rtol=2e-6, atol=1e-12):
            raise FTZError("FULLT ky axis differs from the loaded case")
    if hasattr(data, "length") and (not np.isfinite(data.length) or data.length <= 0 or
            not np.isclose(2*np.pi/data.length, m.dx, rtol=2e-6, atol=1e-12)):
        raise FTZError("FULLT kx spacing differs from the loaded case")
    if hasattr(data, "t"):
        _check_time_prefix(reader.times, data.t)


def reader_for_case(case):
    path = selected_path(case)
    if path is None or path.suffix != ".ftz":
        return None
    reader = FTZReader(path)
    hints = input_values(case)
    for name, actual, expected in (
        ("FULL_T_REAL_ONLY", reader.layout.nc, 1 if hints.get("FULL_T_REAL_ONLY")==1 else 2),
        ("FULL_T_KX0", reader.layout.ns, 1 if hints.get("FULL_T_KX0")==1 else reader.layout.nx),
        ("HIPREC_FLAG", reader.layout.word, 8 if hints.get("HIPREC_FLAG")==1 else 4),
    ):
        if name in hints and actual != expected:
            raise FTZError("FTZ metadata conflicts with case input "+name)
    grids = Path(case)/"out.cgyro.grids"
    if grids.exists():
        nx, ny, dx, dy = case_grid(case)
        if ((nx, ny) != (reader.layout.nx, reader.layout.ny) or
                not np.allclose((dx, dy), (reader.layout.dx, reader.layout.dy), rtol=2e-6, atol=1e-12)):
            raise FTZError("FTZ file does not match the case grid")
    if (Path(case)/"out.cgyro.time").exists():
        _check_time_prefix(reader.times, case_times(case))
    return reader


def fullt_reader_for_case(case, suffix=".cgyro.fullt_asym"):
    """Common native reader. Raw MPI rank/local-Ky order is canonicalized here.

    Absence returns None. Present but incompatible data must raise, never fall
    back to an inferred layout or a synthetic diagnostic.
    """
    if suffix not in (".cgyro.fullt", ".cgyro.fullt_asym"):
        raise FTZError("Unsupported FULLT diagnostic suffix")
    if suffix == ".cgyro.fullt_asym":
        path = selected_path(case)
        if path is None:
            return None
        if path.suffix == ".ftz":
            return reader_for_case(case)
    else:
        paths = [Path(case)/("bin"+suffix), Path(case)/"bin"/("bin"+suffix)]
        if not any(p.is_file() for p in paths):
            return None
    return RawFullT(case, suffix=suffix)


class FTZArrayView:
    """Basic indexing adapter for existing verification tools, not a dense memmap.

    Caches only the last requested (diagnosed Kx, Ky, time window), so repeated
    point-wise pair tests do not repeatedly decompress the same selected map.
    """
    def __init__(self, reader):
        self.reader = reader
        m = reader.layout
        self.shape = m.shape+(m.ny, reader.nt)
        self.dtype = m.dtype
        self.ndim = 6
        self._key = None
        self._selected = None

    def __getitem__(self, key):
        if not isinstance(key, tuple) or len(key) != 6:
            raise IndexError("FTZ verification view requires six explicit basic indices")
        cache_key = repr(key[3:])
        ids = [np.atleast_1d(np.arange(self.shape[j])[key[j]]) for j in (3,4,5)]
        if self._key != cache_key:
            ks, qs, ts = ids
            m = self.reader.layout
            out = np.empty(m.shape[:3]+(len(ks),len(qs),len(ts)),dtype=m.dtype,order="F")
            for ti,t in enumerate(ts):
                for qi,q in enumerate(qs):
                    out[:,:,:,:,qi,ti] = self.reader.block(int(t),int(q))[...,ks]
            self._selected = out
            self._key = cache_key
        tail = tuple(0 if np.isscalar(k) else slice(None) for k in key[3:])
        return self._selected[key[:3]+tail]


class RawFullT:
    """Read the current real-valued native format using explicit run metadata."""

    def __init__(self, case, toroidals_per_proc=None, suffix=".cgyro.fullt_asym"):
        self.case = Path(case)
        values = input_values(case)
        nx, ny, dx, dy = case_grid(case)
        nc = 1 if values.get("FULL_T_REAL_ONLY", 0) == 1 else 2
        ns = 1 if values.get("FULL_T_KX0", 1) == 1 else nx
        word = 8 if values.get("HIPREC_FLAG", 0) == 1 else 4
        self.layout = Layout(nx, ny, nc, ns, word, dx, dy)
        self.ntloc = int(toroidals_per_proc if toroidals_per_proc is not None else
                         values.get("TOROIDALS_PER_PROC", 1))
        if self.ntloc < 1 or ny % self.ntloc:
            raise FTZError("Invalid original TOROIDALS_PER_PROC; specify the run value")
        candidates = [self.case/("bin"+suffix), self.case/"bin"/("bin"+suffix)]
        candidates = [p for p in candidates if p.is_file()]
        if len(candidates) != 1:
            raise FTZError("Expected exactly one original FULLT file")
        self.path = candidates[0]
        self.times = case_times(self.case)
        self.nt = len(self.times)
        expected = self.nt*ny*self.layout.nbytes
        if self.path.stat().st_size != expected:
            raise FTZError("Raw file size conflicts with precision/channels/Kx/time metadata. "
                           "Legacy complex layouts require explicit conversion, not size guessing.")
        self._map = np.memmap(self.path, dtype=self.layout.dtype, mode="r",
                              shape=(nx, 2*ny-1, nc, self.ntloc, ns, ny//self.ntloc, self.nt),
                              order="F")
        stat = self.path.stat()
        self.signature = (str(self.path.resolve()), "raw", stat.st_size, stat.st_mtime_ns,
                          self.layout, self.ntloc, tuple(self.times))

    def block(self, t, q):
        if not 0 <= t < self.nt or not 0 <= q < self.layout.ny:
            raise IndexError("Raw FULLT time or diagnosed ky outside saved data")
        return self._map[:, :, :, q % self.ntloc, :, q//self.ntloc, t]

    read = FTZReader.read
    legacy = FTZReader.legacy
    __enter__ = FTZReader.__enter__
    __exit__ = FTZReader.__exit__

    def close(self):
        self._map._mmap.close()


def convert_case(case, output=None, toroidals_per_proc=None):
    """Convert stopped-run history; original data are never modified/deleted."""
    raw = RawFullT(case, toroidals_per_proc)
    final = Path(output) if output else Path(str(raw.path)+".ftz")
    if final.suffix != ".ftz":
        raise FTZError("Output must end in .ftz")
    final_index = final.with_suffix(".fti")
    partial = final.with_name(final.stem+".partial.ftz")
    partial_index = partial.with_suffix(".fti")
    if any(p.exists() for p in (final, final_index, partial, partial_index)):
        raw.close()
        raise FileExistsError("Destination/partial files already exist; inspect them first")
    initial = raw.path.stat()
    try:
        with FTZWriter(partial, raw.layout) as writer:
            for t, tv in enumerate(raw.times):
                writer.append((raw.block(t, q) for q in range(raw.layout.ny)), float(tv))
        # Reopen the actual persisted bytes, not just the pre-write buffers.
        checked = FTZReader(partial)
        for t in range(raw.nt):
            for q in range(raw.layout.ny):
                if checked.block(t, q).tobytes(order="F") != raw.block(t, q).tobytes(order="F"):
                    raise FTZError("Persisted conversion differs from original raw bytes")
        now = raw.path.stat()
        if (now.st_size, now.st_mtime_ns) != (initial.st_size, initial.st_mtime_ns):
            raise FTZError("Source changed during conversion; stop the simulation first")
        if final.exists() or final_index.exists():
            raise FileExistsError("Destination appeared during conversion")
        partial_index.rename(final_index)
        partial.rename(final)
        return {"data": str(final), "index": str(final_index), "records": raw.nt,
                "raw_bytes": initial.st_size,
                "compressed_bytes": final.stat().st_size+final_index.stat().st_size,
                "verified": "every persisted block compared byte-for-byte"}
    finally:
        raw.close()


def rebuild_index(path, output=None):
    """Scan self-framed committed data. Never overwrite the existing index."""
    path = Path(path)
    output = Path(output) if output else path.with_suffix(".recovered.fti")
    with path.open("rb") as df, output.open("xb") as ix:
        h = df.read(HS)
        layout = Layout.from_header(h)
        ix.write(layout.header(h[48:64], index=True))
        size = os.fstat(df.fileno()).st_size
        step = 1
        while df.tell() < size:
            start, entries, tv = df.tell(), [], None
            incomplete = False
            for q in range(layout.ny):
                off = df.tell()
                bh = df.read(BS)
                if len(bh) < BS:
                    incomplete = True
                    break
                info = block_info(bh, layout, q, step)
                payload = df.read(info[3])
                if len(payload) != info[3]:
                    incomplete = True
                    break
                decode_block(bh+payload, layout, q, step)
                if tv is not None and info[6] != tv:
                    raise FTZError("Mismatched times within a FTZ frame")
                tv = info[6]
                entries.append((off, BS+len(payload), info[4], info[5]))
            if incomplete:
                break
            footer = df.read(FS)
            if len(footer) < FS:
                break
            if footer != _footer(step, layout.ny, start):
                raise FTZError("Corrupt FTZ frame footer; recovered index stops here")
            ix.write(_index_record(step, tv, start, df.tell(), entries))
            step += 1
        _sync(ix)
    return {"index": str(output), "records": step-1}


def export_raw(path, output, toroidals_per_proc=1):
    """Explicit disk extraction, only on user request; exclusive output creation."""
    reader = FTZReader(path)
    m, nloc = reader.layout, int(toroidals_per_proc)
    if nloc < 1 or m.ny % nloc:
        raise FTZError("Invalid output TOROIDALS_PER_PROC")
    with Path(output).open("xb") as f:
        for t in range(reader.nt):
            for rank in range(m.ny//nloc):
                block = np.empty((m.nx, 2*m.ny-1, m.nc, nloc, m.ns), dtype=m.dtype, order="F")
                for local in range(nloc):
                    block[:, :, :, local, :] = reader.block(t, rank*nloc+local)
                f.write(block.tobytes(order="F"))
        _sync(f)


def export_origin(path, output, time_indices=None):
    """Streaming table; expanding a full diagnostic to text can be VERY large."""
    reader = FTZReader(path) if isinstance(path, (str, os.PathLike)) else path
    m = reader.layout
    tids = range(reader.nt) if time_indices is None else time_indices
    with Path(output).open("x", encoding="utf-8", newline="") as f:
        f.write("# Native FULLT; K is diagnosed, not a unique source; time is the output label, "
                "diagnostic is sampled at the last RHS stage; radial Nyquist labels can alias\n")
        f.write("time_index\ttime\tK_y_index\tK_y\tK_x_index\tK_x\tkx_prime_index\tkx_prime\t"
                "ky_prime_index\tky_prime\tchannel\tvalue\n")
        for t in tids:
            for q in range(m.ny):
                a = reader.block(int(t), q)
                for k in range(m.ns):
                    kp = 0 if m.ns == 1 else k-m.nx//2
                    for c in range(m.nc):
                        i, j = np.indices((m.nx, 2*m.ny-1))
                        i, j = i.ravel(), j.ravel()
                        p, yp = i-m.nx//2, j-(m.ny-1)
                        cols = np.column_stack([np.full(i.size, t), np.full(i.size, reader.times[t]),
                                                np.full(i.size, q), np.full(i.size, q*m.dy),
                                                np.full(i.size, kp), np.full(i.size, kp*m.dx),
                                                p, p*m.dx, yp, yp*m.dy, np.full(i.size, c),
                                                a[:, :, c, k].ravel()])
                        np.savetxt(f, cols, delimiter="\t", fmt="%.17g")
