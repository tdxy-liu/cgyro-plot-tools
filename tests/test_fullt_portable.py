"""Bundled codec tests: no custom pygacode, CGYRO executable or MPI required."""
import builtins
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import cgyro_fullt_reader as codec
from cgyro_fullt_compressed_bridge import open_fullt, open_compressed


class PortableFullT(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.case = Path(self.tmp.name)
        self.path = self.case / 'bin.cgyro.fullt_asym.ftz'

    def test_bitwise_all_save_options_and_special_values(self):
        for word in (4, 8):
            special = ([0, 0x80000000, 1, 0x80000001, 0x7f800000, 0xff800000,
                        0x7fc00001, 0x7f800001] if word == 4 else
                       [0, 0x8000000000000000, 1, 0x8000000000000001,
                        0x7ff0000000000000, 0xfff0000000000000,
                        0x7ff8000000000001, 0x7ff0000000000001])
            for nc in (1, 2):
                for ns in (1, 8):
                    m = codec.Layout(8, 4, nc, ns, word, .03, .1)
                    words = np.resize(np.asarray(special, dtype=m.uint), m.nbytes // word)
                    # Deliberately asymmetric pairs, including boundary slots.
                    words[::7] ^= np.asarray(3, dtype=m.uint)
                    raw = words.tobytes()
                    for q in range(m.ny):
                        with self.subTest(word=word, nc=nc, ns=ns, q=q):
                            b = codec.encode_block(raw, m, q, 1, 1.25)
                            restored = codec.decode_block(b, m, q, 1)
                            self.assertEqual(restored.tobytes(order='F'), raw)

    def test_incompressible_raw_fallback(self):
        m = codec.Layout(2, 1, 1, 1, 4)
        raw = bytes(range(m.nbytes))
        block = codec.encode_block(raw, m, 0, 1, 1.)
        self.assertEqual(codec.block_info(block[:codec.BS], m)[2], 0)
        with patch.object(codec, '_zstd', side_effect=AssertionError('must stay lazy')):
            self.assertEqual(codec.decode_block(block, m).tobytes(order='F'), raw)

    def test_missing_dependency_does_not_look_like_corruption(self):
        m = codec.Layout(8, 4, 1, 1, 4)
        b = codec.encode_block(bytes(m.nbytes), m, 1, 1, 1.)
        self.assertEqual(codec.block_info(b[:codec.BS], m)[2], 1)
        original = builtins.__import__
        def deny_zstd(name, *args, **kwargs):
            if name == 'zstandard':
                raise ModuleNotFoundError('test missing decoder')
            return original(name, *args, **kwargs)
        with patch('builtins.__import__', side_effect=deny_zstd):
            with self.assertRaisesRegex(codec.FTZError, 'requirements-compression.txt'):
                codec.decode_block(b, m)

    def test_bundled_reader_ignores_external_pygacode_codec(self):
        m = codec.Layout(8, 4, 1, 1, 4)
        arrays = [np.full(m.shape, q, dtype=m.dtype) for q in range(m.ny)]
        with codec.FTZWriter(self.path, m) as w:
            w.append(arrays, 1.)
        (self.case / 'input.cgyro').write_text('FULL_T_ASYM_COMPRESSION=1D0\n')
        original = builtins.__import__
        def deny_external(name, *args, **kwargs):
            if name.startswith('pygacode'):
                raise AssertionError('GUI FTZ must not import external pygacode')
            return original(name, *args, **kwargs)
        with patch('builtins.__import__', side_effect=deny_external):
            with open_fullt(self.case) as reader:
                self.assertIsInstance(reader, codec.FTZReader)
                self.assertEqual(reader.block(0, 2).tobytes(), arrays[2].tobytes())
            with open_compressed(self.case) as reader:
                self.assertEqual(reader.nt, 1)

    def test_rewind_and_uncommitted_tail(self):
        m = codec.Layout(8, 4, 1, 1, 4)
        blocks = [np.zeros(m.shape, dtype=m.dtype)] * m.ny
        with codec.FTZWriter(self.path, m) as w:
            w.append(blocks, 1.)
            first = codec.FTZReader(self.path)
            w.append(blocks, 2.)
        last = codec.FTZReader(self.path)
        self.assertNotEqual(first.signature, last.signature)
        rn = 48 + 24 * m.ny
        # Model a restart's rollback on only our temporary fixture.
        with self.path.open('r+b') as df, self.path.with_suffix('.fti').open('r+b') as ix:
            df.truncate(first.committed_end)
            ix.truncate(codec.HS + rn)
        last.refresh()
        self.assertEqual(last.nt, 1)
        with self.path.open('ab') as df, self.path.with_suffix('.fti').open('ab') as ix:
            df.write(b'uncommitted data')
            ix.write(b'incomplete index')
        self.assertEqual(codec.FTZReader(self.path).nt, 1)

    def test_full_commit_corruption_and_identity_fail(self):
        m = codec.Layout(8, 4, 1, 1, 4)
        with codec.FTZWriter(self.path, m) as w:
            w.append([np.zeros(m.shape, dtype=m.dtype)] * m.ny, 1.)
        index = self.path.with_suffix('.fti')
        original = index.read_bytes()
        bad = bytearray(original)
        bad[-1] ^= 1
        index.write_bytes(bad)
        with self.assertRaisesRegex(codec.FTZError, 'committed FTZ index'):
            codec.FTZReader(self.path)
        index.write_bytes(m.header(identity=b'X' * 16, index=True) + original[codec.HS:])
        with self.assertRaisesRegex(codec.FTZError, 'identity mismatch'):
            codec.FTZReader(self.path)

    def test_nondivisible_radial_box_uses_integer_division(self):
        nx, ny, nt, box = 8, 4, 2, 3
        g = [ny, 2, 3, nx, nt, 1, 1, box, 2*np.pi/.03, 0, 1]
        g += list(np.arange(nx)-nx//2) + [0] * (nt+2+nt*(nx//box))
        g += list(np.arange(ny)*.1)
        np.savetxt(self.case/'out.cgyro.grids', g)
        with self.assertWarnsRegex(RuntimeWarning, 'not a multiple'):
            ax = codec.case_grid(self.case)
        np.testing.assert_allclose(ax, (nx, ny, .03, .1))


if __name__ == '__main__':
    unittest.main()
