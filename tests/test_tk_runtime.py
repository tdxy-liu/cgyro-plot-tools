"""Private GUI runtimes must never cross Python ABI or load after Tk starts."""

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch

import cgyro_tk_runtime as runtime


class PrivateTkTests(unittest.TestCase):
    def test_missing_runtime_does_not_load_or_change_environment(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(runtime, '_RUNTIME', Path(folder)), \
                patch.object(runtime, '_LOADED', None), \
                patch.object(runtime, 'sys', SimpleNamespace(platform='linux', modules={}, prefix=sys.prefix)), \
                patch.object(runtime.ctypes, 'CDLL') as load, \
                patch.dict(runtime.os.environ, {}, clear=True):
            self.assertFalse(runtime.prepare_tk_runtime())
            load.assert_not_called()
            self.assertNotIn('TK_LIBRARY', runtime.os.environ)

    def test_abi_and_checksum_guards_prevent_loading_incompatible_library(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            scripts = root / 'lib/tk8.6'
            scripts.mkdir(parents=True)
            (scripts / 'tk.tcl').write_text('fixture')
            extension = root / '_tkinter.so'; extension.write_bytes(b'python ABI')
            library = root / 'lib/libtk8.6.so'; library.write_bytes(b'private Tk')
            good = {'python_prefix': sys.prefix,
                    'tkinter_sha256': hashlib.sha256(extension.read_bytes()).hexdigest(),
                    'library_sha256': hashlib.sha256(library.read_bytes()).hexdigest()}
            with patch.object(runtime, '_RUNTIME', root), patch.object(runtime, '_LOADED', None), \
                    patch.object(runtime, 'sys', SimpleNamespace(platform='linux', modules={}, prefix=sys.prefix)), \
                    patch.object(runtime.importlib.util, 'find_spec', return_value=SimpleNamespace(origin=str(extension))), \
                    patch.object(runtime.ctypes, 'CDLL') as load, patch.dict(runtime.os.environ, {}, clear=True):
                for changes in ({'tkinter_sha256': 'wrong'}, {'python_prefix': str(root)}, {'library_sha256': 'wrong'}):
                    (root / 'runtime.json').write_text(json.dumps(dict(good, **changes)))
                    with self.assertLogs('cgyro_tk_runtime', level='WARNING'):
                        self.assertFalse(runtime.prepare_tk_runtime())
                    load.assert_not_called()
                    self.assertNotIn('TK_LIBRARY', runtime.os.environ)
                (root / 'runtime.json').write_text(json.dumps(good))
                self.assertTrue(runtime.prepare_tk_runtime())
                self.assertTrue(runtime.prepare_tk_runtime())
                load.assert_called_once()
                self.assertEqual(runtime.os.environ['TK_LIBRARY'], str(scripts))

    def test_an_existing_tk_interpreter_is_never_replaced(self):
        with patch.object(runtime, 'sys', SimpleNamespace(platform='linux', modules={'_tkinter': object()})), \
                patch.object(runtime, '_LOADED', None), \
                patch.object(runtime.ctypes, 'CDLL') as load:
            self.assertFalse(runtime.prepare_tk_runtime())
            load.assert_not_called()


if __name__ == '__main__':
    unittest.main()
