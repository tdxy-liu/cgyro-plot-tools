"""Optional private Xft Tk for Linux installations with a core-X11 Tk build.

Call before importing tkinter. No download, compilation, or system modification
occurs at startup. A runtime is enabled only for the Python extension it was
built and tested against; absent/incompatible installations keep their own Tk.
"""

import ctypes
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import sys

_RUNTIME = Path(__file__).resolve().parent / '.runtime' / 'tk-xft'
_LOADED = None


def prepare_tk_runtime():
    global _LOADED
    if _LOADED is not None:
        return True
    if not sys.platform.startswith('linux') or '_tkinter' in sys.modules:
        return False
    manifest = _RUNTIME / 'runtime.json'
    if not manifest.is_file():
        return False
    try:
        info = json.loads(manifest.read_text(encoding='utf-8'))
        spec = importlib.util.find_spec('_tkinter')
        if spec is None or not spec.origin:
            return False
        extension = Path(spec.origin).resolve()
        if (Path(info['python_prefix']).resolve() != Path(sys.prefix).resolve()
                or info['tkinter_sha256'] != hashlib.sha256(extension.read_bytes()).hexdigest()):
            logging.getLogger(__name__).warning('Private Tk belongs to another Python runtime; using its installed Tk.')
            return False
        library = _RUNTIME / 'lib' / 'libtk8.6.so'
        scripts = _RUNTIME / 'lib' / 'tk8.6'
        if not (scripts / 'tk.tcl').is_file():
            raise ValueError('missing Tk scripts')
        if hashlib.sha256(library.read_bytes()).hexdigest() != info['library_sha256']:
            raise ValueError('Tk library checksum mismatch')
        # The private library has SONAME libtk8.6.so so _tkinter reuses it,
        # rather than loading a second Tk library from the Conda RPATH.
        _LOADED = ctypes.CDLL(str(library), mode=ctypes.RTLD_GLOBAL)
        os.environ['TK_LIBRARY'] = str(scripts)
        return True
    except (OSError, ValueError, KeyError, TypeError) as exc:
        logging.getLogger(__name__).warning('Private Tk unavailable; using installed Tk: %s', exc)
        return False


def install_xft_runtime(source):
    """Explicit maintenance command; build from an unpacked official Tk source."""
    import _tkinter
    import re
    import shlex
    import shutil
    import subprocess
    import tempfile

    if not sys.platform.startswith('linux') or _tkinter.TK_VERSION != '8.6':
        raise RuntimeError('This installer supports Linux Python linked to Tk 8.6.')
    source = Path(source).resolve()
    header = (source / 'generic' / 'tk.h').read_text()
    if not re.search(r'#\s*define\s+TK_PATCH_LEVEL\s+"8\.6\.', header):
        raise ValueError('Provide an unpacked Tk 8.6 source tree.')
    prefix = Path(sys.prefix).resolve()
    if not (prefix / 'lib/tclConfig.sh').is_file() or not (prefix / 'include/tcl.h').is_file():
        raise RuntimeError('Matching Tcl headers and tclConfig.sh are required in this Python prefix.')
    if (_RUNTIME / 'runtime.json').exists():
        raise FileExistsError('A private Tk runtime is already installed: ' + str(_RUNTIME))
    compiler = shutil.which('gcc')
    if not compiler:
        raise RuntimeError('gcc is required to build Tk.')
    _RUNTIME.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, CC=compiler, CFLAGS='-O2 -I' + shlex.quote(str(prefix / 'include')), LDFLAGS='')
    with tempfile.TemporaryDirectory(prefix='tk-xft-build-', dir=_RUNTIME.parent) as folder:
        commands = (
            [str(source / 'unix/configure'), '--prefix=' + str(_RUNTIME),
             '--with-tcl=' + str(prefix / 'lib'), '--enable-threads', '--enable-xft'],
            ['make', '-j2', 'binaries', 'SHLIB_LD=' + shlex.quote(compiler) + ' -shared -Wl,-soname,libtk8.6.so'],
            ['make', 'install-binaries', 'install-libraries'],
        )
        for command in commands:
            subprocess.run(command, cwd=folder, env=env, check=True)
    library = _RUNTIME / 'lib/libtk8.6.so'
    info = {
        'python_prefix': str(prefix),
        'tkinter_sha256': hashlib.sha256(Path(_tkinter.__file__).read_bytes()).hexdigest(),
        'library_sha256': hashlib.sha256(library.read_bytes()).hexdigest(),
    }
    manifest = _RUNTIME / 'runtime.json'
    manifest.write_text(json.dumps(info, indent=2) + '\n', encoding='utf-8')
    print('Private Xft Tk installed. Restart the application:', _RUNTIME)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Install a private Xft Tk; does not modify Conda or system libraries.')
    parser.add_argument('--install-from', required=True, type=Path, help='Unpacked official Tk 8.6 source directory')
    install_xft_runtime(parser.parse_args().install_from)
