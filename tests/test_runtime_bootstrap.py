"""Startup routing tests; synthetic modules isolate path selection from physics."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import unittest


class RuntimeBootstrap(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name)
        self.app=self.base/'tool'; self.app.mkdir()
        shutil.copyfile(Path(__file__).resolve().parents[1]/'cgyro_comparison_bootstrap.py',
                        self.app/'cgyro_comparison_bootstrap.py')
        self.old=self.make_root('gacode','old')
        self.new=self.make_root('new-root','new')
        self.config=self.app/'cgyro_runtime.local.json'
        self.config.write_text(json.dumps({'gacode_root':str(self.new)}),encoding='utf-8')

    def make_root(self,name,marker):
        root=self.base/name; package=root/'f2py/pygacode'
        (package/'cgyro').mkdir(parents=True)
        (package/'__init__.py').write_text(f"marker={marker!r}\n")
        (package/'cgyro/__init__.py').write_text('')
        (package/'cgyro/data.py').write_text('class cgyrodata: pass\n')
        (package/'cgyro/data_plot.py').write_text('class cgyrodata_plot: pass\n')
        return root

    def run_import(self,extra='',explicit=None,preload=False):
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
        env.pop('GACODE_ROOT',None)
        env.pop('CGYRO_ALLOW_MOCK_DATA',None)
        env['PYTHONPATH']=os.pathsep.join([str(self.app),str(self.old/'f2py')])
        if explicit is not None:
            env['GACODE_ROOT']=str(explicit)
        code=('import pygacode; ' if preload else '')
        code+='import cgyro_comparison_bootstrap; import pygacode; print(pygacode.marker); '+extra
        return subprocess.run([sys.executable,'-c',code],cwd=self.base,env=env,
                              capture_output=True,text=True,timeout=30)

    def test_local_root_beats_old_pythonpath_without_environment(self):
        r=self.run_import()
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(r.stdout.strip(),'new')

    def test_relative_root_is_resolved_against_tool_not_cwd(self):
        self.config.write_text(json.dumps({'gacode_root':'../new-root'}))
        r=self.run_import()
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(r.stdout.strip(),'new')

    def test_explicit_environment_beats_local_config_even_if_malformed(self):
        self.config.write_text('not JSON')
        r=self.run_import(explicit=self.old)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(r.stdout.strip(),'old')

    def test_no_local_config_preserves_legacy_fallback(self):
        self.config.unlink()
        r=self.run_import()
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(r.stdout.strip(),'old')

    def test_invalid_config_does_not_silently_fallback(self):
        for text in ('not JSON','[]','{}','{"gacode_root": 1}',
                     '{"gacode_root": ""}','{"gacode_root": "missing-root"}'):
            with self.subTest(config=text):
                self.config.write_text(text)
                r=self.run_import()
                self.assertNotEqual(r.returncode,0)
                self.assertNotIn('old',r.stdout)
                self.assertIn('RuntimeError',r.stderr)

    def test_loaded_old_pygacode_is_rejected(self):
        r=self.run_import(preload=True)
        self.assertNotEqual(r.returncode,0)
        self.assertIn('older pygacode is already loaded',r.stderr)

    def test_private_dependencies_only_use_current_abi(self):
        abi=sys.implementation.cache_tag+'-'+sysconfig.get_platform()
        current=self.app/'.runtime/python'/abi/'zstandard'
        wrong=self.app/'.runtime/python/wrong-abi/zstandard'
        for path in (current,wrong):
            path.mkdir(parents=True)
            (path/'__init__.py').write_text('marker='+repr(path.parent.name))
        r=self.run_import("import zstandard; print(zstandard.marker)")
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(r.stdout.splitlines(),['new',abi])
        (current/'__init__.py').unlink()  # Remove only our synthetic package entry.
        r=self.run_import("import sys; print(any('wrong-abi' in p for p in sys.path))")
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(r.stdout.splitlines(),['new','False'])


if __name__=='__main__':
    unittest.main()
