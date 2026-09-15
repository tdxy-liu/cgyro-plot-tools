"""Headless raw/FTZ GUI regression tests; no Tk windows are created."""
import os
import importlib.util
import sys
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
import numpy as np
from cgyro_comparison_plotting_energy import EnergyPlotting
from cgyro_data_export import CgyroDataExportMixin
from cgyro_fullt_compressed_bridge import open_compressed
from cgyro_fullt_reader import Layout, FTZWriter, FTZReader, FTZError, FTZArrayView

# External GACODE integration is optional; ordinary clone tests must not need
# a private/custom GACODE checkout. Fail (rather than skip) when opted in.
GACODE_INTEGRATION = os.environ.get('CGYRO_TEST_GACODE_INTEGRATION') == '1'


class Host(EnergyPlotting, CgyroDataExportMixin):
    def _resolve_case_dir(self, data):
        return data.dir

    def _nearest_source_ky_index(self, axis, target):
        return int(np.argmin(np.abs(np.abs(axis)-abs(target))))

    def _display_source_ky_value(self, axis, index, requested):
        return float(axis[index])


class CompressedGUIRegression(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.case=Path(self.tmp.name)
        self.path=self.case/"bin.cgyro.fullt_asym.ftz"
        self.m=Layout(8,4,2,8,4,.03,.1)
        self.arrays=[]
        self.writer=FTZWriter(self.path,self.m)
        self.addCleanup(self.writer.close)
        for t in range(3):
            blocks=[]
            for q in range(4):
                a=(np.arange(np.prod(self.m.shape),dtype="f4")+t*10000+q*1000).reshape(self.m.shape,order="F")
                blocks.append(a)
            self.arrays.append(blocks)
            self.writer.append(blocks,float(t+1))
        (self.case/"input.cgyro").write_text("FULL_T_ASYM_COMPRESSION=1\n")
        self.data=SimpleNamespace(dir=str(self.case),n_n=4,n_radial=8,
            ky=np.arange(4)*.1,p=np.arange(-4,4),length=2*np.pi/.03,t=np.arange(1,4))
        self.host=Host()

    def test_map_trace_and_selected_time_equivalence(self):
        d=self.data
        self.assertTrue(self.host._load_fullt_if_needed(
            d,"test",asym=True,source_ky_value=.2,source_kx_value=.03,time_indices=[0,2]))
        expected=np.stack([self.arrays[t][2][:,:,:,5] for t in (0,2)],axis=-1)
        np.testing.assert_array_equal(d.fullt_asym[0,0],expected)
        trace=self.host._load_fullt_trace_source_slice(d,"test",True,.2,[0,2])
        expected_trace=np.stack([self.arrays[t][2][:,:,0,:] for t in (0,2)],axis=-1).transpose(2,0,1,3)
        np.testing.assert_array_equal(trace["trace"],expected_trace)
        np.testing.assert_allclose(trace["source_kx_axis"],np.arange(-4,4)*.03)

    def test_append_invalidates_map_cache(self):
        d=self.data
        self.host._load_fullt_if_needed(d,"test",True,.2,0.,[0])
        old=d.fullt_asym_ftz_signature
        self.writer.append(self.arrays[0],4.)
        self.host._load_fullt_if_needed(d,"test",True,.2,0.,[0])
        self.assertNotEqual(old,d.fullt_asym_ftz_signature)

    def test_corruption_is_not_synthetic_fallback(self):
        d=self.data
        self.host._load_fullt_if_needed(d,"test",True,.2,0.,[0])
        r=FTZReader(self.path)
        off=r.records[0][4][2][0]+64
        with self.path.open("r+b") as f:
            f.seek(off); b=f.read(1); f.seek(off); f.write(bytes([b[0]^1]))
        with self.assertRaises(FTZError):
            self.host._load_fullt_if_needed(d,"test",True,.2,0.,[0])

    def test_lazy_verifier_view(self):
        view=FTZArrayView(FTZReader(self.path))
        actual=view[:,:,1,5,2,slice(0,3,2)]
        expected=np.stack([self.arrays[t][2][:,:,1,5] for t in (0,2)],axis=-1)
        np.testing.assert_array_equal(actual,expected)
        np.testing.assert_array_equal(view[3,2,0,5,2,slice(0,3,2)],
            [self.arrays[t][2][3,2,0,5] for t in (0,2)])

    def test_export_skips_index_and_does_not_extract_whole_raw_file(self):
        suffixes=self.host._collect_bin_suffixes(str(self.case))
        self.assertEqual(suffixes,[".cgyro.fullt_asym"])
        output=self.case/"export"
        output.mkdir()
        self.assertTrue(self.host._export_generic_bin_file(self.data,suffixes[0],str(output)))
        path=output/"cgyro.fullt_asym.txt"
        lines=path.read_text().splitlines()
        self.assertIn("K_y_index",lines[1])
        self.assertEqual(len(lines)-2,3*4*np.prod(self.m.shape))
        self.assertFalse((self.case/"bin.cgyro.fullt_asym").exists())

    @unittest.skipUnless(GACODE_INTEGRATION, 'optional updated-GACODE integration')
    def test_pygacode_explicit_extraction(self):
        from pygacode.cgyro.data import cgyrodata
        d=object.__new__(cgyrodata)
        d.dir=str(self.case)+os.sep
        _,fmt,arr=d.extract(".cgyro.fullt_asym")
        self.assertEqual(fmt,"ftz")
        expected=np.empty(self.m.shape+(4,3),dtype=self.m.dtype)
        for t in range(3):
            for q in range(4):
                expected[:,:,:,:,q,t]=self.arrays[t][q]
        self.assertEqual(arr.tobytes(),expected.tobytes(order="F"))

    @unittest.skipUnless(GACODE_INTEGRATION, 'optional updated-GACODE integration')
    def test_existing_verifiers_use_lazy_ftz_data(self):
        import pygacode.cgyro.fullt_compressed as codec
        root=Path(codec.__file__).resolve().parents[3]
        def load(name):
            spec=importlib.util.spec_from_file_location(name,root/"cgyro"/"tools"/(name+".py"))
            module=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.add_pygacode_path()
            return module
        sim=self.data
        sim.n_time=3
        sim.BYTE="float32"
        pair=load("verify_fullt_asym_pair")
        record=pair.open_fullt(str(self.case),sim,".cgyro.fullt_asym")
        np.testing.assert_array_equal(record["data"][:,:,0,4,2,slice(None)],
            np.stack([self.arrays[t][2][:,:,0,4] for t in range(3)],axis=-1))
        triad=load("verify_fullt_triad_sum")
        summary=triad.load_fullt_source_summary(str(self.case),sim,".cgyro.fullt_asym",
                                              2,slice(None),.5,np.arange(-3,4))
        expected=[.5*np.sum(self.arrays[t][2][:,:,0,4],dtype=np.float64) for t in range(3)]
        np.testing.assert_array_equal(summary["full_sum"],expected)

    def test_missing_selected_format_is_an_error(self):
        (self.case/"input.cgyro").write_text("FULL_T_ASYM_COMPRESSION=0\n")
        with self.assertRaises(FTZError):
            open_compressed(self.case)

    def test_raw_mpi_map_trace_export(self):
        # The FTZ copy remains present, but flag=0 must select the raw MPI layout.
        m=self.m
        g=[m.ny,2,3,m.nx,2,1,1,2,2*np.pi/m.dx,0,1]
        g+=list(np.arange(-4,4))+[0]*4+[0]*8+list(np.arange(4)*m.dy)
        np.savetxt(self.case/'out.cgyro.grids',g)
        np.savetxt(self.case/'out.cgyro.time',np.column_stack([np.arange(1,4),np.zeros((3,3))]))
        for nloc in (1,2,4):
            with self.subTest(toroidals_per_proc=nloc):
                (self.case/'input.cgyro').write_text(
                    f'FULL_T_ASYM_COMPRESSION=0\nFULL_T_REAL_ONLY=0\nFULL_T_KX0=0\n'
                    f'HIPREC_FLAG=0\nTOROIDALS_PER_PROC={nloc}\n')
                with (self.case/'bin.cgyro.fullt_asym').open('wb') as f:
                    for blocks in self.arrays:
                        for rank in range(4//nloc):
                            f.write(np.stack(blocks[rank*nloc:(rank+1)*nloc],axis=3).tobytes(order='F'))
                self.test_map_trace_and_selected_time_equivalence()
                if GACODE_INTEGRATION:
                    self.test_existing_verifiers_use_lazy_ftz_data()
                    from pygacode.cgyro.data import cgyrodata
                    d=object.__new__(cgyrodata); d.dir=str(self.case)+os.sep
                    _,fmt,arr=d.extract('.cgyro.fullt_asym')
                    self.assertEqual(fmt,'bin')
                    expected=np.stack([np.stack(blocks,axis=-1) for blocks in self.arrays],axis=-1)
                    self.assertEqual(arr.tobytes(),expected.tobytes(order='F'))
                output=self.case/f'raw_export_{nloc}'; output.mkdir()
                self.assertTrue(self.host._export_generic_bin_file(self.data,'.cgyro.fullt_asym',str(output)))
                table=np.loadtxt(output/'cgyro.fullt_asym.txt',skiprows=2)
                for row in table[::137]:
                    t,q,k,p,y,c=int(row[0]),int(row[2]),int(row[4])+4,int(row[6])+4,int(row[8])+3,int(row[10])
                    self.assertEqual(row[11],self.arrays[t][q][p,y,c,k])

    def test_stale_gui_axes_and_times_are_rejected(self):
        self.data.ky*=2
        with self.assertRaisesRegex(FTZError,'ky axis'):
            self.host._load_fullt_if_needed(self.data,'test',True,.2,0.,[0])
        self.data.ky/=2; self.data.t=np.array([999.,2.,3.])
        with self.assertRaisesRegex(FTZError,'time prefix'):
            self.host._load_fullt_trace_source_slice(self.data,'test',True,.2,[0])

    @unittest.skipUnless(GACODE_INTEGRATION, 'optional updated-GACODE integration')
    def test_explicit_gacode_root_precedence_and_bad_root(self):
        import subprocess
        import pygacode.cgyro.fullt_compressed as codec
        root=Path(codec.__file__).resolve().parents[3]
        env=dict(os.environ,GACODE_ROOT=str(root))
        # It must still win if PYTHONPATH already contained this candidate.
        env['PYTHONPATH']=os.pathsep.join([str(root/'f2py'),str(Path(__file__).resolve().parents[1]),
                                        env.get('PYTHONPATH','')])
        cmd=[sys.executable,'-c','import cgyro_comparison_bootstrap; import pygacode; print(pygacode.__file__)']
        result=subprocess.run(cmd,env=env,capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn(str(root),result.stdout)
        env['GACODE_ROOT']=str(self.case/'not-a-repository')
        result=subprocess.run(cmd,env=env,capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('GACODE_ROOT does not contain',result.stderr)


if __name__=="__main__":
    unittest.main()
