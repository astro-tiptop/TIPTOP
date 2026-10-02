import os
import unittest
from unittest import mock

import numpy as np

import tiptop.tiptopUtils as tu
from tiptop.baseSimulation import baseSimulation

try:
    import cupy as cp
    HAS_GPU = cp.cuda.runtime.getDeviceCount() > 0
except Exception:  # cupy missing or no CUDA device
    cp = None
    HAS_GPU = False


class TestArrayBackendConversion(unittest.TestCase):
    """arrayP3toMastsel / arrayMastseltoP3 with P3 and MASTSEL on different backends."""

    def _flags(self, gpu_p3, gpu_mastsel):
        return mock.patch.multiple(tu, gpuP3=gpu_p3, gpuMastsel=gpu_mastsel)

    def test_same_backend_is_identity(self):
        v = np.arange(4.0)
        with self._flags(False, False):
            self.assertIs(tu.arrayP3toMastsel(v), v)
            self.assertIs(tu.arrayMastseltoP3(v), v)

    def test_mastsel_to_p3_cpu_accepts_host_arrays(self):
        # some MASTSEL outputs are host arrays even when MASTSEL runs on GPU
        v = np.arange(4.0)
        with self._flags(False, True):
            out = tu.arrayMastseltoP3(v)
        self.assertIsInstance(out, np.ndarray)
        np.testing.assert_array_equal(out, v)

    @unittest.skipUnless(HAS_GPU, 'requires cupy and a CUDA device')
    def test_mixed_backends_round_trip(self):
        v = np.arange(4.0)
        with self._flags(False, True):
            on_gpu = tu.arrayP3toMastsel(v)
            self.assertIsInstance(on_gpu, cp.ndarray)
            back = tu.arrayMastseltoP3(on_gpu)
        self.assertIsInstance(back, np.ndarray)
        np.testing.assert_array_equal(back, v)
        with self._flags(True, False):
            self.assertIsInstance(tu.arrayMastseltoP3(v), cp.ndarray)


class _NaNFourierModel:
    """Stand-in for P3's fourierModel returning a non-finite HO PSD."""
    dtype = np.float64

    def __init__(self, *args, **kwargs):
        self.PSD = None

    def initComputations(self):
        self.PSD = np.full((8, 8, 1), np.nan)


class TestHoPsdFiniteCheck(unittest.TestCase):

    def test_non_finite_ho_psd_raises_clear_error(self):
        path = os.path.join(os.path.dirname(__file__), '..', 'tiptop', 'perfTest')
        sim = baseSimulation(path, 'MAVIStest', path, 'unused', doPlot=False)
        with mock.patch('tiptop.baseSimulation.fourierModel', _NaNFourierModel):
            with self.assertRaisesRegex(ValueError, 'non-finite'):
                sim._prepare_static_PSF_state(None)


if __name__ == '__main__':
    unittest.main()
