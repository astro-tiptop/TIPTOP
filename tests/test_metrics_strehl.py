import os
import tempfile
import unittest

import numpy as np

from p3.aoSystem.FourierUtils import getStrehl
from tiptop.baseSimulation import baseSimulation
from tiptop.tiptopUtils import cpuArray

"""
computeMetrics builds the diffraction-limited PSF once per wavelength instead
of calling getStrehl(method='max', psfInOnePix=True) for every PSF. The Strehl
values must stay bit-identical to getStrehl, and physically sane.
"""

RAD2MAS = 3600 * 180 * 1000 / np.pi
CASES = ['SOUL', 'ERIS', 'ERIS_LGS', 'MAVIS']

_cache = {}


def _run_case(case):
    """Run a perfTest case; return (sim, PSF copies taken before computeMetrics).

    Parameters
    ----------
    case : str
        Name of the .ini file in tiptop/perfTest (without extension).

    Returns
    -------
    tuple
        (simulation, list over wavelengths of lists of host PSF arrays)
    """
    if case in _cache:
        return _cache[case]
    with tempfile.TemporaryDirectory() as tmpdir:
        sim = baseSimulation('tiptop/perfTest', case, tmpdir, 'x',
                             doConvolve=True, verbose=False)
        sim.doOverallSimulation()
        if sim.nWvl > 1:
            psfs = [[cpuArray(img.sampling).copy() for img in sl] for sl in sim.results]
        else:
            psfs = [[cpuArray(img.sampling).copy() for img in sim.results]]
        # computeMetrics normalizes the result PSFs in place: copies are taken above
        sim.computeMetrics()
    _cache[case] = (sim, psfs)
    return _cache[case]


def _samp(sim, i):
    return sim.wvl[i] * RAD2MAS / (sim.psInMas * 2 * sim.tel_radius)


def _sr_list(sim):
    """Return the SR values as a list (per wavelength) of lists of floats."""
    if sim.nWvl > 1:
        return [[float(cpuArray(x)) for x in sl] for sl in sim.sr]
    return [[float(cpuArray(x)) for x in sim.sr]]


class TestMetricsStrehl(unittest.TestCase):

    def test_sr_identical_to_getStrehl(self):
        for case in CASES:
            with self.subTest(case=case):
                sim, psfs = _run_case(case)
                sr = _sr_list(sim)
                self.assertEqual(len(sr), len(psfs))
                for i, (psf_l, sr_l) in enumerate(zip(psfs, sr)):
                    self.assertEqual(len(psf_l), len(sr_l))
                    for psf, sr_new in zip(psf_l, sr_l):
                        sr_ref = float(cpuArray(getStrehl(
                            psf.copy(), sim.fao.ao.tel.pupil, _samp(sim, i),
                            method='max', psfInOnePix=True)))
                        self.assertEqual(sr_new, sr_ref)

    def test_sr_physical_range(self):
        for case in CASES:
            with self.subTest(case=case):
                sim, _ = _run_case(case)
                for sr_l in _sr_list(sim):
                    for sr in sr_l:
                        self.assertGreater(sr, 0)
                        self.assertLessEqual(sr, 1)

    def test_diffraction_limited_psf_has_near_unit_strehl(self):
        # psf_dl_array is MASTSEL's DL PSF at the reference wavelength, while the
        # reference peak is P3's analytic otf2psf(telescopeOtf): different pixel
        # integration/decimation gives 1-3 % differences, so 1e-3 is not attainable.
        for case in CASES:
            with self.subTest(case=case):
                sim, _ = _run_case(case)
                samp = sim.wvlRef * RAD2MAS / (sim.psInMas * 2 * sim.tel_radius)
                sr = float(cpuArray(getStrehl(
                    np.array(sim.psf_dl_array, dtype=float), sim.fao.ao.tel.pupil,
                    samp, method='max', psfInOnePix=True)))
                self.assertGreater(sr, 0.95)
                self.assertLessEqual(sr, 1.0 + 1e-3)


if __name__ == '__main__':
    unittest.main()
