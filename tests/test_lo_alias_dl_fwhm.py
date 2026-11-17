"""
Reference diffraction-limited FWHM of the LO aliasing term (addAliasError = True in
[sensor_LO]): MASTSEL uses 4 * (FWHM_NGS - FWHM_DL) with FWHM_DL = baseSimulation.NGS_DL_FWHM_mas.
Regression: FWHM_DL used to come from longExposurePsf on the native N grid (pixel
lambda / grid_diameter) but was measured with the NGS PSF pixel LO_PSFsInMas, giving a
FWHM of up to ~2 lambda/D and a negative FWHM_NGS - FWHM_DL for any NGS narrower than that.
Physical expectations:
- the FWHM of a perfect PSF is ~1 lambda/D (1.029 for an unobstructed disc, ~0.99 for an
  ELT-like obstruction), in any case far from 2;
- a perfect PSF is the narrowest one: FWHM_DL <= FWHM_NGS, so the aliasing term is >= 0.
"""

import os
import tempfile
import unittest
from configparser import ConfigParser

import numpy as np

from tiptop.baseSimulation import baseSimulation
from tiptop.tiptopUtils import cpuArray

PERF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tiptop', 'perfTest')
BASE_INI = 'HarmoniLTAO_2'  # single-aperture LO, D = 39 m, eps = 0.28, modified to LO 1400 nm
RAD2MAS = 180 / np.pi * 3600e3


def _run_sim(simcls=baseSimulation):
    """Run HarmoniLTAO_2 with addAliasError = True; return (DL FWHM, NGS FWHM, lambda/D) in mas."""
    config = ConfigParser()
    config.optionxform = str
    config.read(os.path.join(PERF_DIR, BASE_INI + '.ini'))
    config.set('sensor_LO', 'addAliasError', 'True')
    # ANDES-like: the LO PSF pixel (psInMas * LO_wvl / wvlMax = 5.49 mas) is twice the
    # native grid pixel (2.745 mas) and is kept because the LO sensor pixel (7.5 mas) is
    # larger; with a LO pixel < 5.49 mas the PSF is oversampled to the native pixel and the
    # old estimate was correct
    config.set('sensor_LO', 'PixelScale', '7.5')
    config.set('atmosphere', 'Seeing', '2.5')
    config.set('sources_science', 'Wavelength', '[1.02e-06]')
    config.set('sources_LO', 'Wavelength', '[1400e-9]')
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(os.path.join(tmpdir, 'sp.ini'), 'w') as f:
            config.write(f)
        sim = simcls(tmpdir, 'sp', tmpdir, 'x', doPlot=False, verbose=False)
        sim.doOverallSimulation()
        dl = np.atleast_1d(np.array(cpuArray(sim.NGS_DL_FWHM_mas), dtype=float))
        ngs = np.atleast_1d(np.array(cpuArray(sim.NGS_FWHM_mas_field), dtype=float))
        lam_over_d = sim.LO_wvl / (2 * sim.tel_radius) * RAD2MAS
    return dl, ngs, lam_over_d


class TestLoAliasDlFwhm(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dl, cls.ngs, cls.lam_over_d = _run_sim()

    def test_dl_fwhm_is_about_one_lambda_over_d(self):
        ratio = self.dl / self.lam_over_d
        self.assertTrue(np.all(ratio > 0.95) and np.all(ratio < 1.05), ratio)

    def test_dl_fwhm_not_larger_than_ngs_fwhm(self):
        # aliasing term 4 * (FWHM_NGS - FWHM_DL) must not be negative
        self.assertTrue(np.all(self.ngs >= self.dl), (self.ngs, self.dl))


if __name__ == '__main__':
    unittest.main()
