"""
Pixel of the NGS/Focus PSFs with several science wavelengths.

The NGS/Focus PSFs come from MASTSEL psdSetToPsfSet on the PSD grid: their pixel is
wavelength * PSD step * oversampling. P3 builds the grid so that the science pixel psInMas
is exact at the SHORTEST science wavelength, but TIPTOP declared the LO PSF pixel as
psInMas * LO_wvl / wvlMax (LONGEST). With N wavelengths the declared pixel was too small by
wvlMax/wvlMin, so the NGS FWHM and EE (input of the LO noise model) were too small: with
MUSE_LTAO_10wl (480-930 nm) the NGS FWHM was 38.6 mas instead of ~75 mas. With
exactMultiWavelengthPSD=True the oversampling of the reference grid was also wrong, and the
branch where the PSF is computed at the native grid pixel (PSF pixel coarser than the LO
detector pixel) had the same problem.

Physical expectations:
- the pixel declared to the FWHM/EE estimators equals the actual pixel of the PSF returned
  by psdSetToPsfSet (Field.width / number of pixels), in every configuration;
- the NGS PSF is a property of the LO wavelength and of the LO atmosphere/AO residual only:
  its FWHM does not depend on the science wavelengths (MUSE_LTAO_10wl vs MUSE_LTAO).
"""

import os
import tempfile
import unittest
from configparser import ConfigParser

import numpy as np

import tiptop.baseSimulation as bs_module
from tiptop.baseSimulation import baseSimulation
from tiptop.tiptopUtils import cpuArray

PERF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tiptop', 'perfTest')
RAD2MAS = 180 / np.pi * 3600e3


def run_recorded(ini_name, lo=None, module=bs_module, **sim_kwargs):
    """Run a perfTest configuration recording declared and actual NGS/Focus PSF pixels.

    psdSetToPsfSet and _lo_psf_fwhm of `module` are wrapped (and restored) while the
    simulation runs. Only the calls made inside _compute_ngs_psf are kept, in call order.

    Parameters
    ----------
    ini_name : str
        Name (without extension) of a file in tiptop/perfTest.
    lo : dict, optional
        Keys to set in [sensor_LO].
    module : module
        Module holding baseSimulation (default: the current one).
    **sim_kwargs
        Extra keyword arguments of baseSimulation (e.g. exactMultiWavelengthPSD).

    Returns
    -------
    dict
        'actual': pixels [mas] of the PSFs returned by psdSetToPsfSet;
        'declared': pixel_scale [mas] passed to _lo_psf_fwhm, each paired with the actual
        pixel of the most recent psdSetToPsfSet call; 'pairs': list of (declared, actual);
        'fwhm': NGS_FWHM_mas_field; 'LO_wvl', 'wvl': wavelengths.
    """
    config = ConfigParser()
    config.optionxform = str
    config.read(os.path.join(PERF_DIR, ini_name + '.ini'))
    for k, v in (lo or {}).items():
        config.set('sensor_LO', k, str(v))

    state = {'active': False, 'last': None, 'actual': [], 'pairs': []}
    orig_psf, orig_fwhm = module.psdSetToPsfSet, module._lo_psf_fwhm
    orig_compute = module.baseSimulation._compute_ngs_psf

    def psf_wrap(*args, **kwargs):
        out = orig_psf(*args, **kwargs)
        if state['active']:
            f = out[0]
            state['last'] = float(f.width) / f.sampling.shape[0] * RAD2MAS
            state['actual'].append(state['last'])
        return out

    def fwhm_wrap(psf, pixel_scale, *args, **kwargs):
        if state['active']:
            state['pairs'].append((float(pixel_scale), state['last']))
        return orig_fwhm(psf, pixel_scale, *args, **kwargs)

    def compute_wrap(self, *args, **kwargs):
        state['active'] = True
        try:
            return orig_compute(self, *args, **kwargs)
        finally:
            state['active'] = False

    module.psdSetToPsfSet, module._lo_psf_fwhm = psf_wrap, fwhm_wrap
    module.baseSimulation._compute_ngs_psf = compute_wrap
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, 'sp.ini'), 'w') as f:
                config.write(f)
            sim = module.baseSimulation(tmpdir, 'sp', tmpdir, 'x', doPlot=False, verbose=False,
                                        **sim_kwargs)
            sim.doOverallSimulation()
            fwhm = np.atleast_1d(np.array(cpuArray(sim.NGS_FWHM_mas_field), dtype=float))
            wvl = np.atleast_1d(np.array(cpuArray(sim.wvl), dtype=float))
            lo_wvl = float(sim.LO_wvl)
    finally:
        module.psdSetToPsfSet, module._lo_psf_fwhm = orig_psf, orig_fwhm
        module.baseSimulation._compute_ngs_psf = orig_compute

    return dict(actual=state['actual'], pairs=state['pairs'], fwhm=fwhm, wvl=wvl, LO_wvl=lo_wvl)


class TestLoPsfPixelMultiWvl(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.single = run_recorded('MUSE_LTAO')
        cls.multi = run_recorded('MUSE_LTAO_10wl')
        cls.exact = run_recorded('MUSE_LTAO_10wl', exactMultiWavelengthPSD=True)
        cls.native = run_recorded('MUSE_LTAO_10wl', lo={'PixelScale': 5.0})

    def _assert_declared_is_actual(self, res):
        self.assertGreater(len(res['pairs']), 0)
        for declared, actual in res['pairs']:
            self.assertIsNotNone(actual)
            np.testing.assert_allclose(declared, actual, rtol=1e-9)

    def test_pixel_declared_is_actual_shared_grid(self):
        self._assert_declared_is_actual(self.multi)

    def test_pixel_declared_is_actual_exact_multiwavelength(self):
        self._assert_declared_is_actual(self.exact)

    def test_pixel_declared_is_actual_native_pixel_branch(self):
        self._assert_declared_is_actual(self.native)
        # PSF pixel coarser than the 5 mas LO detector pixel: native grid pixel is used
        np.testing.assert_allclose(self.native['pairs'][0][0], 7.527, atol=5e-3)

    def test_ngs_fwhm_independent_of_science_wavelengths(self):
        for res in (self.multi, self.exact):
            np.testing.assert_allclose(res['fwhm'], self.single['fwhm'], rtol=0.03)

    def test_single_wavelength_unchanged(self):
        np.testing.assert_allclose(self.single['pairs'][0][0], 14.0, atol=0.05)
        np.testing.assert_allclose(self.single['fwhm'], 75.9, rtol=0.02)
        self._assert_declared_is_actual(self.single)


if __name__ == '__main__':
    unittest.main()
