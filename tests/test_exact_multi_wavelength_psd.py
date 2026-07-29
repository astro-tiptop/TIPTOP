import unittest
import tempfile
import os
import numpy as np
from configparser import ConfigParser

from tiptop.baseSimulation import baseSimulation

"""
End-to-end regression tests for exactMultiWavelengthPSD=True: P3 computes one
exact PSD grid per science wavelength (frequencyDomain.wvl_grids,
psdPerWavelength=True) instead of one shared, approximate grid, and MASTSEL's
psdSetToPsfSet decimates each wavelength's PSF exactly instead of resampling
with a lossy affine_transform.

Reuses the same lightweight ERIS config as
TestTiptop.test_multiwavelength_eris_undersampled_sr_fwhm in test_all.py
(8m telescope, wavelengths spanning an undersampled-to-oversampled range),
since that is exactly the "lambda_min << lambda_max" regime this feature
targets.
"""


def _make_config(wvl_list):
    config = ConfigParser()
    config.optionxform = str
    config.read('tiptop/perfTest/ERIS.ini')
    config.set('sources_science', 'Wavelength', str(wvl_list))
    config.set('sensor_science', 'PixelScale', '14')
    config.set('sensor_science', 'FieldOfView', '180')
    config.set('telescope', 'Resolution', '64')
    config.set('sensor_HO', 'NumberLenslets', '[16]')

    fd, path = tempfile.mkstemp(suffix='.ini')
    with os.fdopen(fd, 'w') as f:
        config.write(f)
    return path


def _run(wvl_list, tag, exact_multi_wavelength_psd=False):
    ini_path = _make_config(wvl_list)
    try:
        temp_dir = os.path.dirname(ini_path)
        temp_basename = os.path.splitext(os.path.basename(ini_path))[0]
        simulation = baseSimulation(
            temp_dir, temp_basename, temp_dir, tag,
            doConvolve=False, doPlot=False, verbose=False,
            exactMultiWavelengthPSD=exact_multi_wavelength_psd)
        simulation.doOverallSimulation()
        simulation.computeMetrics()
        return simulation
    finally:
        os.remove(ini_path)


class TestExactMultiWavelengthPSD(unittest.TestCase):

    def test_legacy_default_unaffected(self):
        """exactMultiWavelengthPSD defaults to False: self.PSD stays a plain
        array, multiGridPSD is False, matching pre-feature behaviour."""
        sim = _run([0.865e-6, 2.21e-6], 'testLegacyDefault')
        self.assertFalse(sim.multiGridPSD)
        self.assertFalse(isinstance(sim.PSD, list))

    def test_exact_multi_wavelength_runs_and_produces_valid_metrics(self):
        sim = _run([0.865e-6, 2.21e-6], 'testExactRuns', exact_multi_wavelength_psd=True)
        self.assertTrue(sim.multiGridPSD)
        self.assertTrue(isinstance(sim.PSD, list))
        self.assertEqual(len(sim.PSD), 2)

        self.assertEqual(sim.cubeResultsArray.ndim, 4)
        self.assertEqual(sim.cubeResultsArray.shape[0], 2)

        for i in range(2):
            sr_i = np.array(sim.sr[i]).ravel()
            fwhm_i = np.array(sim.fwhm[i]).ravel()
            self.assertTrue(np.all(np.isfinite(sr_i)))
            self.assertTrue(np.all(sr_i > 0))
            self.assertTrue(np.all(sr_i <= 1.0))
            self.assertTrue(np.all(np.isfinite(fwhm_i)))

    def test_each_wavelength_matches_standalone_single_wavelength_run(self):
        """
        The actual point of this feature: a wavelength's slice from a
        multi-wavelength exactMultiWavelengthPSD=True run must be bit-for-bit
        identical to a standalone run requesting only that one wavelength --
        exercising the full P3 -> MASTSEL -> TIPTOP pipeline, not just the
        PSD (as in P3's own tests) or synthetic inputs (as in MASTSEL's).
        """
        multi = _run([0.865e-6, 2.21e-6], 'testExactMulti', exact_multi_wavelength_psd=True)

        for idx, wvl in enumerate([0.865e-6, 2.21e-6]):
            with self.subTest(wvl_nm=wvl * 1e9):
                mono = _run([wvl], f'testExactMono{int(wvl*1e9)}',
                           exact_multi_wavelength_psd=True)
                multi_slice = np.asarray(multi.cubeResultsArray[idx])
                mono_slice = np.asarray(mono.cubeResultsArray)
                self.assertEqual(multi_slice.shape, mono_slice.shape)
                np.testing.assert_array_equal(
                    multi_slice, mono_slice,
                    err_msg=f"{wvl*1e9:.0f}nm slice differs from standalone run"
                )


if __name__ == '__main__':
    unittest.main()
