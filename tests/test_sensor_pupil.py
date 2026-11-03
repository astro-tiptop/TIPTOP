"""
[sensor_LO] / [sensor_Focus] own pupil, built by P3 aoSystem (aoSystem.pupilLO,
aoSystem.pupilFocus) and used by the NGS PSFs: PathPupil (FITS over TelescopeDiameter, rotated
by PupilAngle) or ObscurationRatio (annulus of diameter TelescopeDiameter), otherwise
[sensor_Focus] -> [sensor_LO] pupil and [sensor_LO] -> telescope pupil. The pupil logic is
tested in P3 (tests/test_sensor_pupils.py); here the outputs. Physically predictable results:
- an LO pupil identical to the telescope one reproduces the default NGS metrics exactly;
- the diffraction-limited FWHM scales as lambda / D_pupil: a filled disc of diameter D/2
  gives twice the FWHM of a filled disc of diameter D.
Single-aperture LO (NumberLenslets = [1]) is used for the output-level tests: with
NumberLenslets >= 4 MASTSEL maskSA returns a full square that ignores the pupil.
"""

import os
import tempfile
import unittest
from configparser import ConfigParser

import numpy as np
from astropy.io import fits

from tiptop.baseSimulation import baseSimulation
from tiptop.tiptopUtils import cpuArray


PERF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tiptop', 'perfTest')
BASE_INI = 'HarmoniLTAO_2'  # single-aperture LO, 320 px pupil, D = 39 m

_cache = {}


def _write_ini(tmpdir, name, lo=None, tel=None, remove_lo=()):
    """Write a modified copy of HarmoniLTAO_2.ini.

    Parameters
    ----------
    lo, tel : dict, optional
        Keys to set in [sensor_LO] and [telescope].
    remove_lo : tuple of str
        Keys to remove from [sensor_LO].
    """
    config = ConfigParser()
    config.optionxform = str
    config.read(os.path.join(PERF_DIR, BASE_INI + '.ini'))
    for k, v in (lo or {}).items():
        config.set('sensor_LO', k, str(v))
    for k, v in (tel or {}).items():
        config.set('telescope', k, str(v))
    for k in remove_lo:
        config.remove_option('sensor_LO', k)
    with open(os.path.join(tmpdir, name + '.ini'), 'w') as f:
        config.write(f)


def _make_sim(tmpdir, lo=None, tel=None, run=True):
    """Build a simulation from a temporary .ini variant; run it unless run=False.

    sim.fao (hence sim.fao.ao.pupilLO / pupilFocus) only exists after the run.
    """
    _write_ini(tmpdir, 'sp', lo, tel)
    sim = baseSimulation(tmpdir, 'sp', tmpdir, 'x', doPlot=False, verbose=False)
    if run:
        sim.doOverallSimulation()
    return sim


def _run(key, lo=None, tel=None):
    """Run once per distinct configuration; return host-side NGS results."""
    if key in _cache:
        return _cache[key]
    with tempfile.TemporaryDirectory() as tmpdir:
        sim = _make_sim(tmpdir, lo, tel, run=True)
        out = dict(
            fwhm=np.array(cpuArray(sim.NGS_FWHM_mas_field), dtype=float),
            ee=np.array(cpuArray(sim.NGS_EE_field), dtype=float),
            lo_res=np.array(cpuArray(sim.LO_res), dtype=float),
            dl_fwhm=None if sim.NGS_DL_FWHM_mas is None else np.atleast_1d(np.array(sim.NGS_DL_FWHM_mas, dtype=float)),
        )
    _cache[key] = out
    return out


def _disc_pupil(n, diameter_fraction):
    """Filled disc of diameter diameter_fraction * n pixels, centred on the n x n map."""
    x = (np.arange(n) - n / 2 + 0.5) / n
    return (np.hypot(*np.meshgrid(x, x)) <= diameter_fraction / 2).astype(float)


def _write_fits(path, arr):
    fits.writeto(path, arr, overwrite=True)


class TestDefaultAndEquivalence(unittest.TestCase):

    def test_default_is_telescope_pupil(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sim = _make_sim(tmpdir)
            self.assertIs(sim.fao.ao.pupilLO, sim.fao.ao.tel.pupil)

    def test_lo_pupil_equal_to_telescope_is_bit_identical(self):
        ref = _run('default')
        # PupilAngle explicitly 0 (the telescope's value): equivalent to the default pupil
        same = _run('same', lo={'PathPupil': "'tiptop/data/EELT480pp0.0803m_obs0.283_spider2023.fits'",
                                'PupilAngle': 0.0})
        print(f"\n[equivalence] FWHM {ref['fwhm']} vs {same['fwhm']}, EE {ref['ee']} vs {same['ee']}, "
              f"LO_res {ref['lo_res']} vs {same['lo_res']}")
        for k in ('fwhm', 'ee', 'lo_res'):
            self.assertTrue(np.array_equal(ref[k], same[k]), k)

    def test_different_lo_pupil_changes_output(self):
        # sanity check of the equivalence test: a smaller LO pupil must change the metrics
        ref = _run('default')
        with tempfile.TemporaryDirectory() as tmpdir:
            _write_fits(os.path.join(tmpdir, 'half.fits'), _disc_pupil(320, 0.5))
            other = _run('half', lo={'PathPupil': repr(os.path.join(tmpdir, 'half.fits'))})
        self.assertFalse(np.array_equal(ref['fwhm'], other['fwhm']))


class TestDiffractionScaling(unittest.TestCase):
    """DL FWHM ~ lambda / D_pupil. LO PSF pixel = science PixelScale x 0.85/2.2 = 1.55 mas
    against lambda/D = 4.5 mas (D = 39 m): the core is 2.9 px wide for the full disc and
    5.8 px for the half one, both resolved; the contour FWHM is limited by pixel
    discretisation at the few-% level, hence the 5 % tolerance."""

    def test_half_diameter_doubles_fwhm(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            _write_fits(os.path.join(tmpdir, 'full.fits'), _disc_pupil(320, 1.0))
            _write_fits(os.path.join(tmpdir, 'half.fits'), _disc_pupil(320, 0.5))
            res = {}
            for name in ('full', 'half'):
                res[name] = _run('dl_' + name, lo={'addAliasError': True,
                                                   'PathPupil': repr(os.path.join(tmpdir, name + '.fits'))})
        full, half = res['full']['dl_fwhm'], res['half']['dl_fwhm']
        print(f"\n[DL scaling] FWHM full {full} mas, half {half} mas, ratio {half / full}")
        self.assertEqual(full.shape, half.shape)
        # unobstructed disc: FWHM = 1.03 lambda / D = 4.64 mas at 850 nm, D = 39 m
        self.assertLess(abs(full[0] / (1.029 * 850e-9 / 39 * 180 / np.pi * 3600e3) - 1.0), 0.1)
        np.testing.assert_allclose(half / full, 2.0, rtol=0.05)


class TestFocusFallbackAndErrors(unittest.TestCase):
    # The pupil logic itself (annulus, focus override, error cases) is tested in
    # P3 tests/test_sensor_pupils.py; here only the TIPTOP-side wiring.

    def test_focus_uses_lo_pupil_when_nothing_set(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sim = _make_sim(tmpdir, lo={'ObscurationRatio': 0.3})
            ao = sim.fao.ao
            self.assertIsNot(ao.pupilLO, ao.tel.pupil)
            self.assertIs(ao.pupilFocus, ao.pupilLO)

    def test_config_error_surfaces_from_simulation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sim = _make_sim(tmpdir, run=False, lo={'ObscurationRatio': 0.3,
                                        'PathPupil': "'tiptop/data/EELT480pp0.0803m_obs0.283_spider2023.fits'"})
            with self.assertRaises(ValueError) as cm:
                sim.doOverallSimulation()
            self.assertIn('not both', str(cm.exception))


if __name__ == '__main__':
    unittest.main()
