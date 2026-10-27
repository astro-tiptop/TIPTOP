"""
[sensor_science] PixelResponse: when True the science, open-loop and
diffraction-limited PSFs are integrated over the square detector pixel
(OTF x sinc). The Strehl reference is then pixel-integrated too, so SR must
stay close to the point-sampled value, while FWHM grows and the DL peak drops.
"""

import os
import tempfile
import unittest
from configparser import ConfigParser

import numpy as np
from astropy.io import fits

from scipy.integrate import dblquad
from scipy.special import j1

from tiptop.baseSimulation import baseSimulation, _diffraction_limited_psf
from tiptop.tiptopUtils import cpuArray


RAD2MAS = 3600 * 180 * 1000 / np.pi
PERF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tiptop', 'perfTest')
BASE_INI = 'dummy'  # smallest perfTest configuration (128 px pupil, 1 science direction)

# science PixelScale [mas] -> sampling at 1.65 um of about 2.1, 1.06 and 0.53
PIXEL_SCALES = (20.0, 40.0, 80.0)

_cache = {}


def _write_ini(tmpdir, name, pixel_scale, pixel_response='absent', wavelength='[1.65e-06]'):
    """Write a modified copy of dummy.ini; pixel_response='absent' leaves the key out."""
    config = ConfigParser()
    config.optionxform = str
    config.read(os.path.join(PERF_DIR, BASE_INI + '.ini'))
    # bright HO guide stars, no noise, good seeing: high Strehl (~0.6-0.7), so that SR is
    # not dominated by the 5-decimal rounding of the metric
    config.set('sources_science', 'Wavelength', wavelength)
    config.set('atmosphere', 'Seeing', '0.6')
    config.set('sensor_HO', 'NumberPhotons', '[5000, 5000]')
    config.set('sensor_HO', 'SigmaRON', '0.0')
    config.set('sensor_HO', 'SpotFWHM', '[[0.0, 0.0, 0.0]]')
    config.set('sensor_science', 'PixelScale', str(pixel_scale))
    if pixel_response != 'absent':
        config.set('sensor_science', 'PixelResponse', str(pixel_response))
    with open(os.path.join(tmpdir, name + '.ini'), 'w') as f:
        config.write(f)


def _run(pixel_scale, pixel_response='absent', save=False):
    """Run a simulation with computeMetrics; return a dict of host results."""
    key = (pixel_scale, pixel_response, save)
    if key in _cache:
        return _cache[key]
    with tempfile.TemporaryDirectory() as tmpdir:
        _write_ini(tmpdir, 'pr', pixel_scale, pixel_response)
        sim = baseSimulation(tmpdir, 'pr', tmpdir, 'out', doConvolve=True, verbose=False)
        sim.doOverallSimulation()
        psfs = [cpuArray(img.sampling).copy() for img in sim.results]
        sim.computeMetrics()
        out = dict(
            pixelResponse=sim.pixelResponse,
            psfs=psfs,
            sr=[float(cpuArray(x)) for x in sim.sr],
            fwhm=[float(cpuArray(x)) for x in sim.fwhm],
            psf_dl=np.array(sim.psf_dl_array, dtype=float),
            samp=sim.wvlRef * RAD2MAS / (sim.psInMas * 2 * sim.tel_radius),
            nPix=sim.nPixPSF,
            hdr=None,
        )
        if save:
            sim.saveResults()
            with fits.open(os.path.join(tmpdir, 'out.fits')) as hdul:
                out['hdr'] = [hdul[i].header.copy() for i in range(4)]
    _cache[key] = out
    return out


class TestPixelResponseInput(unittest.TestCase):

    def test_default_is_false(self):
        self.assertIs(_run(PIXEL_SCALES[0])['pixelResponse'], False)

    def test_explicit_values(self):
        self.assertIs(_run(PIXEL_SCALES[0], False)['pixelResponse'], False)
        self.assertIs(_run(PIXEL_SCALES[0], True)['pixelResponse'], True)

    def test_string_value_raises(self):
        # 'false' (lower case) is not parsed as a literal by the .ini reader and
        # stays a (truthy) string: it must be rejected, not silently enabled.
        with tempfile.TemporaryDirectory() as tmpdir:
            _write_ini(tmpdir, 'bad', PIXEL_SCALES[0], 'false')
            with self.assertRaises(ValueError) as cm:
                baseSimulation(tmpdir, 'bad', tmpdir, 'out', verbose=False)
            self.assertIn('PixelResponse', str(cm.exception))


class TestPixelResponsePhysics(unittest.TestCase):

    def test_pixel_integration_effects(self):
        for ps in PIXEL_SCALES:
            with self.subTest(pixel_scale=ps):
                f = _run(ps, False)
                t = _run(ps, True)
                samp = f['samp']
                print(f"\n[PixelScale={ps} mas, sampling={samp:.3f}] "
                      f"SR F={f['sr']} T={t['sr']} ratio={t['sr'][0] / f['sr'][0]:.4f}; "
                      f"FWHM F={f['fwhm']} T={t['fwhm']}; "
                      f"DL peak ratio={t['psf_dl'].max() / f['psf_dl'].max():.4f}")

                # SR: the reference is pixel-integrated too. Below sampling 1 the
                # point-sampled SR is biased high: the sum of point samples, used for the
                # unit-sum normalization, is (1/p^2) sum_n OTF(n/p), and 1/p < D/lambda puts
                # aliases inside the OTF support, larger for the DL reference (A = 1.98 at
                # sampling 0.53) than for the AO PSF (~SR*A): bias (1+A)/(1+SR*A) = 1.33,
                # measured 0.83/0.63. The pixel sinc vanishes at n/p, so the pixel-integrated
                # SR has no such bias: there it is checked against the sampling-1 one.
                # The pixel-integrated SR still depends on PixelScale through the halo and
                # the finite-FoV normalization (0.644 at sampling 2.1, 0.626 at 1.06), so
                # below sampling 1 it is only checked for consistency, at the 5 % level.
                if samp >= 1.0:
                    for sr_t, sr_f in zip(t['sr'], f['sr']):
                        self.assertLess(abs(sr_t / sr_f - 1.0), 0.01)  # measured 0.2-0.4 %
                else:
                    sr_ref = _run(PIXEL_SCALES[1], True)['sr']
                    for sr_t, sr_f, sr_r in zip(t['sr'], f['sr'], sr_ref):
                        self.assertLess(sr_t, sr_f)  # aliasing bias of the point-sampled SR removed
                        self.assertLess(abs(sr_t / sr_r - 1.0), 0.05)

                # FWHM includes the pixel
                self.assertEqual(len(t['fwhm']), len(f['fwhm']))
                for fw_t, fw_f in zip(t['fwhm'], f['fwhm']):
                    self.assertGreater(fw_t, fw_f)

                # PSFs: unit sum and on-axis peak at the centre
                for res in (f, t):
                    for psf in res['psfs']:
                        self.assertAlmostEqual(float(psf.sum()), 1.0, delta=1e-3)
                        c = res['nPix'] // 2
                        self.assertEqual(np.unravel_index(np.argmax(psf), psf.shape), (c, c))

                # DL PSF peak is lowered by the pixel, by an amount set by the sampling
                ratio = t['psf_dl'].max() / f['psf_dl'].max()
                self.assertLess(ratio, 1.0)
                self.assertGreater(ratio, 0.6)
                self.assertLess(ratio, 0.95)

    def test_dl_peak_drop_is_largest_near_sampling_1(self):
        # pixel-integrated / point-sampled DL peak (unit sum): 0.91 at sampling 2.1,
        # 0.70 at 1.06, 0.88 at 0.53 (the point-sampled sum is aliased below 1)
        r = [_run(ps, True)['psf_dl'].max() / _run(ps, False)['psf_dl'].max() for ps in PIXEL_SCALES]
        self.assertLess(r[1], r[0])
        self.assertLess(r[1], r[2])

    def test_fits_header_flag(self):
        for flag in (False, True):
            with self.subTest(flag=flag):
                hdrs = _run(PIXEL_SCALES[0], flag, save=True)['hdr']
                # written only when set (default files unchanged), in the PSF-cube,
                # open-loop and diffraction-limited extensions
                for h in hdrs[1:]:
                    if flag:
                        self.assertIs(h['PIX_RESP'], True)
                    else:
                        self.assertNotIn('PIX_RESP', h)


class TestPerWavelengthPath(unittest.TestCase):
    """Two wavelengths with exactMultiWavelengthPSD=True: per-wavelength PSD grids
    (multiGridPSD), the psdSetToPsfSet call that the single-wavelength runs do not use."""

    @staticmethod
    def _run_multi(pixel_response):
        with tempfile.TemporaryDirectory() as tmpdir:
            _write_ini(tmpdir, 'mw', 40.0, pixel_response, wavelength='[1.2e-06, 1.65e-06]')
            sim = baseSimulation(tmpdir, 'mw', tmpdir, 'out', doConvolve=True, verbose=False,
                                 exactMultiWavelengthPSD=True)
            sim.doOverallSimulation()
            sim.computeMetrics()
            return sim.multiGridPSD, np.array(sim.sr, dtype=float), np.array(sim.fwhm, dtype=float)

    def test_per_wavelength_pixel_response(self):
        grid_f, sr_f, fw_f = self._run_multi(False)
        grid_t, sr_t, fw_t = self._run_multi(True)
        print(f"\n[per-wavelength, 40 mas] SR F={sr_f.ravel()} T={sr_t.ravel()}; "
              f"FWHM F={fw_f.ravel()} T={fw_t.ravel()}")
        self.assertTrue(grid_f and grid_t)
        # sampling 0.78 (1.2 um) and 1.06 (1.65 um)
        self.assertTrue(np.all(fw_t > fw_f))
        self.assertLess(abs(sr_t[1] / sr_f[1] - 1.0).max(), 0.01)  # sampling >= 1
        self.assertTrue(np.all(sr_t[0] < sr_f[0]))  # sampling < 1: aliasing bias removed


class TestDiffractionLimitedReference(unittest.TestCase):
    """Pixel-integrated DL reference against the analytic result: unit-sum peak =
    fraction of the Airy flux in a square pixel of side a = 1/samp (lambda/D units),
    I(r) = (pi/4) (2 J1(pi r) / (pi r))^2 (unit total flux), unobstructed pupil."""

    @staticmethod
    def _airy_pixel_fraction(a):
        def airy(y, x):
            r = np.hypot(x, y)
            return np.pi / 4 if r < 1e-12 else np.pi / 4 * (2 * j1(np.pi * r) / (np.pi * r)) ** 2
        # one octant of the square pixel, times 8
        val, _ = dblquad(lambda y, x: airy(y, x), 0, a / 2, 0, lambda x: x, epsabs=1e-10)
        return 8 * val

    def test_pixel_integrated_airy(self):
        # pixel centres: diameter = array side = D (linspace(-0.5, 0.5, n) gives (n-1)/n D,
        # i.e. a peak ((n-1)/n)^2 too low)
        n = 128
        x = (np.arange(n) - n / 2 + 0.5) / n
        pupil = (np.hypot(*np.meshgrid(x, x)) <= 0.5).astype(float)
        for samp in (2.1, 1.06, 0.53):
            with self.subTest(samp=samp):
                peak = float(cpuArray(_diffraction_limited_psf(pupil, samp, True)).max())
                exact = self._airy_pixel_fraction(1.0 / samp)
                print(f"\n[DL reference, sampling {samp}] peak {peak:.5f}, analytic {exact:.5f}")
                # pupil discretization and enlargeSupport rounding: a few 1e-3
                self.assertLess(abs(peak / exact - 1.0), 0.01)


if __name__ == '__main__':
    unittest.main()
