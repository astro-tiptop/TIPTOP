"""
[sensor_science] PerfectCoronagraph: when True the output cube holds the image behind
a perfect coronagraph, in units of the stellar flux in the field of view. The SR, FWHM
and EE stay those of the AO PSF. Expected physics: the total flux is (1 - SR), the image
equals the PSF far from the axis (the removed coherent part is SR x Airy, small
there) and is much fainter than the PSF peak near the axis.
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
BASE_INI = 'ERIS'  # SCAO, no [sensor_LO], 1.65 um on an 8 m telescope
PIXEL_SCALE = 21.0  # mas: Nyquist sampling at 1.65 um (lambda/D = 42.5 mas)
FOV = 128
SEEING = 2.0
PHOTONS = '[10.0]'

_cache = {}


def _write_ini(tmpdir, name, coro='absent', base=BASE_INI, edit=None):
    """Write a modified copy of a perfTest ini; coro='absent' leaves the key out.

    edit is an optional callable that receives the ConfigParser before writing."""
    config = ConfigParser()
    config.optionxform = str
    config.read(os.path.join(PERF_DIR, base + '.ini'))
    config.set('sensor_science', 'PixelScale', str(PIXEL_SCALE))
    config.set('sensor_science', 'FieldOfView', str(FOV))
    # poor seeing and a faint guide star (SR ~ 0.56): at SR ~ 0.9 the Airy wings are ~10 %
    # of the PSF at 20-30 lambda/D and the point-sampled SR is off by 0.01 (15 % of 1 - SR)
    config.set('atmosphere', 'Seeing', str(SEEING))
    config.set('sensor_HO', 'NumberPhotons', str(PHOTONS))
    if coro != 'absent':
        config.set('sensor_science', 'PerfectCoronagraph', str(coro))
    if edit is not None:
        edit(config)
    with open(os.path.join(tmpdir, name + '.ini'), 'w') as f:
        config.write(f)


def _run(coro='absent'):
    """Run a simulation with metrics, PSF1D and a saved FITS file; return host results."""
    if coro in _cache:
        return _cache[coro]
    with tempfile.TemporaryDirectory() as tmpdir:
        _write_ini(tmpdir, 'pc', coro)
        sim = baseSimulation(tmpdir, 'pc', tmpdir, 'out', doConvolve=True, verbose=False)
        sim.doOverallSimulation()
        sim.computeMetrics()
        sim.computePSF1D()
        sim.saveResults()
        with fits.open(os.path.join(tmpdir, 'out.fits')) as hdul:
            hdr = hdul[1].header.copy()
        out = dict(
            cube=np.array(cpuArray(sim.cubeResultsArray), dtype=float).squeeze(),
            sr=float(cpuArray(sim.sr[0])),
            fwhm=float(cpuArray(sim.fwhm[0])),
            ee=np.asarray(cpuArray(sim.ee), dtype=float),
            psf1d=np.asarray(sim.psf1d, dtype=float).ravel(),
            samp=sim.wvlRef * 3600 * 180 * 1000 / np.pi / (sim.psInMas * 2 * sim.tel_radius),
            hdr=hdr,
        )
    _cache[coro] = out
    return out


def _radial_mean(img, r_min, r_max, samp):
    """Mean of img over the annulus r_min..r_max (lambda/D) around the grid centre."""
    n = img.shape[0]
    y, x = np.indices(img.shape)
    r = np.hypot(x - n // 2, y - n // 2) / samp
    return img[(r >= r_min) & (r < r_max)].mean()


class TestPerfectCoronagraphInput(unittest.TestCase):

    def _expect_error(self, match, base=BASE_INI, coro=True, edit=None):
        with tempfile.TemporaryDirectory() as tmpdir:
            _write_ini(tmpdir, 'bad', coro, base=base, edit=edit)
            with self.assertRaises(ValueError) as cm:
                # the opdMap is read by P3 during the simulation setup, not at construction
                baseSimulation(tmpdir, 'bad', tmpdir, 'out', verbose=False).doOverallSimulation()
            self.assertIn(match, str(cm.exception))

    def test_non_bool_raises(self):
        self._expect_error('PerfectCoronagraph', coro='true')

    def test_lo_loop_raises(self):
        self._expect_error('LO', base='ERIS_LGS')

    def test_jitter_raises(self):
        self._expect_error('jitter', edit=lambda c: c.set('telescope', 'jitter_FWHM', '[5.0, 5.0, 0.0]'))

    def test_static_aberrations_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, 'opd.fits')
            fits.writeto(path, np.zeros((128, 128)))
            self._expect_error('opdMap', edit=lambda c: c.set('telescope', 'PathStaticOn', repr(path)))


class TestPerfectCoronagraphPhysics(unittest.TestCase):

    def test_metrics_unchanged(self):
        n, c = _run(), _run(True)
        print(f"\n[metrics] SR {n['sr']} vs {c['sr']}, FWHM {n['fwhm']} vs {c['fwhm']}")
        self.assertEqual(n['sr'], c['sr'])
        self.assertEqual(n['fwhm'], c['fwhm'])
        np.testing.assert_array_equal(n['ee'], c['ee'])

    def test_flux_is_one_minus_sr(self):
        c = _run(True)
        total = c['cube'].sum()
        print(f"\n[flux] sum {total:.4f}, 1-SR {1 - c['sr']:.4f}, ratio {total / (1 - c['sr']):.4f}")
        self.assertLess(abs(total / (1 - c['sr']) - 1.0), 0.05)
        self.assertGreater(c['cube'].min(), -1e-9)

    def test_far_and_near_axis(self):
        n, c = _run(), _run(True)
        samp = c['samp']
        far = _radial_mean(c['cube'], 20, 30, samp) / _radial_mean(n['cube'], 20, 30, samp)
        near = _radial_mean(c['cube'], 0, 1, samp) / n['cube'].max()
        print(f"\n[axis] sampling {samp:.2f}, SR {c['sr']:.3f}: far ratio {far:.4f}, "
              f"near / PSF peak {near:.4f}")
        self.assertLess(abs(far - 1.0), 0.05)
        self.assertLess(near, 0.2)

    def test_psf1d_centred_and_not_renormalized(self):
        n, c = _run(), _run(True)
        centre = c['cube'][FOV // 2, FOV // 2]
        print(f"\n[psf1d] first sample {c['psf1d'][0]:.3e}, cube centre {centre:.3e}, "
              f"max {c['psf1d'].max():.3e}; PSF peak {n['cube'].max():.3e}, psf1d[0] {n['psf1d'][0]:.3e}")
        # a profile renormalized to unit sum would be 1 / (1 - SR) ~ 2.3 times higher
        self.assertLess(abs(c['psf1d'][0] / centre - 1.0), 0.1)
        self.assertLess(c['psf1d'][0], 0.01 * n['psf1d'][0])  # not the PSF peak

    def test_explicit_zero_jitter_matches_no_jitter(self):
        # an explicit zero jitter still convolves (and unit-sum normalizes) the PSF; the coronagraphic
        # cube must not come out in raw MASTSEL units
        ref = _run(True)['cube']
        for value in ('0', '[0, 0, 30]'):
            key = 'jit' + value
            if key not in _cache:
                with tempfile.TemporaryDirectory() as tmpdir:
                    _write_ini(tmpdir, 'pc', True, edit=lambda c: c.set('telescope', 'jitter_FWHM', value))
                    sim = baseSimulation(tmpdir, 'pc', tmpdir, 'out', doConvolve=True, verbose=False)
                    sim.doOverallSimulation()
                    _cache[key] = np.array(cpuArray(sim.cubeResultsArray), dtype=float).squeeze()
            np.testing.assert_allclose(_cache[key], ref, rtol=1e-9, atol=1e-12 * np.abs(ref).max(),
                                       err_msg=f"jitter_FWHM = {value}")

    def test_fits_header_flag(self):
        self.assertIs(_run(True)['hdr']['CORONAGR'], True)
        self.assertNotIn('CORONAGR', _run()['hdr'])


if __name__ == '__main__':
    unittest.main()
