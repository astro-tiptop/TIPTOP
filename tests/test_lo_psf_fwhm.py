import unittest

import numpy as np
from scipy.special import j1

from tiptop import baseSimulation as bs
from mastsel.mavisPsf import centeredPixelCoords

"""
FWHM of the NGS/Focus PSFs (input of the MASTSEL LO noise model) through
_lo_psf_fwhm -> P3 getFWHM(method='contour', rebin=...).
Physical expectation: the Airy PSF of an unobstructed circular pupil,
I(r) = [2 J1(x)/x]^2 with x = pi r D / lambda, has FWHM = 1.029 lambda/D.
With the science pixel scaled in wavelength the PSF can be below Nyquist
(lambda/D = 1.37 pixel): without rebin the contour has a single pixel above
half maximum and getFWHM collapses to exactly 1 pixel. _fwhm_rebin picks the
upsampling factor ceil(2/samp) that avoids it.
"""

N = 128
D = 38.5
LAMBDA_OVER_D_MAS = 7.5  # ANDES-like: pixel 5.49 mas gives samp = 1.37
LAMBDA = LAMBDA_OVER_D_MAS / bs.rad2mas * D
AIRY_FWHM = 1.029  # in lambda/D


def _airy_psf(pixel_scale: float) -> np.ndarray:
    """Point-sampled Airy PSF (host array) centred on pixel N//2, in the MASTSEL layout."""
    c = centeredPixelCoords(N)[0]
    x = (np.arange(N) - c) * pixel_scale
    r = np.hypot(x[:, None], x[None, :])
    u = np.pi * r / LAMBDA_OVER_D_MAS
    u[u == 0] = 1e-12
    return (2 * j1(u) / u) ** 2


def _fwhm_in_lod(pixel_scale: float, rebin: int) -> float:
    """sqrt(FWHM_x FWHM_y) of the Airy PSF in units of lambda/D."""
    fx, fy = bs._lo_psf_fwhm(_airy_psf(pixel_scale), pixel_scale, rebin=rebin)
    return float(np.sqrt(fx * fy) / LAMBDA_OVER_D_MAS)


class TestFwhmRebin(unittest.TestCase):

    def test_rebin_values(self):
        lod = LAMBDA_OVER_D_MAS
        self.assertEqual(bs._fwhm_rebin(LAMBDA, lod / 2.5, D), 1)
        self.assertEqual(bs._fwhm_rebin(LAMBDA, lod / 2.2, D), 1)
        self.assertEqual(bs._fwhm_rebin(LAMBDA, 5.49, D), 2)
        self.assertEqual(bs._fwhm_rebin(LAMBDA, lod / 0.8, D), 3)

    def test_rebin_subaperture(self):
        self.assertEqual(bs._fwhm_rebin(LAMBDA, 5.49, D / 8), 1)


class TestLoPsfFwhm(unittest.TestCase):

    def test_airy_undersampled_with_rebin(self):
        ps = 5.49
        rebin = bs._fwhm_rebin(LAMBDA, ps, D)
        self.assertAlmostEqual(LAMBDA_OVER_D_MAS / ps, 1.37, delta=0.01)
        self.assertEqual(rebin, 2)
        fwhm = _fwhm_in_lod(ps, rebin)
        self.assertAlmostEqual(fwhm / AIRY_FWHM, 1.0, delta=0.05)

    def test_airy_undersampled_without_rebin_collapses_to_one_pixel(self):
        ps = 5.49
        fx, fy = bs._lo_psf_fwhm(_airy_psf(ps), ps, rebin=1)
        self.assertAlmostEqual(float(fx), ps, places=6)
        self.assertAlmostEqual(float(fy), ps, places=6)

    def test_airy_well_sampled_unchanged(self):
        ps = LAMBDA_OVER_D_MAS / 2.5
        rebin = bs._fwhm_rebin(LAMBDA, ps, D)
        self.assertEqual(rebin, 1)
        fwhm = _fwhm_in_lod(ps, rebin)
        self.assertAlmostEqual(fwhm / AIRY_FWHM, 1.0, delta=0.05)


if __name__ == '__main__':
    unittest.main()
