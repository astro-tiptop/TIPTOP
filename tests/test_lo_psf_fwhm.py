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
upsampling factor ceil(3/samp) (samp < 3) that avoids it. The factor 3 rather
than 2 is needed because the core of an obstructed or square-aperture PSF is
narrower than an Airy: with rebin=1 the contour collapses to one pixel up to
samp ~2.3.
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


def _fft_psf(eps: float, samp: float, square: bool = False, npup: int = N,
             ovs: int = 4) -> tuple[np.ndarray, float]:
    """Exact FFT PSF of a circular, annular or square pupil, point-sampled.

    The pupil is zero-padded so that lambda/D = samp*ovs pixels, then decimated
    by ovs: the result is a point-sampled (not pixel-integrated) PSF centred on
    pixel N//2, as in the MASTSEL layout.

    Parameters
    ----------
    eps : float
        Central obstruction ratio (ignored for a square pupil).
    samp : float
        Requested lambda/D in pixels. The padding is rounded to an integer, so
        the value actually realised is returned.
    square : bool
        Square aperture of side D instead of a circular/annular one.
    npup : int
        Pupil size in pixels (must be >= N/samp so the decimated grid fits).
    ovs : int
        Oversampling of the FFT before the decimation.

    Returns
    -------
    psf : np.ndarray
        (N, N) PSF.
    samp_eff : float
        Realised lambda/D in pixels of the returned PSF.
    """
    m = int(round(npup * samp * ovs))
    m += m % 2
    y, x = np.mgrid[-npup // 2:npup // 2, -npup // 2:npup // 2] + 0.5
    if square:
        pup = np.ones((npup, npup))
    else:
        r = np.hypot(x, y) / (npup / 2)
        pup = ((r <= 1) & (r >= eps)).astype(float)
    padded = np.zeros((m, m))
    padded[:npup, :npup] = pup
    full = np.abs(np.fft.fftshift(np.fft.fft2(padded))) ** 2
    idx = m // 2 + ovs * (np.arange(N) - N // 2)
    return full[np.ix_(idx, idx)], m / (npup * ovs)


def _fwhm_in_lod(pixel_scale: float, rebin: int) -> float:
    """sqrt(FWHM_x FWHM_y) of the Airy PSF in units of lambda/D."""
    fx, fy = bs._lo_psf_fwhm(_airy_psf(pixel_scale), pixel_scale, rebin=rebin)
    return float(np.sqrt(fx * fy) / LAMBDA_OVER_D_MAS)


class TestFwhmRebin(unittest.TestCase):

    def test_rebin_values(self):
        lod = LAMBDA_OVER_D_MAS
        self.assertEqual(bs._fwhm_rebin(LAMBDA, lod / 3.5, D), 1)
        self.assertEqual(bs._fwhm_rebin(LAMBDA, lod / 2.5, D), 2)
        self.assertEqual(bs._fwhm_rebin(LAMBDA, 5.49, D), 3)
        self.assertEqual(bs._fwhm_rebin(LAMBDA, lod / 0.8, D), 4)

    def test_rebin_subaperture(self):
        self.assertEqual(bs._fwhm_rebin(LAMBDA, 5.49, D / 8), 1)


class TestLoPsfFwhm(unittest.TestCase):

    def test_airy_undersampled_with_rebin(self):
        ps = 5.49
        rebin = bs._fwhm_rebin(LAMBDA, ps, D)
        self.assertAlmostEqual(LAMBDA_OVER_D_MAS / ps, 1.37, delta=0.01)
        self.assertEqual(rebin, 3)
        fwhm = _fwhm_in_lod(ps, rebin)
        self.assertAlmostEqual(fwhm / AIRY_FWHM, 1.0, delta=0.05)

    def test_p3_fallback_without_rebin_returns_one_pixel(self):
        """Documents the P3 getFWHM fallback (not a requirement of TIPTOP):
        with a single pixel above half maximum it returns exactly 1 pixel."""
        ps = 5.49
        fx, fy = bs._lo_psf_fwhm(_airy_psf(ps), ps, rebin=1)
        self.assertAlmostEqual(float(fx), ps, places=6)
        self.assertAlmostEqual(float(fy), ps, places=6)

    def test_airy_well_sampled_unchanged(self):
        ps = LAMBDA_OVER_D_MAS / 2.5
        rebin = bs._fwhm_rebin(LAMBDA, ps, D)
        self.assertEqual(rebin, 2)
        fwhm = _fwhm_in_lod(ps, rebin)
        self.assertAlmostEqual(fwhm / AIRY_FWHM, 1.0, delta=0.05)

    def test_airy_nyquist_limit_biased_high(self):
        """Known limitation, not a requirement: at samp = 1 the point samples
        are too coarse and the FWHM is biased high whatever the upsampling."""
        ps = LAMBDA_OVER_D_MAS
        rebin = bs._fwhm_rebin(LAMBDA, ps, D)
        ratio = _fwhm_in_lod(ps, rebin) / AIRY_FWHM
        self.assertGreater(ratio, 1.05)
        self.assertLess(ratio, 1.3)


class TestObstructedAndSquare(unittest.TestCase):
    PUPILS = [('annular eps=0.28', 0.28, False),
              ('annular eps=0.5', 0.5, False),
              ('square', 0.0, True)]

    @staticmethod
    def _sqrt_fwhm(psf: np.ndarray, ps: float, rebin: int) -> float:
        fx, fy = bs._lo_psf_fwhm(psf, ps, rebin=rebin)
        return float(np.sqrt(fx * fy))

    def _reference(self, eps: float, square: bool) -> float:
        """FWHM (mas) from the same PSF sampled at 8 pixels per lambda/D."""
        psf, samp = _fft_psf(eps, 8.0, square)
        return self._sqrt_fwhm(psf, LAMBDA_OVER_D_MAS / samp, 1)

    def test_square_reference_is_analytic(self):
        ref = self._reference(0.0, True) / LAMBDA_OVER_D_MAS
        self.assertAlmostEqual(ref / 0.886, 1.0, delta=0.02)  # +0.9 % point-sampling bias of the reference

    def test_fwhm_with_rebin(self):
        for name, eps, square in self.PUPILS:
            ref = self._reference(eps, square)
            for samp in (1.5, 2.05, 2.2):
                with self.subTest(pupil=name, samp=samp):
                    psf, s_eff = _fft_psf(eps, samp, square)
                    ps = LAMBDA_OVER_D_MAS / s_eff
                    rebin = bs._fwhm_rebin(LAMBDA, ps, D)
                    ratio = self._sqrt_fwhm(psf, ps, rebin) / ref
                    self.assertAlmostEqual(ratio, 1.0, delta=0.03)

    def test_without_rebin_collapses_at_samp_2(self):
        """Why samp < 3 needs upsampling: rebin=1 collapses for some pupils."""
        ratios = []
        for name, eps, square in self.PUPILS:
            psf, s_eff = _fft_psf(eps, 2.05, square)
            ps = LAMBDA_OVER_D_MAS / s_eff
            ratios.append(self._sqrt_fwhm(psf, ps, 1) / self._reference(eps, square))
        self.assertLess(min(ratios), 0.6)


if __name__ == '__main__':
    unittest.main()
