import unittest

import numpy as np
from scipy.interpolate import interp1d
from scipy.special import erf

from p3.aoSystem.FourierUtils import getEncircledEnergy
from tiptop import baseSimulation as bs
from mastsel.mavisPsf import centeredPixelCoords

"""
Encircled energy of a pixel-integrated circular Gaussian (sigma = 32 mas) through
the computeMetrics pipeline: P3 getEncircledEnergy -> _ee_equivalent_area_radius
-> cubic interp1d. Analytic result: EE(R) = 1 - exp(-R^2 / (2 sigma^2)).
"""

N = 256
SIGMA_MAS = 32.0


def _gaussian_image(pixel_scale: float) -> np.ndarray:
    """Pixel-integrated circular Gaussian centred on pixel N//2 (host array, unit sum)."""
    x = (np.arange(N) - centeredPixelCoords(N)[0]) * pixel_scale
    s = np.sqrt(2) * SIGMA_MAS
    g = 0.5 * (erf((x + pixel_scale / 2) / s) - erf((x - pixel_scale / 2) / s))
    return np.outer(g, g)


def _ee_pipeline(pixel_scale: float, radii_mas: np.ndarray) -> np.ndarray:
    """EE at radii_mas as computed in computeMetrics."""
    img = _gaussian_image(pixel_scale)
    ee, _ = getEncircledEnergy(img, pixelscale=pixel_scale, center=centeredPixelCoords(N), nargout=2)
    e, r = bs._ee_equivalent_area_radius(ee, N, pixel_scale)
    return interp1d(r, e, kind='cubic', bounds_error=False)(radii_mas)


class TestEncircledEnergy(unittest.TestCase):

    def _max_rel_err(self, pixel_scale):
        radii = np.arange(1, 11) * 0.5 * pixel_scale  # 0.5 .. 5 px
        ee = _ee_pipeline(pixel_scale, radii)
        ref = 1 - np.exp(-radii**2 / (2 * SIGMA_MAS**2))
        return np.max(np.abs(ee / ref - 1))

    def test_gaussian_ee_10mas(self):
        self.assertLess(self._max_rel_err(10.0), 0.006)

    def test_gaussian_ee_25mas(self):
        self.assertLess(self._max_rel_err(25.0), 0.015)

    def test_equivalent_radii_match_pixel_counts(self):
        ps = 10.0
        _, r = bs._ee_equivalent_area_radius(np.zeros(4), N, ps)
        n_k = np.rint(r[1:]**2 * np.pi / ps**2).astype(int)
        np.testing.assert_allclose(r[1:]**2 * np.pi / ps**2, n_k, atol=1e-9)
        x = np.arange(N) - centeredPixelCoords(N)[0]
        rr = np.round(np.hypot(x[:, None], x[None, :]))
        expected = [int((rr <= k).sum()) for k in range(4)]
        self.assertEqual(expected, [1, 9, 21, 37])
        self.assertEqual(list(n_k), expected)
        self.assertEqual(r[0], 0.0)


if __name__ == '__main__':
    unittest.main()
