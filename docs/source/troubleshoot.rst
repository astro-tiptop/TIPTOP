Troubleshooting
===============

We have some custom error messages that will appear if some variables are set the wrong way or in way that not physically possible.
Here is a list of errors that you may encounter.

ValueError: Error : the PSF field of view is too small to simulate the AO correction area
-----------------------------------------------------------------------------------------

This error appears if ``sensor_science.Fieldofview`` is inferior to ``min(source_science.Wavelength-sensor_science.SpectralBandWidth)/DM.DmPitchs/sensor_science.PixelScale*180*3600*10**3/np.pi`` .
If not set ``sensor_science.SpectralBandWidth = 0``.

Either make ``sensor_science.Fieldofview`` bigger or play with ``min(source_science.Wavelength-sensor_science.SpectralBandWidth)`` or ``DM.DmPitchs`` or ``sensor_science.PixelScale`` to make this term smaller.




