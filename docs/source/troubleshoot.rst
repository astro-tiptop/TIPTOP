Troubleshooting
===============

We have some custom error messages that will appear if some variables are set the wrong way or in way that not physically possible.
Here is a list of errors that you may encounter.

ValueError: Error : the PSF field of view is too small to simulate the AO correction area
-----------------------------------------------------------------------------------------

This error appears if ``sensor_science.Fieldofview`` is inferior to ``min(source_science.Wavelength-sensor_science.SpectralBandWidth)/DM.DmPitchs/sensor_science.PixelScale*180*3600*10**3/np.pi`` .
If not set ``sensor_science.SpectralBandWidth = 0``.

Either make ``sensor_science.Fieldofview`` bigger or play with ``min(source_science.Wavelength-sensor_science.SpectralBandWidth)`` or ``DM.DmPitchs`` or ``sensor_science.PixelScale`` to make this term smaller.


Warning: No information about the tip-tilt star can be retrieved
-------------------------------------------------------------------
/!\Statement needs verification.
This warning arises if source_HO.Height is defined and different from zero and sources_LO is not defined.
This means that for the AO correction only a laser guide star will be considered.


ValueError: ValueError: '[X1,X2,...]' in section 'Y' must have the same length
------------------------------------------------------------------------------
As mentionned in the parameter file description there are some prameter that should be lists and should have the same length as other parameters.
If you have defined at least one of them, and it contains more than one element, you need to define all the other associated parameter, and it must contain a list of the same length.

For example in the 'atmosphere' section If you set 'Cn2Weights'=[0.8,0.2], then you must set 'Cn2Heights', 'WindSpeed', 'WindDirection' and they ned to contain a list with length 2.

OutOfMemoryError: Out of memory allocating Y bytes (allocated so far: X bytes).
-------------------------------------------------------------------------------
If this error arise it means that due to some parameters create behind the scenes arrays that are too big to handle for your GPU.
A simple check is to work out the size of intermediate planes with the following operation::

    np.ceil(2/(scWvl/(scPs*teldia)*180*3600*10**3*np.pi)*scDet
