Parameter files explained
=========================

Introduction
------------

To run TIPTOP you need two things: to execute the function and to create a parameter file. This section explaines
what should your parameter file contain and the various parameters you can set. You can already find example parameter 
files in the `github <https://github.com/FabioRossiArcetri/TIPTOP/tree/main/perfTest>`_ .


The parameter files are divided in sections and they can contain multiple parameter. It is very important that each 
parameter be placed in the appropriate section. The section starts with its name between '[]' and ends either with 
the end of file or with the next section. The order they are placed in the file does not matter.

Overview
--------
The following table resume the what the parameter file should contain and what is mandatory.

+---------------+--------------------------------+
|section        | required?                      |
+===============+================================+
|telescope      | Yes                            |
+---------------+--------------------------------+
|atmosphere     | Yes                            |
+---------------+--------------------------------+
|sources_science| Yes                            |
+---------------+--------------------------------+
|source_HO      | Yes                            |
+---------------+--------------------------------+
|sources_LO     | No/Yes if sensor_LO defined    |
+---------------+--------------------------------+
|sources_Focus  | No/Yes if sensor_Focus defined |
+---------------+--------------------------------+
|sensor_science | Yes                            |
+---------------+--------------------------------+
|sensor_HO      | Yes                            |
+---------------+--------------------------------+
|sensor_LO      | No/Yes if sources_LO defined   |
+---------------+--------------------------------+
|sensor_Focus   | No/Yes if sources_Focus defined|
+---------------+--------------------------------+
|DM             | Yes                            |
+---------------+--------------------------------+
|RTC            | No/Yes if sensor_LO defined    |
+---------------+--------------------------------+
|Computation    | No                             |
+---------------+--------------------------------+


We now go more in detail for each section:

.. warning::

   In the current state, if the parameter file is not found, it raises ``NoSectionError: No section: 'telescope'``.

[telescope]
-----------

+--------------------------+----------+-------+--------------------------------------------------------------------------+
| Parameter                | Required | Type  | Description                                                              |
+==========================+==========+=======+==========================================================================+
|TelescopeDiameter         |Yes       |float  |Set the outer diameter of the telescope pupil in unit of meters.          |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|Resolution                |Yes       |integer|Number of pixels across the pupil diameter                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|ObscurationRatio          |No/Yes if |float  |*Default : 0.0*, Defines the central obstruction                          |
|                          |LO        |       |due to the secondary as a ratio of the TelescopeDiameter                  |
|                          |          |       |                                                                          |
|                          |          |       |*Warning* : MavisLO.py does not have a default value for this parameter   |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|ZenithAngle               |No/Yes if |float  |*Default: 0.0*, Set the pointing direction of the telescope in degree     |
|                          |LO        |       |with respect to the zenith. Used to compute airmass, to scale atmospheric |
|                          |          |       |layers and stars altitude.                                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PupilAngle                |No        |float  |*default : 0.0*, unknown effect                                           |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathPupil                 |No        |string |*default: ''*, path to the pupil model in .fits file (if provided,        |
|                          |          |       |the pupil model is interpolated). if absent or '', not used.              |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathStaticOn              |No        |string |*default: None*, path to a map of static aberrations (nm) in              |
|                          |          |       |.fits file. if absent or '', not used.                                    |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathStaticOff             |No        |string |*default: None*, No clue what this does. if absent or '', not used.       |
|                          |          |       |From P3, not supported in TIPTOP.                                         |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathStaticPos             |No        |string |*default: None*, No clue                                                  |
|                          |          |       |From P3, not supported in TIPTOP.                                         |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathApodizer              |No        |string |*default: ''*, Path to a fits file that contain a binary map corresponding|
|                          |          |       |to a pupil apodizer (TBC). if absent or '', not used.                     |
|                          |          |       |From P3, not supported in TIPTOP.                                         |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathStatModes             |No        |string |*default: ''*, path to a .fits file that contain a cube of map of mode    |
|                          |          |       |in amplitude which lead to a rms of 1 in nanometer of static aberation.   |
|                          |          |       |if absent or '', not used. Unsure how this works.                         |
|                          |          |       |From P3, not supported in TIPTOP.                                         |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|coefficientOfTheStaticMode|not used  |string |*default: ''*, place holder                                               |
|                          |          |       |(TBC) need to find how does the pathStatModes fits file work.             |
|                          |          |       |From P3, not supported in TIPTOP.                                         |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorNm              |No        |float  |*default: 0*, nm RMS of the additional error to be added (an error that   |
|                          |          |       |is not otherwise considered)                                              |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorExp             |No        |float  |*default: -2*, exponent of the power of spatial frequencies used to       |
|                          |          |       |generate the PSD associated with extraErrorNm                             |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorMin             |No        |float  |*default: 0*, minimum spatial frequency for which PSD associated with     |
|                          |          |       |extraErrorNm is > 0                                                       |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorMax             |No        |float  |*default: 0*, maximum spatial frequency for which PSD associated with     |
|                          |          |       |extraErrorNm is > 0                                                       |
|                          |          |       |                                                                          |
|                          |          |       |Note: 0 means maximum frequency is the one present in the spatial         |
|                          |          |       |frequency array of the PSDs.                                              |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorLoNm            |No        |float  |*default: 0*, nm RMS of the additional error to be added (an error that   |
|                          |          |       |is not otherwise considered)                                              |
|                          |          |       |                                                                          |
|                          |          |       |Note: (1) only makes sense if [sensor_LO] is present (2) if not present   |
|                          |          |       |extraErrorNm is used on LO directions.                                    |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorLoExp           |No        |float  |*default: -2*, exponent of the power of spatial frequencies used to       |
|                          |          |       |generate the PSD associated with extraErrorLoNm                           |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorLoMin           |No        |float  |*default: 0*, minimum spatial frequency for which PSD associated with     |
|                          |          |       |extraErrorLoNm is > 0                                                     |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorLoMax           |No        |float  |*default: 0*, maximum spatial frequency for which PSD associated with     |
|                          |          |       |extraErrorLoNm is > 0                                                     |
|                          |          |       |                                                                          |
|                          |          |       |Note: 0 means maximum frequency is the one present in the spatial         |
|                          |          |       |frequency array of the PSDs.                                              |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|jitter_FWHM               |No        |float  |*default: None*, additional kernel to be convolved with PSF, it could be  |
|                          |          |       |a scalar (FWHM in mas) for a round kernel or a list of three values       |
|                          |          |       |[FWHM_mas_max, FWHM_mas_min, angle_rad].                                  |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|glFocusOnNGS              |No        |string |*default: False*, global focus control with natural guide stars.          |
|                          |          |       |Multi-conjugate systems only. Requires NumberLenslets >= 2 in sensor_LO or|
|                          |          |       |a specific global focus sensor (``[sources_Focus]`` and ``[sensor_Focus]``|
|                          |          |       |sections).                                                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|TechnicalFoV              |No/Yes if |float  |*default: ??*, set the size of the technical field of view (diameter) is  |
|                          |LO        |       |Used in laser and multi-conjugate AO systems.                             |
|                          |          |       |                                                                          |
|                          |          |       |*Warning* : This is not optional in MavisLO.py                            |
+--------------------------+----------+-------+--------------------------------------------------------------------------+


[atmosphere]
------------

+-------------------------+---------+-------+--------------------------------------------------------------------------+
| Parameter               | Required| Type  | Description                                                              |
+=========================+=========+=======+==========================================================================+
|Seeing                   |Yes      |float  |Set the seeing at Zenith in arcsec. If not set TIPTOP uses ``r0_value`` . |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|Wavelength               |No/Yes if|float  |*Default : 500e-9*, Wavelength of definition of the atmosphere statistics |
|                         |LO       |       |                                                                          |
|                         |         |       |*Warning*: not optional in MavisLO.py                                     |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|L0                       |No/Yes if|float  |*Default : 25.0*, Outer Scale of the atmosphere  in meters                |
|                         |LO       |       |                                                                          |
|                         |         |       |*Warning*: not optional in MavisLO.py                                     |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|Cn2Weights               |No/Yes   |list of|*Default : [1.0]*, Relative contribution of each layer. The sum of all the|
|                         |if LO    |float  |list element must be 1. Must have the same length as ``Cn2Heights``,      |
|                         |         |       |``WindSpeed`` and ``WindDirection``.                                      |
|                         |         |       |                                                                          |
|                         |         |       |*Warning : required if ``Cn2Heights``, ``WindSpeed`` or ``WindDirection`` |
|                         |         |       |are defined                                                               |
|                         |         |       |*Warning* : extremely confusing error message if absent when it must be   |
|                         |         |       |defined                                                                   |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|Cn2Heights               |No/Yes   |list of|*Default : [0.0]*, altitude of layers in meters.                          |
|                         |if LO    |float  |Must have the same length as ``Cn2Weights``, ``WindSpeed`` and            |
|                         |         |       |``WindDirection``.                                                        |
|                         |         |       |                                                                          |
|                         |         |       |*Warning* : required if ``Cn2Weights``, ``WindSpeed`` or ``WindDirection``|
|                         |         |       |are defined                                                               |
|                         |         |       |*Warning* : extremely confusing error message if absent when it must be   |
|                         |         |       |defined                                                                   |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|WindSpeed                |No/Yes   |list of|*Default : [10.0]*, Wind speed values for each layer in m/s.              |
|                         |if LO    |float  |Must have the same length as ``Cn2Weights``, ``Cn2Heights`` and           |
|                         |         |       |``WindDirection``.                                                        |
|                         |         |       |                                                                          |
|                         |         |       |*Warning* : required if ``Cn2Weights``, ``Cn2Heights`` or                 |
|                         |         |       |``WindDirection`` are defined                                             |
|                         |         |       |*Warning* : extremely confusing error message if absent when it must be   |
|                         |         |       |defined                                                                   |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|WindDirection            |No/Yes   |list of|*Default : [0.0]*, wind direction for each layer in degrees. 0 degree is  |
|                         |if LO    |float  |?? then anticlockwise.                                                    |
|                         |         |       |Must have the same length as ``Cn2Weights``, ``Cn2Heights`` and           |
|                         |         |       |``WindSpeed``.                                                            |
|                         |         |       |                                                                          |
|                         |         |       |*Warning* : required if ``Cn2Weights``, ``Cn2Heights`` or ``WindSpeed``   |
|                         |         |       |are defined                                                               |
|                         |         |       |*Warning* : extremely confusing error message if absent when it must be   |
|                         |         |       |defined                                                                   |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|r0_Value                 |No       |float  |Set the atmospere Fried parameter. If not set TIPTOP uses ``seeing`` .    |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|testWindspeed            |No       |float  |Used only for tests                                                       |
+-------------------------+---------+-------+--------------------------------------------------------------------------+

[sources_science]
-----------------

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Wavelength               |Yes      |list of |list of central wavelengths for each frame in meters. you can have more   |
|                         |         |float   |than one science target. needs explaining why the science sources can be  |
|                         |         |or float|multiple. (polychromatic? several targets? you can set many taget of the  |
|                         |         |        |same wavelength by only setting more than one Zenith and Azimuth but      |
|                         |         |        |leaving the wavelength as a float. It produces one PSF per target. The    |
|                         |         |        |number of PSF is the number of wavelength times the number of             |
|                         |         |        |Azimuth/Zenith couple.                                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Zenith                   |Yes      |list of |Zenithal coordinate in arcsec (distance from axis) of science sources     |
|                         |         |float   |given in ``Wavelength``. Must be the same length as ``Azimuth``           |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|Azimuth                  |Yes      |list of |Azimuthal coordinate in degree (angle from the ref. direction: polar axis |
|                         |         |float   |is x-axis) of science sources given in ``Wavelength``. Must be the same   |
|                         |         |        |length as ``Zenith``                                                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

[sources_HO]
------------

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Wavelength               |Yes      |float   |Sensing wavelength for Hight Order modes in meters,                       |
|                         |         |        |*Warning* : gives a confusing error message if absent                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Zenith                   |No       |list of |*Default : [0.0]*, Zenithal coordinate of each guide stars in arcsec      |
|                         |         |float   |(distance from axis). Must be the same length as ``Azimuth``, Even if     |
|                         |         |        |``Azimutal`` is defined, this is optional.                                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Azimuth                  |No       |list of |*Default : [0.0]*, Azimuthal coordinate in degree (angle from the ref.    |
|                         |         |float   |direction: polar axis is x-axis) of each guide stars.                     |
|                         |         |        |Must be the same length as ``Zenith``, even if ``Zenith`` is defined,     |
|                         |         |        |this is optional.                                                         |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Height                   |No       |float   |*Default : 0.0*, altitude of the guide stars (0 if infinite). Consider    |
|                         |         |        |that all guide star are at the same height.                               |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

[sources_LO]
------------
.. note::

   This section is completely optional (``[sensor_LO]`` section is required to have the LO part simulated)

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Wavelength               |Yes      |float   |Sensing wavelength for Low Order modes in meters                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Zenith                   |Yes      |list of |Zenithal coordinate of each guide stars in arcsec (distance from axis).   |
|                         |         |float   |Must be the same length as ``Azimuth``                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Azimuth                  |Yes      |list of |Azimuthal coordinate in degree (angle from the reference direction: polar |
|                         |         |float   |axis is x-axis) of each guide stars.                                      |
|                         |         |        |Must be the same length as ``Zenith``                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   

[sources_Focus]
---------------
.. note::

   This section is completely optional.
   The ``[sources_Focus]`` section is required to have the global focus part simulated considering specific focus sensors and not the LO sensors.
   This happens when the key ``glFocusOnNGS`` in the ``[telescope]`` section is True and multiple DMs are present.

   Note that the coordinates (``Zenith`` and ``Azimuth``) of the NGSs are the same of the ``[sources_LO]`` section.

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Wavelength               |Yes      |float   |Sensing wavelength for global focus modes in meters                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

[sensor_science]
----------------

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|PixelScale               |Yes      |float   |Pixel/spaxel scale in milliarcsec.                                        |
|                         |         |        |                                                                          |
|                         |         |        |*Warning*: confusing error message if missing                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |float   |Field of view of the camera in pixel/spaxel.                              |
|                         |         |        |                                                                          |
|                         |         |        |*Warning*: confusing error massage if missing                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

.. note::

    Following parameters were added to uniformise all the sensor (HO and LO), but they are not used.

    Binning, NumberPhotons, SpotFWHM, SpectralBandwidth, Transmittance, Dispersion, SigmaRON, Dark, SkyBackground, Gain, ExcessNoiseFactor, Wavelength, FieldOfView

[sensor_HO]
-----------

The High Order WaveFront Sensor can be a pyramid WFS or a Shack-Hartmann. Regardless of the WFS, the following parameters can de defined.

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|NumberLenslets           |No       |list of |*Default : [20]*, Number of WFS lenslets. Used the same way in            |
|                         |         |int     |Shack-Hartmann wavefront sensor and Pyramid. Also used for noise          |
|                         |         |        |computation if `NoiseVariance` is not set.                                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SizeLenslets             |No       |list of |*Default: [Telescope] TelescopeDiameter/[sensor_HO] NumberLenslet*        |
|                         |         |float   |Size of WFS lenslets in meters. used, why a list of float? This overrides |
|                         |         |        |the ratio between telescope size and Number of lenslet used to compute the|
|                         |         |        |matrix size.                                                              |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|PixelScale               |Yes      |integer |High Order WFS pixel scale in [mas],  Not used when a pyramid wavefront   |
|                         |         |        |sensor has been selected.                                                 |
|                         |         |        |                                                                          |
|                         |         |        |*Warning*: gives a confusing error message if missing                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |integer |Number of pixels per subaperture. Not used when a pyramid wavefront sensor|
|                         |         |        |has been selected (4 pixels are used in this case).                       |
|                         |         |        |                                                                          |
|                         |         |        |*Warning*: gives a confusing error message if missing                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WfsType                  |No       |string  |*default : 'Shack-Hartmann'*, type of wavefront sensor used for the High  |
|                         |         |        |Order sensing. Other available option: 'Pyramid'                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberPhotons            |No       |list of |*default : [Inf]*, Flux return in [nph/frame/subaperture]                 |
|                         |         |integer |                                                                          |
|                         |         |        |It can be computed as:                                                    |
|                         |         |        |                                                                          |
|                         |         |        |``(0-magn-flux [ph/s/m2]) * (size of sub-aperture [m])^2                  |
|                         |         |        |* (1/SensorFrameRate_HO) * (total throughput)                             |
|                         |         |        |* (10^(-0.4*magn_source_HO))``                                            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SpotFWHM                 |No       |list of |*defaut: [[0.0, 0.0, 0.0]]*, High Order spot parameters: two axes scale   |
|                         |         |list of |values in milliarcsec (only max value is used) and angle (angle is not    |
|                         |         |float   |used). Why list?                                                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|SpectralBandwidth        |No       |float   |*default: 0.0*, Not used, spectral bandwidth of the filter (imaging mode)?|
|                         |         |        |why specific to the imaging mode? what is the effect?                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Transmittance            |No       |list of |*default: [1.0]*, Used for PSF computation and flux scaling but not with  |
|                         |         |float   |noise computation. Transmittance at the considered wavelengths for        |
|                         |         |        |polychromatic mode. How do you set polychromatic mode? Each element can   |
|                         |         |        |not have a value superior to 1?                                           |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|Dispersion               |No       |list of |*default: [[0.0,0.0]]*, Dispersion x/y at the considered wavelength in    |
|                         |         |list of |pixel. Must be the same size than ``Transmittance``. Chromatic dispertion |
|                         |         |float?  |for PSF computation only. In HarmoniSCAO_1 first the default and the thing|
|                         |         |        |given are not even the same shape but on top the default breaks the must  |
|                         |         |        |be the same size as the transmitance... Also sorry for my ignorance:      |
|                         |         |        |dispersion of what? Isn't this maybe redundant with `SpotFWHM` ?          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Gain                     |No       |float   |*default : 1.0*, Pixel gain. do you mean camera gain or loop goin?        |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ExcessNoiseFactor        |No       |float   |*default : 2.0*, excess noise factor. TODO: default should be 1           |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NoiseVariance            |No       |unknown |*Default : None*?, Noise Variance in rad2. If not empty, this value       |
|                         |         |        |overwrites the analytical noise variance calculation.                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SigmaRON                 |No       |float   |*Default : 0.0*, read-out noise std in [e-], used only if the             |
|                         |         |        |`NoiseVariance` is not set.                                               |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|addMcaoWFsensConeError   |No       |string  |*Default : False*, additional error to consider the reduced sensing volume|
|                         |         |        |due to the cone effect. Multi-conjugate systems only.                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

In the two following section we list the parameters that are specific to each wavefront sensor. If you define a parameter 
for one WFS while another WFS is defined The parameter will be ignired. For example, if you define the parameter SigmaRON,
while WfsType is 'Pyramid', SigmaRON is ignored.

Shack-Hartmann requirement
^^^^^^^^^^^^^^^^^^^^^^^^^^

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Algorithm                |not used |string  |*defaut:'wcog'*, other options: 'cog' (simple center-of-gravity), 'tcog'  |
|                         |         |        |(center-of-gravity with threshold), 'qc' (quad-cell)                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WindowRadiusWCoG         |not used |int     |*default: 2*, FWHM in pixel of the gaussian weighting function            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

Pyramid requirement
^^^^^^^^^^^^^^^^^^^

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Modulation               |Yes      |float   |*default : None*, If the chosen wavefront sensor is the ``'Pyramid'``,    |
|                         |         |        |Spot modulation radius in lambda/D units. This is ignored if the WFS is   |
|                         |         |        |`'Shack-Hartmann'`                                                        |
|                         |         |        |                                                                          |
|                         |         |        |*Warning* : gives a confusing message if missing when required            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Binning                  |No       |integer |*default: 1*, Binning factor of the detector, only used in the pyramid    |
|                         |         |        |case, optional for pyramid                                                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

Can be set but not used
^^^^^^^^^^^^^^^^^^^^^^^

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Dark                     |not used |float   |*default: 0.0*, dark current in [e-/s/pix]                                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SkyBackground            |not used |float   |*default: 0.0*, Sky background [e-/s/pix]                                 |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ThresholdWCoG            |not used |float?  |*default: 0.0*, Threshold Number of pixels for windowing the low order WFS| 
|                         |         |        |pixels                                                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NewValueThrPix           |not used |float   |*default: 0.0*, New value for pixels lower than `ThresholdWCoG`. Is there |
|                         |         |        |a reason to want to force these values to something else?                 |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

[sensor_LO]
-----------

.. note::

   This section is optional, if this section is not present only the HO part will be used (for ex. to simulate a SCAO NGS).

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|PixelScale               |Yes      |float   |LO WFS pixel scale in [mas],                                              |
|                         |         |        |*Warning*: gives a confusing error message if missing                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |integer |not used. Number of pixels per subaperture,                               |
|                         |         |        |*Warning*: gives a confusing error message if missing                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberPhotons            |Yes      |list of |Detected flux in [nph/frame/subaperture], Must be the same length as      |
|                         |         |integer |NumberLenslet                                                             |
|                         |         |        |                                                                          |
|                         |         |        |It can be computed as:                                                    |
|                         |         |        |                                                                          |
|                         |         |        |``(0-magn-flux [ph/s/m2]) * (size of subaperture [m])**2                  |
|                         |         |        |* (1/SensorFrameRate_LO) * (total throughput)                             |
|                         |         |        |* (10**(-0.4*magn_source_LO))``                                           |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberLenslets           |Yes      |list of |*Default : [1]*, number of WFS lenslets, Must be the same length as       |
|                         |         |integer |NumberPhotons                                                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SigmaRON                 |Yes      |float   |*default: 0.0*, read out noise in [e-]                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Dark                     |Yes      |float   |*default: 0.0*, dark current[e-/s/pix]                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SkyBackground            |Yes      |float   |*default: 0.0*, Sky background [e-/s/pix]                                 |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ExcessNoiseFactor        |Yes      |float   |*default: 2.0*, excess noise factor                                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WindowRadiusWCoG         |Yes      |integer |*default: 1*, Radius in pixel of the HWHM of the weights map of the       |
|                         |         |        |weighted CoG the low order WFS pixels                                     |
|                         |         |        |                                                                          |
|                         |         |        |*Warning* : if set to 'optimize', gain is automatically optimized by      |
|                         |         |        |TIPTOP (closest int to half of PSF FWHM), otherwise the float value set is|
|                         |         |        |used.                                                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|ThresholdWCoG            |Yes      |float   |*default: 0.0*, Threshold Number of pixels for windowing the low order WFS|
|                         |         |        |pixels                                                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NewValueThrPix           |Yes      |float   |*default: 0.0*, New value for pixels lower than threshold.                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|filtZernikeCov           |No       |string  |*Default : False*, Filter for the zernike covariance. The zernike cov. is |
|                         |         |        |used to quantify for the TT tomographic (anisoplanatic) error. This filter|
|                         |         |        |accounts for the HO correction of an MCAO system. Multi-conjugate systems |
|                         |         |        |only.                                                                     |
|                         |         |        |                                                                          |
|                         |         |        |*Warning*: Do not use in systems with a single DM.                        |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

Can be set but not used
^^^^^^^^^^^^^^^^^^^^^^^

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Binning                  |not used |integer |*default: 1*, binning factor of the detector                              |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SpotFWHM                 |not used |list of |*default: [[0.0, 0.0, 0.0]]*, Low Order spot scale in [mas]               |
|                         |         |list of |                                                                          |
|                         |         |integer |                                                                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|Gain                     |not used |float   |*default: 1.0*, Camera gain                                               |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Algorithm                |not used |string  |*default: 'wcog'*, CoG computation algorithm                              |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

[sensor_Focus]
-----------

.. note::

   This section is completely optional.
   The ``[sensor_Focus]`` section is required to have the global focus part simulated considering specific focus sensors and not the LO sensors.
   This happens when the key ``glFocusOnNGS`` in the ``[telescope]`` section is True and multiple DMs are present.

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|PixelScale               |Yes      |float   |Focus WFS pixel scale in [mas],                                           |
|                         |         |        |*Warning*: gives a confusing error message if missing                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |integer |not used. Number of pixels per subaperture,                               |
|                         |         |        |*Warning*: gives a confusing error message if missing                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberPhotons            |Yes      |list of |Detected flux in [nph/frame/subaperture], Must be the same length as      |
|                         |         |integer |NumberLenslet                                                             |
|                         |         |        |                                                                          |
|                         |         |        |It can be computed as:                                                    |
|                         |         |        |                                                                          |
|                         |         |        |``(0-magn-flux [ph/s/m2]) * (size of subaperture [m])**2                  |
|                         |         |        |* (1/SensorFrameRate_Focus) * (total throughput)                          |
|                         |         |        |* (10**(-0.4*magn_source_Focus))``                                        |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberLenslets           |Yes      |list of |*Default : [1]*, number of WFS lenslets, Must be the same length as       |
|                         |         |integer |NumberPhotons                                                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SigmaRON                 |Yes      |float   |*default: 0.0*, read out noise in [e-]                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Dark                     |Yes      |float   |*default: 0.0*, dark current[e-/s/pix]                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SkyBackground            |Yes      |float   |*default: 0.0*, Sky background [e-/s/pix]                                 |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ExcessNoiseFactor        |Yes      |float   |*default: 2.0*, excess noise factor                                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WindowRadiusWCoG         |Yes      |integer |*default: 1*, Radius in pixel of the HWHM of the weights map of the       |
|                         |         |        |weighted CoG the global focus WFS pixels                                  |
|                         |         |        |                                                                          |
|                         |         |        |*Warning* : if set to 'optimize', gain is automatically optimized by      |
|                         |         |        |TIPTOP (closest int to half of PSF FWHM), otherwise the float value set is|
|                         |         |        |used.                                                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|ThresholdWCoG            |Yes      |float   |*default: 0.0*, Threshold Number of pixels for windowing the low order WFS|
|                         |         |        |pixels                                                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NewValueThrPix           |Yes      |float   |*default: 0.0*, New value for pixels lower than threshold.                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+


[DM]
----

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|NumberActuators          |Yes      |list of |Number of actuator on the pupil diameter. why a list of int? Must be the  |
|                         |         |integer |same length as DmPitchs. *Warning*: gives a confusing error message if    |
|                         |         |        |missing. *Warning*: not used in TIPTOP!                                   |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|DmPitchs                 |Yes      |list of |DM actuators pitch in meters, on the meta pupil at the conjugasion        |
|                         |         |float   |altitude, used for fitting error computation. Must be the same length as  |
|                         |         |        |NumberActuators? *Warning*: gives a confusing error message if missing    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|InfModel                 |No       |string  |*default: 'gaussian'*, DM influence function model. Not used in TIPTOP but| 
|                         |         |        |used in the psf reconstruction. What are the other possible one?          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|InfCoupling              |No       |list of |*default: [0.2]*, DM influence function model mechanical coupling. Not    | 
|                         |         |float   |used in TIPTOP but used in the psf reconstruction. Unclear what this does.|
|                         |         |        |Must be the same length as NumberActuators?                               |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|DmHeights                |No/Yes if|list of |*default: [0.0]*, DM altitude in meters, Must be the same length as       |
|                         |LO       |float   |NumberActuators and DmPitchs                                              |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|OptimizationZenith       |No       |float   |*default: [0.0]*, Zenith position in arcsec (distance from axis) of the   |
|                         |         |        |direction in which the AO correction is optimized. Must be the same length|
|                         |         |        |as OptimisationAzimuth  and OptimizationWeight. These are for wide field  |
|                         |         |        |AO system, should be a requirement for MCAO and GLAO                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|OptimizationAzimuth      |No       |list of |*default: [0.0]*, Azimuth in degrees (angle from the ref. direction: polar|
|                         |         |float   |axis is x-axis) of the direction in which the AO correction is optimized. |
|                         |         |        |Must be the same length as OptimizationZenith and OptimizationWeight.     |
|                         |         |        |These are for wide field AO system, should be a requirement for MCAO and  |
|                         |         |        |GLAO                                                                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|OptimizationWeight       |No       |float   |*default: [1.0]*, Weights of the optimisation directions. Must be the same|
|                         |         |        |length as OptimizationZenith and OptimizationAzimuth. These are for wide  |
|                         |         |        |field AO system, should be a requirement for MCAO and GLAO.               |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|OptimizationConditioning |No       |float   |*default: 1.0e2*, Matrix Conditioning threshold in the truncated SVD      |
|                         |         |        |inversion.                                                                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberReconstructedLayers|No       |integer |*default: 10*, Only used for wide field AO system, (meaning more than one |
|                         |         |        |guide star is defined). Number of reconstructed layers for tomographic    |
|                         |         |        |systems. Shouldn't this be defaulted to 1 for SCAO sakes?                 |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|AoArea                   |No       |string  |*default: 'circle'*, Shape of the AO-corrected area. Any other options are| 
|                         |         |        |not defined and will give a squarre correction area.                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

[RTC]
-----

.. note::

   This section is optional, if this section is not present the defaul values are used.

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|LoopGain_HO              |No       |float   |*Default : 0.5*, High Order Loop gain. *Warning*: if system to be         |
|                         |         |        |simulated is a multi-conjugate system this parameter is not used.         |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SensorFrameRate_HO       |No       |float   |*Default : 500.0*, High Order loop frequency in [Hz]                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|LoopDelaySteps_HO        |No       |integer |*Default : 2*, High Order loop delay in [frame]                           |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|LoopGain_LO              |No/Yes if|float or|*default: None*, Low Order loop gain, *Warning*: if set to 'optimize',    |
|                         |LO       |string  |gain is automatically optimized by TIPTOP, otherwise the float value set  |
|                         |         |        |is used.                                                                  |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SensorFrameRate_LO       |No/Yes if|float   |*default: None*, Loop frequency in [Hz]. If ``[sensor_LO]`` section is    |
|                         |LO       |        |present it must be set.                                                   |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|LoopDelaySteps_LO        |No/Yes if|integer |*default: None*, Low Order loop delays in [frames]. If ``[sensor_LO]``    |
|                         |LO       |        |section is present it must be set.                                        |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|LoopGain_Focus           |No/Yes if|float or|*default: None*, Global focus loop gain, *Warning*: if set to 'optimize', |
|                         |Focus    |string  |gain is automatically optimized by TIPTOP, otherwise the float value set  |
|                         |         |        |is used.                                                                  |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|SensorFrameRate_Focus    |No/Yes if|float   |*default: None*, Global focus loop frequency in [Hz]. If                  |
|                         |Focus    |        |``[sensor_Focus]`` section is present it must be set.                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|LoopDelaySteps_Focus     |No/Yes if|integer |*default: None*, Global focus loop delays in [frames]. If                 |
|                         |Focus    |        |``[sensor_Focus]`` section is present it must be set.                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ResidualError            |No       |?       |*Default: None*, ?                                                        |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

[COMPUTATION]
-------------

.. note::

   This section is optional, if this section is not present the defaul values are used.

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|simpleVarianceComputation|No       |string  |Set to it to False to activate the more complex and slower MASTSEL LO     |
|                         |         |        |noise computation, that is more accurate in low flux regimes.             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|platform                 |No       |string  |*default: 'GPU'* Set to it to 'CPU' to forcy the library to use numpy     |
|                         |         |        |instead of cupy.                                                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|integralDiscretization1  |No       |float   |*default: 1000.0*, Discretization used in the integrals                   |
|                         |         |        |(astro-tiptop/SEEING library).                                            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|integralDiscretization2  |No       |float   |*default: 4000*, Discretization used in the integrals                     |
|                         |         |        |(astro-tiptop/SEEING library).                                            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
