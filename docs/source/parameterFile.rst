Parameter files explained
=========================

Introduction
------------

To run TIPTOP you need two things: to execute the function and to create a parameter file. This section explaines
what should your parameter file contain and the various parameters you can set. You can already find example parameter 
files in the `github <https://github.com/astro-tiptop/TIPTOP/tree/main/tiptop/perfTest>`_ .


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
|Resolution                |No        |integer|*Default : 256*, Number of pixels across the pupil diameter.              |
|                          |          |       |This value is used in computation of the telescope OTF.                   |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|ObscurationRatio          |No        |float  |*Default : 0.0*, Defines the central obstruction                          |
|                          |          |       |due to the secondary as a ratio of the ``TelescopeDiameter``.             |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|ZenithAngle               |No/Yes if |float  |*Default: 0.0*, Set the pointing direction of the telescope in degree     |
|                          |LO        |       |with respect to the zenith. Used to compute airmass, to scale atmospheric |
|                          |          |       |layers and stars altitude.                                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PupilAngle                |No        |float  |*default : 0.0*, Rotation angle of the telescope pupil in degrees.        |
|                          |          |       |Applied to pupil mask and static aberration maps to match instrument      |
|                          |          |       |orientation.                                                              |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathPupil                 |No        |string |*default: ''*, Path to the pupil model in .fits file (if provided,        |
|                          |          |       |the pupil model is interpolated). if absent or '', not used.              |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathStaticOn              |No        |string |*default: None*, Path to a fits file containing an on-axis static         |
|                          |          |       |aberration map ([nm]). This parameter can be used to add any kind of      |
|                          |          |       |static aberrations. Example: the static aberration of the `ELT M1`_.      |
|                          |          |       |If absent or '', not used.                                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|zCoefStaticOn             |No        |list of|*default: None*, Combination of zernike modes that models an on-axis      |
|                          |          |float  |static aberration. Coefficients are in [nm RMS].                          |
|                          |          |       |Examples: Focus ``[0,0,100]``; Astigmatism ``[0,0,0,100]``;               |
|                          |          |       |Trefoil ``[0,0,0,0,0,0,0,100]``; Spherical ``[0,0,0,0,0,0,0,0,0,100]``.   |
|                          |          |       |If absent not used.                                                       |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathStaticOff             |No        |string |*default: None*, Path to a .fits file that contains field-dependent       |
|                          |          |       |(off-axis) static aberration maps. Must be provided together with         |
|                          |          |       |``PathStaticPos`` specifying the corresponding positions.                 |
|                          |          |       |If absent or '', not used.                                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathStaticPos             |No        |string |*default: None*, Required if ``PathStaticOff``. Path to a fits file that  |
|                          |          |       |contains the field positions [zenith in arcsec, azimuth in rad]           |
|                          |          |       |corresponding to each off-axis static aberration map in ``PathStaticOff``.|
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathApodizer              |No        |string |*default: ''*, Path to a fits file that contains an amplitude apodizer    |
|                          |          |       |map. Used to apply pupil amplitude weighting (transmission mask) in the   |
|                          |          |       |simulation. if absent or '', not used.                                    |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|PathStatModes             |No        |string |*default: ''*, path to a .fits file that contains a cube of static        |
|                          |          |       |aberration modes. Each mode is normalized to have 1 nm RMS amplitude.     |
|                          |          |       |If absent or '', not used.                                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|windPsdFile               |No        |string |*default: ''*, file name of a fits file with a 2D array with a frequency  |
|                          |          |       |vector and PSD of tip and tilt windshake.                                 |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorNm              |No        |float  |*default: 0*, nm RMS of the additional wavefront error to be added (an    |
|                          |          |       |error that is not otherwise considered). This parameter is used to define |
|                          |          |       |a PSD that is summed to the AO PSD. The default power law is f^(-2), but  |
|                          |          |       |is can be modified using the ``extraErrorExp`` parameter. It models a     |
|                          |          |       |generic static aberration.                                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorExp             |No        |float  |*default: -2*, Exponent of the power of spatial frequencies used to       |
|                          |          |       |generate the PSD associated with extraErrorNm                             |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorMin             |No        |float  |*default: 0*, Minimum spatial frequency ([m^(-1)]) for which PSD          |
|                          |          |       |associated with extraErrorNm is > 0                                       |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorMax             |No        |float  |*default: 0*, Maximum spatial frequency ([m^(-1)]) for which PSD          |
|                          |          |       |associated with extraErrorNm is > 0                                       |
|                          |          |       |                                                                          |
|                          |          |       |Note: 0 means maximum frequency is the one present in the spatial         |
|                          |          |       |frequency array of the PSDs.                                              |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|extraErrorLoNm            |No        |float  |*default: 0*, nm RMS of the additional wavefront error to be added (an    |
|                          |          |       |error that is not otherwise considered).                                  |
|                          |          |       |                                                                          |
|                          |          |       |It can be a list of two values, the on-axis error and the error at the    |
|                          |          |       |edge of the technical field ([telescope]TechnicalFoV)                     |
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
|                          |          |       |[FWHM_mas_max, FWHM_mas_min, angle_rad]. It models an additional tip/tilt |
|                          |          |       |jitter (e.g. vibrations, PSF drifts, ...).                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|glFocusOnNGS              |No        |bool   |*default: False*, global focus control with natural guide stars.          |
|                          |          |       |Multi-conjugate systems only. Requires NumberLenslets >= 2 in sensor_LO or|
|                          |          |       |a specific global focus sensor (``[sources_Focus]`` and ``[sensor_Focus]``|
|                          |          |       |sections).                                                                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+
|TechnicalFoV              |No/Yes if |float  |*default: 0.0*, Diameter of the technical field of view in [arcsec]. In   |
|                          |LO        |       |MCAO/LGS configurations, used when ``NumberActuators`` from ``[DM]``      |
|                          |          |       |section is not set: scales the projected DM size with altitude and        |
|                          |          |       |derives the actuator count from ``DmPitchs``. In LO, also sets the        |
|                          |          |       |angular range for interpolating additional low-order error terms          |
|                          |          |       |(``extraErrorLoNm``).                                                     |
|                          |          |       |*Warning*: Mandatory and no default if LO section is used.                |
+--------------------------+----------+-------+--------------------------------------------------------------------------+

.. _ELT M1: https://github.com/astro-tiptop/TIPTOP/blob/main/tiptop/data/ELT_M1_MORFEO_DMs_static_wfe_480px.fits


[atmosphere]
------------

+-------------------------+------------+-------+--------------------------------------------------------------------------+
| Parameter               |  Required  | Type  | Description                                                              |
+=========================+============+=======+==========================================================================+
|Seeing                   |Yes, unless |float  |Set the seeing at Zenith in [arcsec]. Used to compute ``r0`` as           |
|                         |``r0_value``|       |``r0 = 0.976 × λ / Seeing(rad)``. If not set, TipTop uses ``r0_value``.   |
|                         |given       |       |                                                                          |
+-------------------------+------------+-------+--------------------------------------------------------------------------+
|r0_Value                 |Yes, unless |float  |Set the atmosphere Fried parameter ``r0`` in [meters]. Used directly      |
|                         |``Seeing``  |       |if ``Seeing`` is not provided.                                            |
|                         |given       |       |                                                                          |
+-------------------------+------------+-------+--------------------------------------------------------------------------+
|Wavelength               |No/Yes if LO|float  |*Default : 500e-9*, Wavelength at which the atmospheric statistics are    |
|                         |            |       |defined (in meters).                                                      |
|                         |            |       |                                                                          |
|                         |            |       |*Warning*: Mandatory and no default if LO section is used.                |
+-------------------------+------------+-------+--------------------------------------------------------------------------+
|L0                       |No/Yes if LO|float  |*Default : 25.0*, Outer Scale of the atmosphere  in meters                |
|                         |            |       |                                                                          |
|                         |            |       |*Warning*: not optional in MavisLO.py                                     |
+-------------------------+------------+-------+--------------------------------------------------------------------------+
|Cn2Weights               |No/Yes      |list of|*Default : [1.0]*, Relative contribution of each layer. The sum of all the|
|                         |if LO       |float  |list element must be 1. Must have the same length as ``Cn2Heights``,      |
|                         |            |       |``WindSpeed`` and ``WindDirection``.                                      |
|                         |            |       |                                                                          |
|                         |            |       |*Warning : required if ``Cn2Heights``, ``WindSpeed`` or ``WindDirection`` |
|                         |            |       |are defined                                                               |
|                         |            |       |*Warning* : extremely confusing error message if absent when it must be   |
|                         |            |       |defined                                                                   |
+-------------------------+------------+-------+--------------------------------------------------------------------------+
|Cn2Heights               |No/Yes      |list of|*Default : [0.0]*, altitude of layers in [meters].                        |
|                         |if LO       |float  |Must have the same length as ``Cn2Weights``, ``WindSpeed`` and            |
|                         |            |       |``WindDirection``.                                                        |
|                         |            |       |                                                                          |
|                         |            |       |*Warning* : required if ``Cn2Weights``, ``WindSpeed`` or ``WindDirection``|
|                         |            |       |are defined                                                               |
|                         |            |       |*Warning* : extremely confusing error message if absent when it must be   |
|                         |            |       |defined                                                                   |
+-------------------------+------------+-------+--------------------------------------------------------------------------+
|WindSpeed                |No/Yes      |list of|*Default : [10.0]*, Wind speed values for each layer in [m/s].            |
|                         |if LO       |float  |Must have the same length as ``Cn2Weights``, ``Cn2Heights`` and           |
|                         |            |       |``WindDirection``.                                                        |
|                         |            |       |                                                                          |
|                         |            |       |*Warning* : required if ``Cn2Weights``, ``Cn2Heights`` or                 |
|                         |            |       |``WindDirection`` are defined                                             |
|                         |            |       |*Warning* : extremely confusing error message if absent when it must be   |
|                         |            |       |defined                                                                   |
+-------------------------+------------+-------+--------------------------------------------------------------------------+
|WindDirection            |No          |list of|*Default : a list of 0 of the length of WindSpeed*, wind direction for    |
|                         |            |float  |each layer in [degrees]. 0 degree is alogn the x axis then anticlockwise. |
|                         |            |       |Must have the same length as ``Cn2Weights``, ``Cn2Heights`` and           |
|                         |            |       |``WindSpeed``.                                                            |
+-------------------------+------------+-------+--------------------------------------------------------------------------+
|testWindspeed            |No          |float  |Used only for tests                                                       |
+-------------------------+------------+-------+--------------------------------------------------------------------------+

[sources_science]
-----------------

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Wavelength               |Yes      |list of |list of wavelengths in meters.                                            |
|                         |         |float   |                                                                          |
|                         |         |or float|When more than one elements is present the output PSF saved in the fits   |
|                         |         |        |file is a 4D array with dimension (Nw, Ns, Npix, Npix), where Nw is the   |
|                         |         |        |number of wavelengths required ([sources_science] Wavelength), Ns is the  |
|                         |         |        |number of directions required ([sources_science] Zenith and Azimuth) and  |
|                         |         |        |Npix is the size required for the PSFs ([sensor_science] FieldOfView).    |
|                         |         |        |If a single elements is present the fits file is a 3D array with          |
|                         |         |        |dimension (Ns, Npix, Npix).                                               |
|                         |         |        |Instead the profiles will be a 3D array (fourth fits file extension) with |
|                         |         |        |dimensions (2*Nw, Ns, Npix/2). The first Nw elements contain the radius   |
|                         |         |        |and the second Nw elements the profile values (the first radius and       |
|                         |         |        |profile pair is radius=data[0,0,:] profile=data[Nw,0,:], the second is    |
|                         |         |        |radius=data[1,0,:] profile=data[Nw+1,0,:], ...)                           |
|                         |         |        |json file: two lists, radius and psf with dimensions (Nw, Ns, Npix/2).    |
|                         |         |        |                                                                          |
|                         |         |        |In this case more memory is required and small differences with respect   |
|                         |         |        |to monochromatic PSF will be present because: (1) errors Differential     | 
|                         |         |        |refractive anisoplanatism and Chromatism from P3 are computed for a       |
|                         |         |        |single wavelength (the shortest one) (2) effective field-of-view of the   |
|                         |         |        |PSF is typically larger to guarantee that the PSF at the shortest         |
|                         |         |        |wavelength has the required field-of-view (3) The PSF is typically        |
|                         |         |        |computed with a higher sampling to guarantee that the longest wavelength  |
|                         |         |        |has the required sampling and then the PSFs at the shorter wavelengths    |
|                         |         |        |are rebinned.                                                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Zenith                   |Yes      |list of |Zenithal coordinate in arcsec (distance from axis) of science sources.    |
|                         |         |float   |Must be the same length as ``Azimuth``                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|Azimuth                  |Yes      |list of |Azimuthal coordinate in degree (angle from the ref. direction: polar axis |
|                         |         |float   |is x-axis) of science sources. Must be the same length as ``Zenith``      |
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
|Wavelength               |Yes      |float   |Sensing wavelength for Low Order modes in meters.                         |
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
|FieldOfView              |Yes      |integer |Field of view of the camera in pixel/spaxel.                              |
|                         |         |        |                                                                          |
|                         |         |        |*Warning*: confusing error massage if missing                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Super_Sampling           |No       |float   |Desired radial interpolation sampling step in milliarcsec.                |
|                         |         |        |                                                                          |
|                         |         |        |If provided, TipTop performs a 2D polar interpolation of the PSF to       |
|                         |         |        |produce a radial profile resampled at the requested scale.                |
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
|WfsType                  |No       |string  |*default : 'Shack-Hartmann'*, Type of wavefront sensor used for the High  |
|                         |         |        |Order sensing. Other available option: 'Pyramid'                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberLenslets           |No       |list of |*Default : [20]*, Number of WFS lenslets.  Used for both                  |
|                         |         |int     |Shack-Hartmann and Pyramid sensors. Also used for noise computation if    |
|                         |         |        |``NoiseVariance`` is not set.                                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SizeLenslets             |No       |list of |*Default: ``[Telescope] TelescopeDiameter/[sensor_HO] NumberLenslet``*    |
|                         |         |float   |Size of WFS lenslets in meters. Overrides the ratio between telescope     |
|                         |         |        |size and Number of lenslet used to compute the matrix size.               |
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
|NumberPhotons            |No       |list of |*default : [Inf]*, Flux return in [nph/frame/subaperture]                 |
|                         |         |integer |                                                                          |
|                         |         |        |It can be computed as:                                                    |
|                         |         |        |                                                                          |
|                         |         |        |``(0-magn-flux [ph/s/m2]) * (size of sub-aperture [m])^2                  |
|                         |         |        |* (1/SensorFrameRate_HO) * (total throughput)                             |
|                         |         |        |* (10^(-0.4*magn_source_HO))``                                            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SpotFWHM                 |No       |list of |*defaut: [[0.0, 0.0]]*, Represents the instrumental broadening of         |
|                         |         |list of |Shack–Hartmann spot size (FWHM) along x and y, in [milliarcseconds]       |
|                         |         |float   |without turbulence. If set to [[0.0, 0.0]], only atmospheric broadening is|
|                         |         |        |considered. Not used with a Pyramid WFS.                                  |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|SpectralBandwidth        |No       |float   |*default: 0.0*, Spectral fullwidth around each central wavelength         |
|                         |         |        |(in [meters]). If 0, monochromatic simulation.                            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Transmittance            |No       |list of |*default: [1.0]*, Transmission factors at the WFS plane. Expected in the  |
|                         |         |float   |range [0,1].                                                              |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|Dispersion               |No       |list of |*default: [[0.0],[0.0]]*, Chromatic shift of the image on the detector,   |
|                         |         |list of |in pixels. The first sub-list corresponds to x-offsets, the second to     |
|                         |         |float?  |y-offsets. Must have the same number of elements as ``Transmittance``.    |
|                         |         |        |Used only in PSF computation to account for wavelength-dependent shifts   |
|                         |         |        |(e.g. due to residual atmospheric dispersion).                            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Gain                     |No       |float   |*default : 1.0*, Detector pixel gain.                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ExcessNoiseFactor        |No       |float   |*default : 1.0*, excess noise factor.                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NoiseVariance            |No       |list of |*Default : [None]*, Noise Variance in rad2. If set, this value            |
|                         |         |float   |overrides the analytical noise variance calculation.                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SigmaRON                 |No       |float   |*Default : 0.0*, read-out noise std in [e-], used only if the             |
|                         |         |        |`NoiseVariance` is not set.                                               |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|addMcaoWFsensConeError   |No       |bool    |*Default : False*, additional error to consider the reduced sensing volume|
|                         |         |        |due to the cone effect. Multi-conjugate systems only.                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

In the two following section we list the parameters that are specific to each wavefront sensor. If you define a parameter 
for one WFS while another WFS is defined The parameter will be ignored. For example, if you define the parameter SigmaRON,
while WfsType is 'Pyramid', SigmaRON is ignored.

Shack-Hartmann requirement
^^^^^^^^^^^^^^^^^^^^^^^^^^

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Algorithm                |No       |string  |*defaut:'wcog'*, other options: 'cog' (simple center-of-gravity), 'tcog'  |
|                         |         |        |(center-of-gravity with threshold), 'qc' (quad-cell)                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WindowRadiusWCoG         |No       |integer |*default: 5*, FWHM in pixel of the gaussian weighting function            |
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
|Dark                     |No       |float   |*default: 0.0*, dark current in [e-/s/pix]                                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SkyBackground            |No       |float   |*default: 0.0*, Sky background [e-/s/pix]                                 |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ThresholdWCoG            |No       |float?  |*default: 0.0*, Threshold Number of pixels for windowing the low order WFS| 
|                         |         |        |pixels                                                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NewValueThrPix           |No       |float   |*default: 0.0*, New value for pixels lower than `ThresholdWCoG`. Is there |
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
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |integer |Number of pixels per subaperture.  Not used when a Pyramid wavefront      |
|                         |         |        |sensor has been selected (4 pixels are used in this case).                |
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
|ExcessNoiseFactor        |Yes      |float   |*default: 1.0*, excess noise factor                                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WindowRadiusWCoG         |Yes      |integer |*default: 1*, Radius in pixel of the FWHM of the weights map of the       |
|                         |         |or      |weighted CoG the low order WFS pixels                                     |
|                         |         |string  |                                                                          |
|                         |         |        |*Warning* : if set to 'optimize', gain is automatically optimized by      |
|                         |         |        |TIPTOP (closest int to half of PSF FWHM), otherwise the float value set is|
|                         |         |        |used.                                                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|ThresholdWCoG            |Yes      |float   |*default: 0.0*, Threshold Number of pixels for windowing the low order WFS|
|                         |         |        |pixels                                                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NewValueThrPix           |Yes      |float   |*default: 0.0*, New value for pixels lower than threshold.                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|filtZernikeCov           |No       |bool    |*Default : False*, Filter for the zernike covariance. The zernike cov. is |
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
|Binning                  |No       |integer |*default: 1*, binning factor of the detector                              |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SpotFWHM                 |No       |list of |*default: [[0.0, 0.0]]*, Low Order spot scale in [mas]                    |
|                         |         |list of |                                                                          |
|                         |         |integer |                                                                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|Gain                     |No       |float   |*default: 1.0*, Camera gain                                               |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Algorithm                |No       |string  |*default: 'wcog'*, CoG computation algorithm                              |
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
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |integer |not used. Number of pixels per subaperture,                               |
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
|ExcessNoiseFactor        |Yes      |float   |Excess noise factor                                                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WindowRadiusWCoG         |Yes      |integer |*default: 1*, Radius in pixel of the HWHM of the weights map of the       |
|                         |         |or      |weighted CoG the global focus WFS pixels                                  |
|                         |         |string  |                                                                          |
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
|DmPitchs                 |Yes      |list of |DM actuators pitch in meters, on the meta pupil at the conjugation        |
|                         |         |float   |altitude, used for fitting error computation.                             |
|                         |         |        |                                                                          |
|                         |         |        |*Warning*: if it is smaller than [sensor_HO] SizeLenslets                 |
|                         |         |        |(= [Telescope] TelescopeDiameter/[sensor_HO] NumberLenslet ) aliasing     |
|                         |         |        |error will be significant.                                                |
|                         |         |        |                                                                          |
|                         |         |        |Must be the same length as NumberActuators                                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberActuators          |No       |list of |*default: computed from ``TelescopeDiameter``, ``TechnicalFoV``,          |
|                         |         |integer |``DMHeights`` and ``DMPitchs``.* Number of actuator on the pupil diameter.|
|                         |         |        |Must be the same length as ``DmPitchs``.                                  |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|InfModel                 |No       |string  |*default: 'gaussian'*, DM influence function model. Supported values:     |
|                         |         |        |``'gaussian'`` or ``'xinetics'``.                                         |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|InfCoupling              |No       |list of |*default: [0.2]*, Mechanical coupling coefficient (0–1) between DM        |
|                         |         |float   |actuators. Controls the width of the influence function. Must have the    |
|                         |         |        |same length as ``NumberActuators`` (one value per DM).                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|DmHeights                |No/Yes if|list of |*default: [0.0]*, DM altitude in meters, Must be the same length as       |
|                         |LO or    |float   |NumberActuators and DmPitchs                                              |
|                         |multi DMs|        |                                                                          |
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
|NumberReconstructedLayers|No       |integer |*default: 10*, Number of reconstructed atmospheric layers for tomographic |
|                         |         |        |AO systems (multi-guide-star).                                            |
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
|LoopGain_HO              |No       |float   |*Default : 0.5*, High Order Loop gain.                                    |
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

[COMPUTATION]
-------------

.. note::

   This section is optional, if this section is not present the defaul values are used.

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|platform                 |No       |string  |*default: 'GPU'* Set to it to 'CPU' to forcy the library to use numpy     |
|                         |         |        |instead of cupy.                                                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|integralDiscretization1  |No       |float   |*default: 1000.0*, Discretization used in the integrals                   |
|                         |         |        |(astro-tiptop/SEEING library).                                            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|integralDiscretization2  |No       |float   |*default: 4000*, Discretization used in the integrals                     |
|                         |         |        |(astro-tiptop/SEEING library).                                            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
