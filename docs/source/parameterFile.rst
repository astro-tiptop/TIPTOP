Parameter files explained
=========================

Introduction
------------

To run tiptop you need two things: to execute the function and to create a parameter file. This section explaines
what should your parameter file contain and the various parameters you can set. You can already find example parameter 
files in the `github <https://github.com/FabioRossiArcetri/TIPTOP/tree/main/perfTest>`_ .


The parameter files are divided in sections and they can contain multiple parameter. It is very important that each 
parameter be placed in the appropriate section. The section starts with its name between '[]' and ends either with 
the end of file or with the next section. The order they are placed in the file does not matter.

Overview
--------
The following table resume the what the parameter file should contain and what is mandatory.

+---------------+------------+
|section        | required?  |
+===============+============+
|telescope      | Yes        |
+---------------+------------+
|atmosphere     | Yes        |
+---------------+------------+
|sources_science| Yes        |
+---------------+------------+
|source_HO      | Yes        |
+---------------+------------+
|sources_LO     | No         |
+---------------+------------+
|sensor_science | Yes        |
+---------------+------------+
|sensor_HO      | Yes        |
+---------------+------------+
|sensor_LO      | No         |
+---------------+------------+
|DM             | Yes        |
+---------------+------------+
|RTC            | No         |
+---------------+------------+
|Computation    | No         |
+---------------+------------+


We now go more in detail for each section:

.. warning::

   In the current state, if the parameter file is not found, it raises ``NoSectionError: No section: 'telescope'``.

[telescope]
-----------
+-------------------------+---------+-------+--------------------------------------------------------------------------+
| Parameter               | Required| Type  | Description                                                              |
+=========================+=========+=======+==========================================================================+
|TelescopeDiameter        |Yes      |float  |Set the outer diameter of the telescope pupil in unit of meters.          |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|Resolution               |Yes      |integer|Number of pixels across the pupil diameter                                |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|ObscurationRatio         |No       |float  |*Default : 0.0*, Defines the central obstruction                          |
|                         |         |       |due to the secondary as a ratio of the TelescopeDiameter                  |
|                         |         |       |*Warning* : MavisLO.py does not have a default value for this parameter   |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|ZenithAngle              |No       |float  |*Default: 0.0*, Set the pointing direction of the telescope in degree     |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|PupilAngle               |No       |float  |*default : 0.0*, unknown effect                                           |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|PathPupil                |No       |string |*default: ''*, path to the pupil model in .fits file (if provided,        |
|                         |         |       |the pupil model is interpolated).if absent or '', not used                |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|PathStaticOn             |No       |string |*default: None*, path to a map of static aberrations (nm) in              |
|                         |         |       | .fits file. if absent or '', not used                                    |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|PathStaticOff            |No       |string |*default: None*, No clue what this does. if absent or '', not used        |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|PathStaticPos            |No       |string |*default: None*, No clue                                                  |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|PathApodizer             |No       |string |*default: ''*, Path to a fits file that contain a binary map corresponding|
|                         |         |       |to a pupil apodizer (TBC). if absent or '', not used                      |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|PathStatModes            |No       |string |*default: ''*, path to a .fits file that contain a cube of map of mode    |
|                         |         |       |in amplitude which lead to a rms of 1 in nanometer of static aberation.   |
|                         |         |       |if absent or '', not used. Unsure how this works.                         |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|coeficientOfTheStaticMode|not used |string |*default: ''*, place holder                                               |
|                         |         |       |(TBC) need to find how does the pathStatModes fits file work.             |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|extraErrorNm             |No       |float  |*default: 0*, nm RMS of the additional error to be added (an error that   |
|                         |         |       |is not otherwise considered)                                              |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|extraErrorExp            |No       |float  |*default: -2*, exponent of the power of spatial frequencies used to       |
|                         |         |       | generate the PSD associated with extraErrorNm                            |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|extraErrorMin            |No       |float  |*default: 0*, minimum spatial frequency for which PSD associated with     |
|                         |         |       | extraErrorNm is > 0                                                      |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|TechnicalFoV             |Yes if LO|float  |*default: ??*, set the size of the technical field of view                |
|                         |         |       |   *Warning* : This is not optional in MavisLO.py                         |
+-------------------------+---------+-------+--------------------------------------------------------------------------+



[atmosphere]
------------

+-------------------------+---------+-------+--------------------------------------------------------------------------+
| Parameter               | Required| Type  | Description                                                              |
+=========================+=========+=======+==========================================================================+
|Seeing                   |Yes      |float  |Set the seeing at Zenith in arcsec. If not set TIPTOP uses ``r0_value`` . |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|Wavelength               |Yes if LO|float  |*Default : 500e-9*, Wavelength of definition of the atmosphere statistics |
|                         |         |       |   Warning: not optional in MavisLO.py                                    |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|L0                       |Yes if LO|float  |*Default : 25.0*, Outer Scale of the atmosphere  in meters                |
|                         |         |       |Warning: not optionnal in MavisLO.py                                      |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|Cn2Weights               |No\*     |list of|*Default : [1.0]*, Relative contribution of each layer. The sum of all the|
|                         |         |float  |list element must be 1. Must have the same length as ``Cn2Heights``,      |
|                         |         |       |``WindSpeed`` and ``WindDirection``.                                      |
|                         |         |       |\*Warning : required if ``Cn2Heights``, ``WindSpeed`` or ``WindDirection``|
|                         |         |       |are defined                                                               |
|                         |         |       |Warning : extremly confusing error message if absent when it must be      |
|                         |         |       |defined                                                                   |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|Cn2Heights               |No\*     |list of|*Default : [0.0]*, altitude of layers in meters.                          |
|                         |         |float  |Must have the same length as ``Cn2Weights``, ``WindSpeed`` and            |
|                         |         |       |``WindDirection``.                                                        |
|                         |         |       |\*Warning : required if ``Cn2Weights``, ``WindSpeed`` or ``WindDirection``|
|                         |         |       |are defined                                                               |
|                         |         |       |Warning : extremly confusing error message if absent when it must be      |
|                         |         |       |defined                                                                   |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|WindSpeed                |No\*     |list of|*Default : [10.0]*, Wind speed values for each layer in m/s.              |
|                         |         |float  |Must have the same length as ``Cn2Weights``, ``Cn2Heights`` and           |
|                         |         |       |``WindDirection``.                                                        |
|                         |         |       |\*Warning : required if ``Cn2Weights``, ``Cn2Heights`` or                 |
|                         |         |       |``WindDirection`` are defined                                             |
|                         |         |       |Warning : extremly confusing error message if absent when it must be      |
|                         |         |       |defined                                                                   |
+-------------------------+---------+-------+--------------------------------------------------------------------------+
|WindDirection            |No\*     |list of|*Default : [0.0]*, wind direction for each layer in degrees. 0 degree is  |
|                         |         |float  |?? then anticlockwise.                                                    |
|                         |         |       |Must have the same length as ``Cn2Weights``, ``Cn2Heights`` and           |
|                         |         |       |``WindSpeed``.                                                            |
|                         |         |       |\*Warning : required if ``Cn2Weights``, ``Cn2Heights`` or ``WindSpeed``   |
|                         |         |       |are defined                                                               |
|                         |         |       |Warning : extremly confusing error message if absent when it must be      |
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
|Zenith                   |Yes      |list of |Zenithal coordinate in arcsec of Wavelength sources given in              |
|                         |         |float   |``Wavelength``. Must be the same length as ``Azimuth``                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|Azimuth                  |Yes      |list of |Azimuthal coordinate in degree of Wavelength sources given in             |
|                         |         |float   |``Wavelength``. Must be the same length as ``Zenith``                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

[sources_HO]
------------

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|Wavelength               |Yes      |float   |Sensing wavelength for Hight Order modes in meters,                       |
|                         |         |        |Warning : gives a confusing error message if absent                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Zenith                   |No       |list of |*Default : [0.0]*, Zenithal coordinate of each guide stars in arcsec.     |
|                         |         |float   |Must be the same length as ``Azimuth``, Even if ``Azimutal`` is defined,  |
|                         |         |        |this is optionnal.                                                        |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Azimuth                  |No       |list of |*Default : [0.0]*, Azimuthal coordinate in degree of each guide stars.    |
|                         |         |float   |Must be the same length as ``Zenith``, even if ``Zenith`` is defined,     |
|                         |         |        |this is optionnal.                                                        |
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
|Zenith                   |Yes      |list of |Zenithal coordinate of each guide stars in arcsec.                        |
|                         |         |float   |Must be the same length as ``Azimuth``                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Azimuth                  |Yes      |list of |Azimuthal coordinate in degree of each guide stars.                       |
|                         |         |float   |Must be the same length as ``Zenith``                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   

[sensor_science]
----------------

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|PixelScale               |Yes      |float   |Pixel/spaxel scale in mili arcsec.                                        |
|                         |         |        |Warning: confusing error message if missing                               |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |float   |Field of view of the camera in pixel/spaxel.                              |
|                         |         |        |Warning: confusing error massage if missing                               |
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
|PixelScale               |Yes      |integer |High Order WFS pixel scale in [mas], unclear what are the units if we     |
|                         |         |        |chose a pyramid wavefront sensor.                                         |
|                         |         |        |Warning: gives a confusing error message if missing                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |integer |Number of pixels per subaperture.                                         |
|                         |         |        |Warning: gives a confusing error message if missing                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WfsType                  |No       |string  |*default : 'Shack-Hartmann'*, type of wavefront sensor used for the High  |
|                         |         |        |Order sensing. Other available option: 'Pyramid'                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberPhotons            |Yes      |list of |Flux return in [nph/frame/subaperture]                                    |
|                         |         |integer |Warning: confusing error message if missing                               |
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
|Gain                     |No       |float   |*default:1.0*, Pixel gain. do you mean camera gain or loop goin?          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ExcessNoiseFactor        |No       |float   |*default: 2.0*, excess noise factor. TODO: default should be 1            |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NoiseVariance            |No       |unknown |*Default : None*?, Noise Variance in rad2. If not empty, this value       |
|                         |         |        |overwrites the analytical noise variance calculation.                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+

In the two following section we list the parameters that are specific to each wavefront sensor. If you define a parameter 
for one WFS while another WFS is defined The parameter will be ignired. For example, if you define the parameter SigmaRON,
while WfsType is 'Pyramid', SigmaRON is ignored.

Shack-Hartmann requirement
^^^^^^^^^^^^^^^^^^^^^^^^^^

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|SigmaRON                 |Yes      |float   |read-out noise std in [e-], used only if the `NoiseVariance` is not set.  |
|                         |         |        |Note: this is optional if the ``WfsType`` == ``'Pyramid'``                |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
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
|                         |         |        |Warning : gives a confusing message if missing when required              |
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
|                         |         |        |Warning: gives a confusing error message if missing                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|FieldOfView              |Yes      |integer |not used. Number of pixels per subaperture,                               |
|                         |         |        |Warning: gives a confusing error message if missing                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberPhotons            |Yes      |list of |Detected flux in [nph/frame/subaperture], Must be the same length as      |
|                         |         |integer |NumberLenslet                                                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NumberLenslets           |No       |list of |*Default : [1]*, number of WFS lenslets, Must be the same length as       |
|                         |         |integer |NumberPhotons                                                             |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SigmaRON                 |No       |float   |*default: 0.0*, read out noise in [e-]                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|Dark                     |No       |float   |*default: 0.0*, dark current[e-/s/pix]                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SkyBackground            |No       |float   |*default: 0.0*, Sky background [e-/s/pix]                                 |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|ExcessNoiseFactor        |No       |float   |*default: 2.0*, excess noise factor                                       |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|WindowRadiusWCoG         |No       |integer |*default: 1*,2 Used instead of field of view, Number of pixels for        |
|                         |         |        |windiwing the low order WFS pixels                                        |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|ThresholdWCoG            |No       |float   |*default: 0.0*, Threshold Number of pixels for windowing the low order WFS|
|                         |         |        |pixels                                                                    |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|NewValueThrPix           |No       |float   |*default: 0.0*, New value for pixels lower than threshold.                |
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

[DM]
----

+-------------------------+---------+--------+--------------------------------------------------------------------------+
| Parameter               | Required| Type   | Description                                                              |
+=========================+=========+========+==========================================================================+
|NumberActuators          |Yes      |list of |Number of actuator on the pupil diameter. why a list of int? Must be the  |
|                         |         |integer |same length as DmPitchs. Warning: gives a confusing error message if      |
|                         |         |        |missing. Warning: not used in TIPTOP!                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|DmPitchs                 |Yes      |list of |DM actuators pitch in meters, on the meta pupil at the conjugasion        |
|                         |         |float   |altitude, used for fitting error computation. Must be the same length as  |
|                         |         |        |NumberActuators? Warning: gives a confusing error message if missing      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+    
|InfModel                 |No       |string  |*default: 'gaussian'*, DM influence function model. Not used in tiptop but| 
|                         |         |        |used in the psf reconstruction. What are the other possible one?          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|InfCoupling              |No       |list of |*default: [0.2]*, DM influence function model mechanical coupling. used in| 
|                         |         |float   |the psf reconstruction. Unclear what this does. Must be the same length as|
|                         |         |        |NumberActuators?                                                          |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|DmHeights                |No       |list of |*default: [0.0]*, DM altitude in meters, Must be the same length as       |
|                         |         |float   |NumberActuators and DmPitchs                                              |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|OptimizationZenith       |No       |float   |*default: [0.0]*, Zenith position in arcsec of the direction in which the |
|                         |         |        |AO correction is optimized.   Must be the same length as                  |
|                         |         |        |OptimisationAzimuth  and OptimizationWeight. These are for wide field AO  |
|                         |         |        |system, should be a requirement for MCAO and GLAO                         |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|OptimizationAzimuth      |No       |list of |*default: [0.0]*, Azimuth in degrees  of the direction in which the AO    |
|                         |         |float   |correction is optimized. Must be the same length as OptimizationZenith    |
|                         |         |        |and OptimizationWeight. These are for wide field AO system, should be a   |
|                         |         |        |requirement for MCAO and GLAO                                             |
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
|LoopGain_HO              |No       |float   |*Default : 0.5*, High Order Loop gain. Warning: if system to be simulated |
|                         |         |        |is a multi-conjugate system this parameter is not used.                   |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|SensorFrameRate_HO       |No       |float   |*Default : 500.0*, High Order loop frequency in [Hz]                      |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|LoopDelaySteps_HO        |No       |integer |*Default : 2*, High Order loop delay in [frame]                           |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|LoopGain_LO              |No       |float or|*default: None*, Low Order loop gain, Warning: if set to 'optimize', gain |
|                         |         |string  |is automatically optimized by tiptop, otherwise the float value set is    |
|                         |         |        |used.                                                                     |
+-------------------------+---------+--------+--------------------------------------------------------------------------+   
|SensorFrameRate_LO       |Yes      |float   |*default: None*, Loop frequency in [Hz]. If ``[sensor_LO]`` section is    |
|                         |         |        |present it must be set.                                                   |
+-------------------------+---------+--------+--------------------------------------------------------------------------+
|LoopDelaySteps_LO        |No       |integer |*default: None*, Low Order loop delays in [frames]. If ``[sensor_LO]``    |
|                         |         |        |section is present it must be set.                                        |
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
|                         |         |        |noise computation.                                                        |
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