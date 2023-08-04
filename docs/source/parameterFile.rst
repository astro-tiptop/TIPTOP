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


The mandatory sections and their content are:

.. warning::

   In the current state, if the parameter file is not found, it raises ``NoSectionError: No section: 'telescope'``.

[telescope]
^^^^^^^^^^^

.. option:: TelescopeDiameter

   **Required**, 
   *type : float*, 
   Set the outer diameter of the telescope pupil in unit of meters.

.. option:: ObscurationRatio

   **Optionnal**, 
   *type : float*, 
   *Default : 0.0*,
   Defines the central obstruction due to the secondary as a ratio of the TelescopeDiameter
   *Warning* : MavisLO.py does not have a defualt value for this parameter 

.. option:: Resolution

   **Required**, 
   *type : int*, 
   Number of pixels across the pupil diameter

.. option:: ZenithAngle

   **Optional**, 
   *type : float*, 
   *Default: 0.0*, 
   Set the pointing direction of the telescope in degree


.. option:: TechnicalFoV

   **Optional**, 
   *type : float*, 
   *default: ??*, 
   set the size of the technical field of view
   *Warning* : This is not optional in MavisLO.py

.. option:: PathPupil

   **Optional**, 
   *type : str*, 
   *default: ''*, 
   path to the pupil model in .fits file (if provided, the pupil model is interpolated).if absent or '', not used

.. option:: PathStaticOn

   **Optional**, 
   *type : str*, 
   *default: None*, 
   path to a map of static aberrations (nm) in .fits file. if absent or '', not used

.. option:: PathStaticOff

   **Optional**, 
   *type : str*, 
   *default: None*, 
   No clue what this does. if absent or '', not used

.. option:: PathStaticPos

   **Optional**, 
   *type : str*, 
   *default: None*, 
   No clue

.. option::  PathApodizer

   **Optional**, 
   *type : str*, 
   *default: ''*, 
   Path to a fits file that contain a binary map corresponding to a pupil apodizer (TBC). if absent or '', not used

.. option:: PathStatModes
   
   **Optional**, 
   *type : str*, 
   *default: ''*, 
   path to a .fits file that contain a cube of map of mode in amplitude which lead to a rms of 1 in nanometer of static aberation. if absent or '', not used. Unsure how this works.

.. option:: coeficientOfTheStaticMode
   
   **not used**, 
   *type : str*, 
   *default: ''*, 
   place holder 
   (TBC) need to find how does the pathStatModes fits file work.

.. option:: extraErrorNm
   
   **Optional**, 
   *type : float*, 
   *default: 0*, 
   nm RMS of the additional error to be added (an error that is not otherwise considered)

.. option:: extraErrorExp
   
   **Optional**, 
   *type : float*, 
   *default: -2*, 
   exponent of the power of spatial frequencies used to generate the PSD associated with extraErrorNm

.. option:: extraErrorMin
   
   **Optional**, 
   *type : float*, 
   *default: 0*, 
   minimum spatial frequency for which PSD associated with extraErrorNm is > 0

[atmosphere]
^^^^^^^^^^^^

.. option:: Seeing

   **Required**, 
   *type : float*, 
   Set the seeing at Zenith in arcsec. 
   If not set TIPTOP uses ``r0_value`` .

.. option:: Wavelength

   **Optional**, 
   *type : float*, 
   *Default : 500e-9*, 
   Wavelength of definition of the atmosphere statistics
   Warning: not optional in MavisLO.py

.. option:: L0

   **Optional**, 
   *type : float*, 
   *Default : 25.0*, 
   Outer Scale of the atmosphere  in meters
   Warning: not optionnal in MavisLO.py

.. option:: Cn2Weights

   **Optional**, 
   *type : list of float*, 
   *Default : [1.0]*, 
   Relative contribution of each layer. The sum of all the list element must be 1. 
   Must have the same length as ``Cn2Heights``, ``WindSpeed`` and ``WindDirection``.
   Warning : required if ``Cn2Heights``, ``WindSpeed`` or ``WindDirection`` are defined
   Warning : extremly confusing error message if absent when it must be defined

.. option:: Cn2Heights

   **Optional**, 
   *type : list of float*, 
   *Default : [0.0]*, 
   altitude of layers in meters.
   Must have the same length as ``Cn2Weights``, ``WindSpeed`` and ``WindDirection``.
   Warning : required if ``Cn2Weights``, ``WindSpeed`` or ``WindDirection`` are defined
   Warning : extremly confusing error message if absent when it must be defined

.. option:: WindSpeed

   **Optional**, 
   *Type : list of float*, 
   *Default : [10.0]*, 
   Wind speed values for each layer in m/s. 
   Must have the same length as ``Cn2Weights``, ``Cn2Heights`` and ``WindDirection``.
   Warning : required if ``Cn2Weights``, ``Cn2Heights`` or ``WindDirection`` are defined
   Warning : extremly confusing error message if absent when it must be defined

.. option:: WindDirection

   **Optional**, 
   *Type : list of float*, 
   *Default : [0.0]*, 
   wind direction for each layer in degrees. 0 degree is ?? then anticlockwise.
   Must have the same length as ``Cn2Weights``, ``Cn2Heights`` and ``WindSpeed``.
   Warning : required if ``Cn2Weights``, ``Cn2Heights`` or ``WindSpeed`` are defined
   Warning : extremly confusing error message if absent when it must be defined

.. option:: r0_Value
   
   **Optionnal**, 
   *type : float*, 
   set the atmospere Fried parameter.
   If not set TIPTOP uses ``seeing`` .

.. option:: testWindspeed

   **Optionnal**, 
   *type : float*, 
   used only for tests

[sources_science]
^^^^^^^^^^^^^^^^^

.. option:: Wavelength

   **Required**, 
   *Type : list of float or float*, 
   list of central wavelengths for each frame in meters. you can have more than one science target. needs explaining why the science sources can be multiple. (polychromatic? several targets? you can set many taget of the same wavelength by only setting more than one Zenith and Azimuth but leaving the wavelength as a float. It produces one PSF per target. The number of PSF is the number of wavelength times the number of Azimuth/Zenith couple.

.. option:: Zenith

   **Required**, 
   *Type : list of float*, 
   Zenithal coordinate in arcsec of Wavelength sources given in ``Wavelength``.
   Must be the same length as ``Azimuth``

.. option:: Azimuth

   **Required**, 
   *Type : list of float*, 
   Azimuthal coordinate in degree of Wavelength sources given in ``Wavelength``.
   Must be the same length as ``Zenith``

[sources_HO]
^^^^^^^^^^^^

Typically the wavelength is the same for all guide star (at least in Laser guide star)

.. option:: Wavelength

   **Required**, 
   *type : float*, 
   Sensing wavelength for Hight Order modes in meters
   Warning : gives a confusing error message if absent

.. option:: Zenith

   **Optional**, 
   *Type : list of float*, 
   *Default : [0.0]*
   Zenithal coordinate of each guide stars in arcsec.
   Must be the same length as ``Azimuth``
   Even if ``Azimutal`` is defined, this is optionnal.
   
.. option:: Azimuth

   **Optional**, 
   *Type : list of float*, 
   *Default : [0.0]*
   Azimuthal coordinate in degree of each guide stars.
   Must be the same length as ``Zenith``
   Even if ``Zenith`` is defined, this is optionnal.

.. option:: Height
   
   **Optional**, 
   *Type : float*, 
   *Default : 0.0*, 
   altitude of the guide stars (0 if infinite). Consider that all guide star are at the same height.

[sources_LO]
^^^^^^^^^^^^
.. note::

   This section is completely optional (``[sensor_LO]`` section is required to have the LO part simulated)

.. option:: Wavelength

   **Required**, 
   *type : float*, 
   Sensing wavelength for Low Order modes in meters

.. option:: Zenith

   **Required**, 
   *Type : list of float*, 
   Zenithal coordinate of each guide stars in arcsec.
   Must be the same length as ``Azimuth``
   
.. option:: Azimuth

   **Required**, 
   *Type : list of float*, 
   Azimuthal coordinate in degree of each guide stars.
   Must be the same length as ``Zenith``

[sensor_science]
^^^^^^^^^^^^^^^^

.. option:: PixelScale

   **Required**, 
   *type : float*, 
   pixel/spaxel scale in mili arcsec
   Warning: confusing error message if missing


.. option:: FieldOfView

   **Required**, 
   *type : float*, 
   Field of view of the camera in pixel/spaxel. need confirmation on the optionality of this paramiter. 
   Warning: confusing error massage if missing

.. note::

    Following parameters was added to uniform all the sensor (HO and LO), but they are not used.

    Binning, NumberPhotons, SpotFWHM, SpectralBandwidth, Transmittance, Dispersion, SigmaRON, Dark, SkyBackground, Gain, ExcessNoiseFactor, Wavelength, FieldOfView

[sensor_HO]
^^^^^^^^^^^

Used regardless of the WFS, desired behaviour, 

.. option:: NumberLenslets

   **Optional**, 
   *type: list of int*, 
   *Default : [20]*
   Number of WFS lenslets. Used the same way in Shack-Hartmann wavefront sensor and Pyramid. Also used for noise computation if `NoiseVariance` is not set. 

.. option:: SizeLenslets                                                   
   
   **Optionnal**,
   *type: list of float*, 
   *Default: [Telescope] TelescopeDiameter/[sensor_HO] NumberLenslet*
   Size of WFS lenslets in meters. used, why a list of float? This overrides the ratio between telescope size and Number of lenslet used to compute the matrix size.

.. option:: PixelScale

   **Required**, 
   *type: int*, 
   High Order WFS pixel scale in [mas], unclear what are the units if we chose a pyramid wavefront sensor
   Warning: gives a confusing error message if missing 

.. option:: FieldOfView

   **Required**, 
   *type: int*, 
   Number of pixels per subaperture. 
   Warning: gives a confusing error message if missing 

.. option:: Binning
   
   **Optional**, 
   *type: int*, 
   *default: 1*, 
   Binning factor of the detector, only used in the pyramid case, optional for pyramid

.. option:: WfsType
   
   **Optional**, 
   *type: str*, 
   *default : 'Shack-Hartmann'*, 
   type of wavefront sensor used for the High Order sensing.
   Other available option: 'Pyramid'

.. option:: NumberPhotons  

   **Required**, 
   *type: list of int*, 
   Flux return in [nph/frame/subaperture]
   Warning: extremly confusing error message if missing

.. option:: SpotFWHM    
   
   **Optional**, 
   *type: list of list of float*, 
   *defaut: [[0.0, 0.0, 0.0]]*, 
   High Order spot parameters: two axes scale values in milliarcsec (only max value is used) and angle (angle is not used). Why list?

.. option:: SpectralBandwidth
   
   **Optional**, 
   *type: float*, 
   *default: 0.0*, 
   Not used
   Spectral bandwidth of the filter (imaging mode)? why specific to the imaging mode? what is the effect?

.. option:: Transmittance
   
   **Optional**, 
   *type: list of float*, 
   *default: [1.0]*, 
   Used for PSF computation and flux scaling but not with noise computation
   Transmittance at the considered wavelengths for polychromatic mode. How do you set polychromatic mode? Each element can not have a value superior to 1?

.. option:: Dispersion
   
   **Optional**, 
   *type: apparently list of list of float?*, 
   *default: [[0.0,0.0]]*, 
   Dispersion x/y at the considered wavelength in pixel. Must be the same size than ``Transmittance``. Chromatic dispertion for PSF computation only.
   In HarmoniSCAO_1 first the default and the thing given are not even the same shape but on top the default breaks the must be the same size as the transmitance...
   Also sorry for my ignorance: dispersion of what? Isn't this maybe redundant with `SpotFWHM` ?

.. option:: Gain 
   
   **Optional**, 
   *type: float*, 
   *default:1.0*, 
   Pixel gain. do you mean camera gain?

.. option:: ExcessNoiseFactor
   
   **Optional**, 
   *type: float*, 
   *default: 2.0*,
   excess noise factor.
   TODO: default should be 1

.. option:: NoiseVariance

   **Optional**, 
   *type: unknown*, 
   *Default : None*?, 
   Noise Variance in rad2. If not empty, this value overwrites the analytical noise variance calculation.

Shack-Hartmann requirement
""""""""""""""""""""""""""

.. option:: SigmaRON 

   **Required?**, 
   *type: float*, 
   read-out noise std in [e-], used only if the `NoiseVariance` is not set. 
   Note: this is optional if the ``WfsType`` == ``'Pyramid'``

.. option:: Algorithm
   
   **not used**, 
   *type: str*, 
   *defaut:'wcog'*, 
   other options: 'cog' (simple center-of-gravity), 'tcog' (center-of-gravity with threshold), 'qc' (quad-cell)
    
.. option:: WindowRadiusWCoG 
   
   **not used**, 
   *type: int?*, 
   *default: 2*, 
   FWHM in pixel of the gaussian weighting function

Pyramid requirement
"""""""""""""""""""

.. option:: Modulation
   
   **Required if WfsType == 'Pyramid'**, 
   *type: float*, 
   *default : None*, 
   If the chosen wavefront sensor is the ``'Pyramid'``, Spot modulation radius in lambda/D units. This is ignored if the WFS is `'Shack-Hartmann'`
   Warning : gives really confusing message if missing when required

Can be set but not used
"""""""""""""""""""""""

.. option:: Dark
   
   **not used**, 
   *type: float*, 
   *default: 0.0*, 
    
   dark current in [e-/s/pix]

.. option:: SkyBackground
   
   **not used**, 
   *type: float*, 
   *default: 0.0*, 
   Sky background [e-/s/pix]

.. option:: ThresholdWCoG
   
   **not used**, 
   *type: float?*, 
   *default: 0.0*, 
   Threshold Number of pixels for windowing the low order WFS pixels

.. option:: NewValueThrPix 
   
   **not used**, 
   *type: float*, 
   *default: 0.0*, 
   New value for pixels lower than `ThresholdWCoG`. Is there a reason to want to force these values to something else?

[sensor_LO]
^^^^^^^^^^^

.. note::

   This section is optional, if this section is not present only the HO part will be used (for ex. to simulate a SCAO NGS).

.. option:: PixelScale

   **Required**, 
   *type: float*, 
   LO WFS pixel scale in [mas]
   Warning: gives a confusing error message if missing

.. option:: FieldOfView 

   **Required**, 
   *type: int*, 
   not used. 
   Number of pixels per subaperture
   Warning: gives a confusing error message if missing

.. option:: NumberPhotons 

   **Required**, 
   *type: list of int*, 
   detected flux in [nph/frame/subaperture]
   Must be the same length as NumberLenslet

.. option:: NumberLenslets

   **Optional**, 
   *type: list of int*, 
   *Default : [1]*
   number of WFS lenslets
   Must be the same length as NumberPhotons

.. option:: SigmaRON   

   **Optional**, 
   *type: float*, 
   *default: 0.0*,
   read out noise in [e-]

.. option:: Dark

   **Optional**, 
   *type: float*, 
   *default: 0.0*,
   dark current[e-/s/pix]

.. option:: SkyBackground

   **Optional**, 
   *type: float*, 
   *default: 0.0*,
   sky background [e-/s/pix]

.. option:: ExcessNoiseFactor

   **Optional**, 
   *type: float*, 
   *default: 2.0*,
   excess noise factor

.. option:: WindowRadiusWCoG

   **Optional**, 
   *type: int*, 
   *default: 1*,2
   used instead of field of view
   Number of pixels for windiwing the low order WFS pixels
    
.. option:: ThresholdWCoG

   **Optional**, 
   *type: float*, 
   *default: 0.0*,
   Threshold Number of pixels for windowing the low order WFS pixels

.. option:: NewValueThrPix

   **Optional**, 
   *type: float*, 
   *default: 0.0*,
   New value for pixels lower than threshold.

Can be set but not used
"""""""""""""""""""""""

.. option:: Binning   

   **not used**, 
   *type: int*, 
   *default: 1*, 
   binning factor of the detector

.. option:: SpotFWHM   

   **not used**, 
   *type: list of list of int*, 
   *default: [[0.0, 0.0, 0.0]]*,
   Low Order spot scale in [mas]

.. option:: Gain

   **not used**, 
   *type: float*, 
   *default: 1.0*,
   camera gain

.. option:: Algorithm

   **not used**, 
   *type: str*, 
   *default: 'wcog'*,
   CoG computation algorithm

[DM]
^^^^

.. option:: NumberActuators

   **Required**, 
   *type: list of int*, 
   Number of actuator on the pupil diameter. why a list of int?
   Must be the same length as DmPitchs
   Warning: gives a confusing error message if missing
   Warning: not used in TIPTOP!

.. option:: DmPitchs

   **Required**, 
   *type: list of float*, 
   DM actuators pitch in meters, on the meta pupil at the conjugasion altitude, used for fitting error computation.
   Must be the same length as NumberActuators?
   Warning: gives a confusing error message if missing

.. option:: InfModel

   **Optional**, 
   *type: str*, 
   *default: 'gaussian'*,
   DM influence function model. Not used in tiptop but used in the psf reconstruction. What are the other possible one?

.. option:: InfCoupling

   **Optional**, 
   *type: list of float*, 
   *default: [0.2]*,
   DM influence function model mechanical coupling. used in the psf reconstruction. Unclear to me what this does.
   Must be the same length as NumberActuators?

.. option:: DmHeights 

   **Optional**, 
   *type: list of float*, 
   *default: [0.0]*,
   DM altitude in meters 
   Must be the same length as NumberActuators and DmPitchs

.. option:: OptimizationZenith

   **Optional**, 
   *type: float*, 
   *default: [0.0]*,
   Zenith position in arcsec of the direction in which the AO correction is optimized.
   Must be the same length as OptimisationAzimuth  and OptimizationWeight
   These are for wide field AO system, should be a requirement for MCAO and GLAO 

.. option:: OptimizationAzimuth

   **Optional**, 
   *type: list of float*, 
   *default: [0.0]*,
   Azimuth in degrees  of the direction in which the AO correction is optimized
   Must be the same length as OptimizationZenith and OptimizationWeight
   These are for wide field AO system, should be a requirement for MCAO and GLAO 

.. option:: OptimizationWeight

   **Optional**, 
   *type: float*, 
   *default: [1.0]*,
   Weights of the optimisation directions 
   Must be the same length as OptimizationZenith and OptimizationAzimuth
   These are for wide field AO system, should be a requirement for MCAO and GLAO 

.. option:: OptimizationConditioning

   **Optional**, 
   *type: float*, 
   *default: 1.0e2*,
   Matrix Conditioning threshold in the truncated SVD inversion. 

.. option:: NumberReconstructedLayers

   **Optional**, 
   *type: int*, 
   *default: 10*,
   Only used for wide field AO system, (meaning more than one guide star is defined)
   Number of reconstructed layers for tomographic systems. Shouldn't this be defaulted to 1 for SCAO sakes?

.. option:: AoArea

   **Optional**, 
   *type: str*, 
   *default: 'circle'*,
   Shape of the AO-corrected area. Any other options are not defined and will give a squarre correction area.  

[RTC]
^^^^^

.. option:: LoopGain_HO

   **Optional**, 
   *Type : float*, 
   *Default : 0.5*, 
   High Order Loop gain.
   Warning: if system to be simulated is a multi-conjugate system this parameter is not used.

.. option:: SensorFrameRate_HO

   **Optional**, 
   *type: float*, 
   *Default : 500.0*,
   High Order loop frequency in [Hz]

.. option:: LoopDelaySteps_HO

   **Optional**, 
   *type: int*, 
   *Default : 2*, 
   High Order loop delay in [frame]

.. option:: LoopGain_LO

   **Optional**, 
   *type: float or string*, 
   *default: None*,
   Low Order loop gain
   Warning: if set to 'optimize', gain is automatically optimized by tiptop, otherwise the float value set is used.

.. option:: SensorFrameRate_LO

   **Required**, 
   *type: float*, 
   *default: None*,
   Loop frequency in [Hz]
   If ``[sensor_LO]`` section is present it must be set.  

.. option:: LoopDelaySteps_LO

   **Optional**, 
   *type: int*, 
   *default: None*,
   Low Order loop delays in [frames]
   If ``[sensor_LO]`` section is present it must be set.

.. option:: ResidualError

   **Optional**
   *Type : ?*
   *Default: None*
   ?

[COMPUTATION]
^^^^^^^^^^^

.. note::

   This section is optional, if this section is not present the defaul values are used.

.. option:: simpleVarianceComputation

   **Optional**, 
   *type : str*, 
   Set to it to False to activate the more complex and slower MASTSEL LO noise computation.


.. option:: platform

   **Optional**, 
   *type : str*, 
   *default: 'GPU'*
   Set to it to 'CPU' to forcy the library to use numpy instead of cupy.

.. option:: integralDiscretization1

   **Optional**, 
   *type : float*, 
   *default: 1000*
   Discretization used in the integrals (astro-tiptop/SEEING library).

.. option:: integralDiscretization2

   **Optional**, 
   *type : float*, 
   *default: 4000*
   Discretization used in the integrals (astro-tiptop/SEEING library).
