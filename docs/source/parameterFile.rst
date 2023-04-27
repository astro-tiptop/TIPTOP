Parameter files explained
=========================

Introduction
------------

To run tiptop you need two things: to execute the function and to create a parameter file. This section explaines
what should your parameter file contain and the various parameters you can set. You can already find example parameter 
files in the `github <https://github.com/FabioRossiArcetri/TIPTOP/tree/main/perfTest>`_ .

There is support for two different parameter file format : the '.ini' and the '.yml'. Either format is structured the same way, only the syntax changes. For reference on the syntaxe refert to the `configparser <https://docs.python.org/3/library/configparser.html>`_ documentation for the .ini. For the .yml refert to the `pyYAML <https://pyyaml.org/wiki/PyYAMLDocumentation>`_ documentation.

The parameter files are divided in sections and they can contain multiple parameter. It is very important that each 
parameter be placed in the appropriate section. The order of the sections in the file does not matter.

Inles specified otherwise the following list of section must appear in the parameter file.

telescope
^^^^^^^^^

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

atmosphere
^^^^^^^^^^

.. option:: Seeing

   **Required**, 
   *type : float*, 
   Set the seeing at Zenith in arcsec. Might need fixing. seeing looks like it is converted to R0
   note : has a confusing error message if absent



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
   set the atmospere Fried parameter. Might need fixing. 
   temporary use: ``R0_Value`` =0 and the seeing to the wanted value.
   If not set it is calculated with ``seeing`` .

.. option:: testWindspeed

   **Optionnal**, 
   *type : float*, 
   test variable that should not be here?
   Required at the moment, for LO (which requires it) because the LO requires the average windspeed. If set to zero takes the average of WindSpeed weighted by the Cn2Weights. 
   Update 19/09/2022: no longer required? it does not output a message when removed. Completely ignored if not set.


sources_science
^^^^^^^^^^^^^^^

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


sources_HO
^^^^^^^^^^

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
   altitude of the guide stars (0 if infinite). Consider that all guide star are at the same heigh.


sources_LO
^^^^^^^^^^
.. note::

   This section is completely optional


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

sensor_science
^^^^^^^^^^^^^^

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

    We have examples of the following parameters being set but we do not understand if they are used in the science sensor case. (none of these are used for the science detector, initialised for one specific class but has no effect.). 
    Binning = 1
    NumberPhotons = [1500]
    SpotFWHM = [[0.0, 0.0, 0.0]]
    SpectralBandwidth = 0
    Transmittance = [1.0]
    Dispersion = [[0.0],[0.0]]
    SigmaRON = [0.1]
    Dark = 0.0
    SkyBackground = 0.0
    Gain = 1.0
    ExcessNoiseFactor = 1.0
    Wavelength = [0.55e-06] if the wavelength is specified here it is ignored. 
    FieldOfView = 1024


sensor_HO
^^^^^^^^^

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
   Number of pixels on the detector. 
   TODOI : change this behaviour as it makes no sense. Guido found that this is divided by `NumberLenslet`. Used for the noise. 
   Warning: gives a confusing error message if missing 

.. option:: Binning
   
   **Optional**, 
   *type: int*, 
   *default: 1*, 
   Binning factor of the detector, only used in the pyramid case, optionnal for pyramid

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
   Not used
   High Order spot scale in mili arcsec. Why list of list and why three?

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
   Dispersion x/y at the considered wavelength in pixel. Must be the same size than ``Transmittance``. Chromatic dispertion for PSF computation only.  In HarmoniSCAO_1 first the default and the thing given are not even the same shape but on top the default breaks the must be the same size as the transmitance... Also sorry for my ignorance: dispersion of what? Isn't this maybe redundant with `SpotFWHM` ?



.. option:: Gain 
   
   **Optional**, 
   *type: float*, 
   *default:1.0*, 
   Pixel gain. do you mean camera gain?

.. option:: ExcessNoiseFactor
   
   **Optional**, 
   *type: float*, 
   *default: 2.0*,
   excess noise factor. Why the hell would you by default have excess noise?

.. option:: NoiseVariance

   **Optional**, 
   *type: unknown*, 
   *Default : None*?, 
   Noise Variance in rd2. If not empty, this value overwrites the analytical noise variance calculation.







Shack-Hartmann requirement
""""""""""""""""""""""""""

.. option:: SigmaRON 

   **Required?**, 
   *type: float*, 
   read-out noise std in [e-], used only if the `NoiseVariance` is not set. 
   Note: this is optionnal if the ``WfsType`` == ``'Pyramid'``

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

.. option:: Algorithm
   
   **not used**, 
   *type: str*, 
   *defaut:'wcog'*, 
   not used, even if it actually has a default value. The default seems to be normal center of gravity.
   CoG computation algorithm, What are the other options?

    
.. option:: WindowRadiusWCoG 
   
   **not used**, 
   *type: int?*, 
   *default: 2*, 
   Number of pixels for ?windiwing? the low order WFS pixels

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


sensor_LO
^^^^^^^^^

.. note::

   This section is optional


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


DM
^^

.. option:: NumberActuators

   **Required**, 
   *type: list of int*, 
   Number of actuator on the pupil diameter. why a list of int?
   Must be the same length as DmPitchs
   Warning: gives a confusing error message if missing

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

RTC
^^^

.. option:: LoopGain_HO

   **Optional**, 
   *Type : float*, 
   *Default : 0.5*, 
   High Order Loop gain

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
   *type: float*?, 
   *default: None*,
   not used, auto matically optimized by tiptop.
   Low Order loop gain

.. option:: SensorFrameRate_LO

   **Required**, 
   *type: float*, 
   *default: None*,
   Loop frequency in [Hz]
   This is confusing : this is not optional if the ``[sensor_LO]`` is set.  

.. option:: LoopDelaySteps_LO

   **Optional**, 
   *type: int*, 
   *default: None*,
   Low Order loop delays in [frames]

.. option:: ResidualError

   **Optional**
   *Type : ?*
   *Default: None*
   ?



