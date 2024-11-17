How to set up ...
=================

The telescope
-------------

This section contains the parameters related to the telescopen::

   [telescope]
   TelescopeDiameter = 38.5
   ZenithAngle = 30.0
   ObscurationRatio = 0.28
   Resolution = 480

to the pupil::

   # path to the pupil model in .fits file - optional (if provided, the pupil model is interpolated) - default: ''
   PathPupil = 'data/EELT480pp0.0803m_obs0.283_spider2023.fits'
   PupilAngle = 0.0

to additional aberrations (in the main path or NCP)::

   PathStaticOn = '../P3/aoSystem/data/ELT_CALIBRATION/CombinedError_Wavefront_nm.fits'
   PathApodizer = ''
   PathStatModes = ''

or::

   # extra error in the science FoV (error not included in TIPTOP like NCPA, optics quality, ...)
   extraErrorNm = 160
   extraErrorExp = -2
   extraErrorMin = 0
   # extra error in the technical FoV  (error not included in TIPTOP like NCPA, optics quality, ...)
   extraErrorLoNm = 132
   extraErrorLoExp = -2
   extraErrorLoMin = 0

to windshake and additional tilt jitter::

   # ELT tip & tilt wind shake when wind speed on M2 is 8 m/s
   windPsdFile = 'data/morfeo_windshake8ms_psd_2022_1k.fits'
   # jitter_FWHM --> 10 nm RMS tip error is 0.505arcesc
   jitter_FWHM = 3.5

to the size of the technical field::

   TechnicalFoV = 160.0

and to the global focus control::

   # ground layer focus is controlled with NGS WFS
   glFocusOnNGS = True


The atmosphere
--------------

This section contains the atmospheric parameters, seeing (or r\ :sub:`0`\), L\ :sub:`0`\, C\ :sub:`n`\ :sup:`2`\  profile and wind profile::

   [atmosphere]
   Wavelength = 500e-9
   Seeing = 0.8
   L0 = 22.0
   Cn2Weights = [0.59, 0.02, 0.04, 0.06, 0.01, 0.05, 0.09, 0.04, 0.05, 0.05]
   Cn2Heights = [30, 140, 281, 562, 1125, 2250, 4500, 7750, 11000, 14000]
   WindSpeed = [6.6, 5.9, 5.1, 4.5, 5.1, 8.3, 16.3, 10.2, 14.3, 17.5]
   WindDirection = [0., 0., 0., 0., 90., -90., -90., 90., 0., 0.]

The adaptive optics system
--------------------------
Before diving into the different sections you need to figure what kind of AO system you want among the following:


* :ref:`Single Conjugate Adaptive Optics <SCAO>`

* :ref:`Multi Conjugate Adaptive Optics <MCAO>`

* :ref:`Laser Tomography Adaptive Optics <LTAO>`

* :ref:`Ground Layer Adaptive Optics <GLAO>`


.. _SCAO:

Single Conjugate Adaptive Optics
--------------------------------

The Single Conjugate Adaptive Optics system is described here: `ESO - AO MODES - SCAO <https://www.eso.org/sci/facilities/develop/ao/ao_modes/.html#scao>`_ 

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

The Wavefront Sensor in a SCAO system (natural guide star only) can be a Pyramid sensor::

   [sensor_HO]
   # WFS type - optional - default : Shack-Hartmann
   WfsType = 'Pyramid'
   # Spot modulation radius in lambda/D units for pyramid WFS - optional - default : None
   Modulation = 3
   # HO WFS pixel scale in [mas] - required
   PixelScale = 220      
   # Number of pixels per subaperture - required
   FieldOfView = 600         
   # Flux return in [nph/frame/subaperture] - required
   NumberPhotons = [500]                  
   # read-out noise std in [e-] - required
   SigmaRON = 1.0               
   # dark current[e-/s/pix] - optional - default: 0.0
   Dark = 0.2
   # Sky background [e-/s/pix] - optional - default: 0.0           
   SkyBackground = 0.6
   # excess noise factor - optional - default: 2.0                     
   ExcessNoiseFactor = 1.0 
   # Number of WFS lenslets - required
   NumberLenslets = [100]

or a Shack-Hartmann sensor::

   [sensor_HO]
   WfsType = 'Shack-Hartmann'
   Modulation = None
   PixelScale = 832
   FieldOfView = 6
   Binning = 1
   NumberPhotons = [100.0]
   SigmaRON = 0.2
   ExcessNoiseFactor = 2.0
   # CoG computation algorithm - optional  -defaut:'wcog'
   Algorithm = 'wcog' 
   # Number of pixels for windiwing the low order WFS pixels - optional - default: 2      
   WindowRadiusWCoG = 2
   NumberLenslets = [40]

Then a section with the guide star must be present::

   [sources_HO]
   # Sensing wavelength for HO modes in meters - required
   Wavelength = 950e-9
   # list of polar coordinates of the guide stars sources; zenith in arcsec and azimuth in degrees - optional - default [0.0]
   Zenith = [0.0]
   Azimuth = [0.0]
   # altitude of the guide stars (0 if infinite) - optional - default: 0.0
   Height = 0.0

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The deformable mirror is used to achieve wavefront control and correction and this section contains the following parameters::

   [DM]
   # DM actuators pitch in meters - required
   NumberActuators = [80]
   # DM actuators pitch in meters - required
   DmPitchs = [0.38]
   # DM influence function model - optional - default: 'gaussian'
   InfModel = 'gaussian'
   # DM influence function model  mechanical coupling- optional - default: [0.2]
   InfCoupling = [0.2]
   # DM altitude in m - optional - default: [0.0]
   DmHeights = [600.0] 
   # Zenith position in arcsec - optional - default: [0.0]
   OptimizationZenith = [0] 
   # Azimuth in degrees - optional - default: [0.0]
   OptimizationAzimuth = [0] 
   # Weights - optional - default: [1.0]
   OptimizationWeight  = [1]   
   # Matrix Conditioning - optional - default: 1e2
   OptimizationConditioning = 1.0e2 
   # Number of reconstructed layers for tomographic systems - optional - default: 10
   NumberReconstructedLayers= 10
   # Shape of the AO-corrected area - optional - default: 'circle'
   AoArea = 'circle'

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~

This section contains the details of the control, the framerate and the delay::

   [RTC]
   # HO Loop gain - required
   LoopGain_HO = 0.5                             
   # HO loop frequency in [Hz] - required
   SensorFrameRate_HO = 500.0
   # HO loop frame delay - required
   LoopDelaySteps_HO = 1

.. _MCAO:

Multi Conjugate Adaptive Optics
-------------------------------

The Multi Conjugate Adaptive Optics system is described here: `ESO - AO MODES - MCAO <https://www.eso.org/sci/facilities/develop/ao/ao_modes/.html#mcao>`_ 

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

For an MCAO system with multiple NAtural Guide Stars (NGSs), the difference here is that ``NumberPhotons`` and ``NumberLenslets`` are vectors as are ``Zenith`` and ``Azimuth`` in the ``soruces_HO`` section.

If the system have Laser Guide Stars (LGSs) and NGSs these sections are present::

   [sources_LO]
   Wavelength = [1650e-9]
   Zenith = [66.6, 79.3, 69.0]
   Azimuth = [221.7, 321.2, 106.6]
   
   [sensor_LO]
   PixelScale = 16.0
   FieldOfView = 100
   Binning = 1
   # zero magnitude flux 8.17e11ph/s (H band)
   # magnitudes 10.7, 16.3, 14.5
   # 2x2 sub-apertures and 250 Hz framerate
   # --> 8.17e11*10**(-[10.7,16.3,14.5]/2.5)/4/250.
   NumberPhotons = [42900,247,1300]
   SpotFWHM = [[0.0,0.0,0.0]]
   SigmaRON = 0.5
   Dark = 40.0
   SkyBackground = 120.0
   Gain = 1.0
   ExcessNoiseFactor = 1.3
   # note 2x2 is required to provide focus control
   # (see glFocusOnNGS in telescope section)
   NumberLenslets = [2, 2, 2]
   Algorithm = 'wcog'
   WindowRadiusWCoG = 'optimize'
   ThresholdWCoG = 0.0
   NewValueThrPix = 0.0
   noNoise = False
   filtZernikeCov = True

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

For an MCAO system with multiple DM, the difference here is that ``NumberActuators``, ``DmPitchs``, ``InfCoupling`` and ``DmHeights`` are vectors.

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~

If the system have Laser Guide Stars (LGSs) and NGSs this section has the following parameters::

   LoopGain_LO = 'optimize'
   SensorFrameRate_LO = 250.0
   LoopDelaySteps_LO = 1

.. _LTAO:

Laser Tomography Adaptive Optics
--------------------------------

The Laser Tomography Adaptive Optics system is described here: `ESO - AO MODES - LTAO <https://www.eso.org/sci/facilities/develop/ao/ao_modes/.html#ltao>`_ 

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

As for the MCAO system when LGSs and NGSs are present the sections ``[sources_LO]`` and ``[sensor_LO]`` must be added.

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

Like for a SCAO system.

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~

Like for a MCAO system with LGSs and NGSs.

.. _GLAO:

Gound Layer Adaptive Optics
---------------------------

The Ground Layery Adaptive Optics system is described here: `ESO - AO MODES - GLAO <https://www.eso.org/sci/facilities/develop/ao/ao_modes/.html#glao>`_ 

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

Like for a MCAO system.

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

Like for a SCAO system.

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~

Like for a MCAO system.


