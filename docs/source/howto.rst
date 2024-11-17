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

   ; path to the pupil model in .fits file - optional (if provided, the pupil model is interpolated) - default: ''
   PathPupil = 'data/EELT480pp0.0803m_obs0.283_spider2023.fits'
   PupilAngle = 0.0

to additional aberrations (in the main path or NCP) and other stuff::

   PathStaticOn = '/home/frossi/dev/TIPTOP/P3/aoSystem/data/ELT_CALIBRATION/CombinedError_Wavefront_nm.fits'
   PathApodizer = ''
   PathStatModes = ''
   # extra error in the science FoV (error not included in TIPTOP like NCPA, optics quality, ...)
   extraErrorNm = 160
   extraErrorExp = -2
   extraErrorMin = 0
   # extra error in the technical FoV  (error not included in TIPTOP like NCPA, optics quality, ...)
   # Total is 160 nm without 90 nm of contingency because contingency should be applied only once
   extraErrorLoNm = 132
   extraErrorLoExp = -2
   extraErrorLoMin = 0

to windshake and attional tilt jitter::

   # ELT tip & tilt wind shake when wind speed on M2 is 8 m/s
   windPsdFile = 'data/morfeo_windshake8ms_psd_2022_1k.fits'
   # jitter_FWHM --> 10 nm RMS tip error is 0.505arcesc
   # extra error on tip/tilt 70 nm (3.5) to consider tilt error due to aliasing not included by TIPTOP 
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
   ;WFS type - optional - default : Shack-Hartmann
   WfsType = 'Pyramid'
   ;Spot modulation radius in lambda/D units for pyramid WFS - optional - default : None
   Modulation = 3
   ;HO WFS pixel scale in [mas] - required
   PixelScale = 220      
   ;Number of pixels per subaperture - required
   FieldOfView = 600         
   ;Flux return in [nph/frame/subaperture] - required
   NumberPhotons = [500]                  
   ;read-out noise std in [e-] - required
   SigmaRON = 1.0               
   ; dark current[e-/s/pix] - optional - default: 0.0
   Dark = 0.2
   ;Sky background [e-/s/pix] - optional - default: 0.0           
   SkyBackground = 0.6
   ;excess noise factor - optional - default: 2.0                     
   ExcessNoiseFactor = 1.0 
   ;Number of WFS lenslets - required
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
   ;CoG computation algorithm - optional  -defaut:'wcog'
   Algorithm = 'wcog' 
   ;Number of pixels for windiwing the low order WFS pixels - optional - default: 2      
   WindowRadiusWCoG = 2
   NumberLenslets = [40]


The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~

.. _MCAO:

Multi Conjugate Adaptive Optics
-------------------------------

The Multi Conjugate Adaptive Optics system is described here: `ESO - AO MODES - MCAO <https://www.eso.org/sci/facilities/develop/ao/ao_modes/.html#mcao>`_ 

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~


.. _LTAO:

Laser Tomography Adaptive Optics
--------------------------------

The Laser Tomography Adaptive Optics system is described here: `ESO - AO MODES - LTAO <https://www.eso.org/sci/facilities/develop/ao/ao_modes/.html#ltao>`_ 

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~

.. _GLAO:

Gound Layer Adaptive Optics
---------------------------

The Ground Layery Adaptive Optics system is described here: `ESO - AO MODES - GLAO <https://www.eso.org/sci/facilities/develop/ao/ao_modes/.html#glao>`_ 

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~



