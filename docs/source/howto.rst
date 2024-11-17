How to set up ...
=================

The telescope
-------------

This section contains the parameters related to the telescopen::

   [telescope]
   TelescopeDiameter = 8.
   ZenithAngle = 30.0
   ObscurationRatio = 0.16
   Resolution = 128
   TechnicalFoV = 120

and to some general aspects of the simulatio::

   PathPupil = 'data/EELT480pp0.0803m_obs0.283_spider2023.fits'
   PathApodizer = ''
   PathStatModes = ''
   PupilAngle = 0.0
   TechnicalFoV = 160.0
   # extra error in the science FoV (error not included in TIPTOP like NCPA, optics quality, ...)
   extraErrorNm = 160
   extraErrorExp = -2
   extraErrorMin = 0
   # extra error in the technical FoV  (error not included in TIPTOP like NCPA, optics quality, ...)
   extraErrorLoNm = 132
   extraErrorLoExp = -2
   extraErrorLoMin = 0
   # ELT tip & tilt wind shake when wind speed on M2 is 8 m/s
   windPsdFile = 'data/morfeo_windshake8ms_psd_2022_1k.fits'
   # jitter_FWHM --> 10 nm RMS tip error is 0.505arcesc
   # extra error on tip/tilt 70 nm (3.5) to consider tilt error due to aliasing not included by TIPTOP 
   jitter_FWHM = 3.5
   # ground layer focus is controlled with NGS WFS
   glFocusOnNGS = True

The atmosphere
--------------

This section contains the atmospheric parameters::

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

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~

.. _MCAO:

Multi Conjugate Adaptive Optics
-------------------------------

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~


.. _LTAO:

Laser Tomography Adaptive Optics
--------------------------------

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~

.. _GLAO:

Gound Layer Adaptive Optics
---------------------------

The wavefront sensor
~~~~~~~~~~~~~~~~~~~~

The deformable mirror
~~~~~~~~~~~~~~~~~~~~~

The real time controler
~~~~~~~~~~~~~~~~~~~~~~~



