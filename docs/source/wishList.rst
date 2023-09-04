Wish List
=========

* Fix the behaviour of the pixelScale/FieldOfView/numberoflenslet/sizeoflenslet
* adding the noise Variance computation for the pyramid in High Order
* Computation for low Orders for the pyramid. 
* (Done) Adding the center of gravity algorithms
* (Done) for sensor_LO change the NumberPhotons to nph/frame/subaperture to be consistent with HO.
* Fix the FieldOfView/PixelScale/NumberOfLenslet behaviour to something that makes sense.
* Discuss with Cedric the use of FieldOfView VS  WindowRadiusWCoG for LO sensor. Also the center of gravity algorithms
* (Done: fixed now optional) NewValueThrPix: is this a really needed parameter or could we not parametrize this? 
* minimal parameter file example for configuring each of the type of AO
* AoArea other shapes than 'circle' and 'square'
* Zenith and Azimuth: in ``[sources_science]`` and ``[sources_LO]`` they are required but optional in ``[sources_LO]``. Make this consistent : either all required or all optional. 
* If an option in the parameter file is missing please choose between the custom message and the basic Python ``NoOptionError``. ALC would prefert the basic ``NoOptionError`` as it lets python deal with error management and lead to less confusing error messages (see below)
* Is it normal that in ``[sensor_HO]`` , ``NumberLenslet`` ``SizeLenslet`` is optionnal? (similar question for ``[sensor_LO]``, and for the entire ``[RTC]`` section)
* in aoSystem.py line 707 to 720, The code figures out what is the AO mode (SCAO, LTAO... ) Why is this not simply a parameter?





The matlab code works in spatial frequency domain.
pyramid from christoph verinaud, pyr responce goes from 0 to 1 
find 
umod line 279 fourrierModel.py
try to plot Sx in fourrier model. 
for the pyramid only use the high order part

test high order system 
extract the error breakdown from tiptop





https://arxiv.org/pdf/2003.06011.pdf

New 
rod conan 
julia shatokina

Comment on some behaviour
=========================

Confusing section missing error messages
----------------------------------------

There is not verification that the section atmosphere exist in the params.ini. If this section is absent this is the error message:

.. code-block:: python

   %%%%%%%% ERROR %%%%%%%%
   You must provide a value for the atmosphere seeing
   
   Traceback (most recent call last):
   
     File ~\simulations\python\TIPTOP\TIPTOP-EXAMPLE.py:25 in <module>
       overallSimulation(fullInputPath, parametersFile, baseOutpuPath, outputFile,
   
     File ~\simulations\python\TIPTOP\tiptop\tiptop.py:95 in overallSimulation
       fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False,
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:159 in __init__
       self.bounds = self.defineBounds()
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:183 in defineBounds
       bounds_down += list(np.zeros(self.ao.src.nSrc))
   
   AttributeError: 'aoSystem' object has no attribute 'src'


Consider changing this behaviour?

There is no verification that source_HO section is present in the params.ini. In case of absence this is the error message:

.. code-block:: python

   %%%%%%%% ERROR %%%%%%%%
   You must provide a value for the wavelength of the HO source (sources_HO)
   
   Traceback (most recent call last):
   
     File ~\simulations\python\TIPTOP\TIPTOP-EXAMPLE.py:25 in <module>
       overallSimulation(fullInputPath, parametersFile, baseOutpuPath, outputFile,
   
     File ~\simulations\python\TIPTOP\tiptop\tiptop.py:95 in overallSimulation
       fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False,
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:74 in __init__
       self.ao = aoSystem(path_ini,path_root=path_root,getPSDatNGSpositions=getPSDatNGSpositions)
   
   TypeError: __init__() should return None, not 'int'

Same for sensor science:

.. code-block:: python

   %%%%%%%% ERROR %%%%%%%%
   You must provide a value for the science detector (sensor_science) pixel scale
   
   Traceback (most recent call last):
   
     File ~\simulations\python\TIPTOP\TIPTOP-EXAMPLE.py:25 in <module>
       overallSimulation(fullInputPath, parametersFile, baseOutpuPath, outputFile,
   
     File ~\simulations\python\TIPTOP\tiptop\tiptop.py:95 in overallSimulation
       fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False,
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:159 in __init__
       self.bounds = self.defineBounds()
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:186 in defineBounds
       bounds_down += list(-self.freq.nPix//2 * np.ones(2*self.ao.src.nSrc))
   
   AttributeError: 'fourierModel' object has no attribute 'freq'

Same for sensor_HO:

.. code-block:: python

   %%%%%%%% ERROR %%%%%%%%
   You must provide a value for the HO detector (sensor_HO) pixel scale
   
   Traceback (most recent call last):
   
     File ~\simulations\python\TIPTOP\TIPTOP-EXAMPLE.py:25 in <module>
       overallSimulation(fullInputPath, parametersFile, baseOutpuPath, outputFile,
   
     File ~\simulations\python\TIPTOP\tiptop\tiptop.py:95 in overallSimulation
       fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False,
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:159 in __init__
       self.bounds = self.defineBounds()
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:186 in defineBounds
       bounds_down += list(-self.freq.nPix//2 * np.ones(2*self.ao.src.nSrc))
   
   AttributeError: 'fourierModel' object has no attribute 'freq'

And for the DM section:

.. code-block:: python 

   %%%%%%%% ERROR %%%%%%%%
   You must provide a value for the Dm number of actuators (NumberActuators)
   
   Traceback (most recent call last):
   
     File ~\simulations\python\TIPTOP\TIPTOP-EXAMPLE.py:25 in <module>
       overallSimulation(fullInputPath, parametersFile, baseOutpuPath, outputFile,
   
     File ~\simulations\python\TIPTOP\tiptop\tiptop.py:95 in overallSimulation
       fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False,
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:159 in __init__
       self.bounds = self.defineBounds()
   
     File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:186 in defineBounds
       bounds_down += list(-self.freq.nPix//2 * np.ones(2*self.ao.src.nSrc))
   
   AttributeError: 'fourierModel' object has no attribute 'freq'

confusing messages when parameters are missing
----------------------------------------------

If ``resolution`` is missing in section ``[telescope]`` , the error message is 
confusing.

.. code-block:: python

    %%%%%%%% ERROR %%%%%%%%
    You must provide a value for the pupil (telescope) resolution
    
    Traceback (most recent call last):
    
      File ~\simulations\python\TIPTOP\TIPTOP-EXAMPLE.py:25 in <module>
        overallSimulation(fullInputPath, parametersFile, baseOutpuPath, outputFile,
    
      File ~\simulations\python\TIPTOP\tiptop\tiptop.py:95 in overallSimulation
        fao = fourierModel(fullPathFilename, calcPSF=False, verbose=False,
    
      File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:159 in __init__
        self.bounds = self.defineBounds()
    
      File ~\simulations\python\TIPTOP\P3\aoSystem\fourierModel.py:183 in defineBounds
        bounds_down += list(np.zeros(self.ao.src.nSrc))
    
    AttributeError: 'aoSystem' object has no attribute 'src'

The following have the same confusing error message:

* [atmosphere] seeing
* [atmosphere] Cn2Weights
* [sources_HO] Wavelength
* [sensor_science] PixelScale
* [sensor_science] FieldOfView
* [sensor_HO] PixelScale
* [sensor_HO] FieldOfView
* [sensor_HO] NumberPhotons
* [sensor_LO] PixelScale
* [sensor_LO] FieldOfView




