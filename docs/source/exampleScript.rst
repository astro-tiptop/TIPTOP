Example Scripts
===============

We list here a example scripts or tips to acess informations

Minimum Working example
-----------------------

To get started with TIPTOP you need a parameter file that contains at minimum the following arguments


.. tabs::

    .. tab:: params.ini

        .. code-block:: python

            [telescope]
            TelescopeDiameter =8.
            Resolution = 128

            [atmosphere]
            Seeing = 0.6

            [sources_science]
            Wavelength = [1.6e-6]
            Zenith = [0.]
            Azimuth = [0.]

            [sources_HO]
            Wavelength = [750e-9]

            [sensor_science]
            PixelScale = 40
            FieldOfView = 256

            [sensor_HO]
            PixelScale = 832
            FieldOfView = 10

            [DM]
            NumberActuators = [20]
            DmPitchs = [0.25]


    .. tab:: params.yml

        .. code-block:: yaml

            telescope:
              TelescopeDiameter =8.
              Resolution = 128

            atmosphere:
              Seeing = 0.6

            sources_science:
              Wavelength = [1.6e-6]
              Zenith = [0.]
              Azimuth = [0.]

            sources_HO:
              Wavelength = [750e-9]

            sensor_science:
              PixelScale = 40
              FieldOfView = 256

            sensor_HO:
              PixelScale = 832
              FieldOfView = 10

            DM:
              NumberActuators = [20]
              DmPitchs = [0.25]

Then you need to create a .py script that contain as a minimum the following:

.. code-block:: python

    from tiptop import *
    #overallSimulation(<path to parameter file>,<parameter file name WITHOUT extension>,<path to save the result>,<name of the file to save>)
    overallSimulation("./", "params", './', 'test')

This script will run TIPTOP with the parameters in the parameter file while saving the instrument image(s) in the local file under the name test.fits
