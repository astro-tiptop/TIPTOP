Usage
=====

.. _installation:

Installation
------------

We recommend conda.

Create a Anaconda environment for TIPTOP with python>=3.9 (below for example
3.11)::

   conda create --name tiptop python=3.11
   conda activate tiptop

We only support python>=3.9.

if you want to beneficiate from GPU acceleration, and you have a GPU, install cupy::

   conda install -c conda-forge cupy

All other dependencies should be installed automatically wheather you choose to install through the pip packaging or by dowanloading the repository.
If you require a more sofisticated user interface than a command prompt, any IDE that uses iPython or python command prompt will work.

Example of IDE you could use the team has tested:

Jupyter and Jupyter-lab (many files are Jupyter Notebooks)::

   conda install jupyter
   conda install -c conda-forge jupyterlab

If you do not wish to use the provided jupyter notebooks you can convert them using the jupyter library. ::
   
   jupyter nbconvert --to python targetNotebook.ipynb

Spyder is another possibility::
   
   conda install spyder

Install from pypi
^^^^^^^^^^^^^^^^^

To install the last release of TIPTOP with its dependencies::

   pip install astro-tiptop

Install from git repo
^^^^^^^^^^^^^^^^^^^^^

First clone the repository::

   git clone https://github.com/astro-tiptop/TIPTOP.git

Navigate to the folder in which you have cloned TIPTOP and install it (remove
``--user`` to install to all users)::

   pip install -e --user .

Should you want to do your own developments, or fix bugs for us, you will need to download and install the following libraries.

- https://github.com/astro-tiptop/MASTSEL
- https://github.com/astro-tiptop/SYMAO
- https://github.com/astro-tiptop/SEEING
- https://github.com/astro-tiptop/P3/

Quickstart
----------
We recommand you first try your instalation of tiptop with the provided example TIPTOP-EXAMPLE.py or TIPTOP-EXAMPLE.ipynb if you use jupyter
If you installed from the repository, open a command prompt (you need the anaconda power shell in windows) and navigate to the repository ::
   
   conda activate tiptop
   python TIPTOP-EXAMPLE.py

If it executes without an error message your installation was sucessfull.
If you installed from pypi either you download the example script from the repository with the example parameter files or you use the following minimal working example.
First you need a minimum working parameter file.
Here is some example code to use it:

First the parameter file if used in the .ini format should look like this::

    [telescope]
    TelescopeDiameter=8.
    Resolution = 128
    
    [atmosphere]
    Seeing = 0.6
    
    [sources_science]
    Wavelength = [1.6e-6]
    Zenith = [0.]
    Azimuth = [0.]
    
    [sources_HO]
    Wavelength = 750e-9
    
    [sensor_science]
    PixelScale = 40
    FieldOfView = 256 
    
    [sensor_HO]
    PixelScale = 832
    FieldOfView = 6
    NumberPhotons = [200.]
    SigmaRON = 0.
    
    [DM]
    NumberActuators = [20]
    DmPitchs = [0.25]

You can simply copy paste the above in a text file which you rename minimalPar.ini.
To run a simulation with these parameter run the following script in the same folder as your parameter file.

.. code-block::

    from tiptop.tiptop import *
    plt.ion()
    
    overallSimulation("./", "minimalPar", './', 'test')

