Usage
=====


.. _installation:

Installation
------------

We recomend conda.

Create a Anaconda environment for TIPTOP with python 3.X (below for example 3.8):

.. code-block:: console

   conda create -–name tiptop python=3.8
   conda activate tiptop

Install cupy:

.. code-block:: console

   conda install -c conda-forge cupy

Install other libraries:

.. code-block:: console

   conda install ipython matplotlib scipy astropy sympy


Optionally you can get jupyter and jupiter-lab (many files are Jupyter Notebooks):

.. code-block:: console

   conda install jupyter
   conda install conda-forge jupyterlab


Get the library
^^^^^^^^^^^^^^^

There are two ways clone TIPTOP:

1. with sub-modules:

.. code-block:: console

   git clone --recurse-submodules --remote-submodules https://github.com/FabioRossiArcetri/TIPTOP.git

2. without sub-modules (cloned separately):

.. code-block:: console

   git clone https://github.com/FabioRossiArcetri/TIPTOP.git

In the second case please clone also the repos:

- https://github.com/FabioRossiArcetri/MASTSEL
- https://github.com/FabioRossiArcetri/SYMAO
- https://github.com/FabioRossiArcetri/SEEING
- https://github.com/oliviermartin-lam/P3/

as reported in the setup.py file.

Install with pip
^^^^^^^^^^^^^^^^

Navigate to the folder in which you have cloned TIPTOP, then navigate to the folder P3 and install P3 (remove --user to install to all users)

.. code-block:: console

   pip install -e --user .

Navigate out of P3, navigate inside MASTSEL and install MASTSEL:

.. code-block:: console

   pip install -e --user .

Navigate out of MASTSEL, navigate inside SEEING and install SEEING:

.. code-block:: console

   pip install -e --user .

Navigate out of SEEING, navigate inside of SYMAO and install SYMAO:

.. code-block:: console

   pip install -e --user .

Navigate out of SEEING and install TIPTOP:

.. code-block:: console

   pip install -e --user .


Quickstart
----------

To try execute the project you can use ``tiptop.overallSimulation()``

.. py:function:: tiptop.overallSimulation(path2param,paramFileName)

   return nothing but create a fits file containing the PSF
   
   :param path2param: required path to the parameter file
   :type path2param: str
   :param paramFileName: required name of the parameter file to be used without the extention
   :type paramFileName: str
   :return: nothing
   :rtype: None