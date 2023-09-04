Usage
=====

.. _installation:

Installation
------------

We recommend conda.

Create a Anaconda environment for TIPTOP with python>=3.9 (below for example
3.11)::

   conda create -–name tiptop python=3.11
   conda activate tiptop

Install cupy if you have a GPU::

   conda install -c conda-forge cupy

Install other libraries (this step can be skipped because it is managed by pip following libraries dependencies)::

   conda install ipython matplotlib scipy astropy sympy

Optionally you can get Jupyter and Jupyter-lab (many files are Jupyter
Notebooks)::

   conda install jupyter
   conda install -c conda-forge jupyterlab

Install with pip
^^^^^^^^^^^^^^^^

To install the last release of TIPTOP with its dependencies::

   pip install astro-tiptop

Install with git
^^^^^^^^^^^^^^^^

First clone the repository::

   git clone https://github.com/astro-tiptop/TIPTOP.git

Navigate to the folder in which you have cloned TIPTOP and install it (remove
``--user`` to install to all users)::

   pip install -e --user .

You may also want to clone and install the other repositories, otherwise their
last release is installed as dependency of TIPTOP:

- https://github.com/astro-tiptop/MASTSEL
- https://github.com/astro-tiptop/SYMAO
- https://github.com/astro-tiptop/SEEING
- https://github.com/astro-tiptop/P3/

Quickstart
----------

To try execute the project you can use ``tiptop.overallSimulation()``
Here is some example code to use it:

.. code-block::

   from tiptop import overallSimulation
   overallSimulation("perfTest", HarmoniSCAO_1, "perfTest",
                     "testPyramid", doPlot=False, doConvolve=False)
