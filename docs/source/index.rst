.. TIPTOP documentation master file, created by
   sphinx-quickstart on Wed May 18 15:52:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Welcome to TIPTOP's documentation!
==================================

This documentation is still under development. If you find confusing elements please submit an issue on `Github`_. and use the label documentation.
Also Some heavy modification is to be expected in TIPTOP. the documentation will try to follow.

In order to be able to easily predict the AO performance, we have developed this fast algorithm producing the expected Adaptive Optics (AO see https://en.wikipedia.org/wiki/Adaptive_optics) Point Spread Function (PSF) for any of the existing AO observing modes (Single-Conjugate-AO, Laser-Tomographic-AO, Multi-Conjugate-AO, Ground-Layer-AO), and any atmospheric conditions. This TIPTOP tool takes its roots in an analytical approach, where the simulations are done in the Fourier domain. This allows to reach a very fast computation time (few seconds per PSF with GPU acceleration), and efficiently explore the wide parameter space.

References
----------

Reference: "TIPTOP: a new tool to efficiently predict your favorite AO PSF" SPIE 2020 (ARXIV: https://doi.org/10.48550/arXiv.2101.06486).


Check out the :doc:`usage` section for further information, including how to
:ref:`install <installation>` the project.

.. note::

   This documentation is under development





.. toctree::
   
   usage
   parameterFile
   wishList
   howto



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _Github: https://github.com/FabioRossiArcetri/TIPTOP/issues
