Usage
=====


.. _installation:

Installation
------------

.. code-block:: console

     pip install -e .

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