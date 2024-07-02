Troubleshooting
===============

We have some custom error messages that will appear if some variables are set the wrong way or in way that not physically possible.
Here is a list of errors that you may encounter.

ValueError: Error : the PSF field of view is too small to simulate the AO correction area
-----------------------------------------------------------------------------------------
This error appears if ``sensor_science.Fieldofview`` is inferior to ``min(source_science.Wavelength-sensor_science.SpectralBandWidth)/DM.DmPitchs/sensor_science.PixelScale*180*3600*10**3/np.pi`` .
If not set ``sensor_science.SpectralBandWidth = 0``.

Either make ``sensor_science.Fieldofview`` bigger or play with ``min(source_science.Wavelength-sensor_science.SpectralBandWidth)`` or ``DM.DmPitchs`` or ``sensor_science.PixelScale`` to make this term smaller.


Warning: No information about the tip-tilt star can be retrieved
----------------------------------------------------------------
/!\Statement needs verification.
This warning arises if source_HO.Height is defined and different from zero and sources_LO is not defined.
This means that for the AO correction only a laser guide star will be considered.


ValueError: ValueError: '[X1,X2,...]' in section 'Y' must have the same length
------------------------------------------------------------------------------
As mentionned in the parameter file description there are some prameter that should be lists and should have the same length as other parameters.
If you have defined at least one of them, and it contains more than one element, you need to define all the other associated parameter, and it must contain a list of the same length.

For example in the 'atmosphere' section If you set 'Cn2Weights'=[0.8,0.2], then you must set 'Cn2Heights', 'WindSpeed', 'WindDirection' and they need to contain a list with length 2.

OutOfMemoryError: Out of memory allocating A bytes (allocated so far: B bytes).
-------------------------------------------------------------------------------
This is a GPU error, you gave parameters that creates arrays that are too big for your GPU memmory. 
If you do not care about GPU we recommend that you start your tiptop script as follow::

   import os
   os.environ.get('MASTSEL_DISABLE_GPU')
   os.environ['MASTSEL_DISABLE_GPU'] = 'FALSE'
   os.environ.get('P3_DISABLE_GPU')
   os.environ['P3_DISABLE_GPU'] = 'FALSE'
   os.environ.get('SEEING_DISABLE_GPU')
   os.environ['SEEING_DISABLE_GPU'] = 'FALSE'
   #do not import TIPTOP before setting the previous global environment variables
   from tiptop.tiptop import * 

If you want to beneficiate from your GPU acceleration you need to find in the following list the responsible.

.. note::

   If your case is not listed, report your issue on [github](https://github.com/astro-tiptop/TIPTOP/issues/new), provide the full traceback of the error message and add the label documentation to your issue.
   Or if you feel up to the task of Documenting this issue yourself, following is a a couple lines of code to help you find which allocations take the most space::

      import cupy
      import traceback
      import sys
      mempool = cupy.get_default_memory_pool()
      
      def myAllocator(size):
          global mempool
          print('\n{} GB'.format(size/10**9))
          traceback.print_stack(limit = 10,file=sys.stdout)
          return mempool.malloc(size)
      
      cupy.cuda.set_allocator(myAllocator)
      
      from tiptop.tiptop import *
      #here your tiptop code

psdSetToPsfSet 
^^^^^^^^^^^^^^
If somewhere in the traceback of the error you find something as follow::

   
  File ~\TIPTOP\tiptop\baseSimulation.py:353 in doOverallSimulation
    pointings_SR, psdPointingsArray, psfLongExpPointingsArr, pointings_FWHM_mas = self.***psdSetToPsfSet***(self.N,

If you rerun the script add the option verbose=True in::

   overallSimulation({your other arguments},***verbose = True***)

In the printed message somwhere it should have printed::

   fao.PSD.shape: (X, Y, Z)

With X,Y,Z being integers, you can then compute (X*Y)**2*Z/A (A being the memory that TIPTOP tried to allocate given in the OutOfMemoryError message). If the result of the division is exactly 64 or 32, the following should help you change the parameters to make it possible for your GPU to compute what you want.
The PSD array size is directly impacted by the following parameters:

- sensor_science.FieldofView
- telescope.TelescopeDiameter
- source_science.Wavelength
- sensor_science.PixelScale

With the following relationship::

    np.ceil(2/(source_science.Wavelength*(180*3600*10**3/np.pi)/(sensor_science.PixelScale*telescope.TelescopeDiameter)))*sensor_science.FieldofView

The result of this operation should be the X dimension of the fao.PSD.shape. You should now find the parameter that you can change to reduce the size of the array until it has a size that holds on you GPU.
