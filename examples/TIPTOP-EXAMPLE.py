#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# get_ipython().run_line_magic('matplotlib', 'inline')

from tiptop.tiptop import *

from matplotlib import rc
rc("text", usetex=False)

from pathlib import Path

base_path = Path(__file__).resolve().parent.parent


# In[ ]:


overallSimulation(str(base_path / "tiptop/perfTest"), "SOUL", str(base_path / "tiptop/perfTest"), 'testSOUL', doPlot=True, doConvolve=True)


# In[ ]:


overallSimulation(str(base_path / "tiptop/perfTest"), "MAVIS", str(base_path / "tiptop/perfTest"), 'testMAVIS', doPlot=True, doConvolve=True)


# In[ ]:


overallSimulation(str(base_path / "tiptop/perfTest"), "SPHERE", str(base_path / "tiptop/perfTest"), 'testSPHERE', doPlot=True, doConvolve=True)