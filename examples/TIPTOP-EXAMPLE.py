#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# get_ipython().run_line_magic('matplotlib', 'inline')

from tiptop.tiptop import *

from matplotlib import rc
rc("text", usetex=False)

import os
from pathlib import Path
import tiptop

base_path = Path(tiptop.__file__).resolve().parents[1]
os.chdir(base_path)


# In[ ]:


overallSimulation("tiptop/perfTest", "SOUL", 'tiptop/perfTest', 'testSOUL', doPlot=True, doConvolve=True)


# In[ ]:


overallSimulation("tiptop/perfTest", "MAVIS", 'tiptop/perfTest', 'testMAVIS', doPlot=True, doConvolve=True)


# In[ ]:


overallSimulation("tiptop/perfTest", "SPHERE", 'tiptop/perfTest', 'testSPHERE', doPlot=True, doConvolve=True)