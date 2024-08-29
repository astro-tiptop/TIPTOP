#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# get_ipython().run_line_magic('matplotlib', 'inline')

from tiptop.tiptop import *
rc("text", usetex=False)


# In[ ]:


overallSimulation("perfTest", "SOUL", 'perfTest', 'testSOUL', doPlot=True, doConvolve=True)


# In[ ]:


overallSimulation("perfTest", "MAVIS", 'perfTest', 'testMAVIS', doPlot=True, doConvolve=True)


# In[ ]:


overallSimulation("perfTest", "SPHERE", 'perfTest', 'testSPHERE', doPlot=True, doConvolve=True)

