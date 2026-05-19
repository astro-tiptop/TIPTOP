#%%
%reload_ext autoreload
%autoreload 2

import os
from tiptop.TipTorchSimulation import baseSimulation

# TODO: override TipTorch paths
# TODO: automatic 

#%%
root = os.path.dirname(os.path.abspath(__file__))

simulation = baseSimulation(
    path = os.path.join(root, "../tiptop/perfTest"),
    parametersFile = "muse_ltao.ini",
    outputDir = os.path.join(root, '../tiptop/perfTest'),
    outputFile = 'testNFM',
    doConvolve = True,
    doPlot = True,
)

# %%
simulation.doOverallSimulation()