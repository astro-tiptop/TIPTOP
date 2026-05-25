from tiptop.tiptop import *

from matplotlib import rc
rc("text", usetex=False)

import os
from pathlib import Path
import tiptop

base_path = Path(tiptop.__file__).resolve().parents[1]
os.chdir(base_path)

sr, fw, ee, covs, simul = asterismSelection("Test1000", "tiptop/astTest", "MAVISast", 'tiptop/astTest', 'testMAVIS', doPlot=False)