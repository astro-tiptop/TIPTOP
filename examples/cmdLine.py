from tiptop.tiptop import *

from matplotlib import rc
rc("text", usetex=False)

from pathlib import Path

base_path = Path(__file__).resolve().parent.parent

sr, fw, ee, covs, simul = asterismSelection("Test1000", str(base_path / "tiptop/astTest"), "MAVISast", str(base_path / "tiptop/astTest"), 'testMAVIS', doPlot=False)