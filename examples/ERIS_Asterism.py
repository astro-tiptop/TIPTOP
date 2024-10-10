import time
import datetime
#%matplotlib inline
from tiptop.tiptop import *
rc("text", usetex=False)

# 1
# uncomment this to preform the computation of the 1000 example fields
# can be done once for all, as long as the system parameters are not changed
# can be run multiple times (i.e. if a parameter is changed)
# the first time will also generate the input data, to regenerate the data delete the files savedRandom*.npy
sr, fw, ee, covs, simul = asterismSelection("TestERISRandom", "astTest", "ERISastRandom", 'astTest', 'testERISRandom', doPlot=False)

# 2
# reloads the data, assumes 1 was done before 
sr, fw, ee, covs, simul = reloadAsterismSelection("TestERISRandom", "astTest", "ERISastRandom", 'astTest', 
                                                  'testERISRandom', doPlot=False)

# 3
# use the first 800 fields to fit the heuristic model
simul.fitHeuristicModel(0, 799, 'randomData0_799')

# 4
# test the heuristic model perfomances on unseen data
simul.testHeuristicModel(799, 999, 'randomData0_799', [])

# 5
# example code to use the heuristic model
# first all the stars are evaluated using the heuristic model, and sorted in order of increasing penalty
# then the actual jitters, strehls and fwhm are computed, one at a time 

simulation = asterismSimulation("TestERISSingles", "astTest", "ERISastSingles1", 'astTest', 
                                'testERISSingles', doPlot=False, addSrAndFwhm=True, verbose=False)

print('*************************************************************')
approxResult = simulation.runHeuristicModel()
print(approxResult)
print('*************************************************************')
fullSimulationResultSingle = simulation.computeAsterisms(50, 2, doConvolve=True)
print(fullSimulationResultSingle)
print('*************************************************************')
fullSimulationResult = simulation.computeAsterisms(50, doConvolve=False)
print(fullSimulationResult)
print('*************************************************************')
fullSimulationResult = simulation.computeAsterisms(50, doConvolve=True)
print(fullSimulationResult)

