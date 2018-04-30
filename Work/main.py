# main.py
# All of the analysis is functionally called upon within this file, compiling
# all of the results into a neat and readable manner.

import numpy as np

# Import associated files.
import conditions as data
import simulation as sim
import plotting as plot

# Materials List:
materialList = \
[
    'air',
    'brick',
    'wood',
    'copper'
]

# Material Enumeration:
class material:
    air, brick, wood, copper = range(len(materialList))

# Season Enumeration:
class season:
    summer, winter = range(2)

save = False
time = 10*24*3600

#T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.air, 0.15, 25, dt=100)
#plot.timeEvolution(season.summer, material.wood, T_plot, dx, dt, steps, save)

airData, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.air, 0.15, 25, time=time, dt=100)
brickData, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.brick, 0.15, 25, time=time, dt=100)
woodData, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.wood, 0.15, 25, time=time, dt=100)
copperData, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.copper, 0.15, 25, time=time, dt=100)

allData = np.zeros((4, data.spatialCells))
allData[0] = airData[-1, :]
allData[1] = brickData[-1, :]
allData[2] = woodData[-1, :]
allData[3] = copperData[-1, :]

plot.materialAnalysis(season.summer, None, allData, data.x, dx, save)

#T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(season.winter, material.iron, 1, 50, dt = 10)
#plot.timeEvolution(season.summer, material.iron, T_plot, dx, dt, steps)

#save = True
#for i in range(2):
#    for j in range(len(materialList)):
#        T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(i, j, 0.12, 50)
#        plot.timeEvolution(i, j, T_plot, dx, dt, steps, save)
