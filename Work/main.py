# main.py
# All of the analysis is functionally called upon within this file, compiling
# all of the results into a neat and readable manner.

import numpy as np

# Import associated files.
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

time = 2*24*3600
save = False
T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.wood, 0.12, 50, time=time, dt=100)
plot.timeEvolution(season.summer, material.wood, T_plot, dx, dt, steps, save)

#T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(season.winter, material.iron, 1, 50, dt = 10)
#plot.timeEvolution(season.summer, material.iron, T_plot, dx, dt, steps)

#save = True
#for i in range(2):
#    for j in range(len(materialList)):
#        T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(i, j, 0.12, 50)
#        plot.timeEvolution(i, j, T_plot, dx, dt, steps, save)
