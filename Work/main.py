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
    'iron'
]

# Material Enumeration:
class material:
    air, iron = range(len(materialList))

# Season Enumeration:
class season:
    summer, winter = range(2)

T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.iron, 1, 50, dt = 100)
plot.timeEvolution(season.summer, material.iron, T_plot, dx, dt, steps)

#T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(season.winter, material.iron, 1, 50, dt = 10)
#plot.timeEvolution(season.summer, material.iron, T_plot, dx, dt, steps)
