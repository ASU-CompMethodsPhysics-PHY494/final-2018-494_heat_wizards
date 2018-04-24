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

# Enumerations:
class material:
    air, iron = range(len(materialList))

class season:
    summer, winter = range(2)

T_plot, (dx, dt, step) = sim.CrankNicolson_1D(season.summer, material.iron, 1, 50)
plot.plot_T(T_plot, dx, dt, step)
