# plotting.py
# All of the plotting code used in the analysis to render the various house
# conditions is located within this file.

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# This plotting function renders the temperature distribution in 1D as a func-
# tion of time, resulting in a 3D plot depicting the temperature T at a time and
# position pair (t, x).
def timeEvolution(season, insulation, data, dx, dt, steps):
    # Define variables to easily access the data components.
    t, x = np.meshgrid(range(data.shape[0]), range(data.shape[1]))
    T = data[t, x]

    # Define a figure as a base on which the final plot will be displayed.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection = "3d")

    # Set the wireframe to be active and based on the data's discretization.
    ax.plot_wireframe(t*dt*steps, x*dx, T)

    # Set the plot's labels.
    ax.set_xlabel(r"Time $t$ (s)")
    ax.set_ylabel(r"Position $x$ (m)")
    ax.set_zlabel(r"Temperature $T$ (K)")

    # Set the plot to be compactly displayed.
    fig.tight_layout()

    plt.show()
    return ax

# This plotting function renders the temperature distribution in 1D as a func-
# tion of insulation material, resulting in a 2D plot depicting the final dis-
# tribution after the 48-hour simulation for the various materials being inves-
# tigated. The function is generalized so that the material list can be altered
# to include more or less materials.
def materialAnalysis(season, materialList, data, dx):
    #
    pass

# This plotting function renders the temperature distribution in 1D as a func-
# tion of wall thickness, resulting in a 2D plot depicting the final distribu-
# tion after the 48-hour simulation for the various wall thicknesses being in-
# vestigated. The function is generalized so that the thickness list can be al-
# tered to include more or less thicknesses.
def wallAnalysis(season, thicknessList, data, dx):
    #
    pass

# This plotting function renders the temperature distribution in 2D, resulting
# in a 3D plot depicting the temperature T at an (x,y) coordinate in space at
# the end of the 48-hour simulation period.
def temperatureDistribution(season, insulation, data, dx, dy, dt, steps):
    #
    pass
