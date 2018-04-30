# plotting.py
# All of the plotting code used in the analysis to render the various house
# conditions is located within this file.

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define title's and filenames for the various plot types.

#=============#
# Plot Titles #
#=============#

title_timeEvo = \
[
    [
        "Time Evolving 1D Temperature Diffusion\nSeason: Summer; Material: Air",
        "Time Evolving 1D Temperature Diffusion\nSeason: Summer; Material: Brick",
        "Time Evolving 1D Temperature Diffusion\nSeason: Summer; Material: Wood",
        "Time Evolving 1D Temperature Diffusion\nSeason: Summer; Material: Copper"
    ],
    [
        "Time Evolving 1D Temperature Diffusion\nSeason: Winter; Material: Air",
        "Time Evolving 1D Temperature Diffusion\nSeason: Winter; Material: Brick",
        "Time Evolving 1D Temperature Diffusion\nSeason: Winter; Material: Wood",
        "Time Evolving 1D Temperature Diffusion\nSeason: Winter; Material: Copper"
    ]
]

title_materialAnalysis = \
[
    "Material Analysis\nSeason: Summer",
    "Material Analysis\nSeason: Winter"
]

title_wallAnalysis = \
[
    "Wall Thickness Analysis\nSeason: Summer",
    "Wall Thickness Analysis\nSeason: Winter"
]

title_tempDistribution = \
[
    [
        "2D Temperature Diffusion\nSeason: Summer; Material: Air",
        "2D Temperature Diffusion\nSeason: Summer; Material: Brick",
        "2D Temperature Diffusion\nSeason: Summer; Material: Wood",
        "2D Temperature Diffusion\nSeason: Summer; Material: Copper",
    ],
    [
        "2D Temperature Diffusion\nSeason: Winter; Material: Air",
        "2D Temperature Diffusion\nSeason: Winter; Material: Brick",
        "2D Temperature Diffusion\nSeason: Winter; Material: Wood",
        "2D Temperature Diffusion\nSeason: Winter; Material: Copper",
    ]
]

#================#
# Plot Filenames #
#================#

filename_timeEvo = \
[
    [
        "./TimeEvoPlots/TimeEvo_Summer_Air.png",
        "./TimeEvoPlots/TimeEvo_Summer_Brick.png",
        "./TimeEvoPlots/TimeEvo_Summer_Wood.png",
        "./TimeEvoPlots/TimeEvo_Summer_Copper.png"
    ],
    [
        "./TimeEvoPlots/TimeEvo_Winter_Air.png",
        "./TimeEvoPlots/TimeEvo_Winter_Brick.png",
        "./TimeEvoPlots/TimeEvo_Winter_Wood.png",
        "./TimeEvoPlots/TimeEvo_Winter_Copper.png"
    ]
]

filename_materialAnalysis = \
[
    "./MaterialAnalysisPlots/MaterialAnalysis_Summer.png",
    "./MaterialAnalysisPlots/MaterialAnalysis_Winter.png"
]

filename_wallAnalysis = \
[
    "./WallThicknessAnalysisPlots/WallThicknessAnalysis_Summer.png",
    "./WallThicknessAnalysisPlots/WallThicknessAnalysis_Winter.png"
]

filename_tempDistribution = \
[
    [
        "./TempDistributionPlots/TempDistribution_Summer_Air.png",
        "./TempDistributionPlots/TempDistribution_Summer_Brick.png",
        "./TempDistributionPlots/TempDistribution_Summer_Wood.png",
        "./TempDistributionPlots/TempDistribution_Summer_Copper.png"
    ],
    [
        "./TempDistributionPlots/TempDistribution_Winter_Air.png",
        "./TempDistributionPlots/TempDistribution_Winter_Brick.png",
        "./TempDistributionPlots/TempDistribution_Winter_Wood.png",
        "./TempDistributionPlots/TempDistribution_Winter_Copper.png"
    ]
]

#====================#
# Plotting Functions #
#====================#

# This plotting function renders the temperature distribution in 1D as a func-
# tion of time, resulting in a 3D plot depicting the temperature T at a time and
# position pair (t, x).
def timeEvolution(season, insulation, data, dx, dt, steps, save):
    # Define variables to easily access the data components.
    t, x = np.meshgrid(range(data.shape[0]), range(data.shape[1]))
    T = data[t, x]

    # Alter position data to account for varying wall discretization.
    #x = np.transpose(x)
    #x[:, :] = xData[:]
    #x = np.transpose(x)

    # Define a figure as a base on which the final plot will be displayed.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection = "3d")

    # Set the plot type to be a wireframe plot.
    ax.plot_wireframe(t*dt*steps, x*dx, T)

    # Set the plot's title.
    ax.set_title(title_timeEvo[season][insulation])

    # Set the plot's labels.
    ax.set_xlabel(r"Time $t$ (s)")
    ax.set_ylabel(r"Position $x$ (m)")
    ax.set_zlabel(r"Temperature $T$ (K)")

    # Set the plot to be compactly displayed.
    fig.tight_layout()

    if save:
        # Save plot to a file.
        fig.savefig(filename_timeEvo[season][insulation])

    # Render plot.
    plt.show()

# This plotting function renders the temperature distribution in 1D as a func-
# tion of insulation material, resulting in a 2D plot depicting the final dis-
# tribution after the 48-hour simulation for the various materials being inves-
# tigated.
def materialAnalysis(season, materialList, data, xData, save):
    # Define a figure as a base on which the final plot will be displayed.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # Add data sets to the figure.
    ax.plot(xData, data[0, :], color = 'skyblue', linestyle = 'solid',
            label = materialList[0])
    ax.plot(xData, data[1, :], color = 'black', linestyle = 'dotted',
            label = materialList[1])
    ax.plot(xData, data[2, :], color = 'forestgreen', linestyle = 'dashdot',
            label = materialList[2])
    ax.plot(xData, data[3, :], color = 'slategrey', linestyle = 'dashed',
            label = materialList[3])

    # Set the plot's title.
    ax.set_title(title_materialAnalysis[season])

    # Set the plot's labels, legend, and layout.
    ax.set_xlabel(r"Position $x$ (m)")
    ax.set_ylabel(r"Temperature $T$ (K)")

    ax.legend(loc = "best")
    ax.grid()

    fig.tight_layout()

    if save:
        # Save plot to a file.
        fig.savefig(filename_materialAnalysis[season])

    # Render plot.
    plt.show()

# This plotting function renders the temperature distribution in 1D as a func-
# tion of wall thickness, resulting in a 2D plot depicting the final distribu-
# tion after the 48-hour simulation for the various wall thicknesses being in-
# vestigated.
def wallAnalysis(season, thicknessList, data, xData, save):
    # Define a figure as a base on which the final plot will be displayed.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # Add data sets to the figure.
    ax.plot(xData[0], data[0], color = 'skyblue', linestyle = 'solid',
            label = thicknessList[0])
    ax.plot(xData[1], data[1], color = 'black', linestyle = 'dotted',
            label = thicknessList[1])
    ax.plot(xData[2], data[2], color = 'forestgreen', linestyle = 'dashdot',
            label = thicknessList[2])
    ax.plot(xData[3], data[3], color = 'slategrey', linestyle = 'dashed',
            label = thicknessList[3])

    # Set the plot's title.
    ax.set_title(title_wallAnalysis[season])

    # Set the plot's labels, legend, and layout.
    ax.set_xlabel(r"Position $x$ (m)")
    ax.set_ylabel(r"Temperature $T$ (K)")

    ax.legend(loc = "best")
    ax.grid()

    fig.tight_layout()

    if save:
        # Save plot to a file.
        fig.savefig(filename_wallAnalysis[season])

    # Render plot.
    plt.show()

# This plotting function renders the temperature distribution in 2D, resulting
# in a 3D plot depicting the temperature T at an (x,y) coordinate in space at
# the end of the 48-hour simulation period.
def temperatureDistribution(season, insulation, data, dx, dy, dt, steps, save):
    #
    pass
