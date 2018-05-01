# main.py
# All of the analysis is functionally called upon within this file, compiling
# all of the results into a neat and readable manner.

import numpy as np

# Import associated files.
import conditions as data
import simulation as sim
import plotting as plot

# Parameter determining whether plots are saved to files or not.
save = True

# Parameters to define how fine to sample in time.
time = data.time
dt = data.dt

#===============#
# Seasonal Data #
#===============#

# Season Enumeration:
class season:
    summer, winter = range(2)

#===============#
# Material Data #
#===============#

# Materials List:
materialList = \
[
    'Air',
    'Brick',
    'Wood',
    'Copper'
]

# Material Enumeration:
class material:
    air, brick, wood, copper = range(len(materialList))

#=====================#
# Wall Thickness Data #
#=====================#

# Thickness List:
thicknessList = \
[
    '0.10 m',
    '0.15 m',
    '0.20 m',
    '1.00 m'
]

# Wall Thickness Enumeration:
class wall:
    thin, average, thick, outlier = range(len(thicknessList))

#======================#
# Default Square House #
#      Dimensions      #
#======================#

# Physical Distances (m)
length = 25
wall_thickness = \
[
    0.10,   # Thin
    0.15,   # Average
    0.20,   # Thick
    1.00    # Outlier
]

interiorCells = int(length // data.dx)
wallCells = \
[
    int(wall_thickness[0] // data.dw),
    int(wall_thickness[1] // data.dw),
    int(wall_thickness[2] // data.dw),
    int(wall_thickness[3] // data.dw)
]

spatialCells = \
[
    wallCells[0] + interiorCells + wallCells[0],
    wallCells[1] + interiorCells + wallCells[1],
    wallCells[2] + interiorCells + wallCells[2],
    wallCells[3] + interiorCells + wallCells[3]
]

# Create position arrays accurately reflecting discretizations.
xData = []
wData = []

for index in range(len(spatialCells)):
    xTemp = np.zeros(spatialCells[index])
    indexArray = np.array(range(spatialCells[index]))

    # Left Wall:
    xTemp[:wallCells[index]] = data.dw * indexArray[:wallCells[index]]

    # Interior:
    xTemp[wallCells[index]:-wallCells[index]] = wall_thickness[index] + \
                                                data.dx * indexArray[:interiorCells]

    # Right Wall:
    xTemp[-wallCells[index]:] = wall_thickness[index] + length + \
                                data.dw * indexArray[:wallCells[index]]

    xData.append(xTemp)

    # Well-Defined Wall
    wTemp = data.dw_zoom * np.array(range(int(wall_thickness[index] // data.dw_zoom)))
    wData.append(wTemp)

#==============================================================================#
#---------------------------------- Analysis ----------------------------------#
#==============================================================================#

# Cycle simulations over seasons.
for seasons in range(2):
    #===================#
    # 1D Time Evolution #
    #      Analysis     #
    #===================#

    # Cycle simulations over materials.
    for materials in range(len(materialList)):
        # Simulate:
        T_plot = sim.CrankNicolson_1D(seasons, materials,
                                      wall_thickness[wall.average], length,
                                      time = time, dt = dt)

        # Plot Results:
        plot.timeEvolution(seasons, materials, T_plot,
                           data.dx, data.dt, data.steps, save)

    #======================#
    # 1D Material Analysis #
    #======================#

    # Parameter determining the domain of the results.
    # Initial analysis is set to cover the entire house.
    zoom = False

    # Simulate:
    airData = sim.CrankNicolson_1D(seasons, material.air, wall_thickness[wall.average], length, time = time, dt = dt)
    brickData = sim.CrankNicolson_1D(seasons, material.brick, wall_thickness[wall.average], length, time = time, dt = dt)
    woodData = sim.CrankNicolson_1D(seasons, material.wood, wall_thickness[wall.average], length, time = time, dt = dt)
    copperData = sim.CrankNicolson_1D(seasons, material.copper, wall_thickness[wall.average], length, time = time, dt = dt)

    materialData = []
    materialData.append(airData[-1, :])
    materialData.append(brickData[-1, :])
    materialData.append(woodData[-1, :])
    materialData.append(copperData[-1, :])

    # Plot Results:
    plot.materialAnalysis(seasons, materialList, materialData,
                          xData[wall.average], wData[wall.average], zoom, save)

    # Following analysis is set to cover a single wall.
    zoom = True

    # Simulate:
    airData = sim.CrankNicolson_1D(seasons, material.air, wall_thickness[wall.average], length, time = time, dw = data.dw_zoom, dt = dt)
    brickData = sim.CrankNicolson_1D(seasons, material.brick, wall_thickness[wall.average], length, time = time, dw = data.dw_zoom, dt = dt)
    woodData = sim.CrankNicolson_1D(seasons, material.wood, wall_thickness[wall.average], length, time = time, dw = data.dw_zoom, dt = dt)
    copperData = sim.CrankNicolson_1D(seasons, material.copper, wall_thickness[wall.average], length, time = time, dw = data.dw_zoom, dt = dt)

    materialData = []
    materialData.append(airData[-1, :])
    materialData.append(brickData[-1, :])
    materialData.append(woodData[-1, :])
    materialData.append(copperData[-1, :])

    # Plot Results:
    plot.materialAnalysis(seasons, materialList, materialData,
                          xData[wall.average], wData[wall.average], zoom, save)

    #===================#
    # 1D Wall Thickness #
    #      Analysis     #
    #===================#

    # Note: Wall thickness analysis is done using copper as the constant insul-
    #       ation material since it has the most apparent effect in the results
    #       due to its physical properties.

    # Simulate:
    thinWallData = sim.CrankNicolson_1D(seasons, material.copper, wall_thickness[wall.thin], length, time = time, dt = dt)
    averageWallData = sim.CrankNicolson_1D(seasons, material.copper, wall_thickness[wall.average], length, time = time, dt = dt)
    thickWallData = sim.CrankNicolson_1D(seasons, material.copper, wall_thickness[wall.thick], length, time = time, dt = dt)
    outlierWallData = sim.CrankNicolson_1D(seasons, material.copper, wall_thickness[wall.outlier], length, time = time, dt = dt)

    wallData = []
    wallData.append(thinWallData[-1, :])
    wallData.append(averageWallData[-1, :])
    wallData.append(thickWallData[-1, :])
    wallData.append(outlierWallData[-1, :])

    # Plot Results:
    plot.wallAnalysis(seasons, thicknessList, wallData, xData, save)
