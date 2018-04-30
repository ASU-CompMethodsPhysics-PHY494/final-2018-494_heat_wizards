# main.py
# All of the analysis is functionally called upon within this file, compiling
# all of the results into a neat and readable manner.

import numpy as np

# Import associated files.
import conditions as data
import simulation as sim
import plotting as plot

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

#===========#
# Main Code #
#===========#

save = False
time = 10*24*3600

T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.air, wall_thickness[wall.average], length, time=time, dt=100)
plot.timeEvolution(season.summer, material.wood, T_plot, dx, dt, steps, save)



airData, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.air, wall_thickness[wall.average], length, time=time, dt=100)
brickData, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.brick, wall_thickness[wall.average], length, time=time, dt=100)
woodData, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.wood, wall_thickness[wall.average], length, time=time, dt=100)
copperData, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.copper, wall_thickness[wall.average], length, time=time, dt=100)

materialData = np.zeros((4, spatialCells[wall.average]))
materialData[0] = airData[-1, :]
materialData[1] = brickData[-1, :]
materialData[2] = woodData[-1, :]
materialData[3] = copperData[-1, :]

plot.materialAnalysis(season.summer, materialList, materialData, xData[wall.average], save)



wall0Data, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.copper, wall_thickness[wall.thin], length, time=time, dt=100)
wall1Data, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.copper, wall_thickness[wall.average], length, time=time, dt=100)
wall2Data, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.copper, wall_thickness[wall.thick], length, time=time, dt=100)
wall3Data, (dx, dt, steps) = sim.CrankNicolson_1D(season.summer, material.copper, wall_thickness[wall.outlier], length, time=time, dt=100)

wallData = []
wallData.append(wall0Data[-1, :])
wallData.append(wall1Data[-1, :])
wallData.append(wall2Data[-1, :])
wallData.append(wall3Data[-1, :])

plot.wallAnalysis(season.summer, thicknessList, wallData, xData, save)



#T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(season.winter, material.iron, 1, 50, dt = 10)
#plot.timeEvolution(season.summer, material.iron, T_plot, dx, dt, steps)

#save = True
#for i in range(2):
#    for j in range(len(materialList)):
#        T_plot, (dx, dt, steps) = sim.CrankNicolson_1D(i, j, 0.12, 50)
#        plot.timeEvolution(i, j, T_plot, dx, dt, steps, save)
