# simulation.py
# This file contains all of the relevant functions necessary to conduct simu-
# lations for various house conditions, whether the variable is insulation
# material, wall thickness, or house size/shape.

import numpy as np

# Import associated files.
import conditions

def CrankNicolson_1D(season, insulation, wall_thickness, length,
                     time = conditions.time, dx = conditions.dx,
                     dt = conditions.dt, steps = conditions.steps):
    # Quantize Space and Time:
    temporalCells = int(time // dt)             # Time
    interiorCells = int(length // dx)           # Interior
    wallCells = int(wall_thickness // dx)       # Wall

    spatialCells = wallCells + interiorCells + wallCells

    # Material Properties:
    kappa = np.zeros(spatialCells)
    c_heat = np.zeros(spatialCells)
    rho = np.zeros(spatialCells)

    # Define material properties for each spatial cell to account for the diff-
    # erence between the insulation within the wall and the air within the in-
    # terior of the house.

    # Wall Cells (Left) - Insulation Material
    kappa[:wallCells] = conditions.kappa[insulation]
    c_heat[:wallCells] = conditions.c_heat[insulation]
    rho[:wallCells] = conditions.rho[insulation]

    # Interior Cells - Air
    kappa[wallCells:-wallCells] = conditions.kappa[0]
    c_heat[wallCells:-wallCells] = conditions.c_heat[0]
    rho[wallCells:-wallCells] = conditions.rho[0]

    # Wall Cells (Right) - Insulation Material
    kappa[-wallCells:] = conditions.kappa[insulation]
    c_heat[-wallCells:] = conditions.c_heat[insulation]
    rho[-wallCells:] = conditions.rho[insulation]

    """ Inefficient Method:
    for cellIndex in range(spatialCells):
        # Wall Cells (Left)
        if cellIndex < wallCells:
            # Insulation Material
            kappa[cellIndex] = 237
            c_heat[cellIndex] = 900
            rho[cellIndex] = 2700
        # Interior Cells
        elif cellIndex < (wallCells + interiorCells):
            # Air
            kappa[cellIndex] = 2.624
            c_heat[cellIndex] = 1000
            rho[cellIndex] = 1.204
        # Wall Cells (Right)
        else:
            # Insulation Material
            kappa[cellIndex] = 237
            c_heat[cellIndex] = 900
            rho[cellIndex] = 2700
    """

    # Define the thermal diffusivity constant at each spatial cell.
    eta = kappa * dt / (c_heat * rho * dx**2)

    # Create an array to store the temperature at each spatial cell and create
    # an additional 2D array to store the temperature values in space and time
    # for plotting purposes.
    # Note: The steps parameter represents the number of time steps being record-
    #       ed for the final plot. By default, the plot will be rendered based
    #       on 20 steps in time while the iterations are based on dt.
    T = np.zeros(spatialCells)
    timeEvolvingT = np.zeros((int(np.ceil(temporalCells / steps)) + 1,
                              spatialCells))

    # Initial Condition(s):
    T[1:-1] = conditions.T0

    # Boundary Condition(s):
    T[0] = T[-1] = conditions.boundary(season, 0)

    #--------------------------------------------------------------------------#

    # Define the tridiagonal matrix representing the system of equations solving
    # the problem. The matrix has the following dimensions,
    #
    #   (spatialCells - 2) x (spatialCells - 2)
    #
    # and solves the equation,
    #
    #   M_eta * T[1:-1, t+1] = bT

    # Define the dimensionality of the square tridiagonal matrix, M_eta.
    dimension = spatialCells - 2

    # Define relevant constants when solving the matrix equation.
    alpha = 2/eta + 2
    beta = 2/eta - 2

    # Explicitly define the tridiagonal matrix, M_eta.
    M_eta = np.diagflat(-np.ones(dimension - 1), 1) + \
            np.diagflat(alpha[1:-1], 0) + \
            np.diagflat(-np.ones(dimension - 1), -1)

    # Initialize the terms at the current timestep.
    bT = np.zeros(dimension)

    # Define the inverse to the tridiagonal matrix, M_eta, once in order to im-
    # prove the efficiency of the algorithm since this will create the inverse
    # once rather than during each iteration using np.linalg.solve().
    inv_M_eta = np.linalg.inv(M_eta)

    # Define a token to represent the steps in time being stored for plotting
    # and initialize the initial time step to the initial conditions of the
    # system.
    t_index = 0
    timeEvolvingT[t_index, :] = T

    # Simulate!
    # The simulation is conducted for all timesteps not including the initial
    # time, t = 0, since that has been defined by the conditions.
    for t in range(1, temporalCells):
        # Define the temperature values at the previous timestep that contribute
        # to the next timestep's temperatures.
        bT[:] = T[:-2] + beta[1:-1] * T[1:-1] + T[2:]

        # Update the boundary values to account for the boundary at the next
        # timestep.
        bT[0] += conditions.boundary(season, t * dt)
        bT[-1] += conditions.boundary(season, t * dt)

        # Solve the matrix equation using the predefined inverse rather than
        # using np.linalg.solve().
        T[1:-1] = np.dot(inv_M_eta, bT)

        # Update the boundary temperatures.
        T[0] = T[-1] = conditions.boundary(season, t * dt)

        # In the case the timestep has reached an incremental value for the
        # timesteps being plotted or the final timestep in the simulation, the
        # temperature data is stored within the plotting data.
        if t % steps == 0 or t == temporalCells - 1:
            t_index += 1
            timeEvolvingT[t_index, :] = T

    parameters = (dx, dt, steps)
    return timeEvolvingT, parameters

def CrankNicolson_2D(season, insulation, wall_thickness, length, width,
                     time = conditions.time, dx = conditions.dx,
                     dy = conditions.dy, dt = conditions.dt,
                     steps = conditions.steps):
    # Quantize Space and Time:
    temporalCells = int(time // dt)             # Time
    xInteriorCells = int(length // dx)          # Interior in x
    yInteriorCells = int(width // dy)           # Interior in y
    wallCells = int(wall_thickness // dx)       # Wall

    xSpatialCells = wallCells + xInteriorCells + wallCells
    ySpatialCells = wallCells + yInteriorCells + wallCells

    # Material Properties:
    kappa = np.zeros((xSpatialCells, ySpatialCells))
    c_heat = np.zeros((xSpatialCells, ySpatialCells))
    rho = np.zeros((xSpatialCells, ySpatialCells))

    # Define material properties for each spatial cell to account for the diff-
    # erence between the insulation within the wall and the air within the in-
    # terior of the house.

    # Wall Cells (Bottom) - Insulation Material
    kappa[:wallCells, :] = conditions.kappa[insulation]
    c_heat[:wallCells, :] = conditions.c_heat[insulation]
    rho[:wallCells, :] = conditions.rho[insulation]

    # Wall Cells (Left) - Insulation Material
    kappa[:, :wallCells] = conditions.kappa[insulation]
    c_heat[:, :wallCells] = conditions.c_heat[insulation]
    rho[:, :wallCells] = conditions.rho[insulation]

    # Interior Cells - Air
    kappa[wallCells:-wallCells, wallCells:-wallCells] = conditions.kappa[0]
    c_heat[wallCells:-wallCells, wallCells:-wallCells] = conditions.c_heat[0]
    rho[wallCells:-wallCells, wallCells:-wallCells] = conditions.rho[0]

    # Wall Cells (Top) - Insulation Material
    kappa[-wallCells:, :] = conditions.kappa[insulation]
    c_heat[-wallCells:, :] = conditions.c_heat[insulation]
    rho[-wallCells:, :] = conditions.rho[insulation]

    # Wall Cells (Right) - Insulation Material
    kappa[:, -wallCells:] = conditions.kappa[insulation]
    c_heat[: -wallCells:] = conditions.c_heat[insulation]
    rho[: -wallCells:] = conditions.rho[insulation]

    # Define the thermal diffusivity constant at each spatial cell.
    eta = kappa * dt / (c_heat * rho * dx**2)

    # Create an array to store the temperature at each spatial cell and create
    # an additional 3D array to store the temperature values in space and time
    # for plotting purposes.
    # Note: The steps parameter represents the number of time steps being record-
    #       ed for the final plot. By default, the plot will be rendered based
    #       on 20 steps in time while the iterations are based on dt.
    T = np.zeros((xSpatialCells, ySpatialCells))
    timeEvolvingT = np.zeros((int(np.ceil(temporalCells / steps)) + 1,
                              xSpatialCells, ySpatialCells))

    # Initial Condition(s):
    T[1:-1, 1:-1] = conditions.T0

    # Boundary Condition(s):
    T[0, :] = T[:, 0] = T[-1, :] = T[:, -1] = conditions.boundary(season, 0)

    #--------------------------------------------------------------------------#

    # Define the tridiagonal matrix representing the system of equations solving
    # the problem. The matrix has the following dimensions,
    #
    #   (spatialCells - 2) x (spatialCells - 2)
    #
    # and solves the equation,
    #
    #   M_eta * T[1:-1, t+1] = bT

    # Define the dimensionality of the square tridiagonal matrix, M_eta.
    xDimension = xSpatialCells - 2
    yDimension = ySpatialCells - 2

    # Define relevant constants when solving the matrix equation.
    alpha = 2/eta + 2
    beta = 2/eta - 2

    # Explicitly define the tridiagonal matrix, M_eta.
    M_eta = np.diagflat(-np.ones(dimension - 1), 1) + \
            np.diagflat(alpha[1:-1], 0) + \
            np.diagflat(-np.ones(dimension - 1), -1)

    # Initialize the terms at the current timestep.
    bT = np.zeros(dimension)

    # Define the inverse to the tridiagonal matrix, M_eta, once in order to im-
    # prove the efficiency of the algorithm since this will create the inverse
    # once rather than during each iteration using np.linalg.solve().
    inv_M_eta = np.linalg.inv(M_eta)

    # Define a token to represent the steps in time being stored for plotting
    # and initialize the initial time step to the initial conditions of the
    # system.
    t_index = 0
    timeEvolvingT[t_index, :] = T

    # Simulate!
    # The simulation is conducted for all timesteps not including the initial
    # time, t = 0, since that has been defined by the conditions.
    for t in range(1, temporalCells):
        # Define the temperature values at the previous timestep that contribute
        # to the next timestep's temperatures.
        bT[:] = T[:-2] + beta[1:-1] * T[1:-1] + T[2:]

        # Update the boundary values to account for the boundary at the next
        # timestep.
        bT[0] += conditions.boundary(season, t * dt)
        bT[-1] += conditions.boundary(season, t * dt)

        # Solve the matrix equation using the predefined inverse rather than
        # using np.linalg.solve().
        T[1:-1] = np.dot(inv_M_eta, bT)

        # Update the boundary temperatures.
        T[0] = T[-1] = conditions.boundary(season, t * dt)

        # In the case the timestep has reached an incremental value for the
        # timesteps being plotted or the final timestep in the simulation, the
        # temperature data is stored within the plotting data.
        if t % steps == 0 or t == temporalCells - 1:
            t_index += 1
            timeEvolvingT[t_index, :] = T

    parameters = (dx, dy, dt, steps)
    return timeEvolvingT, parameters
