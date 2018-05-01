# simulation.py
# This file contains all of the relevant functions necessary to conduct simu-
# lations for various house conditions, whether the variable is insulation
# material, wall thickness, or house size/shape.

import numpy as np

# Import associated files.
import conditions

#==============================================================================#
#-------------------------------- 1D Simulator --------------------------------#
#==============================================================================#

def CrankNicolson_1D(season, insulation, wall_thickness, length,
                     time = conditions.time, dx = conditions.dx,
                     dw = conditions.dw, dt = conditions.dt,
                     steps = conditions.steps):
    # Quantize Space and Time:
    temporalCells = int(time // dt)             # Time
    interiorCells = int(length // dx)           # Interior
    wallCells = int(wall_thickness // dw)       # Wall

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
        current_t = t * dt

        # Define the temperature values at the previous timestep that contribute
        # to the next timestep's temperatures.
        bT[:] = T[:-2] + beta[1:-1] * T[1:-1] + T[2:]

        # Update the boundary values to account for the boundary at the next
        # timestep.
        bT[0] += conditions.boundary(season, current_t)
        bT[-1] += conditions.boundary(season, current_t)

        # Solve the matrix equation using the predefined inverse rather than
        # using np.linalg.solve().
        T[1:-1] = np.dot(inv_M_eta, bT)

        # Update the boundary temperatures.
        T[0] = T[-1] = conditions.boundary(season, current_t)

        # In the case the timestep has reached an incremental value for the
        # timesteps being plotted or the final timestep in the simulation, the
        # temperature data is stored within the plotting data.
        if t % steps == 0 or t == temporalCells - 1:
            t_index += 1
            timeEvolvingT[t_index, :] = T

    # parameters = (dx, dt, steps)
    return timeEvolvingT

#==============================================================================#
#-------------------------------- 2D Simulator --------------------------------#
#==============================================================================#

def CrankNicolson_2D(season, insulation, wall_thickness, length, width,
                     time = conditions.time, dx = conditions.dx,
                     dy = conditions.dy, dw = conditions.dw, dt = conditions.dt,
                     steps = conditions.steps):
    # Quantize Space and Time:
    temporalCells = int(time // dt)             # Time
    xInteriorCells = int(length // dx)          # Interior in x
    yInteriorCells = int(width // dy)           # Interior in y
    wallCells = int(wall_thickness // dw)       # Wall

    xSpatialCells = wallCells + xInteriorCells + wallCells
    ySpatialCells = wallCells + yInteriorCells + wallCells

    # Material Properties:
    xKappa = np.zeros(xSpatialCells)
    yKappa = np.zeros(ySpatialCells)

    xC_heat = np.zeros(xSpatialCells)
    yC_heat = np.zeros(ySpatialCells)

    xRho = np.zeros(xSpatialCells)
    yRho = np.zeros(ySpatialCells)

    # Define material properties for each spatial cell to account for the diff-
    # erence between the insulation within the wall and the air within the in-
    # terior of the house.

    # x Dimension:
    # Wall Cells (Left in x) - Insulation Material
    xKappa[:wallCells] = conditions.kappa[insulation]
    xC_heat[:wallCells] = conditions.c_heat[insulation]
    xRho[:wallCells] = conditions.rho[insulation]

    # Interior Cells - Air
    xKappa[wallCells:-wallCells] = conditions.kappa[0]
    xC_heat[wallCells:-wallCells] = conditions.c_heat[0]
    xRho[wallCells:-wallCells] = conditions.rho[0]

    # Wall Cells (Right in x) - Insulation Material
    xKappa[-wallCells:] = conditions.kappa[insulation]
    xC_heat[-wallCells:] = conditions.c_heat[insulation]
    xRho[-wallCells:] = conditions.rho[insulation]

    # y Dimension:
    # Wall Cells (Left in y) - Insulation Material
    yKappa[:wallCells] = conditions.kappa[insulation]
    yC_heat[:wallCells] = conditions.c_heat[insulation]
    yRho[:wallCells] = conditions.rho[insulation]

    # Interior Cells - Air
    yKappa[wallCells:-wallCells] = conditions.kappa[0]
    yC_heat[wallCells:-wallCells] = conditions.c_heat[0]
    yRho[wallCells:-wallCells] = conditions.rho[0]

    # Wall Cells (Right in y) - Insulation Material
    yKappa[-wallCells:] = conditions.kappa[insulation]
    yC_heat[-wallCells:] = conditions.c_heat[insulation]
    yRho[-wallCells:] = conditions.rho[insulation]

    # Define the thermal diffusivity constant at each spatial cell.
    xEta = xKappa * dt / (xC_heat * xRho * dx**2)
    yEta = yKappa * dt / (yC_heat * yRho * dx**2)

    """ Original Plan:
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
    """

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

    # Using the alternating direction implicit (ADI) method, the 2D problem be-
    # comes a case in which the finite difference problem is split into two sim-
    # plified problems: one in which the x-derivative is taken implicitly while
    # the other takes the y-derivative implicitly. The result is a system of
    # equations that are symmetric and tridiagonal, allowing the 1D method to
    # be applied twice: once in each spatial direction for a half-timestep.

    # Define the dimensionality of the square tridiagonal matrices. These two
    # values should be identical.
    xDimension = xSpatialCells - 2
    yDimension = ySpatialCells - 2

    # Define relevant constants when solving the matrix equation.
    xAlpha = 2/xEta + 2
    yAlpha = 2/yEta + 2

    xBeta = 2/xEta - 2
    yBeta = 2/yEta - 2

    # Explicitly define the tridiagonal matrices, xM_eta and yM_eta.
    xM_eta = np.diagflat(-np.ones(xDimension - 1), 1) + \
             np.diagflat(xAlpha[1:-1], 0) + \
             np.diagflat(-np.ones(xDimension - 1), -1)

    yM_eta = np.diagflat(-np.ones(yDimension - 1), 1) + \
             np.diagflat(yAlpha[1:-1], 0) + \
             np.diagflat(-np.ones(yDimension - 1), -1)

    # Initialize the terms at the intermediary timestep.
    xbT = np.zeros(xDimension)
    ybT = np.zeros(yDimension)

    # Define the inverse to the tridiagonal matrices once in order to improve
    # the efficiency of the algorithm since this will create the inverse once
    # rather than during each iteration using np.linalg.solve().
    inv_xM_eta = np.linalg.inv(xM_eta)
    inv_yM_eta = np.linalg.inv(yM_eta)

    # Define a token to represent the steps in time being stored for plotting
    # and initialize the initial time step to the initial conditions of the
    # system.
    t_index = 0
    timeEvolvingT[t_index, :, :] = T

    # Simulate!
    # The simulation is conducted for all timesteps not including the initial
    # time, t = 0, since that has been defined by the conditions.
    for t in range(1, temporalCells):
        # Apply the Crank-Nicolson algorithm in 1D over a half-timestep in the
        # y direction.
        for index in range(ySpatialCells):
            # Define the temperature values at the previous timestep that con-
            # tribute to the next half-timestep's temperatures.
            ybT[:] = T[:-2, index] + yBeta[1:-1] * T[1:-1, index] + T[2:, index]

            # Update the boundary values to account for the boundary at the next
            # half-timestep.
            ybT[0] += conditions.boundary(season, t * dt)
            ybT[-1] += conditions.boundary(season, t * dt)

            # Solve the matrix equation using the predefined inverse rather than
            # using np.linalg.solve().
            T[1:-1, index] = np.dot(inv_xM_eta, ybT)

            # Update the boundary temperatures.
            T[0, index] = T[-1, index] = conditions.boundary(season, t * dt)

        # Apply the Crank-Nicolson algorithm in 1D over the remaining half time-
        # step in the remaining x direction.
        for index in range(xSpatialCells):
            # Define the temperature values at the previous half-timestep that
            # contribute to the next half-timestep's temperatures.
            xbT[:] = T[index, :-2] + xBeta[1:-1] * T[index, 1:-1] + T[index, 2:]

            # Update the boundary values to account for the boundary at the next
            # half-timestep.
            xbT[0] += conditions.boundary(season, t * dt)
            xbT[-1] += conditions.boundary(season, t * dt)

            # Solve the matrix equation using the predefined inverse rather than
            # using np.linalg.solve().
            T[index, 1:-1] = np.dot(inv_yM_eta, xbT)

            # Update the boundary temperatures.
            T[index, 0] = T[index, -1] = conditions.boundary(season, t * dt)

        # In the case the timestep has reached an incremental value for the
        # timesteps being plotted or the final timestep in the simulation, the
        # temperature data is stored within the plotting data.
        if t % steps == 0 or t == temporalCells - 1:
            t_index += 1
            timeEvolvingT[t_index, :, :] = T

    # parameters = (dx, dy, dt, steps)
    return timeEvolvingT
