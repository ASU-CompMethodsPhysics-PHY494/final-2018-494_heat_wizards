# simulation.py
# This file contains all of the relevant functions necessary to conduct simu-
# lations for various house conditions, whether the variable is insulation
# material, wall thickness, or house size/shape.

import numpy as np

# Import associated files.
import conditions

def CrankNicolson_T(L_rod=1, t_max=3000, Dx=0.02, Dt=2, T0=293, Tb=311,
                     step=20, verbose=True):
    Nx = int(L_rod // Dx)
    Nt = int(t_max // Dt)

    #Kappa = 237 # W/(m K)
    #CHeat = 900 # J/K
    #rho = 2700  # kg/m^3

    Kappa = 2.624 # W/(m K)
    CHeat = 1000  # J/K
    rho = 1.204   # kg/m^3

    eta = Kappa * Dt / (CHeat * rho * Dx**2)

    if verbose:
        print("Nx = {0}, Nt = {1}".format(Nx, Nt))
        print("eta = {0}".format(eta))

    T = np.zeros(Nx)
    T_plot = np.zeros((int(np.ceil(Nt/step)) + 1, Nx))

    # initial conditions
    T[1:-1] = T0
    # boundary conditions
    T[0] = T[-1] = Tb

    #---------------------
    # M_eta * T[1:-1, j+1] = bT
    # M_eta * xT = bT
    # Nx-2 x Nx-2 matrix: tridiagonal
    NM = Nx - 2
    alpha = 2/eta + 2
    beta = 2/eta - 2
    M_eta = np.diagflat(-np.ones(NM-1), 1) \
            + np.diagflat(alpha * np.ones(NM), 0) \
            + np.diagflat(-np.ones(NM-1), -1)
    bT = np.zeros(NM)

    t_index = 0
    T_plot[t_index, :] = T
    for jt in range(1, Nt):
        bT[:] = T[:-2] + beta*T[1:-1] + T[2:]
        # boundaries are special cases
        bT[0] += T[0]  #  + T[0,j+1]
        bT[-1] += T[-1] #  + T[-1,j+1]

        T[1:-1] = np.linalg.solve(M_eta, bT)

        if jt % step == 0 or jt == Nt-1:
            t_index += 1
            T_plot[t_index, :] = T
            if verbose:
                print("Iteration {0:5d}".format(jt), end="\r")
    else:
        if verbose:
            print("Completed {0:5d} iterations: t={1} s".format(jt, jt*Dt))

    parameters = (Dx, Dt, step)
    return T_plot, parameters

def CrankNicolson_1D(length, time, dx, dt, insulation, wall_thickness, step = 20):
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
    for cellIndex in range(spatialCells):
        # Wall Cells (Left)
        if cellIndex < wallCells:
            kappa[cellIndex] = 237
            c_heat[cellIndex] = 900
            rho[cellIndex] = 2700
        # Interior Cells
        elif cellIndex < (wallCells + interiorCells):
            kappa[cellIndex] = 2.624
            c_heat[cellIndex] = 1000
            rho[cellIndex] = 1.204
        # Wall Cells (Right)
        else:
            kappa[cellIndex] = 237
            c_heat[cellIndex] = 900
            rho[cellIndex] = 2700

    # Define the thermal diffusivity constant at each spatial cell.
    eta = kappa * dt / (c_heat * rho * dx**2)

    # Create an array to store the temperature at each spatial cell and create
    # an additional 2D array to store the temperature values in space and time
    # for plotting purposes.
    T = np.zeros(spatialCells)
    timeEvolvingT = np.zeros((int(np.ceil(temporalCells / step)) + 1, spatialCells))

    # Initial Condition(s):
    T[1:-1] = 293 # About room temp

    # Boundary Condition(s):
    T[0] = T[-1] = 311 # 100 degrees F

    #--------------------------------------------------------------------------#

    # Define the tridiagonal matrix representing the system of equations solving
    # the problem. The matrix has the following dimensions,
    #
    #   (spatialCells - 2) x (spatialCells - 2)
    #
    # and solves the equation,
    #
    #   M_eta * T[1:-1, t+1] = bT

    dimension = spatialCells - 2

    alpha = 2/eta + 2
    beta = 2/eta - 2

    M_eta = np.diagflat(-np.ones(dimension - 1), 1) + \
            np.diagflat(-np.ones(dimension - 1), -1)
    for i in range(dimension - 1):
        M_eta[i, i] = alpha[i + 1]

    bT = np.zeros(dimension)

    inv_M_eta = np.linalg.inv(M_eta)

    t_index = 0
    timeEvolvingT[t_index, :] = T

    for jt in range(1, temporalCells):
        bT[:] = T[:-2] + T[2:]
        for i in range(dimension - 1):
            bT[i] = beta[i + 1] * T[i + 1]

        bT[0] += 311 + np.sin(jt + 1)
        bT[-1] += 311 + np.sin(jt + 1)

        T[1:-1] = np.dot(inv_M_eta, bT)

        if jt % step == 0 or jt == temporalCells - 1:
            t_index += 1
            timeEvolvingT[t_index, :] = T

    parameters = (dx, dt, step)
    return timeEvolvingT, parameters
