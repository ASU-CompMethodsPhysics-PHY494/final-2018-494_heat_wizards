# simulation.py
# This file contains all of the relevant functions necessary to conduct simu-
# lations for various house conditions, whether the variable is insulation
# material, wall thickness, or house size/shape.
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
%matplotlib notebook

%matplotlib inline

L_rod = 1.    # m
t_max = 3000. # s

Dx = 0.02   # m
Dt = 2    # s

Nx = int(L_rod // Dx)
Nt = int(t_max // Dt)

Kappa =  # W/(m K)
CHeat =  # J/K
rho =   # kg/m^3

T0 =  # K
Tb =  # K


eta = Kappa * Dt/(CHeat * rho * Dx*Dx)


step = 20  # plot solution every n steps

print("Nx = {0}, Nt = {1}".format(Nx, Nt))
print("eta = {0}".format(eta))

T = np.zeros(Nx)
T_new = np.zeros_like(T)
T_plot = np.zeros((Nt//step + 1, Nx))

# initial conditions
T[1:-1] = T0 
# boundary conditions
