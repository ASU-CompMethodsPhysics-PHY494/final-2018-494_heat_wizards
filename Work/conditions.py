# conditions.py
# This file contains all of the boundary and initial conditions that will be
# used to define the various problems being simulated.

import numpy as np

# Constants:
# Note: All values documented below and used within the simulations are defined
#       in terms of SI base units.
#
#           Distance    = Meters (m)
#           Time        = Seconds (s)
#           Temperature = Kelvin (K)

#====================#
# Default Parameters #
#====================#
time = 48 * 3600    # 48-Hour Simulation Period

# Discretization
dx = 0.1
dy = 0.1

dw = 0.01

dt = 2

# Plotting Steps
steps = 20

#=====================#
# Material Properties #
#=====================#

# Materials List:
materialList = \
[
    'air',
    'brick',
    'wood',
    'copper'
]

# Enumeration:
class material:
    air, brick, wood, copper = range(len(materialList))

# Initialize Data Containers:
kappa = np.zeros(len(materialList))
c_heat = np.zeros(len(materialList))
rho = np.zeros(len(materialList))

# Thermal Conductivity:
kappa[material.air] = 0.02624
kappa[material.brick] = 0.8
kappa[material.wood] = 0.17
kappa[material.copper] = 401

# Specific Heat:
c_heat[material.air] = 1000
c_heat[material.brick] = 900
c_heat[material.wood] = 2000
c_heat[material.copper] = 390

# Density:
rho[material.air] = 1.177
rho[material.brick] = 1900
rho[material.wood] = 750
rho[material.copper] = 8790

#======================#
# Initial and Boundary #
#      Conditions      #
#======================#

# Initial Interior Temperature:
T0 = 75 # Room Temperature (degrees F)

# The average high and low temperatures are taken from a collection of data sta-
# tistically averaged from 1981 until 2010 in Phoenix, AZ. In particular, the
# summer temperatures are taken from an average over June to September while the
# winter temperatures are taken from an average over November to February based
# on the seasons in the Northern Hemisphere.
summer_high = 103.5 # Summer High Temperature (degrees F)
summer_low = 79.25  # Summer Low Temperature (degrees F)

winter_high = 71    # Winter High Temperature (degrees F)
winter_low = 48     # Winter Low Temperature (degrees F)

# Convert from Fahrenheit to Kelvin:
def FahrenheitToKelvin(temperatureF):
    temperatureK = (5/9) * (temperatureF + 459.67)
    return temperatureK

# Initial Interior Temperature:
T0 = FahrenheitToKelvin(T0)

# Seasonal Temperature Extremes:
summer_high = FahrenheitToKelvin(summer_high)
summer_low = FahrenheitToKelvin(summer_low)

winter_high = FahrenheitToKelvin(winter_high)
winter_low = FahrenheitToKelvin(winter_low)

# Temperature Averages:
summer_average = (summer_high + summer_low) / 2
winter_average = (winter_high + winter_low) / 2

average = np.array([summer_average, winter_average])

# Temperature Variations/Amplitudes:
summer_amplitude = (summer_high - summer_low) / 2
winter_amplitude = (winter_high - winter_low) / 2

amplitude = np.array([summer_amplitude, winter_amplitude])

# Temperature Offset:
# The high temperatures occur near 5 pm and the low temperatures occur near 3
# am; however for convenience, a 24-hour period is used, setting the high temp-
# erature to occur at 4 pm while the low temperature is set to occur at 4 am.
delta = (2/3) * np.pi

# Boundary conditions:
def boundary(season, t):
    return average[season] + amplitude[season] * np.cos((2*np.pi*t)/(24*3600) + delta)
