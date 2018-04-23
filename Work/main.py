# main.py
# All of the analysis is functionally called upon within this file, compiling
# all of the results into a neat and readable manner.

import simulation as sim
import plotting as plot

#T_plot, (dx, dt, step) = sim.CrankNicolson_T(L_rod = 10, t_max = 7200, Dx = 0.1)
T_plot, (dx, dt, step) = sim.CrankNicolson_1D(50, 172800, 0.1, 2, None, 1)
plot.plot_T(T_plot, dx, dt, step)
