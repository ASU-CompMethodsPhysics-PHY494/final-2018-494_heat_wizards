# plotting.py
# All of the plotting code used in the analysis to render the various house
# conditions is located within this file.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_T(T_plot, Dx, Dt, step):
    X, Y = np.meshgrid(range(T_plot.shape[0]), range(T_plot.shape[1]))
    Z = T_plot[X, Y]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(X*Dt*step, Y*Dx, Z)
    ax.set_xlabel(r"time $t$ (s)")
    ax.set_ylabel(r"position $x$ (m)")
    ax.set_zlabel(r"temperature $T$ (K)")
    fig.tight_layout()
    plt.show()
    return ax
