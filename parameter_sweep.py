# %% Import Libraries
import numpy as np
from IS_A_T_Calc_2 import IS_A_T_Calc
import matplotlib.pyplot as plt
from matplotlib import cm

# %% System Parameters
K = .01
tau = 2

ome_vec = np.pi * np.arange(1,100,10)
zet_vec = np.arange(0, 0.7, .1)
ome, zet = np.meshgrid(ome_vec, zet_vec)
A = np.zeros_like(ome)
T = np.zeros_like(ome)

calc = IS_A_T_Calc(tau, K)

# %% Simulate
for j in range(len(zet_vec)): # zet
       for i in range(len(ome_vec)): # ome
              A_i, T_i = calc.calculate_zet(ome[j,i], zet[j,i])
              A[j,i] = A_i
              T[j,i] = T_i

# %% Visualization
def style_3d_axes(ax):

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.line.set_color("black")
        axis.line.set_linewidth(2)

    ax.xaxis._axinfo["grid"].update({"linewidth":0.7, "linestyle":"--", "color":"lightgray"})
    ax.yaxis._axinfo["grid"].update({"linewidth":0.7, "linestyle":"--", "color":"lightgray"})
    ax.zaxis._axinfo["grid"].update({"linewidth":0.7, "linestyle":"--", "color":"lightgray"})

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.pane.set_edgecolor("black")
        axis.pane.set_alpha(0.0)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(projection='3d')

surf = ax.plot_surface(
    ome, zet, A,
    cmap=cm.plasma,
    linewidth=0,
    antialiased=True,
    alpha=0.7
)

ax.set(xlabel=r'$\omega$', ylabel=r'$\zeta$', zlabel=r'$A$')
fig.colorbar(surf, shrink=0.6, aspect=10, pad=0.1)
ax.view_init(elev=25, azim=135)
style_3d_axes(ax)   # apply styling

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(projection='3d')

surf = ax.plot_surface(
    ome, zet, T,
    cmap=cm.plasma,
    linewidth=0,
    antialiased=True,
    alpha=0.7
)

ax.set(xlabel=r'$\omega$', ylabel=r'$\zeta$', zlabel=r'$T$')
fig.colorbar(surf, shrink=0.6, aspect=10, pad=0.1)
ax.view_init(elev=25, azim=135)
style_3d_axes(ax)   # apply styling

# plt.show()

# %%
