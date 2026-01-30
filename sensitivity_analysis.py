# %% Import Libraries
import numpy as np
import matplotlib.pyplot as plt

import os
os.chdir("...")
print(os.getcwd())
print(os.listdir())

# Local
from Libraries.IS_A_T_Calc import IS_A_T_Calc

# %% System Parameters
K = .01
tau = 2
calc = IS_A_T_Calc(tau, K)
fixed_ome = 3 * np.pi  # Fixed omega 
fixed_zeta = 0.707  # Fixed zeta

ome_vec1 = np.pi * np.linspace(3,100,50)
zet_vec1 = np.arange(0, 1, 0.1)
A_ome1, A_zeta1 = np.zeros_like(ome_vec1), np.zeros_like(zet_vec1)
T_ome1, T_zeta1 = np.zeros_like(ome_vec1), np.zeros_like(zet_vec1)

for j in range(len(zet_vec1)): # zet
     A_zeta1[j], T_zeta1[j] = calc.calculate_zet(fixed_ome, zet_vec1[j])

for i in range(len(ome_vec1)): # ome
     A_ome1[i], T_ome1[i] = calc.calculate_zet(ome_vec1[i], fixed_zeta)

# %% Visualization
fig, ax1 = plt.subplots(1, 2, figsize=(12, 4))
ax1[0].plot(zet_vec1, A_zeta1, 'm-.', linewidth=2, label=r'$A(\zeta)\; | \; \omega_\mathrm{n} = 3\pi$')
ax1[0].plot(zet_vec1, T_zeta1, 'r-.', linewidth=2, label=r'$T(\zeta)\; | \; \omega_\mathrm{n} = 3\pi$')
ax1[0].set_xlim([0, 0.9])
# ax1[0].set_ylim([-2.5, 2.5])
ax1[0].set_xlabel(r"$\zeta \ \mathrm{[-]}$", fontsize=14)
ax1[0].set_ylabel(r"$A \ \mathrm{[-]}, \; T \ \mathrm{[s]}$", fontsize=14)
ax1[0].tick_params(axis='both', which='major', labelsize=14)
for spine in ax1[0].spines.values():
    spine.set_linewidth(2)
ax1[0].grid(True, which='both', linewidth=0.7, color='gray', alpha=0.5)
ax1[0].legend(fontsize=14)

ax1[1].plot(ome_vec1, A_ome1, 'b-.', linewidth=2, label=r'$A(\omega_\mathrm{n})\; | \; \zeta = 0.707$')
ax1[1].plot(ome_vec1, T_ome1, 'c-.', linewidth=2, label=r'$T(\omega_\mathrm{n})\; | \; \zeta = 0.707$')
ax1[1].set_xlim([3*np.pi, 100*np.pi])
ax1[1].set_xlabel(r"$\omega_\mathrm{n} \ \mathrm{[rad/s]}$", fontsize=14)
ax1[1].set_ylabel(r"$A \ \mathrm{[-]}, \; T \ \mathrm{[s]}$", fontsize=14)
ax1[1].tick_params(axis='both', which='major', labelsize=14)
for spine in ax1[1].spines.values():
    spine.set_linewidth(2)
ax1[1].grid(True, which='both', linewidth=0.7, color='gray', alpha=0.5)
ax1[1].legend(fontsize=14)

# fig.savefig("/Users/nyinyia/Documents/09_LSU_GIT/A_TDF/Figures/sen_A_T.pdf",
#             dpi=300,
#             bbox_inches="tight")

# %%