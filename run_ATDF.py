# %% Import Libraries
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import os
os.chdir("...")
print(os.getcwd())
print(os.listdir())

# Local
from Estimator import Estimator
from IS_A_T_Calc import IS_A_T_Calc
from SecondOrderSystem import SecondOrderSystem

# %% System Parameters
# Zeta and Omega_n values to be tested
zeta_test = np.array([.1, .5, .707, 1.0])
# zeta_test = np.array([0.9])
# zeta_test = np.linspace(0.0,1.0,24)
omega_test = np.pi * np.array([1])

duration = 10
identification_duration = 2
dt = .0001
ref = 1

iden_ref = .01 * ref
time = np.arange(0,duration,dt)
iden_time = np.arange(0, identification_duration, dt)

# Storage Arrays
est_zeta_store = np.ndarray((len(zeta_test),len(omega_test)))
est_omega_n_store = np.ndarray((len(zeta_test),len(omega_test)))
A_store = np.ndarray((len(zeta_test),len(omega_test)))
T_store = np.ndarray((len(zeta_test),len(omega_test)))
err_omega_n_store = np.ndarray((len(zeta_test)))

x = np.ndarray((len(zeta_test), len(time)))
dx = np.ndarray((len(zeta_test), len(time)))
ref_vec = np.ndarray((len(zeta_test), len(time)))

# %% Simulate
for k, w in enumerate(omega_test):

    sys_omega_n = w
    # Loop thru zeta
    for j, z in enumerate(zeta_test):
        # State reset
        x[j,0] = 0
        dx[j,0] = 0
        state = np.array([x[j,0],dx[j,0]])

        # Initialize System
        sys_zeta = z
        sys = SecondOrderSystem(sys_omega_n, sys_zeta)
        print('Omega: ',w)
        print('Zeta: ',z)

        # Simulation Time Loop
        for i, t in enumerate(time[:-1]):
            # During Identification Time
            if t <= identification_duration:
                ref_vec[j,i] = iden_ref

            # End of Identification Time
            if t == iden_time[-1]:
                # print('Reached iden_time')
                iden_x = x[j,:len(iden_time)]
                estimate = Estimator(sys, iden_time, iden_ref*np.ones_like(iden_time))
                est_zeta, est_omega_n = estimate.estimate(iden_x)

                est_zeta_store[j,k] = est_zeta
                est_omega_n_store[j,k] = est_omega_n

                err_omega_n = np.abs((sys_omega_n-est_omega_n)/sys_omega_n)*100
                err_omega_n_store[j] = err_omega_n

                tau = identification_duration
                K = iden_ref
                A_T_Calc = IS_A_T_Calc(tau= tau, K= K)
                # Initialize IS for Zeta = 0 
                if est_zeta <= 0.02:
                    A, T = A_T_Calc.calculate_zet(est_omega_n, est_zeta)
                    # print('Zeta = 0')

                elif est_zeta <= .51:
                    # print('Zeta < 0.5')
                    A, T = A_T_Calc.calculate_zet(est_omega_n, est_zeta)
                    
                elif est_zeta >= .51 and est_zeta != 1:
                    # print('Zeta > 0.5')
                    A, T = A_T_Calc.calculate_zet(est_omega_n, est_zeta)
                else:
                    A = 0
                    T = 0
                    print('IS Disabled')
                    # if est_zeta == 1: print('Zeta = 1')

                A_store[j,k] = A
                T_store[j,k] = T
                    
            if t > time[len(iden_time)]:
                ref_vec[j,i] = (K + A*np.heaviside(t - (tau+T), 1/2) + (1-K-A)*np.heaviside(t - (tau+2*T), 1/2))*ref

            # Solve ODE for time step
            out = sp.integrate.solve_ivp(sys.ODE, [time[i], time[i+1]], state, args = (ref_vec[j,i],), rtol = 1e-8, atol = 1e-10)
            x[j,i+1] = out.y[0,-1]
            dx[j,i+1] = out.y[1,-1]
            state = [out.y[0,-1],out.y[1,-1]]

            if t == time[:-1][-1]:
                ref_vec[j,i+1] = ref_vec[j,i]

    # Plotting (9 Plots for zeta values for each omega value)
    fig, ax = plt.subplots(3, 3, figsize=(16,12))
    plt.suptitle(r'$\omega_{n}$ = ' + str(w/np.pi) + r'$\pi$', fontsize = 20)
    for ii in range(len(zeta_test)):
        row = ii // 3
        col = ii % 3
        ax[row, col].plot(time, ref_vec[ii,:], 'r-',linewidth = 2, label = r'$x^{*}$')
        ax[row, col].plot(time, x[ii,:], 'b-',linewidth = 2, label = r'$x$')
        ax[row, col].set_xlim([0, duration])
        ax[row, col].grid()
        ax[row, col].set_title(r'Zeta = '+ str(zeta_test[ii]), fontsize = 14)
        ax[row, col].set_xlabel(r"$t \ \mathrm{[s]}$", fontsize = 14)
        ax[row, col].set_ylabel(r'$x$ [m]', fontsize = 14)
        ax[row, col].legend()
    plt.tight_layout()

# %%