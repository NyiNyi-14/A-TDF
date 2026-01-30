# %% Import Libraries
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import os
os.chdir("...")
print(os.getcwd())
print(os.listdir())

# Local
from Libraries.Estimator import Estimator
from Libraries.IS_A_T_Calc import IS_A_T_Calc
from Libraries.SecondOrderSystem import SecondOrderSystem

# %% System Parameters
# Zeta and Omega_n values to be tested
# zeta_test = np.array([.1, .5, .707, 1.0])
zeta_test = np.array([0.707])
# zeta_test = np.linspace(0.0,1.0,24)
omega_test = np.pi * np.array([30])

# Time Inputs
duration = 12
identification_duration = 2
dt = .01
# ref = 1

# iden_ref = .01 * ref
time = np.arange(0,duration,dt)
iden_time = np.arange(0, identification_duration, dt)

R = np.ones(time.shape[0])
R[int(4/dt):int(6/dt)] = 2
R[int(6/dt):int(8/dt)] = 0
R[int(8/dt):int(10/dt)] = -2
R[int(10/dt):] = -1
ref = R
iden_ref = .01 * ref

# Storage Arrays
est_zeta_store = np.ndarray((len(zeta_test),len(omega_test)))
est_omega_n_store = np.ndarray((len(zeta_test),len(omega_test)))
A_store = np.ndarray((len(zeta_test),len(omega_test)))
T_store = np.ndarray((len(zeta_test),len(omega_test)))
err_omega_n_store = np.ndarray((len(zeta_test)))

x = np.ndarray((len(zeta_test), len(time)))
dx = np.ndarray((len(zeta_test), len(time)))
ref_vec = np.ndarray((len(zeta_test), len(time)))

RR = np.zeros((1, len(time)))

R1 = np.zeros((1, len(time)))
R2 = np.zeros((1, len(time)))
R3 = np.zeros((1, len(time)))
R4 = np.zeros((1, len(time)))
R5 = np.zeros((1, len(time)))

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
                ref_vec[j,i] = iden_ref[0]

            # End of Identification Time
            if t == iden_time[-1]:
                # print('Reached iden_time')
                iden_x = x[j,:len(iden_time)]
                estimate = Estimator(sys, iden_time, iden_ref[0]*np.ones_like(iden_time))
                est_zeta, est_omega_n = estimate.estimate(iden_x)

                est_zeta_store[j,k] = est_zeta
                est_omega_n_store[j,k] = est_omega_n

                err_omega_n = np.abs((sys_omega_n-est_omega_n)/sys_omega_n)*100
                err_omega_n_store[j] = err_omega_n

                tau = identification_duration
                K = iden_ref[0]
                A_T_Calc = IS_A_T_Calc(tau= tau, K= K)
                # Initialize IS for Zeta = 0 

                A, T = A_T_Calc.calculate_zet(est_omega_n, est_zeta)

                A_store[j,k] = A
                T_store[j,k] = T
                # print(str(A))
                # print(str(T))
                # print('==========')
                    
            if t > time[len(iden_time)]:
                # ref_vec[j,i] = (K + A*np.heaviside(t - (tau+T), 1/2) + (1-K-A)*np.heaviside(t - (tau+2*T), 1/2))*ref[i]

                step_wise = np.diff(ref)
                indices = np.where(step_wise)[0] + 1
                time_ch = time[indices]
                value_ch = ref[indices]
                # print("Change times:", time_ch)
                # print("Change values:", value_ch)

                ref_vec[j,i] = (K + A*np.heaviside(t - (tau+T), 1/2) + (1-K-A)*np.heaviside(t - (tau+2*T), 1/2))*ref[i]
                R1[:, i] = (K + A*np.heaviside(t - (tau+T), 1/2) + (1-K-A)*np.heaviside(t - (tau+2*T), 1/2))*ref[0]
                R2[:,i] = R1[:, i] + (A*np.heaviside(t - (time_ch[0]+T), 1/2) + (1-A)*np.heaviside(t - (time_ch[0]+2*T), 1/2))*R1[:, i]
                R3[:,i] = R2[:, i] + (ref[indices[1]] - ref[indices[0]]) * (A*np.heaviside(t - (time_ch[1]+T), 1/2) + (1-A)*np.heaviside(t - (time_ch[1]+2*T), 1/2))
                R4[:,i] = R3[:, i] + (ref[indices[2]] - ref[indices[1]]) * (A*np.heaviside(t - (time_ch[2]+T), 1/2) + (1-A)*np.heaviside(t - (time_ch[2]+2*T), 1/2))
                R5[:,i] = R4[:,i] + (ref[indices[3]] - ref[indices[2]]) * (A*np.heaviside(t - (time_ch[3]+T), 1/2) + (1-A)*np.heaviside(t - (time_ch[3]+2*T), 1/2))

                ref_vec = R5

            # Solve ODE for time step
            out = sp.integrate.solve_ivp(sys.ODE, [time[i], time[i+1]], state, args = (ref_vec[j,i],), rtol = 1e-8, atol = 1e-10)
            x[j,i+1] = out.y[0,-1]
            dx[j,i+1] = out.y[1,-1]
            state = [out.y[0,-1],out.y[1,-1]]

            if t == time[:-1][-1]:
                ref_vec[j,i+1] = ref_vec[j,i]

# plt.plot(time, R)
# plt.plot(time, R5[0,:], 'g--')

# %% Visualization
fig, ax1 = plt.subplots(1, 1, figsize=(10, 4))
ax1.plot(time, R, 'm-.', linewidth=2, label=r'$x^{*}$')
ax1.plot(time, ref_vec[0, :], 'r--', linewidth=2, label=r'$x_{\mathrm{TDF}}^{*}$')
ax1.plot(time, x[0, :], 'b-', linewidth=2, label=r'$x$')
ax1.axvline(
    x=identification_duration,
    color='lightblue',
    linestyle='-',
    linewidth=2,
)
ax1.axvspan(0, identification_duration,
            color='lightblue', alpha=0.3, label="Identification Phase")
ax1.set_xlim([0, duration])
ax1.set_ylim([-2.5, 2.5])
ax1.set_xlabel(r"$t \ \mathrm{[s]}$", fontsize=14)
ax1.set_ylabel(r"$x \ \mathrm{[m]}$", fontsize=14)
ax1.tick_params(axis='both', which='major', labelsize=14)

for spine in ax1.spines.values():
    spine.set_linewidth(2)

ax1.grid(True, which='both', linewidth=0.7, color='gray', alpha=0.5)
ax1.legend(fontsize=14)

# %%