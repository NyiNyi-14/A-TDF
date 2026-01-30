# %% Import Libraries
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from Libraries.PerformanceMetrics import PerformanceMetrics

# %%
class Estimator:

    def __init__(self,System,time_vec,ref_vec,):
        self.System = System
        self.duration = 10
        self.dt = time_vec[1] - time_vec[0]
        self.time = time_vec
        self.ref = ref_vec

    def estimate(self,out):        
        self.out = out
        self.peak_indices, _ = find_peaks(self.out)
        self.min_indicies, _ = find_peaks(-self.out)
        self.peak_values = self.out[self.peak_indices]
        self.min_values = self.out[self.min_indicies]
        self.peak_times = self.time[self.peak_indices]

        if np.size(self.peak_values) == 0:
            self.zet_est = 1
            eval = PerformanceMetrics(self.time, self.ref)
            ts = eval.settling_time(self.out)
            self.settling_time = ts
            # how to find omega_n
            self.ome_est = 4/ (self.zet_est * self.settling_time)

        elif np.size(self.peak_values) >= 2 and np.abs((self.peak_values[1] - self.peak_values[0])/self.peak_values[0]) <= .0001:
            self.zet_est = 0
            self.min_times = self.time[self.min_indicies]

            self.t_delta_vec = np.zeros(len(self.peak_values)-1)
            for i in range(len(self.peak_values)-1):
                self.t_delta_vec[i] = self.peak_times[i+1] - self.peak_times[i]
            pos_t_delta = np.mean(self.t_delta_vec)

            self.t_delta_vec = np.zeros(len(self.min_values)-1)
            for i in range(len(self.min_values)-1):
                self.t_delta_vec[i] = self.min_times[i+1] - self.min_times[i]
            neg_t_delta = np.mean(self.t_delta_vec)

            self.t_delta = (pos_t_delta + neg_t_delta)/2
            self.ome_est = 2*np.pi/self.t_delta

        else:
            Mp = (self.peak_values[0] - self.ref[-1])/ self.ref[-1]
            self.zet_est = -np.log(Mp) / np.sqrt(np.pi**2 + np.log(Mp)**2)
            eval = PerformanceMetrics(self.time, self.ref)
            ts = eval.settling_time(self.out)
            self.settling_time = ts
            # how to find omega_n
            self.ome_est = np.pi / (self.peak_times[0] * np.sqrt(1-self.zet_est**2))

        return [self.zet_est, self.ome_est]
    

    def run_test(self):
        state = np.array([0.0, 0.0])  # y(0), y_dot(0)
        states = []

        for i, t in enumerate(self.time):
            u = self.ref[i]  # Step input
            sol = solve_ivp(self.System.ODE, [t, t+self.dt], state, args=(u,))
            state = sol.y[:, -1]
            states.append(state)

        out = np.array(states)
        self.out = out

        return out
    
    def graph(self):
        plt.figure(figsize=(8, 6))
        plt.plot(self.time, self.out[:,0],'b-', label=r"$x$",linewidth=2)

        plt.xlim(0,self.time[-1])
        plt.xlabel("Time [s]",fontsize = 20)
        plt.ylabel(r"$x$",fontsize = 20)
        # plt.legend(fontsize = 20)
        plt.grid()
        plt.tight_layout()
        plt.tick_params(axis='both', labelsize=20)

# %%