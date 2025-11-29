# %% Import Libraries
import numpy as np
import sympy as sy

# %%
class IS_A_T_Calc:
    def __init__(self, tau, K):
        self.tau = tau
        self.K = K

    def calculate_zet0(self,omgn):
        self.omgn = omgn
        phi = self.omgn * self.tau
        A = (self.K + 3*self.K*np.cos(phi) - 1)/(3*self.K*np.sin(phi))
        B = (self.K + 3*self.K*np.cos(phi) - 1)**2/(9*self.K**2*np.sin(phi)**2) + 1
        C = (self.K + 3*self.K*np.cos(phi) - 1)/(2*self.K*np.sin(phi)) - \
            (self.K*np.cos(phi) - self.K + 1)/(2*self.K*np.sin(phi)) + \
            (self.K + 3*self.K*np.cos(phi) - 1)**3/(27*self.K**3*np.sin(phi)**3)
        D = (C**2 - B**3)**(1/2)
        U = np.cbrt(C + D)

        zsol_vec = np.zeros(4, dtype=complex)
        zsol_vec[0] = 0
        zsol_vec[1] = A + B/U + U
        zsol_vec[2] = A - U/2 - B/(2*U) + (3**(1/2)*(U - B/U)*1j)/2
        zsol_vec[3] = A - U/2 - B/(2*U) - (3**(1/2)*(U - B/U)*1j)/2

        constraint_vec = np.ones(len(zsol_vec))
        for ind in range(len(zsol_vec)):
            if np.real(zsol_vec[ind]) == 0 and np.imag(zsol_vec[ind]) == 0:
                constraint_vec[ind] = 0
            elif np.real(zsol_vec[ind]) == 0:
                constraint_vec[ind] = 1        
            else:
                constraint_vec[ind] = np.abs(np.imag(zsol_vec[ind]))/np.abs(np.real(zsol_vec[ind]))

        self.zsol_real = np.zeros(len(zsol_vec))
        for i in range(len(constraint_vec)):
            if constraint_vec[i] < 1e-6:
                self.zsol_real[i] = np.real(zsol_vec[i])

        def fct_T(k):
            return 2*np.atan(self.zsol_real)/self.omgn + 2*np.pi*k/self.omgn

        k_vec = np.linspace(-3,3,7)
        bit_save = 1

        for i in range(len(k_vec)):
            self.T_temp = fct_T(k_vec[i])

            # Select smallest positive T
            if sum(self.T_temp > 0) > 0 and bit_save == 1:
                if sum(self.T_temp > 0) > 1:
                    T = min(self.T_temp[self.T_temp > 0])
                else:
                    T = np.float64(self.T_temp[self.T_temp > 0][0])

                if (np.cos(self.omgn*(T + self.tau)) - np.cos(self.omgn*(2*T + self.tau))) == 0:
                    A = 0.5
                else:
                    A = -(self.K - np.cos(self.omgn*(2*T + self.tau))*(self.K - 1))/(np.cos(self.omgn*(T + self.tau)) \
                        - np.cos(self.omgn*(2*T + self.tau)))
                    
                bit_save = 0

        return A, T
    
    def calculate_zet(self, omgn,zeta):
        self.omgn = omgn
        self.zet = zeta
        A,T = sy.symbols('A,T')
        # print('Start Solving')
        
        real = self.K + \
            A*sy.exp(self.zet*self.omgn*(self.tau+T))*sy.cos(self.omgn*(self.tau+T)*sy.sqrt(1-self.zet**2)) + \
            (1-self.K-A)*sy.exp(self.zet*self.omgn*(self.tau+2*T))*sy.cos(self.omgn*(self.tau+2*T)*sy.sqrt(1-self.zet**2))

        imag = A*sy.exp(self.zet*self.omgn*(self.tau+T))*sy.sin(self.omgn*(self.tau+T)*sy.sqrt(1-self.zet**2)) + \
                (1-self.K-A)*sy.exp(self.zet*self.omgn*(self.tau+2*T))*sy.sin(self.omgn*(self.tau+2*T)*sy.sqrt(1-self.zet**2))

        # print('T_approx')
        T_approx = np.pi/self.omgn

        for iii in range(100):    
            # print('Test loop')
            A_sol, T_sol = sy.nsolve([real,imag], [A, T], [0,(1.2 + .1*iii)*T_approx], verify = False)
            A_ret = np.float64(A_sol)
            T_ret = np.float64(T_sol)
            
            if  A_ret > 0 and A_ret < (1-self.K) and T_ret > 0.4*T_approx:
                break

        return A_ret, T_ret
    
# %%
