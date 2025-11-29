# %%
class SecondOrderSystem:
    def __init__(self, omega_n, zeta):
        self.omega_n = omega_n
        self.zeta = zeta
 
    def ODE(self, t, state, u):
        y, y_dot = state
        y_ddot = -2 * self.zeta * self.omega_n * y_dot - self.omega_n**2 * y + self.omega_n**2 * u
        return [y_dot, y_ddot]
    
    # %%