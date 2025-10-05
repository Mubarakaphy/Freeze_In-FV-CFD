
import numpy as np
from scipy.integrate import solve_ivp
from .physics import RHS_of_x

def solve_reference(x_start=1e-3, x_end=1e3, n_plot=800):
    sol = solve_ivp(lambda x,y: RHS_of_x(x), [x_start, x_end], [0.0],
                    rtol=1e-8, atol=1e-12, dense_output=True)
    x_plot = np.logspace(np.log10(x_start), np.log10(x_end), n_plot)
    Y = sol.sol(x_plot)[0]
    return {"x_plot": x_plot, "Y": Y, "sol": sol}
