
import numpy as np

Mpl   = 2.435e18
mA    = 100.0
mChi  = 10.0
GammaA= 1e-22
gA    = 2.0
g_star   = 100.0
g_star_s = 100.0
gchi  = 2.0

def s_of_T(T):
    return (2*np.pi**2/45.0) * g_star_s * T**3

def H_of_T(T):
    return np.sqrt((4*np.pi**3*g_star/45.0)) * T**2 / Mpl

def n_eq_MB(m, T, g):
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        val = g * (m*T/(2*np.pi))**1.5 * np.exp(-m/T)
        val = np.where(np.isfinite(val), val, 0.0)
    return val

def C_of_T(T):
    return 2.0 * n_eq_MB(mA, T, gA) * GammaA

def RHS_of_x(x, mchi=None, gchi=None):
    if mchi is None:
        mchi = globals().get("mChi", 10.0)
    T = mchi / x
    T = np.maximum(T, 1e-300)
    return (C_of_T(T) / s_of_T(T)) / (H_of_T(T) * x)
