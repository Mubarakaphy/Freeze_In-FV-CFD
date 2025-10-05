
import numpy as np
from .physics import s_of_T, H_of_T, C_of_T, mA, RHS_of_x

def build_p_grid(pmin=1e-6, pmax=1e3, Np=600, log=True):
    edges = np.logspace(np.log10(pmin), np.log10(pmax), Np+1) if log else np.linspace(pmin, pmax, Np+1)
    centers = 0.5*(edges[:-1] + edges[1:])
    dp = edges[1:] - edges[:-1]
    return edges, centers, dp

def build_x_marks(xmin=1e-3, xmax=1e3, Nx=3000):
    return np.logspace(np.log10(xmin), np.log10(xmax), Nx)

def build_x_grid(xmin=1e-3, xmax=1e3, Nx=3000):
    return build_x_marks(xmin, xmax, Nx)

def narrow_gaussian(p, p0, width_frac=0.05):
    sigma = max(width_frac*p0, 1e-30)
    return np.exp(-0.5*((p-p0)/sigma)**2) / (np.sqrt(2*np.pi)*sigma)

def source_S(x, p, mchi, gchi=2.0, width_frac=0.05):
    T = mchi/x
    s = s_of_T(T); H = H_of_T(T); C = C_of_T(T)
    rhs_int = C/(H*x)
    pstar = 0.5*mA
    phi = narrow_gaussian(p, pstar, width_frac=width_frac)
    norm = np.trapz(p**2 * phi, p)
    S = (2*np.pi**2/gchi) * rhs_int * (phi / (norm + 1e-300))
    return S, pstar, s

def minmod(a, b):
    s = np.sign(a) + np.sign(b)
    return 0.5 * s * np.minimum(np.abs(a), np.abs(b))

def muscl_reconstruct(f, v_iface):
    Np = f.size
    df = np.diff(f)
    slope = np.zeros_like(f)
    if Np >= 3:
        slope[1:-1] = minmod(df[:-1], df[1:])
    slope[0] = 0.0; slope[-1] = 0.0
    left_face  = f - 0.5 * slope
    right_face = f + 0.5 * slope
    left_state  = np.empty(Np + 1)
    right_state = np.empty(Np + 1)
    left_state[0]  = left_face[0]; right_state[0] = left_face[0]
    left_state[1:Np]  = right_face[:-1]
    right_state[1:Np] = left_face[1:]
    left_state[Np]  = right_face[-1]; right_state[Np] = right_face[-1]
    assert v_iface.shape == (Np+1,)
    return np.where(v_iface >= 0.0, left_state, right_state)

def march_pde(mchi=10.0, gchi=2.0, Np=600, Nx=3000,
              pmin=1e-6, pmax=1e3, xmin=1e-3, xmax=1e3,
              inj_width=0.02,
              scheme="muscl", time_integrator="rk2",
              debug=False, x_switch=1e-2,
              snapshots_x=None):
    pedges, p, dp = build_p_grid(pmin, pmax, Np, log=True)
    x_marks = build_x_marks(xmin, xmax, Nx)

    # ODE pre-integration to x_switch
    Y0 = 0.0
    for i in range(len(x_marks)-1):
        x_i, x_ip1 = x_marks[i], x_marks[i+1]
        if x_ip1 > x_switch:
            if x_i < x_switch:
                dx = x_switch - x_i
                Y0 += RHS_of_x(x_i, mchi=mchi) * dx
            break
        dx = x_ip1 - x_i
        Y0 += 0.5*(RHS_of_x(x_i, mchi=mchi) + RHS_of_x(x_ip1, mchi=mchi)) * dx

    f = np.zeros_like(p)
    x_hist = [x_switch]; Y_hist = [Y0]; Y_src_hist = [Y0]

    # snapshots bookkeeping
    snap_targets = []
    if snapshots_x:
        sx = np.array(sorted(set(float(v) for v in snapshots_x)))
        snap_idx = np.searchsorted(x_marks, sx, side="left")
        snap_idx = np.clip(snap_idx, 0, len(x_marks)-1)
        snap_targets = list(zip(sx, snap_idx))
    snapshots = {}

    def compute_Y(T, f_arr):
        return (gchi/(2*np.pi**2*s_of_T(T))) * np.trapz(p**2 * f_arr, p)

    def rhs_operator(f_in, x_now):
        v_iface = -pedges / x_now
        if scheme.lower() == "muscl":
            f_up = muscl_reconstruct(f_in, v_iface)
        else:
            f_left  = np.concatenate(([f_in[0]], f_in))
            f_right = np.concatenate((f_in, [f_in[-1]]))
            f_up = np.where(v_iface >= 0.0, f_left, f_right)
        F = v_iface * f_up
        divF = (F[1:] - F[:-1]) / dp
        S, _, s_now = source_S(x_now, p, mchi, gchi=gchi, width_frac=inj_width)
        df_dx = -divF + S
        dY_src = (gchi/(2*np.pi**2*s_now)) * np.trapz(p**2 * S, p)
        return df_dx, dY_src

    start_idx = np.searchsorted(x_marks, x_switch, side="left")
    Y_src_total = Y0
    for i in range(start_idx, len(x_marks)-1):
        if snap_targets:
            to_keep = []
            for xv, idx in snap_targets:
                if i >= idx and xv not in snapshots:
                    snapshots[xv] = {"x": x_marks[i], "p": p.copy(), "f": f.copy()}
                else:
                    to_keep.append((xv, idx))
            snap_targets = to_keep

        x_now = x_marks[i]; x_next = x_marks[i+1]; dx = x_next - x_now
        if time_integrator.lower() == "rk2":
            df1, dY1 = rhs_operator(f, x_now)
            f_star = f + dx * df1
            df2, dY2 = rhs_operator(f_star, x_next)
            f = 0.5*(f + f_star + dx*df2)
            Y_src_total += 0.5*(dY1 + dY2) * dx
        else:
            df, dY = rhs_operator(f, x_now)
            f = f + dx * df
            Y_src_total += dY * dx
        Y_val = compute_Y(mchi/x_next, f)
        x_hist.append(x_next); Y_hist.append(Y_val); Y_src_hist.append(Y_src_total)

    return {"p": p, "x": np.array(x_hist), "f": f,
            "Y_hist": np.array(Y_hist), "Y_src_hist": np.array(Y_src_hist),
            "Y_final": (Y_hist[-1] if Y_hist else 0.0),
            "snapshots": snapshots}
