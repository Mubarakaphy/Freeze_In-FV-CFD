import os, sys
import numpy as np
import matplotlib.pyplot as plt

# Make local src importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from freezein_fv.pde_solver import march_pde
from freezein_fv.ivp_reference import solve_reference

def main():
    # plotting + solver options
    x_switch = 1e-2

    # Run PDE (fixed-step, no CFL)
    out = march_pde(
        Np=800, Nx=5000,
        inj_width=0.02,
        scheme="muscl", time_integrator="rk2",
        x_switch=x_switch
    )
    x_pde   = out["x"]
    Y_src   = out["Y_src_hist"]   # source-only diagnostic (should match ODE)
    Y_fromf = out["Y_hist"]       # from integrating f (may drift)

    # ODE reference
    ref = solve_reference()
    x_ref, Y_ref = ref["x_plot"], ref["Y"]

    # ------------------ Plot results (publication style) -----------------------------
    fig = plt.figure(figsize=(8.2, 4.8))
    ax = plt.gca()

    # main curves
    ax.loglog(x_ref, Y_ref, lw=2.0, label="ODE (solve_ivp)", color="C0")
    ax.loglog(x_pde, Y_src, "--", lw=2.0, label="PDE (source-only diagnostic)", color="C1")
    ax.loglog(x_pde, Y_fromf, ":", lw=2.0, label="PDE (integrated f)", color="C2")

    # transition marker
    ax.axvline(x_switch, color="C7", ls="--", lw=1.0)
    # annotate a bit above the curve where we switch
    idx = np.searchsorted(x_ref, x_switch, side="left")
    y_here = Y_ref[min(idx, len(Y_ref)-1)]
    ax.text(x_switch*1.05, y_here*1.3, "ODE → PDE\ntransition", fontsize=9, va="center")

    # labels / title
    ax.set_xlabel(r"$x = m_\chi/T$", fontsize=12)
    ax.set_ylabel(r"$Y = n_\chi/s$", fontsize=12)
    ax.set_title("Yield comparison", fontsize=13)

    # legend + style tweaks
    ax.legend(frameon=False, fontsize=10, loc="lower right")
    ax.grid(False)                         # remove background grid entirely
    ax.set_facecolor("white")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    fig.tight_layout()

    # save
    out_dir = os.path.join(os.path.dirname(__file__), "..", "docs", "figures")
    os.makedirs(out_dir, exist_ok=True)
    fig_path = os.path.join(out_dir, "pde_yield_comparison.png")
    plt.savefig(fig_path, dpi=150)
    print(f"Saved -> {os.path.relpath(fig_path)}")

    # summary
    print("Final Y (ODE)          :", Y_ref[-1])
    print("Final Y (PDE src-only) :", Y_src[-1])
    print("Final Y (PDE integrate):", Y_fromf[-1])

    # show if running interactively
    # plt.show()

if __name__ == "__main__":
    main()
