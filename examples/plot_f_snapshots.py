
import os, sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 120
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from freezein_fv.pde_solver import march_pde
from freezein_fv.physics import mA

def ensure_dir(path):
    os.makedirs(path, exist_ok=True); return path

def add_pstar_marker(ax, y_fraction=0.05, text=r"$p^\star \simeq m_A/2$"):
    pstar = 0.5 * mA
    ax.axvline(pstar, linestyle="--", linewidth=1)
    ymin, ymax = ax.get_ylim()
    ytext = 10**(np.log10(ymin) + y_fraction*(np.log10(ymax) - np.log10(ymin)))
    ax.text(pstar*1.05, ytext, text, fontsize=9, rotation=90, va="bottom")

def main():
    snaps = [1e-2, 3e-2, 1e-1, 3e-1, 1.0]
    out = march_pde(Np=600, Nx=2500, inj_width=0.02, scheme="muscl", time_integrator="rk2", x_switch=1e-2, snapshots_x=snaps)
    S = out.get("snapshots", {})
    if not S:
        print("No snapshots returned."); return
    items = sorted([(k, v) for k, v in S.items()], key=lambda kv: kv[1]["x"])
    out_dir = ensure_dir(os.path.join(os.path.dirname(__file__), "..", "docs", "figures"))

    # Fig 1: absolute p^2 f
    fig1, ax1 = plt.subplots(figsize=(9,5.2))
    for req_x, data in items:
        p, f, x_act = data["p"], data["f"], data["x"]
        ax1.loglog(p, p**2 * f, label=f"x≈{x_act:.2e} (req {req_x:.0e})")
    ax1.set_xlabel("p  [GeV]"); ax1.set_ylabel(r"$p^2 f_\chi(p)$  (arb.)")
    ax1.set_title("Momentum-space snapshots of $f_\chi(p)$")
    ax1.grid(True, which="both", ls=":"); ax1.legend()
    add_pstar_marker(ax1)
    fig1.tight_layout()
    f1_png = os.path.join(out_dir, "pde_f_snapshots.png"); f1_pdf = os.path.join(out_dir, "pde_f_snapshots.pdf")
    fig1.savefig(f1_png, dpi=150); fig1.savefig(f1_pdf)

    # Fig 2: normalized shapes
    fig2, ax2 = plt.subplots(figsize=(9,5.2))
    for req_x, data in items:
        p, f, x_act = data["p"], data["f"], data["x"]
        area = np.trapz(p**2 * f, p)
        ax2.loglog(p, (p**2 * f)/(area + 1e-300), label=f"x≈{x_act:.2e}")
    ax2.set_xlabel("p  [GeV]"); ax2.set_ylabel(r"$\frac{p^2 f_\chi(p)}{\int p^2 f_\chi\,dp}$")
    ax2.set_title("Snapshot shapes (normalized in $p^2$ measure)")
    ax2.grid(True, which="both", ls=":"); ax2.legend()
    add_pstar_marker(ax2)
    fig2.tight_layout()
    f2_png = os.path.join(out_dir, "pde_f_snapshots_shape.png"); f2_pdf = os.path.join(out_dir, "pde_f_snapshots_shape.pdf")
    fig2.savefig(f2_png, dpi=150); fig2.savefig(f2_pdf)

    # Fig 3: energy-density integrand p^3 f
    fig3, ax3 = plt.subplots(figsize=(9,5.2))
    for req_x, data in items:
        p, f, x_act = data["p"], data["f"], data["x"]
        ax3.loglog(p, p**3 * f, label=f"x≈{x_act:.2e}")
    ax3.set_xlabel("p  [GeV]"); ax3.set_ylabel(r"$p^3 f_\chi(p)$  (arb.)")
    ax3.set_title("Energy–density integrand snapshots")
    ax3.grid(True, which="both", ls=":"); ax3.legend()
    add_pstar_marker(ax3)
    fig3.tight_layout()
    f3_png = os.path.join(out_dir, "pde_f_snapshots_energy.png"); f3_pdf = os.path.join(out_dir, "pde_f_snapshots_energy.pdf")
    fig3.savefig(f3_png, dpi=150); fig3.savefig(f3_pdf)

    print("Saved ->", os.path.relpath(f1_png))
    print("Saved ->", os.path.relpath(f1_pdf))
    print("Saved ->", os.path.relpath(f2_png))
    print("Saved ->", os.path.relpath(f2_pdf))
    print("Saved ->", os.path.relpath(f3_png))
    print("Saved ->", os.path.relpath(f3_pdf))

if __name__ == "__main__":
    main()
