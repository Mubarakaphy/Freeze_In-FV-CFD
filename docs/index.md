## Results

### 1. Freeze-in yield comparison

This plot compares the **ODE benchmark** (solved with `scipy.solve_ivp`)  
against the **CFD-style finite-volume PDE solver**.  
The two curves agree extremely well for the “source-only” diagnostic,  
confirming that the finite-volume discretization conserves the freeze-in yield.

<p align="center">
  <img src="figures/pde_yield_comparison.png" width="700" alt="Freeze-in yield comparison plot">
</p>

---

### 2. Momentum-space snapshots

Snapshots of the DM phase-space distribution $f_\chi(p)$ at different $x = m_\chi/T$.  
As the Universe cools, the injection band around $p \simeq m_A/2$  
gradually redshifts and the distribution broadens.

<p align="center">
  <img src="figures/pde_f_snapshots.png" width="700" alt="Momentum-space snapshots f(p)">
</p>

---

### 3. Normalized distribution shapes

Here, $p^2 f_\chi(p)$ is normalized to unity, highlighting how the  
distribution evolves in shape rather than amplitude.

<p align="center">
  <img src="figures/pde_f_snapshots_shape.png" width="700" alt="Normalized momentum distribution">
</p>

---

### 4. Energy-weighted spectra

The $p^3 f_\chi(p)$ weighting emphasizes contributions to the comoving  
energy density, useful for comparing kinetic vs thermal regimes.

<p align="center">
  <img src="figures/pde_f_snapshots_energy.png" width="700" alt="Energy-weighted momentum spectra">
</p>

