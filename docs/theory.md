# Theory Note — Finite-Volume Boltzmann Solver for Freeze-in Dark Matter

## 1. Physical Background

In the early Universe, the number density of a feebly interacting dark matter species χ evolves according to the Boltzmann equation

$$
\frac{dn_\chi}{dt} + 3H n_\chi = C(T),
$$

where $H(T)$ is the Hubble expansion rate and $C(T)$ is the collision (source) term, typically arising from the decay or annihilation of heavier particles in equilibrium.

It is customary to express the abundance in terms of the **yield**

$$
Y \equiv \frac{n_\chi}{s},
$$

where $s(T)$ is the entropy density. Using $x \equiv m_\chi / T$ as an evolution variable, one obtains the ordinary differential equation

$$
\frac{dY}{dx} = \frac{1}{s(T)} \frac{C(T)}{H(T) \, x}.
$$

This equation describes **freeze-in** production when the initial abundance is negligible and $Y$ gradually increases as the Universe cools.

---

## 2. Finite-Volume PDE Formulation

Rather than solving for $Y$ directly, we evolve the full phase-space distribution $f_\chi(p, x)$ governed by the comoving Boltzmann equation

$$
\frac{\partial f_\chi}{\partial x}+ \frac{\partial}{\partial p} \big[ (-p/x) f_\chi \big] = S(x, p),
$$

where $S(x,p)$ is the source term describing injection of χ particles with momentum $p$ at temperature $T = m_\chi/x$.

Integrating over momentum gives back the ODE for $Y$:

$$
Y(x) = \frac{g_\chi}{2\pi^2 s(T)} \int_0^\infty p^2 f_\chi(p, x) \, dp,
$$

and consistency requires

$$
\int_0^\infty p^2 S(x,p)\, dp = \frac{2\pi^2}{g_\chi} \frac{C(T)}{H(T)\,x}.
$$

---

## 3. Numerical Discretization

We discretize $(p, x)$ space using a **finite-volume (FV)** scheme inspired by Computational Fluid Dynamics (CFD):

- Logarithmic grids in $p$ and $x$.
- Conservative update of cell-averaged quantities.
- Upwind fluxes reconstructed using the **MUSCL (minmod)** limiter.
- Explicit marching in $x$ analogous to time-integration in advection problems.

This approach ensures numerical conservation of the integrated yield $Y$ and allows the solver to capture redshifting and injection features in momentum space.

---

## 4. Analogy to CFD

The PDE form of the Boltzmann equation is structurally identical to the **advection equation**

$$
\frac{\partial f}{\partial t} + \nabla \cdot (\mathbf{v} f) = S,
$$

where here:
- $x$ acts as a pseudo-time variable (temperature evolution),
- $p$ plays the role of spatial coordinate,
- $v_p = -p/x$ is the “velocity field” in momentum space.

This analogy enables direct reuse of CFD algorithms and opens the path to implementing cosmological Boltzmann solvers in existing CFD frameworks such as **OpenFOAM** or **COMSOL**.

---

## 5. References and Further Reading

1. J. McDonald, *Thermally generated gauge singlet scalars as self-interacting dark matter*, Phys. Rev. Lett. **88**, 091304 (2002).  
2. L. Hall, K. Jedamzik, J. March-Russell, and S. West, *Freeze-In Production of FIMP Dark Matter*, JHEP **03** (2010) 080.  
3. C. Patrignani et al. (Particle Data Group), *Cosmological Parameters*, *Prog. Theor. Exp. Phys.* **2020**, 083C01.  
4. Hirsch et al., *Computational Methods for Fluid Dynamics*, Springer (2015).  
5. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, Springer (2009).  

---

> **Summary:**  
> This note outlines the bridge between cosmological kinetic theory and CFD.  
> By treating $f_\chi(p, x)$ as a transported scalar field in momentum space,  
> we demonstrate that early-Universe dark matter dynamics can be evolved using  
> the same finite-volume principles used in classical fluid mechanics.
