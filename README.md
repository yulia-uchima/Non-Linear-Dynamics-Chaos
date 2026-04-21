# Non-Linear Dynamics & Chaos: Numerical Solvers

This repository contains a collection of numerical simulations and solvers for non-linear dynamical systems, focusing on chaotic behavior and complex physical oscillations. The core of this project is a generalized implementation of numerical integration methods applied to classical physics problems.

## Featured Project: Lorenz System Solver (OOP)

The highlight of this repository is a robust solver for the **Lorenz Equations**, implemented using **Object-Oriented Programming (OOP)** in C++. 

### Mathematical Background
The Lorenz system is a set of three ordinary differential equations (ODEs) that describe the simplified dynamics of atmospheric convection:

$$\frac{dx}{dt} = \sigma(y - x)$$
$$\frac{dy}{dt} = x(\rho - z) - y$$
$$\frac{dz}{dt} = xy - \beta z$$

For specific parameters ($\sigma=10, \rho=28, \beta=8/3$), the system exhibits chaotic behavior, known as the **Lorenz Attractor**.

### Technical Implementation
* **Generalized RK4 Class:** I developed a reusable C++ class that implements the **4th Order Runge-Kutta (RK4)** method. This class is designed to handle generic systems of differential equations, demonstrating high-level abstraction and code reusability.
* **Numerical Stability:** The solver is optimized for long-term integration to capture the butterfly effect characteristic of chaotic systems.
* **Data Analysis:** Output data is processed to visualize the phase space trajectories.

### Visualizing Chaos
Below is a simulation result showing the evolution of the attractor in phase space:

![Lorenz System Trajectory](Lorez-Equations/SolutionCase4.png)

## Other Included Simulations
Beyond the Lorenz system, this repository includes:
* **Non-linear Oscillators:** Numerical solutions for damped and driven oscillations.
* **Chaotic Systems:** Exploration of sensitivity to initial conditions.
* **Rigid Body Dynamics:** Modeling of complex rotations and movements.

## Technologies Used
* **C++:** High-performance numerical computing and OOP.
* **Python (NumPy/Matplotlib):** Used for data cleaning, statistical validation, and high-quality plotting.

## Structure
```text
.
├── Lorez-Equations/
│   ├── SolutionCase4.png     # Visual result of the chaotic attractor
│   └── src/                  # C++ source code (RK4 Class and Main)
├── Non-Linear-Oscillators/   # Scripts for various oscillation models
└── docs/                     # Academic reports and derivations
