# GMRF-wind

A lightweight, training-free, physics-informed **2D indoor airflow field estimation** using Gaussian Markov Random Fields (GMRF).

[![Paper](https://img.shields.io/badge/Paper-Building%20and%20Environment-blue)](https://doi.org/10.1016/j.buildenv.2026.114957)
[![ROS 2](https://img.shields.io/badge/ROS%202-Humble%20%7C%20Jazzy-brightgreen)](https://docs.ros.org/)

---

## Overview

`GMRF-wind` provides a real-time, training-free spatial estimation framework to reconstruct continuous 2D wind velocity maps ($\mathbf{W}$) from a set of sparse, noisy wind vector observations ($\mathbf{Z}$) and occupancy grid geometries. 

Instead of relying on heavy Computational Fluid Dynamics (CFD) solvers or data-intensive machine learning models, this package formulates wind field estimation as a **Maximum A Posteriori (MAP)** inference problem on a Gaussian Markov Random Field graph. By minimizing a total energy function $E(\mathbf{W}, \mathbf{Z})$, the algorithm efficiently solves a linear sparse system to deliver macro-scale wind maps in milliseconds.

This repository provides the official ROS 2 wrapper and C++ core implementation of the methodology presented in:

> **A Physics-Informed Gaussian Markov Random Field Framework for Indoor Airflow Field Estimation**  
> *Javier Monroy, Pepe Ojeda, and Javier Gonzalez-Jimenez*  
> *Building and Environment*, Vol. 288, 2026.  
> 🔗 [Read Full Paper (Open Access)](https://doi.org/10.1016/j.buildenv.2026.114957)

## The Code

The repository is organized as a ROS 2 workspace-style package collection: one small message package for custom interfaces (**gmrf_msgs**) and one main mapping package that contains the C++ solver, launch/config files, scripts, and data assets for estimating 2D indoor airflow fields. The code is essentially divided into the core library (in the subpackage **gmrf_wind_core**) and a ROS-wrapper node (**gmrf_wind_mapping**).

Whitin it, there are two ROS executables, the main “nodes” you should care about:

1) **gmrf_wind_mapping_node** implemented in gmrf_node.cpp:  this is the main runtime node for online wind mapping.
- accepts an environment occupancy map (from MapServer)
- receives sparse wind observations from anemometers (via topic subscription)
- converts sensor measurements into map-frame coordinates using TF
- inserts those observations into a CGMRF_map (the core class, ROS-independent)
- solves the MAP estimation repeatedly
- publishes the result as RViz markers (for visualization and debug)
- exposes two services: **WindEstimation**: returns U/V velocity components, can query either the whole grid or specific points. **AddWindObservation**: lets clients add additional sensor observations directly, useful for external sources or debugging.

This is the node used in the standard launch file **gmrf_wind_launch.py**.

2) **gmrf_validation** implemented in gmrf_validation.cpp: this is an evaluation/benchmark node designed to operate offline for numerical testing.
- it loads an occupancy map and CFD ground-truth dataset (from file)
- runs the estimator after "simulating" measurements
- compares prediction quality using metrics like MAE/RMSE or similar optimization criteria
- it is designed for validation experiments rather than real-time deployment

This is the node launched by **gmrf_validation_launch.py**.



---

## Key Features

- ⚡ **Real-Time & Computational Efficiency:** Solves spatial estimates in milliseconds on CPU, suitable for low-power onboard robotic processors.
- 🚫 **Training-Free & Zero Setup:** Requires no offline datasets, neural network training, or parameter tuning per environment.
- 🛡️ **Physics-Grounded Constraints:** Embeds fundamental fluid transport mechanics directly into the graph precision matrix:
  - **Incompressibility / Continuity Constraint:** $\nabla \cdot \mathbf{w} = 0$
  - **Advection-Diffusion Momentum Proxy:** Smooths flow direction along streamlines while preserving spatial gradients.
  - **No-Penetration Boundary Conditions:** Prevents unphysical airflow through walls and solid obstacles ($\mathbf{w} \cdot \mathbf{n}_{\text{obs}} = 0$).
- 🤖 **Robotics Ready:** Native ROS 2 integration accepting live OccupancyGrid maps and point wind observations (e.g., from anemometers mounted on mobile robots).

---

## How It Works

The framework embeds simplified physical constraints derived from the Navier–Stokes equations—including mass conservation, advection, and viscous diffusion—into the estimation process, ensuring physically consistent flow reconstruction at a fraction of the computational cost of CFD. It defines specialized energy factors encoding geometric spatial priors, physical conservation laws, and sensor observation models:

$$E(\mathbf{W}, \mathbf{Z}) = E_{\text{z}}(\mathbf{W}, \mathbf{Z}) + E_{\text{o}}(\mathbf{W}) + E_{\text{physics}}(\mathbf{W})$$

1. **Observation Factor ($E_{\text{z}}$):** Pulls the local velocity vector towards incoming anemometer measurements.
2. **Spatial Prior ($E_{\text{o}}$):** enforces a boundary condition without penetration at obstacles.
3. **Physics Factors ($E_{\text{physics}}$):** Constrains neighbor-to-neighbor transitions according to indoor mass continuity, advection and diffusion.

The resulting sparse Gaussian precision matrix ($\mathbf{Q}$) allows solving the linear system $\mathbf{Q}\mathbf{W} = \mathbf{b}$ continuously, as new measurements arrive.



---

## Intended Applications

`GMRF-wind` is designed for operational robotics and environmental monitoring applications where fast, macro-scale flow awareness is required:

- **Robot-Assisted Gas Source Localization (GSL):** Guiding mobile inspection robots toward hazardous gas leaks by tracking active wind corridors.
- **Indoor Air Quality (IAQ) & HVAC Optimization:** Rapid mapping of ventilation patterns, stagnation zones, and pollutant dispersion routes.

---

## Quick Start

### Prerequisites

- ROS 2 (Humble / Iron / Jazzy)
- Eigen3
- OpenCV / PCL (for map handling)

Although GMRF-W is a self-contained pkg, the implementation considers anemometer sensor readings which depends on an external pkg defining some "olfaction" related msgs. This pkg is available in a different repository named olfaction_msgs (https://github.com/MAPIRlab/olfaction_msgs).

### Citation

If it is relevant to your research, you can cite the paper with the following BibTex: 
```
@ARTICLE{monroy_bae_2026,
    author = {Monroy, Javier and Ojeda, Pepe and Gonzalez-Jimenez, Javier},
     title = {A Physics-Informed Gaussian Markov Random Field Framework for Indoor Airflow Field Estimation},
   journal = {Building and Environment},
      year = {2026},
       url = {https://doi.org/10.1016/j.buildenv.2026.114957},
       doi = {10.1016/j.buildenv.2026.114957}
}

```


