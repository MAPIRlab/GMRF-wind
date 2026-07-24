# GMRF-wind

A Gaussian Markov Random Field (GMRF) is a specific type of Markov Random Field (MRF) where the set of random variables follows a multivariate Gaussian (Normal) distribution. It is a statistical model used to describe the dependencies among a collection of random variables, often representing spatial or temporal data. This repository applies this framwork to the estimation of 2D wind maps (W) from a set of wind vector observations (Z) and prior knowledge encapsulating physical constraints, providing a ROS2 wrapped implementation of the algorithm presented in this paper: https://www.sciencedirect.com/science/article/pii/S0360132326007614?via%3Dihub

The core of the Gaussian Markov Random Field (GMRF) framework is the definition of energy terms (or factors) that encode the relationships between adjacent cells and observations, ultimately leading to the Maximum a Posteriori (MAP) estimation through the minimization of the total energy function, E(W,Z). This allows the GMRF-W approach to function as a real-time 2D approximation of more complex Computational Fluid Dynamics (CFD) techniques, ideal for robotics that adquire new observations as they inspect the environment.

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
Although GMRF-W is a self-contained pkg, the implementation considers anemometer sensor readings which depends on an external pkg defining some "olfaction" related msgs. This pkg is available in a different repository named olfaction_msgs (https://github.com/MAPIRlab/olfaction_msgs).
