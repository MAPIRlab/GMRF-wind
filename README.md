# GMRF-wind

A Gaussian Markov Random Field (GMRF) is a specific type of Markov Random Field (MRF) where the set of random variables follows a multivariate Gaussian (Normal) distribution. It is a statistical model used to describe the dependencies among a collection of random variables, often representing spatial or temporal data. This repository applies this framwork to the estimation of 2D wind maps (W) from a set of wind vector observations (Z) and prior knowledge encapsulating physical constraints, providing a ROS2 wrapped implementation of the algorithm presented in this paper: https://ieeexplore.ieee.org/document/7968883

The core of the Gaussian Markov Random Field (GMRF) framework is the definition of energy terms (or factors) that encode the relationships between adjacent cells and observations, ultimately leading to the Maximum a Posteriori (MAP) estimation through the minimization of the total energy function, E(W,Z). For GMRF-W (used for Online Estimation of 2D Wind Maps), four primary ”energies”
are considered. All four terms have direct or conceptual physical relevance. The overall energy function E(W,Z) is the sum of these four factors:
     E(W,Z) = Ez(W,Z) + Em(W) + Eo(W) + Er(W)

By combining these four energy terms, the GMRF-W framework is able to estimate a 2D wind map that is consistent with observations, respects the presence of obstacles, and adheres to the law of mass conservation for incompressible flow. This allows the GMRF-W approach to function as a real-time approximation of more complex Computational Fluid Dynamics (CFD) techniques.

If it is relevant to your research, you can cite the paper with the following BibTex: 

@INPROCEEDINGS{jmonroy_isoen_2017,
     author = {Monroy, Javier and Jaimez, Mariano and Gonzalez-Jimenez, Javier},
      title = {Online Estimation of 2D Wind Maps for Olfactory Robots},
  booktitle = {International Symposium on Olfaction and Electronic Nose (ISOEN)},
       year = {2017},
   location = {Montreal (Canada)},
        doi = {10.1109/ISOEN.2017.7968883},
      pages = {1--3}
}
