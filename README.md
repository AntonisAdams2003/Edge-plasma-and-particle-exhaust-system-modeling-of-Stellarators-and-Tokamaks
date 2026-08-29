# Edge-plasma-and-particle-exhaust-system-modeling-of-Stellarators-and-Tokamaks

## **OVERVIEW**
This repository contains the codes that supported my Master's thesis: "Simplified Modeling of the Edge Plasma Physics and Neutral Gas Exhaust in the Thermonuclear Reactor Wendelstein 7-X"**.** These codes contain solvers and models for edge plasma physics and neutral gas exhaust.

## **REPOSITORY STRUCTURE**
**The project is organized as follows:**
- **'**Two_Point_Model_Solver.py**':** This script contains solvers for 3 different simplified 0D edge plasma transport models.
- **'**Divertor_Subdivertor_Solver.py**':** This script **contains analytical, conductance, and Dirk solvers for the divertor and subdivertor.**
- **'**Post_Processing.py**':** This script combines the aforementioned two (Two_Point_Model_Solver.py & Divertor_Subdivertor_Solver.py) and provides a parametric post processing script that solves the problem from the upstream to the subdivertor.
- **'**Subdivertor_Models_Comparison.py**':** This script uses some analytical expressions that estimate the subdivertor pressure (NOT the divertor one though) and compares them with experimental measurements provided in [Haak's paper].
- **'**Parallel_vs_Perpendicular_Conduction.py**':** This script compares the heat conduction terms of parallel and perpendicular (to the magnetic field lines) directions.
- **'**Fitted_Closed_Form_nt.py**':** This script tries to find a closed form expression for the target particle density given by the Extended Stellarator Two-Point model.
- **'**Bayesian_MAP_and_MCMC_Solver.py**':** This script solves a toy-problem of free-parameter estimation using a simplified analytical Tokamak model and pseudo-data.

## **KEY DESIGN DECISIONS**
**Two_Point_Model_Solver.py:**
- **Import via:** `from Two_Point_Model_Solver import (...Solvers name...)`
- **TTPM:** Solves the Tokamak Two-Point model, introduced by Stangeby, using Stangeby's approximation to derive a closed-form solution.
- **TTPM_numerical:** Solves the Tokamak Two-Point model of Stangeby numerically.
- **Extended_STPM:** Solves the Extended Stellarator Two-Point model, introduced by Feng, numerically. It contains 3 loss parameters that account for the neglected physics.

**Subdivertor_Models_Comparison.py:**
- **To use these models separately, import them as:** `from Subdivertor_Models_Comparison import (...Solvers name...)`
- **Model comparisons:** A comparison of the models 'Litovoli-Haak_subdivertor' and 'ParticleBalance_subdivertor' is made with [Varoutis 2024 Fig11].
- **Varoutis_subdivertor:** Uses the regression expressions introduced in [Varoutis 2025].
- **Dirk_subdivertor:** Uses **[Incomplete in original draft]**.
- **Litovoli_Haak_subdivertor:** Uses the conductance model, for 2 reservoirs, introduced in [Haak-litovoli].
- **ParticleBalance_subdivertor:** Uses the particle balance model introduced in my thesis.

**Parallel_vs_Perpendicular_Conduction.py:**
- **The figure produced shows the curves where q_parallel_conduction = q_perpendicular_conduction, for both Tokamaks and Stellarators and for both electrons and ions.**

**Fitted_Closed_Form_nt.py:**
- **Specifically, it fits the numerical results of the Two-Point model in a logarithmic function, in order to find the log(n_target) with respect to n_upstream, loss_alpha, loss_fcool, loss_fconv : `log(nt) = g(nu,loss_alpha,loss_fcool, loss_fconv)`**
- **CAUTION:** This expression may only be useful for initial guesses in a numerical solver and should never be used to estimate the target density. The fit is only good at the logarithmic scale. When we exponentiate to get the actual n_t expression the error grows significantly.

**Bayesian_MAP_and_MCMC_Solver.py:**
- **Specifically, it solves the Maximum A Posteriori problem and then applies the Markov Chain Monte Carlo algorithms to estimate the posterior distribution of the free parameter. Two extreme cases are examined: strong and weak prior.**

## **LIBRARIES REQUIRED**
**To run the scripts, you will need to install the following Python libraries:**
- **[List dependencies here]**

## **HOW TO RUN THE CODE**
**To reproduce the analysis, please follow these steps:**
- **[Provide execution steps here]**

## **NOTES**
- **[Provide any additional notes here]**

## **AUTHOR**
**[Your Name], Email: [Your Email]**
