# Edge-plasma-and-particle-exhaust-system-modeling-of-Stellarators-and-Tokamaks

## **OVERVIEW**
This repository contains the codes that supported my Master's thesis: "Simplified Modeling of the Edge Plasma Physics and Neutral Gas Exhaust in the Thermonuclear Reactor Wendelstein 7-X". These scripts contain solvers and models for the analysis of the edge plasma physics and neutral gas exhaust system of Stellarator and Tokamak reactors.
Specifically, the edge plasma transport is modeled through the 0D Two-Point model (where extended parameters are considered) while the neutral gas transport in the divertor module is modeled with particle balance or conductance models.

## **REPOSITORY STRUCTURE**
**The project is organized as follows:**
- **'**Two_Point_Model_Solver.py**':** Contains solvers for 3 different simplified 0D edge plasma transport models.
- **'**Divertor_Subdivertor_Solver.py**':** Contains analytical, conductance, and Dirk solvers for the divertor and subdivertor.
- **'**Post_Processing.py**':** This script combines the aforementioned two (Two_Point_Model_Solver.py & Divertor_Subdivertor_Solver.py) and provides a parametric post processing script that solves the problem from the upstream to the subdivertor.
- **'**Subdivertor_Models_Comparison.py**':** This script uses some analytical expressions that estimate the subdivertor pressure (NOT the divertor one though) and compares them with experimental measurements provided in [Haak's paper].
- **'**Parallel_vs_Perpendicular_Conduction.py**':** This script compares the heat conduction terms of parallel and perpendicular (to the magnetic field lines) directions.
- **'**Fitted_Closed_Form_nt.py**':** This script tries to find a closed form expression for the target particle density given by the Extended Stellarator Two-Point model.
- **'**Bayesian_MAP_and_MCMC_Solver.py**':** This script solves a toy-problem of free-parameter estimation using a simplified analytical Tokamak model and pseudo-data.

## 1) Two_Point_Model_Solver.py

This script containts solvers for 3 different simplfied 0D edge plasma transport model. Select the desired solver and import it as `from Two_Point_Model_Solver import (...Solvers name...)`

### TTPM
Uses the analytical solutions of the Stangeby's Tokamak Two-Point model, derived in [1]. 

* **Inputs:** upstream plasma density (n_upst [m^-3]), connection length (Lc [m]), perpendicular  diffusion term (x [m^2/s), and parallel heat flux on the SOL (q_SOL [MW/m^2)
* **Outputs:** target particle density (n_targ [m^-3]), upstream plasma temperature (T_upst [eV]), target plasma temperature (T_targ [eV])

### TTPM_numerical
Solves the Tokamak Two-Point model of Stangeby numerically

* **Inputs:** same as above
* **Outputs:** same as above

### Extended_STPM
Solves the Extended Stellarator Two-Point model, inroduced by Feng, numerically. It containes 3 loss parameters that account for the neglected physics.

* **Inputs:** upstream plasma density (n_upst [m^-3]), Momentum loss factor (loss_alpha [eV^1/2]), Cooling losses (loss_fcool [-]), Convection fraction (loss_fconv [-]) connection length (Lc [m]), perpendicular diffusion term (x [m^2/s), Theta pitch (theta [-]), and parallel heat flux on the SOL (q_SOL [MW/m^2)
* **Outputs:** same as above

## 2) Divertor_Subdivertor_Solver.py

This script 

### Divertor_Subdivertor_Analytical
* **Inputs:**
* **Outputs:**

### Divertor_Subdivertor_Conductance
* **Inputs:**
* **Outputs:**

### Divertor_Subdivertor_Dirk
* **Inputs:**
* **Outputs:**

## 3) Post_Processing.py

This scripts combines the aforementioned two (Two_Point_Model_Solver.py & Divertor_Subdivertor_Solver.py) and provides a paramtric post processing script that solves the problem from the upsream to the subdivertor.

* **Inputs:** 
* **Outputs:**

## 4) Subdivertor_Models_Comparison.py

This script uses some analytical expressions that estimate the subdivertor pressure (NOT the divertor one though) and compares them with experimental meassurements provided in [Haak's paper]. Adtionally, a comparison of the models 'Litovoli-Haak_subdivertor' and 'ParticleBalance_subdivertor' is made with [Varoutis 2024 Fig11]. To use these models seperatelly, one can import them in another script as follows: `from Subdivertor_Models_Comparison import (...Solvers name...)`.
The individual models are:

### Varoutis_subdivertor
Uses the regression expressions introduced in [Varoutis 2025]

* **Inputs:** N_in_AEH, N_in_AEP
* **Outputs:** p_sub_AEH,p_sub_AEP

### Dirk_subdivertor
Uses 

* **Inputs:** N_in_AEH, N_in_AEP
* **Outputs:** p_sub_AEH,p_sub_AEP , h_AEH,h_AEP,w_AEH,w_AEP

### Litovoli_Haak_subdivertor
Uses the conductance model, for 2 reservoirs, introduced in [Haak-litovoli]

* **Inputs:** N_in_AEH, N_in_AEP
* **Outputs:** p_sub_AEH, p_sub_AEP, N_out_in_AEH, N_out_in_AEP, N_leak_in_AEH, N_leak_in_AEP, N_pump_in_AEH, N_pump_in_AEP

### ParticleBalance_subdivertor
Uses the particle balance model introduced in my thesis

* **Inputs:** N_in_AEH, N_in_AEP
* **Outputs:** p_sub_AEH, p_sub_AEP, N_out_in_AEH, N_out_in_AEP, N_leak_in_AEH, N_leak_in_AEP, N_pump_in_AEH, N_pump_in_AEP

---

## 5) Parallel_vs_Perpenicular_Conduction.py **[Note: Consider moving to an Archive folder]**

This script compares the heat conduction terms of parallel and perpendicular (to the magnetic fied lines) directions. The figure produced shows the curves where q_parallel_conduction = q_perpaendicular_conduction, for both Tokamaks and Stellarators and for both electrons and ions.                                                     

* **Inputs:** None
* **Outputs:** figure of curves

## 6) Fitted_Closed_Form_nt.py

This script tries to find a closed form expression for the target particle density given by the Extended Stellarator Two-Point model. Specifically, it fits the numerical results of the Two-Point model in a logarithmic function, in order to find the log(n_target) with respect to n_upstream, loss_alpha, loss_fcool, loss_fconv : `log(nt) = g(nu,loss_alpha,loss_fcool, loss_fconv)`

* **Inputs:** None
* **Outputs:** the weights of the regresin expression + figure comparing the expression to the data

> **CAUTION:** This expression may only be useful for initial guesses in a numerical solver and should never be used to estimate the target density. The fit is onle good at the logarithmic scale. When we exponate to get the actual n_t expression the erro grows significantly.

## 7) Bayesian_MAP_and_MCMC_Solver.py

This script solves a toy-problem of free-parameter estimation using a simplified analyical Tokamak model and pseudo-data. Specifically, it solves the Maximum A Poseriori problem and then applies the Marcov Chain Monte Carlo algorithms to estimate the posterior ditribution of the free parameter. Two extreme cases are examined: strong and week prior. 

* **Inputs:** None
* **Outputs:** Figures of the posterior and prior distributions# Edge-plasma-and-particle-exhaust-system-modeling-of-Stellarators-and-Tokamaks

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

## **REFERENCES**
- [1]: Stangeby

## **AUTHOR**
**[Your Name], Email: [Your Email]**

