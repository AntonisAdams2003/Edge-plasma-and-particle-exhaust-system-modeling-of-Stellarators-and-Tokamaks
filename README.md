# Edge-plasma-and-particle-exhaust-system-modeling-of-Stellarators-and-Tokamaks

## **OVERVIEW**
This repository contains the codes that supported my Master's thesis: "Simplified Modeling of the Edge Plasma Physics and Neutral Gas Exhaust in the Thermonuclear Reactor Wendelstein 7-X". 

Specifically, the edge plasma transport is modeled through the 0D Two-Point model (where extended parameters are considered), while the neutral gas transport in the divertor module is modeled through particle balance or conductance models.

## **REPOSITORY STRUCTURE**
**The project is organized as follows:**


## 1) Two_Point_Model_Solver.py*

This script contains solvers for three different simplified 0D edge plasma transport models. Select the desired solver and import it as `from Two_Point_Model_Solver import (...Solvers name...)`

### TTPM
Uses the analytical solutions of the Stangeby's Tokamak Two-Point model, derived in [1]. 

* **Inputs:** upstream plasma density (n_upst [m^-3]), connection length (Lc [m]), perpendicular  diffusion term (x [m^2/s]), and parallel heat flux on the SOL (q_SOL [MW/m^2])
* **Outputs:** target particle density (n_targ [m^-3]), upstream plasma temperature (T_upst [eV]), target plasma temperature (T_targ [eV])

### TTPM_numerical
Solves the Stangeby's Tokamak Two-Point model numerically

* **Inputs:** same as above
* **Outputs:** same as above

### Extended_STPM
Solves the Extended Stellarator Two-Point model [2], introduced by Feng, numerically. It contains three loss parameters that account for the neglected physics.

* **Inputs:** upstream plasma density (n_upst [m^-3]), momentum loss factor (loss_alpha [eV^1/2]), cooling losses (loss_fcool [-]), convection fraction (loss_fconv [-]) connection length (Lc [m]), perpendicular diffusion term (x [m^2/s]), theta pitch (theta [-]), and parallel heat flux on the SOL (q_SOL [MW/m^2])
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

This script compares the heat conduction terms of parallel and perpendicular (to the magnetic field lines) directions. The figure produced shows the curves where q_parallel_conduction = q_perpandicular_conduction, for both Tokamaks and Stellarators and for both electrons and ions.                                                     

* **Inputs:** None
* **Outputs:** figure of curves


## 6) Fitted_Closed_Form_nt.py

This script tries to find a closed form expression for the target particle density given by the Extended Stellarator Two-Point model. Specifically, it fits the numerical results of the Two-Point model in a logarithmic function, in order to find the log(n_target) with respect to n_upstream, loss_alpha, loss_fcool, loss_fconv : `log(nt) = g(nu,loss_alpha,loss_fcool, loss_fconv)`

* **Inputs:** None
* **Outputs:** the weights of the regression expression + figure comparing the expression to the data

> **CAUTION:** This expression may only be useful for initial guesses in a numerical solver and should never be used to estimate the target density. The fit is only good at the logarithmic scale. When we exponantiate to get the actual n_t expression the error grows significantly.


## 7) Bayesian_MAP_and_MCMC_Solver.py

This script solves a toy-problem of free-parameter estimation using a simplified analytical Tokamak model and pseudo-data. Specifically, it solves the Maximum A Poseriori problem and then applies the Markov Chain Monte Carlo algorithms to estimate the posterior distribution of the free parameter. Two extreme cases are examined: strong and week prior. 


## **LIBRARIES REQUIRED**
**To run the scripts, you will need to install the following Python libraries:**
- numpy
- pandas

## **HOW TO RUN THE CODE**
**To reproduce the analysis, please follow these steps:**
- Never run seperately the scripts 1) and 2). They may only be used by calling their functions
- The Post_Processing.py script is used to solve the upstream to subdivertor flow, given ONLY the upstream conditions
- Models solely for the subdivertor are contained in the Subdivertor_Models_Comparison.py. They may be used only if the particle flow (N_in) that passes the pumping gaps is known

## **NOTES**
- **[Provide any additional notes here]**

## **REFERENCES**
- [1]: https://www.routledge.com/The-Plasma-Boundary-of-Magnetic-Fusion-Devices/Stangeby/p/book/9780750305594
- [2]: https://iopscience.iop.org/article/10.1088/1741-4326/ae855c

## **AUTHOR**
**[Your Name], Email: [Your Email]**

