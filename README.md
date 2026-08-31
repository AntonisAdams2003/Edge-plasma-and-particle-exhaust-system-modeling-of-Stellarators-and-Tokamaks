# **Edge Plasma and Particle Exhaust System Modeling of Stellarators and Tokamaks**

## **Overview**
This repository contains the source code supporting my Master's thesis: *"Simplified Modeling of the Edge Plasma Physics and Neutral Gas Exhaust in the Thermonuclear Reactor Wendelstein 7-X"*. 

Specifically, it models edge plasma transport using an extended 0D Two-Point model, and neutral gas transport in the divertor module via particle balance and conductance models.

## **Prerequisites**
The scripts are written in python and require the installation of the following libraries:
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `math`
- `sys`
- `PyQt5`
- `os` (only for Bayesian_MAP script)
- `emcee` (only for Bayesian_MAP script)

---

## **Repository Structure**

#### `Two_Point_Model_Solver.py`
Contains solvers for three simplified 0D edge plasma transport models.
Import them via: `from Two_Point_Model_Solver import [SolverName]`

*   **`TTPM`**: Analytical solutions of Stangeby's Tokamak Two-Point model [1]. 
    *   **Inputs**: `n_upst` [m⁻³], `Lc` [m], `x` [m²/s], `q_SOL` [MW/m²]
    *   **Outputs**: `n_targ` [m⁻³], `T_upst` [eV], `T_targ` [eV]
*   **`TTPM_numerical`**: Numerical solver for Stangeby's model.
    *   **Inputs / Outputs**: Same as `TTPM`.
*   **`Extended_STPM`**: Numerical solver for Feng's Extended Stellarator Two-Point model [2] incorporating three loss parameters.
    *   **Inputs**: `n_upst`, `loss_alpha` [eV^(1/2)], `loss_fcool` [-], `loss_fconv` [-], `Lc`, `x`, `theta` [-], `q_SOL`
    *   **Outputs**: Same as `TTPM`.




#### `Divertor_Subdivertor_Solver.py`
**[Description needed]**
*   **`Divertor_Subdivertor_Analytical`**
    *   **Inputs:** | **Outputs:**
*   **`Divertor_Subdivertor_Conductance`**
    *   **Inputs:** | **Outputs:**
*   **`Divertor_Subdivertor_Dirk`**
    *   **Inputs:** | **Outputs:**





#### `Post_Processing.py`
Integrates the two core modules above. Provides a parametric post-processing script that solves the problem from the upstream to the subdivertor.
*   **Inputs:** 
*   **Outputs:** 




#### `Subdivertor_Models_Comparison.py`
Compares closed-form expressions estimating subdivertor pressure (NOT divertor pressure) against experimental measurements [3]. Moreover, it provides bar plots for two of the models, to compare them with high-fidelity code results from [4]. 
Import via:** `from Subdivertor_Models_Comparison import [ModelName]`

*   **`Varoutis_subdivertor`**: Regression expressions from [5].
    *   **Inputs**: `N_in_AEH`, `N_in_AEP` | **Outputs**: `p_sub_AEH`, `p_sub_AEP`
*   **`Dirk_subdivertor`**: **[Description needed]**
    *   **Inputs**: `N_in_AEH`, `N_in_AEP` | **Outputs**: `p_sub_AEH`, `p_sub_AEP`, `h_AEH`, `h_AEP`, `w_AEH`, `w_AEP`
*   **`Litovoli_Haak_subdivertor`**: 2-reservoir conductance model from **[]**.
    *   **Inputs**: `N_in_AEH`, `N_in_AEP` | **Outputs**: `p_sub_X`, `N_out_in_X`, `N_leak_in_X`, `N_pump_in_X` **(where X = AEH/AEP)**
*   **`ParticleBalance_subdivertor`**: Particle balance model **introduced in this thesis**.
    *   **Inputs / Outputs: Same as `Litovoli_Haak_subdivertor`.**



#### `Parallel_vs_Perpendicular_Conduction.py`
Compares the parallel and perpendicular heat conduction terms of the heat equation. Outputs a figure showing curves where `q_parallel_conduction = q_perpendicular_conduction` for Tokamaks/Stellarators and electrons/ions.
*   **Inputs**: None
*   **Outputs**: Figure of equilibrium curves.




#### `Bayesian_MAP_and_MCMC_Solver.py`
Solves a toy-problem of parameter estimation (MAP), using a simplified analytical Tokamak model and pseudo-data. Additionally, uses Markov Chain Monte Carlo (MCMC) to estimate the posterior distribution under strong vs. weak priors.

## Filler codes

#### `Fitted_Closed_Form_nt.py`
Fits numerical results of the Extended Stellarator Two-Point model to a logarithmic function to approximate target particle density: `log(nt) = g(nu, loss_alpha, loss_fcool, loss_fconv)`.
*   **Inputs**: None
*   **Outputs**: Regression weights, comparison figure.

**CAUTION:** This expression is strictly for providing initial guesses to a numerical solver. Do not use it to estimate target density directly (exponentiating the logarithmic fit produces significant errors).


---


## **Usage Guidelines**
**Follow these rules to reproduce the analysis:**
- Do not execute `Two_Point_Model_Solver.py` or `Divertor_Subdivertor_Solver.py` directly; they are module libraries meant to be imported.
- Use `Post_Processing.py` to evaluate the upstream-to-subdivertor flow using only upstream conditions.
- Use the isolated subdivertor models in `Subdivertor_Models_Comparison.py` ONLY when the particle flow (`N_in`) through the pumping gaps is given.


## **References**
1. **[Stangeby, P.C. - The Plasma Boundary of Magnetic Fusion Devices]** (https://www.routledge.com/The-Plasma-Boundary-of-Magnetic-Fusion-Devices/Stangeby/p/book/9780750305594)
2. **[N. Maaziz et al 2026 Nucl. Fusion]** (https://iopscience.iop.org/article/10.1088/1741-4326/ae855c)
3. **[V Haak et al 2023 Plasma Phys. Control. Fusion]** (https://iopscience.iop.org/article/10.1088/1361-6587/acc8fb/meta)
4. **[S. Varoutis et al 2024 Nucl. Fusion]** (https://iopscience.iop.org/article/10.1088/1741-4326/ad49b5)
5. **[S. Varoutis et al 2025 Nucl. Fusion]** (https://iopscience.iop.org/article/10.1088/1741-4326/addbf1/meta)

## **Author**
**Antonis Adamopoulos**
**📧 [antonis.adamopoulos2003@gmail.com](mailto:antonis.adamopoulos2003@gmail.com) | [aadamopoulos@ethz.ch](mailto:aadamopoulos@ethz.ch)**
