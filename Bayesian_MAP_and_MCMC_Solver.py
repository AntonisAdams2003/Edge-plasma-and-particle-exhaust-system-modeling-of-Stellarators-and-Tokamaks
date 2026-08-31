import numpy as np
from numpy import pi, sqrt
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import emcee
import os
import scipy.stats as stats

#######################################################################################################
## In this code the Maximum A Posteriori optimization problem is solved for a single free parameter, ##
## for the simple case of the analytical Stangeby's solution. The desired free-parameter is the      ##
## momentum loss fraction (f_mom)                                                                    ##
####################################################################################################### 
  

## PREDICTION MODEL
# Define the known constants and geometrical dimentions
m_H = 1.67262e-27  # [kg] mass hydrogen ion
e = 1.60218e-19    # electron charge
# Define constants by assuming typical values
gamma = 8   # Sheath transmission coefficient
koe = 2000  # [W/m eV^7/2] Spitzer-Härm heat conductivity for electrons
P_SOL = 50.0   # [MW] Heating power in the SOL
q_par = P_SOL*1e6  # The parallel heat flux [W/m^2]
Lc = 1         # [m] Estimated average connection length (=400 for W7-X, =1 for tokamaks)

# Prodiction model
c_const = gamma**2*e**3/(4*m_H) * (7/2*q_par*Lc/koe)**(6/7) / q_par**2
# Define the nondimensionalization constant
nondimen_const = 1e19

# Define the Two point model
def model(nu,fmom):
    nt = (c_const * nondimen_const**2) * fmom**3 * nu**3
    return nt


def run_map_mcmc_case(N, fmom_true, Z_noise_percentage, p_mean_val, p_std_val, case_name):
    print(f"- {case_name}")

    ## TRAINING DATA AND INPUTS ##
    x_data = np.linspace(1,4e1,N)
    y_true = model(x_data, fmom_true)

    epsilon_noise_std = Z_noise_percentage * y_true
    y_data = np.clip(y_true + np.random.normal(0, epsilon_noise_std, N), 1e-10, np.inf)
    error_std = Z_noise_percentage * y_true

    param_mean = np.array([p_mean_val])
    param_std = np.array([p_std_val])
    param_bounds = [(1e-5, 1.0)] 
    initial_guess = param_mean/2

    ## OBJECTIVE FUNCTION & MAP PROBLEM
    def objective_func(parameters, x, y, err_std, p_mean, p_std):
        loss_param = parameters[0]
        y_predict = model(x, loss_param)
        sum_square_error = np.sum((y - y_predict)**2 / (2 * error_std**2))
        sum_square_param = (loss_param - p_mean[0])**2 / (2 * p_std[0]**2)
        return sum_square_error + sum_square_param

    ## MINIMIZATION SOLVER ##
    opt_parameters = minimize(
        objective_func,
        initial_guess,
        args=(x_data, y_data, error_std, param_mean, param_std),
        method='L-BFGS-B',
        bounds=param_bounds,
        options={'disp': True, 'maxiter': 500}
    )

    print(f"Most probable parameter value: {opt_parameters.x}")
    
    ## MARKOV CHAIN MONTE CARLO ##
    def log_posterior(parameters, x, y, err_std, p_mean, p_std):
        lower_bound, upper_bound = param_bounds[0]
        if not (lower_bound <= parameters[0] <= upper_bound):
            return -np.inf
        return -objective_func(parameters, x, y, err_std, p_mean, p_std)

    ndim = 1
    nwalkers = 50
    map_estimate = opt_parameters.x
    initial_pos = np.clip(map_estimate + 1e-4 * np.random.randn(nwalkers, ndim), 1e-5, 1.0 - 1e-5)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(x_data, y_data, error_std, param_mean, param_std))
    N_samples = 10000 
    sampler.run_mcmc(initial_pos, N_samples)
    samples = sampler.get_chain(discard=round(N_samples/10), flat=True)

    ## PLOTTING ##
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    nu_data = x_data * nondimen_const
    nt_data = y_data * nondimen_const
    fmom_opt = opt_parameters.x
    nt_pred = model(x_data, fmom_opt) * nondimen_const

    ax1.loglog(nu_data, nt_data, marker='D', markersize=7, color='#000000', label="Data", linestyle='', markerfacecolor='#000000')
    ax1.loglog(nu_data, nt_pred, label=r'Prediction using $f_{mom,opt}$', color='#5e2c84', linewidth=2.5)
    ax1.set_xlabel(r'$n_u$ [m$^{-3}$]', fontsize=18)
    ax1.set_ylabel(r'$n_t$ [m$^{-3}$]', fontsize=18)
    ax1.set_title(f'Density build-up', fontsize=18)
    ax1.legend(loc='upper left', frameon=True, fontsize=16)
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5, color='#d3d3d3')
    ax1.tick_params(which='both', labelsize=20)
    plt.tight_layout() 

    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.hist(samples, bins=50, density=True, color="#4890c4", edgecolor='#0b5394', alpha=0.7, label='MCMC Posterior')
    ax2.axvline(fmom_opt[0], color='#000000', linestyle='--', linewidth=2.5, label=f'MAP value: {fmom_opt[0]:.3f}')
    ax2.axvline(fmom_true, color="#eb2626", linestyle='-', linewidth=2.5, zorder=5, label=f'True value: {fmom_true:.3f}')

    fmom_space = np.linspace(0.0, 1.0, 500)
    ax2.plot(fmom_space, stats.norm.pdf(fmom_space, param_mean[0], param_std[0]), color="#050E1B", linestyle='-', linewidth=2, label='Prior Distribution')

    ax2.set_xlabel(r'$f_{mom}$', fontsize=18)
    ax2.set_ylabel('Probability Density', fontsize=18)
    ax2.set_xlim([0.3, 0.9]) 
    ax2.legend(loc='upper right', frameon=True, fontsize=16)
    ax2.grid(True, which='major', linestyle='--', linewidth=0.5, color='#d3d3d3')
    ax2.tick_params(which='both', labelsize=20)
    ax2.set_title(f'Parameter Uncertainty', fontsize=16)
    plt.tight_layout()

    # Ensure directory exists and dynamically name saved files
    os.makedirs('Thesis_Figures', exist_ok=True)
    fig1.savefig(f'Thesis_Figures/{case_name}_models_prediction.svg', format='svg', bbox_inches='tight')
    fig2.savefig(f'Thesis_Figures/{case_name}_posterior_PDF.svg', format='svg', bbox_inches='tight')
    
    # Close figures to prevent memory leaks and overlapping plots in the loop
    plt.close(fig1)
    plt.close(fig2)


if __name__ == "__main__":
    # Case 1: Ideal Scenario 
    run_map_mcmc_case(N=20, fmom_true=0.5, Z_noise_percentage=0.1, p_mean_val=0.54, p_std_val=0.05, case_name="Case_1")

    # Case 2: Worst-case Scenario
    run_map_mcmc_case(N=4, fmom_true=0.5, Z_noise_percentage=0.8, p_mean_val=0.8, p_std_val=0.05, case_name="Case_2")

    print("\nAll cases executed successfully. SVGs saved in 'Thesis_Figures/' directory.")