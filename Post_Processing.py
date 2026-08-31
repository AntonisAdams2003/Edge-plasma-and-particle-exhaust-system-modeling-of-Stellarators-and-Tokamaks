import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Qt5Agg')
import numpy as np
from numpy import pi
import pandas as pd
from scipy.interpolate import interp1d
from Two_Point_Model_Solver import Extended_STPM  # STELLARATOR TPM_Solver
from Two_Point_Model_Solver import TTPM, TTPM_numerical  # TOKAMAK TPM_Solver
from Divertor_Subdivertor_Solver import Divertor_Subdivertor_Analytical  # Divertor-Subdivertor Model

 

##################################################################################################
## This code aims to plot the results of the Two-Point model and Divertor_Subdivertor balances  ##
## in order to provide a qualitative results for some of the most important variables           ##
## like (density, temperature,...)                                                              ##
##################################################################################################

### Define the parameters of the Two-Point Model ###

# Define the vector nu
nu_upper = 1.0e21   # [m^-3] upper limit of the nu vector
nu_lower = 1.0e18   # [m^-3] lower limit of the nu vector
nu_elem = 90      # elements of the nu vector
nu = np.linspace(nu_lower, nu_upper, nu_elem)

# Define the loss parameters
alpha = [0.0 , 1.0 , 2.0 , 5.0 , 10.0]   # momentum losse factor
f_cool = [0.0 , 0.25 , 0.5 , 0.8 , 0.9]   # colling losses
f_conv = [0.0 , 0.3 , 0.6 , 0.8 , 0.9]   # convection fraction
# Define the other parameters
Lc = 180          # [m] Estimated average connection length (=400 for W7-X, =1 for tokamaks)
x = 1.5           # [m^2/s] perpendicular thermal diffusivities for electrons and ions
theta = np.array([0.25, 0.5, 1.0, 2.0, 5.0])*0.001      # Field line pitch  (=0.001 for stellarators, =0.1 for tokamaks)
q_SOL = [500.0 , 250.0 , 100.0 , 50.0 , 20.0]    # [MW/m^2] Heating power in the SOL


# Defining the parameters that we want to plot
alpha_elem = np.size(alpha)
fcool_elem = np.size(f_cool)
fconv_elem = np.size(f_conv)
q_elem = np.size(q_SOL)
theta_elem = np.size(theta)



#########################
### Stellarator Model ###
#########################

# - TPM OUTPUT: [0] n_targ, [1] T_upst, [2] T_targ, [3] p_targ
nt = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])     # [0] n_targ
Tu = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])     # [1] T_upst
Tt = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])     # [2] T_targ

#-OUTPUT: [0] Gamma_targ, [1] N_in_AEH, [2] N_in_AEP, [3] p_divertor_AEH, [4] p_divertor_AEP, [5] p_subdivertor_AEH, [6] p_subdivertor_AEP
Gamma_t = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])   # [0] Gamma_targ
N_in_AEH = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])   # [1] N_in_AEH
N_in_AEP = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])   # [2] N_in_AEP
p_div_AEH = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])   # [3] p_divertor_AEH
p_div_AEP = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])   # [4] p_divertor_AEP
p_sub_AEH = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])   # [5] p_subdivertor_AEH
p_sub_AEP = np.zeros([nu_elem,alpha_elem,fcool_elem,fconv_elem,theta_elem,q_elem])   # [6] p_subdivertor_AEP

# Calling the solver
for i1 in range (nu_elem):
    for i2 in range (alpha_elem):
        for i3 in range (fcool_elem):
            for i4 in range (fconv_elem):
                for i5 in range (theta_elem):
                    for i6 in range (q_elem):
 
                        n_upst = nu[i1]
                        loss_alpha = alpha[i2]
                        loss_fcool = f_cool[i3]
                        loss_fconv = f_conv[i4]
                        theta_pitch = theta[i5]
                        power_sol = q_SOL[i6]
                        

                        ## UPSTREAM -> TARGET ##

                        results_STPM = Extended_STPM(n_upst,loss_alpha,loss_fcool,loss_fconv , Lc,x,theta_pitch,power_sol)

                        # Target Particle Density
                        nt[i1,i2,i3,i4,i5,i6] = results_STPM[0]  # [m^-3]
                        # Upstream Temperature
                        Tu[i1,i2,i3,i4,i5,i6] = results_STPM[1]     # [eV]
                        # Target Temperature
                        Tt[i1,i2,i3,i4,i5,i6] = results_STPM[2]    # [eV]


            
                        ## TARGET -> SUBDIVERTOR ##

                        # Choose one of the three models of Divertor_Subdivetor_Solver.py
                        results_div = Divertor_Subdivertor_Analytical(results_STPM[0],results_STPM[2] , theta_pitch)

                        Gamma_t[i1,i2,i3,i4,i5,i6] = results_div[0]   # [particles/m^2/s]
                        # Particle flux in AEH
                        N_in_AEH[i1,i2,i3,i4,i5,i6] = results_div[1]   # [particles/s]
                        # Particle flux in AEP
                        N_in_AEP[i1,i2,i3,i4,i5,i6] = results_div[2]   # [particles/s]
                        # Divertor Pressure at AEH
                        p_div_AEH[i1,i2,i3,i4,i5,i6] = results_div[3]
                        # Divertor Pressure at AEP
                        p_div_AEP[i1,i2,i3,i4,i5,i6] = results_div[4]
                        # Subdivertor Pressure at AEH
                        p_sub_AEH[i1,i2,i3,i4,i5,i6] = results_div[5]
                        # Subdivertor Pressure at AEP
                        p_sub_AEP[i1,i2,i3,i4,i5,i6] = results_div[6]



#####################
### Tokamak Model ###
#####################

## Solve the TTPM numerically to take the exact models solution ##

#- TPM_numerical OUTPUT: [0] n_targ, [1] T_upst, [2] T_targ
nt_tok_numerical = np.zeros([nu_elem])     # [0] n_targ
Tu_tok_numerical = np.zeros([nu_elem])     # [1] T_upst
Tt_tok_numerical = np.zeros([nu_elem])     # [2] T_targ
for i in range (nu_elem):
    n_upst = nu[i]
    results_tomak = TTPM_numerical(n_upst , Lc,x,q_SOL[2])

    # Target Particle Density
    nt_tok_numerical[i] = results_tomak[0]   # [m^-3]
    # Upstream Temperature
    Tu_tok_numerical[i] = results_tomak[1]    # [eV]
    # Target Temperature
    Tt_tok_numerical[i] = results_tomak[2]    # [eV]


# Solve the analytical TTPM to extract the scaling of nt ##

#- TPM OUTPUT: [0] n_targ, [1] T_upst, [2] T_targ
nt_tok = np.zeros([nu_elem])     # [0] n_targ
for i in range (nu_elem):
    n_upst = nu[i]
    results_tomak = TTPM(n_upst , Lc,x,q_SOL[2])

    # Target Particle Density
    nt_tok[i] = results_tomak[0]   # [m^-3]



## DEFINE HAAK'S EXPERIMENTAL DATA ##
# AEH Haaks data
N_test_Haak_AEH = np.array([3.23e20 , 1.15e20])  # Discharge [.008 , .031]
p_div_test_Haak_AEH = np.array([2.6e-4 , 3.5e-5])
p_sub_test_Haak_AEH = np.array([1.4e-4 , 2.4e-5])
# AEP Haaks data
N_test_Haak_AEP = np.array([2.66e19 , 2.78e20])  # Discharge [.008 , .031]
p_div_test_Haak_AEP = np.array([1.5e-4 , 5.9e-4])
p_sub_test_Haak_AEP = np.array([8.2e-5 , 4.1e-4])



############################
### PLOTTING the RESULTS ###
############################

## Plotting Configuration ##
plt.rcParams.update({
    'text.usetex': False,
    'font.family': 'serif',
    'font.size': 26,    ## CHANGE
    'axes.labelsize': 26,  ## CHANGE
    'legend.fontsize': 16,   ## CHANGE
    'font.weight': 'normal',
    'axes.linewidth': 1.5,
    'lines.linewidth': 2.0,
    'axes.grid': True,
    'grid.alpha': 0.6,
    'mathtext.fontset': 'cm',       # Use Computer Modern (LaTeX font) for math
    'axes.formatter.use_mathtext': True,  # Renders 1e-4 as 10^{-4}
    'axes.formatter.limits': (-2, 3),     # Force scientific notation if < 10^-2 or > 10^3
    'axes.titlesize': 25,
    'xtick.labelsize': 23,
    'ytick.labelsize': 23
})

# Labels for plots
labels_alpha = [fr'$\alpha={val}$' for val in alpha]
labels_fcool = [fr'$f_{{cool}}={val}$' for val in f_cool]
labels_fconv = [fr'$f_{{conv}}={val}$' for val in f_conv]
labels_theta = [fr'$\Theta={val}$' for val in theta]
labels_q = [fr'$q_{{\parallel}}={val}$' for val in q_SOL]
# Colors for the plots
colors = ['#000004', '#3b0f70', '#8c2981', '#de4968', '#fe9f6d']

# format Python 
def format_sci(val):
   a, b = f"{val:.1e}".split('e')
   return fr"{a} \times 10^{{{int(b)}}}"

# Baseline parameters for the plots
baseline_parameters = {
    'alpha': fr'$\alpha = {alpha[2]}$',
    'fcool': fr'$f_{{cool}} = {f_cool[0]:.2f}$',
    'fconv': fr'$f_{{conv}} = {f_conv[0]:.2f}$',
    'theta': fr'$\Theta = {format_sci(theta[2])}$', 
    'q_par': fr'$q_{{\parallel}} = {q_SOL[2]}$' 
}

# This function takes the axes that we are plotting and removes the variale that is NOT fixed
def add_fixed_params_box(ax, varying_param_key):
    # Filter out the varying parameter using a list comprehension
    fixed_lines = [
        text for key, text in baseline_parameters.items() 
        if key != varying_param_key
    ]
    # Join the remaining strings with a newline
    param_text = "Fixed Params:\n" + "\n".join(fixed_lines)
    ax.text(0.95, 0.05, param_text, 
            transform=ax.transAxes, 
            fontsize=16,
            verticalalignment='bottom', 
            horizontalalignment='right',
            bbox=dict(boxstyle='square,pad=0.8', facecolor='white', edgecolor='black', linewidth=1)
            )



# Define a function that plots the two-point model results
def plot_TPM(param_array, labels, param_key, a_i, fcool_i, fconv_i, theta_i, q_i):

    if param_key == 'theta':
        labels = [fr'$\Theta = {format_sci(val)}$' for val in param_array]

    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    axs[0].set_prop_cycle(color=colors)
    
    axs[0].loglog(nu, nu/2, '--k' , label='$n_t = n_u/2$')
    axs[0].loglog(nu, nt_tok, '--k' , label=r'$n_t \propto n_u^3$')

    for idx in range(len(param_array)):
        current_color = colors[idx % len(colors)]
        
        # Replaces the variable index with the loop index
        idx_tuple = (
            slice(None), 
            idx if param_key == 'alpha' else a_i,
            idx if param_key == 'fcool' else fcool_i,
            idx if param_key == 'fconv' else fconv_i,
            idx if param_key == 'theta' else theta_i,
            idx if param_key == 'q_par' else q_i
        )
        
        axs[0].loglog(nu, nt[idx_tuple], color=current_color, label=labels[idx])
        axs[1].loglog(nu, Tt[idx_tuple], '-', color=current_color)
        axs[1].loglog(nu, Tu[idx_tuple], '--', color=current_color)

    # TITLE FORMATTING MAPPING
    title_var = param_key
    if param_key == 'q_par': title_var = r'q_{\parallel}'
    elif param_key == 'theta': title_var = r'\Theta'
    elif param_key == 'alpha': title_var = r'\alpha'
    elif param_key == 'fcool': title_var = 'f_{cool}'
    elif param_key == 'fconv': title_var = 'f_{conv}'

    axs[0].set_title(fr'$n_t = f(n_u; {title_var})$')
    axs[0].set_ylabel(r'$n_t$ [m$^{-3}$]')
    axs[0].legend(frameon=False)
    axs[0].axis([1.0e18, 1e21, 5.0e17, 1e23])

    axs[1].set_title(fr'$T = f(n_u; {title_var})$')
    axs[1].set_xlabel(r'$n_u$ [m$^{-3}$]')
    axs[1].set_ylabel(r'$T$ [eV]')
    axs[1].plot([], [], '-', color='black', label='$T_t$')
    axs[1].plot([], [], '--', color='black', label='$T_u$')
    axs[1].legend(frameon=False)
    axs[1].axis([1.0e18, 1e21, 0.1e0, 2.5e2])
    
    add_fixed_params_box(axs[0], param_key)

    for ax in axs:
        ax.grid(False)
        ax.tick_params(axis='both', which='major', direction='in', length=6, width=1, top=True, right=True)
        ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1, top=True, right=True)

    plt.tight_layout()
    return fig


# Function that plots the subdivertor results
def plot_subdivertor(var_AEP, var_AEH, y_label, title_base, param_array, labels, param_key, a_i, fcool_i, fconv_i, theta_i, q_i):

    if param_key == 'theta':
        labels = [fr'$\Theta = {format_sci(val)}$' for val in param_array]

    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    axs[0].set_prop_cycle(color=colors)
    axs[1].set_prop_cycle(color=colors)

    for idx in range(len(param_array)):
        current_color = colors[idx % len(colors)]
        
        idx_tuple = (
            slice(None), 
            idx if param_key == 'alpha' else a_i,
            idx if param_key == 'fcool' else fcool_i,
            idx if param_key == 'fconv' else fconv_i,
            idx if param_key == 'theta' else theta_i,
            idx if param_key == 'q_par' else q_i
        )

        axs[0].loglog(nu, var_AEP[idx_tuple], color=current_color, label=labels[idx])
        axs[1].loglog(nu, var_AEH[idx_tuple], color=current_color, label=labels[idx])

    axs[0].set_title(fr'${{{title_base}}}_{{AEP}} = f(n_u; \alpha)$')
    axs[0].set_ylabel(y_label)
    axs[0].legend(frameon=False)

    axs[1].set_title(fr'${{{title_base}}}_{{AEH}} = f(n_u; \alpha)$')
    axs[1].set_xlabel(r'$n_u$ [m$^{-3}$]')
    axs[1].set_ylabel(y_label)
    axs[1].legend(frameon=False)

    for ax in axs:
        ax.grid(False)
        ax.tick_params(axis='both', which='major', direction='in', length=6, width=1, top=True, right=True)
        ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1, top=True, right=True)

    plt.tight_layout()
    return fig



## Plotting the Figures 1-5

on_figures_15 = 1
if on_figures_15 == 1:
    ## Figures 1-5 ##
    # Plot the two-point model
    fig1 = plot_TPM(alpha, labels_alpha, 'alpha', 0, 0, 0, 2, 2)
    fig2 = plot_TPM(f_cool, labels_fcool, 'fcool', 2, 0, 0, 2, 2)
    fig3 = plot_TPM(f_conv, labels_fconv, 'fconv', 2, 0, 0, 2, 2)
    fig4 = plot_TPM(theta, labels_theta, 'theta', 2, 0, 0, 0, 2)
    fig5 = plot_TPM(q_SOL, labels_q, 'q_par', 2, 0, 0, 2, 0)

    # Save the figures
    fig1.savefig('Thesis_Figures/STPM_parameter_alpha.svg', format='svg', bbox_inches='tight')
    fig2.savefig('Thesis_Figures/STPM_parameter_fcool.svg', format='svg', bbox_inches='tight')
    fig3.savefig('Thesis_Figures/STPM_parameter_fconv.svg', format='svg', bbox_inches='tight')
    fig4.savefig('Thesis_Figures/STPM_parameter_theta.svg', format='svg', bbox_inches='tight')
    fig5.savefig('Thesis_Figures/STPM_parameter_q.svg', format='svg', bbox_inches='tight')



## Plotting Figure 6: TTPM vs STPM

on_figure_6 = 1
if on_figure_6 == 1:

    # Figure 6: Tokamak vs Stellarator
    fig6, axs6 = plt.subplots(2, 1, figsize=(10, 10))

    a_idx = 2  
    f_conv_idx = 0  
    f_cool_idx = 0
    theta_idx = 2
    q_idx = 2

    # Subplot 1
    axs6[0].loglog(nu, nu/2, '--k' , label='$n_t = n_u/2$')
    axs6[0].loglog(nu, nt_tok, '--k' , label=r'$n_t \propto n_u^3$')
    #

    axs6[0].loglog(nu, nt[:, 0, f_cool_idx, f_conv_idx, theta_idx, q_idx], linestyle='-', color='#3b0f70', label='STPM')
    axs6[0].loglog(nu, nt[:, a_idx, f_cool_idx, f_conv_idx, theta_idx, q_idx], linestyle='-', color='#de4968',label=r'Extended STPM, $f_{mom}(\alpha=2)$')
    axs6[0].loglog(nu, nt_tok_numerical, linestyle='-', color='#fe9f6d',label='TTPM')

    # Subplot 2
    axs6[1].loglog(nu, Tu[:, 0, f_cool_idx, f_conv_idx, theta_idx, q_idx], linestyle='--', color='#3b0f70')   # STPM
    axs6[1].loglog(nu, Tu[:, a_idx, f_cool_idx, f_conv_idx, theta_idx, q_idx], linestyle='--', color='#de4968')   # Extended STPM
    axs6[1].loglog(nu, Tu_tok_numerical, linestyle='--', color='#fe9f6d')  # TTPM

    axs6[1].loglog(nu, Tt[:, 0, f_cool_idx, f_conv_idx, theta_idx, q_idx], linestyle='-', color='#3b0f70')
    axs6[1].loglog(nu, Tt[:, a_idx, f_cool_idx, f_conv_idx, theta_idx, q_idx], linestyle='-', color='#de4968')
    axs6[1].loglog(nu, Tt_tok_numerical, linestyle='-', color='#fe9f6d')

    axs6[0].set_title(fr'$n_t = f(n_u)$, for $ \alpha = {alpha[a_idx]}$')
    axs6[0].set_ylabel(r'$n_t$ [m$^{-3}$]')
    axs6[0].legend(frameon=False)
    axs6[0].axis([1.0e18, 1e21, 5.0e17, 1e23])

    axs6[1].set_title(fr'$T = f(n_u)$')
    axs6[1].set_xlabel(r'$n_u$ [m$^{-3}$]')
    axs6[1].set_ylabel(r'$T$ [eV]')
    axs6[1].plot([], [], '-', color='black', label='$T_t$')
    axs6[1].plot([], [], '--', color='black', label='$T_u$')
    axs6[1].legend(frameon=False)
    axs6[1].axis([1.0e18, 1e21, 0.1e0, 2.5e2])

    # Add vertical dashed lines for physical regimes
    axs6[0].axvline(4e18, color='gray', linestyle='--', linewidth=1)
    axs6[0].axvline(3e19, color='gray', linestyle='--', linewidth=1)
    axs6[1].axvline(4e18, color='gray', linestyle='--', linewidth=1)
    axs6[1].axvline(3e19, color='gray', linestyle='--', linewidth=1)

    # Add text annotations for region names
    axs6[1].text(2e18, 5, 'sheath\nlimited', color='gray', fontsize=16, ha='center', va='center')
    axs6[1].text(1.1e19, 5, 'conduction\nlimited', color='gray', fontsize=16, ha='center', va='center')
    axs6[1].text(2e20, 5, 'diffusion\nlimited', color='gray', fontsize=16, ha='center', va='center')

    for ax in axs6:
        ax.grid(False)
        ax.tick_params(axis='both', which='major', direction='in', length=6, width=1, top=True, right=True)
        ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1, top=True, right=True)

    plt.tight_layout()


    # Save figure
    fig6.savefig('Thesis_Figures/TTPM_vs_STPM.svg', format='svg', bbox_inches='tight')





## Plotting figures 7-9

on_figures_79 = 1
if on_figures_79 == 1: 
    # Plot the subdivertor variables
    fig7 = plot_subdivertor(N_in_AEP, N_in_AEH, r'$N_{in}$ [s$^{-1}$]', 'N_{in}', alpha, labels_alpha, 'alpha', 0, 0, 0, 2, 2)
    fig8 = plot_subdivertor(p_div_AEP, p_div_AEH, r'$p_{div}$ [mbar]', 'p_{div}', alpha, labels_alpha, 'alpha', 0, 0, 0, 2, 2)
    fig9 = plot_subdivertor(p_sub_AEP, p_sub_AEH, r'$p_{sub}$ [mbar]', 'p_{sub}', alpha, labels_alpha, 'alpha', 0, 0, 0, 2, 2)

    # Save figures
    fig7.savefig('Thesis_Figures/Nin_parametric.svg', format='svg', bbox_inches='tight')
    fig8.savefig('Thesis_Figures/pdiv_parametric.svg', format='svg', bbox_inches='tight')
    fig9.savefig('Thesis_Figures/psub_parametric.svg', format='svg', bbox_inches='tight')





def generate_comparison_table(N_exp, p_div_exp, p_sub_exp, N_model, p_div_model, p_sub_model, nu_model, config_name):
    inv_N_func = interp1d(N_model, nu_model, kind='linear', bounds_error=False, fill_value="extrapolate")
    p_div_func = interp1d(nu_model, p_div_model, kind='linear', bounds_error=False, fill_value="extrapolate")
    p_sub_func = interp1d(nu_model, p_sub_model, kind='linear', bounds_error=False, fill_value="extrapolate")

    nu_pred = inv_N_func(N_exp)
    p_div_pred = p_div_func(nu_pred)
    p_sub_pred = p_sub_func(nu_pred)

    df = pd.DataFrame({
        'Configuration': config_name,
        'Discharge': ['0.008', '0.031'],
        'N_in_exp': N_exp,
        'p_div_exp': p_div_exp,
        'p_div_pred': p_div_pred,
        'p_sub_exp': p_sub_exp,
        'p_sub_pred': p_sub_pred
        })
    return df



# Compare with experimental results
alpha_indices_to_test = [2]  # Example: Index 2 -> alpha=2.0, Index 4 -> alpha=10.0. Modify this list to test others.
all_dfs = []

for a_idx in alpha_indices_to_test:
    current_idx = (slice(None), a_idx, 0, 0, 2, 2)   # [alpha, fcool, fconv, theta, q_par]

    df_AEH = generate_comparison_table(
        N_test_Haak_AEH, p_div_test_Haak_AEH, p_sub_test_Haak_AEH, N_in_AEH[current_idx], p_div_AEH[current_idx], p_sub_AEH[current_idx], nu, f"AEH (alpha={alpha[a_idx]})"
    )

    df_AEP = generate_comparison_table(
        N_test_Haak_AEP, p_div_test_Haak_AEP, p_sub_test_Haak_AEP, N_in_AEP[current_idx], p_div_AEP[current_idx], p_sub_AEP[current_idx], nu, f"AEP (alpha={alpha[a_idx]})"
    )
    
    all_dfs.extend([df_AEH, df_AEP])

comparison_table = pd.concat(all_dfs, ignore_index=True)

pd.set_option('display.float_format', lambda x: f'{x:.2e}')
print("\n--- EXPERIMENTAL VS MODEL PREDICTION TABLE ---")
print(comparison_table.to_string(index=False))
print("----------------------------------------------\n")


plt.show()