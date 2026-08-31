import matplotlib
matplotlib.use('Qt5Agg')
import numpy as np
from numpy import linspace, log
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from Two_Point_Model_Solver import Extended_STPM

##################################################################################################
### This code aims to find a closed-form logaritmic expression for the target particle density ###
##################################################################################################

# Specifically, it fits the numerical results of the Two-Point model in a logarithmic function, in order to 
# find the log(n_target) with respect to n_upstream, loss_alpha, loss_fcool, loss_fconv

# Define the function's form log(nt) = g(nu,loss_alpha,loss_fcool, loss_fconv)
def log_nt_fit(paramet , a1,a2,a3,a4,a5):
    nu,alpha,f_cool,f_conv = paramet
    return a1 + a2*log(nu) + a3*log(alpha+1) + a4*log(f_cool) + a5*log(f_conv) # form (We write a1 instead of log(a1) to be <0)



# Extract data from the numerical solution of the Two-Point Model
def log_nt_data(n_upst, loss_alpha,loss_fcool,loss_fconv,Lc,x,theta,P_SOL):
    results = Extended_STPM(n_upst, loss_alpha,loss_fcool,loss_fconv,Lc,x,theta,P_SOL)  # Call the solver
    return log(results[0])  # target density


# Generate the training space
loss_alpha = linspace(0.0,10.0,30)
loss_fcool = linspace(0.5,0.9,20)
loss_fconv = linspace(0.5,0.9,20)
n_upst = linspace(1e17,5e20,200)
# Define the other parameters
Lc = 180           # [m] Estimated average connection length (=400 for W7-X, =1 for tokamaks)
x = 1.5
theta = 10e-3      # Field line pitch  (=0.001 for stellarators, =0.1 for tokamaks)
P_SOL = 100.0        # [MW] Heating power in the SO

# Generate the training data
log_nt_STPM_data = np.zeros((len(n_upst),len(loss_alpha),len(loss_fcool),len(loss_fconv)))
for k1 in range (len(n_upst)):
    for k2 in range (len(loss_alpha)):
        for k3 in range (len(loss_fcool)):
                for k4 in range (len(loss_fconv)):
                    log_nt_STPM_data[k1,k2,k3,k4] = log_nt_data(n_upst[k1], loss_alpha[k2], loss_fcool[k3], loss_fconv[k4],Lc,x,theta,P_SOL)


# Prossess the data so the function curve_fit can handle them
nu,alpha,fcool, fconv = np.meshgrid(n_upst, loss_alpha, loss_fcool, loss_fconv, indexing='ij')
x_data = (nu.ravel(), alpha.ravel(), fcool.ravel(), fconv.ravel())
y_data = log_nt_STPM_data.ravel()


# Fitting options
initial_guess = [-22.55, 1.5, -1.69, 1.54, -1.25]  # initial guess (very close to the solution; minimizes the iterrations)
# Fit the data
sol_opt, sol_cov = curve_fit(log_nt_fit, x_data, y_data, p0=initial_guess)
# Extract constants ai
a1_fit, a2_fit, a3_fit, a4_fit, a5_fit = sol_opt


## COMPARING THE RESULTS ##

# Calculate the function
log_nt_fitted = log_nt_fit((nu, alpha, fcool, fconv), *sol_opt)
RMSE = np.sqrt(np.mean((log_nt_fitted - log_nt_STPM_data)**2))
print(RMSE)   


# Print the results
a1_raw, a2_fit, a3_fit, a4_fit, a5_fit = sol_opt
a1_fit = np.exp(a1_raw)
print(f"The parameters of the prediction: nt = a1 * nu^a2 * (alpha+1)^a3 * fcool^a4 * fconv^a5")
print(f"Constants ai: [{a1_fit}, {a2_fit}, {a3_fit}, {a4_fit}, {a5_fit}]")



##########################
## PLOTTING THE RESULTS ##
##########################

# nt : Compare the solutions for nt = g(nu) 
plt.figure(1, figsize=(10, 6))
plt.plot(n_upst , log_nt_STPM_data[:,-1,-1,-1] , label='STPM Data (alpha=10, fcool=0.9, fconv=0.9)')
plt.plot(n_upst , log_nt_fitted[:,-1,-1,-1] , label='Fitted Function')
plt.xlabel('Upstream Density')
plt.ylabel('Logarithmic Target Density')
plt.title('Target Density as a function of nu')
plt.grid()
plt.legend()


# alpha : Compare the solutions for nt = g(alpha) 
plt.figure(2, figsize=(10, 6))
plt.plot(loss_alpha , log_nt_STPM_data[-1,:,-1,-1] , label='STPM Data (nu=5e20, fcool=0.9, fconv=0.9)')
plt.plot(loss_alpha , log_nt_fitted[-1,:,-1,-1] , label='Fitted Function')
plt.xlabel('alpha')
plt.ylabel('Logarithmic Target Density')
plt.title('Target Density as a function of alpha')
plt.grid()
plt.legend()


# fcool : Compare the solutions for nt = g(f_cool) 
plt.figure(3, figsize=(10, 6))
plt.plot(loss_fcool , log_nt_STPM_data[-1,-1,:,-1] , label='STPM Data (nu=5e20, alpha=10, fconv=0.9)')
plt.plot(loss_fcool , log_nt_fitted[-1,-1,:,-1] , label='Fitted Function')
plt.xlabel('f_cool')
plt.ylabel('Logarithmic Target Density')
plt.title('Target Density as a function of f_cool')
plt.grid()
plt.legend()



# fconv : Compare the solutions for nt = g(f_conv) 
plt.figure(4, figsize=(10, 6))
plt.plot(loss_fconv , log_nt_STPM_data[-1,-1,-1,:] , label='STPM Data (nu=5e20, alpha=10, fcool=0.9)')
plt.plot(loss_fconv , log_nt_fitted[-1,-1,-1,:] , label='Fitted Function')
plt.xlabel('f_conv')
plt.ylabel('Logarithmic Target Density')
plt.title('Target Density as a function of f_conv')
plt.grid()
plt.legend()
plt.show()