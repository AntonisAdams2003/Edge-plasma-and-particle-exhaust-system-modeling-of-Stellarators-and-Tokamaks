import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.use('Qt5Agg')
import numpy as np
from numpy import pi, sqrt, cos, exp
import math
from scipy.optimize import fsolve
import pandas as pd
import sys
sys.stdout.reconfigure(encoding='utf-8')



#############################################################################################
## This code assumes that the particle influx passing the pumping gaps (N_in) is known     ##
## and provides some functions that estimate the subdivertor pressure (p_sub) in [mbar].   ##
## Moreover, it compares the models with experimental data and modeling data from          ##
## high-fidelity codes.                                                                    ##
## Finally, it provides with an estimaton for the Knudsen number in the AEH and AEP ports. ##
#############################################################################################


# Define the known constants and geometric dimentions
m_H = 1.67262e-27  # [kg] mass hydrogen ion
m_H2 = 3.34524e-27  # [kg] mass of H2  = 2 x m_H
e = 1.60218e-19    # electron charge
kB = 1.380649*1e-23   # [J/K] Boltzmann's constant
a_r = 0.5   # [m] minor radius of W7-X
R = 5.5   # [m] major radius of W7-X
d = 0.4    # [m]  diameter of vacuum cylinders cross-sections
A_pg_AEH = 0.153   # [m^2] punping gap area for AEH (low div)
A_pg_AEP = 0.045   # [m^2] punping gap area for AEP (high div)

# Leakages
A_leakages_AEH = np.array([1.97,5.79,5.71,5.41,1.1,7.81,7.44,7.14,6.21,10.3,10.9,3.77,3.47])*1e-3
A_leakages_AEP = np.array([3.02,3.83,4.11,3.58,3.39,4.78,14.4,12.5,5.59])*1e-3

A_leak_AEH = np.sum(A_leakages_AEH)   # 7.7e-2 [m^2] leakage area for low divertor
A_leak_AEP = np.sum(A_leakages_AEP)   # 5.52e-2 [m^2] leakage area for high divertor
# A_leak_AEH = 0.3  # [m^2] leakage area for low divertor
# A_leak_AEP = 0.22 # [m^2] leakage area for high divertor


# Define constants by assuming typical values
gamma = 8   # Sheath transmission coefficient
koe = 2000  # [W/m eV^7/2] Spitzer-Härm heat conductivity for electrons
koi = 60    # [W/m eV^7/2] Spitzer-Härm heat conductivity for ions
Tvv = 303   # [K] Temperature of the vacuum vessel
T_div = 300  # [K] Temperature at the Divertor
T_sub = 300  # [K] Temperature at the Sub-divertor
T_leak = 300
T_0 = 600    # [K] Temperature of the H that enters
ji_AEH = 0.06  # probability of particle pumped out for AEH
ji_AEP = 0.0264  # probability of particle pumped out for AEP

# Imposed bulk flow velocity
s = 0.5/sqrt(pi)
th = 0.0
B = exp(-s**2*cos(th)**2) + sqrt(pi)*s*cos(th)*(1 + math.erf(s*cos(th)))
# Pumping speed
S_eff_AEP = (1.0/4)*ji_AEP*(pi*d**2/4)*sqrt(8*kB*Tvv/(pi*m_H2))   # [m^3/s] V.Haak -> 1.18 AEP  
S_eff_AEH = (1.0/4)*ji_AEH*(pi*d**2/4)*sqrt(8*kB*Tvv/(pi*m_H2))   # [m^3/s] V.Haak -> 2.35 AEH

# Thermal velocities of divertor and subdivertor gas
vth_0 = sqrt(8*kB*T_0/(pi*m_H2))   # Divertor
vth_leak = sqrt(8*kB*T_leak/(pi*m_H2))   # Subdivertor



##############
### MODELS ###
##############


## VAROUTIS ##
def Varoutis_subdivertor(N_in_AEH, N_in_AEP):
    p_sub_AEH = 0.01 * 5e-25 * N_in_AEH**1.0747
    p_sub_AEP = 0.01 * 4e-24 * N_in_AEP**1.0779
    return [p_sub_AEH,p_sub_AEP]




## DIRK ##
def Dirk_subdivertor(N_in_AEH,N_in_AEP):
  
    ## HAAK'S PAPER DATA ##
    A_pg_AEH_Dirk = A_pg_AEH #0.15   # [m^2] punping gap area for AEH (low div)
    A_pg_AEP_Dirk = A_pg_AEP #0.06   # [m^2] punping gap area for AEP (high div)

    # Leakage area
    A_leak_AEH_Dirk = A_leak_AEH #0.3
    A_leak_AEP_Dirk = A_leak_AEP #0.22

    # Pressure in the divertor
    h_AEH = sqrt(pi*m_H2*kB*T_div/2.0) / (A_leak_AEH_Dirk + A_pg_AEH_Dirk*S_eff_AEH/(A_pg_AEH_Dirk*sqrt(kB*T_sub/(2.0*pi*m_H2)) + S_eff_AEH))
    h_AEP = sqrt(pi*m_H2*kB*T_div/2.0) / (A_leak_AEP_Dirk + A_pg_AEP_Dirk*S_eff_AEP/(A_pg_AEP_Dirk*sqrt(kB*T_sub/(2.0*pi*m_H2)) + S_eff_AEP))
    p_divertor_AEH = (h_AEH * N_in_AEH) * 0.01  # [mbar]
    p_divertor_AEP = (h_AEP * N_in_AEP) * 0.01  # [mbar]
    
    # Reduce constants of pressure as the gas flows to the sub-divertor
    w_AEH = A_pg_AEH_Dirk*sqrt(kB*T_div) / (sqrt(2*pi*m_H2)*S_eff_AEH + A_leak_AEH_Dirk*sqrt(kB*T_sub))
    w_AEP = A_pg_AEP_Dirk*sqrt(kB*T_div) / (sqrt(2*pi*m_H2)*S_eff_AEP + A_leak_AEP_Dirk*sqrt(kB*T_sub))
    # Estimating the sub-divertor's pressure 
    p_sub_AEH = p_divertor_AEH * w_AEH   # [mbar]
    p_sub_AEP = p_divertor_AEP * w_AEP   # [mbar]
    
    return [p_sub_AEH,p_sub_AEP]




## CONDUCTANCE LITOVOLI-HAAK 2026 ##
def Litovoli_Haak_subdivertor(N_in_AEH, N_in_AEP):

    # We should choose the assumption of in what themperature does the outfluxed gases leave the pumping gap
    T_outflux = T_0   # Choose: T_leak or T_0
    vth_outflux = sqrt(8*kB*T_outflux/(pi*m_H2))  # thermal speed

    # Define the model
    def System_of_Equations(p, N_in_1, N_in_2):
        p1, p2 = p

        ## THROUGHPUTS (Q) ##
        # leakages
        q_leak_1 = B * 0.25 * A_leak_AEH * p1 * vth_leak   # AEH
        q_leak_2 = B * 0.25 * A_leak_AEP * p2 * vth_leak   # AEP
        # pump
        q_pump_1 = S_eff_AEH * p1  # AEH
        q_pump_2 = S_eff_AEP * p2  # AEP
        # outflux
        q_outflux_1 = B * 0.25 * A_pg_AEH * p1 * vth_outflux   # AEH
        q_outflux_2 = B * 0.25 * A_pg_AEP * p2 * vth_outflux   # AEP

        C12 = 0.0
        q_C12 = C12 * (p1 - p2)

        ## PARTICLE BALANCE (N=q/kT) ##
        # We divide the outgoing throughputs by their respective temperatures to balance PARTICLES.
        eq1 = N_in_1*kB*T_0 - q_outflux_1 - q_leak_1 - q_C12 - q_pump_1
        eq2 = N_in_2*kB*T_0 - q_outflux_2 - q_leak_2 + q_C12 - q_pump_2

        return [eq1, eq2]

    # Initial guess for solving the system
    p_initial = [1e-3, 1e-3]

    # Call the solver to solve the system (Passing N_in directly)
    p_optimal = fsolve(System_of_Equations, p_initial, args=(N_in_AEH, N_in_AEP))

    # Pressures in Pascals
    p_sub_Pa_AEH = p_optimal[0]
    p_sub_Pa_AEP = p_optimal[1]

    # Find the neutral gas flows 
    N_pump_AEH = S_eff_AEH * p_sub_Pa_AEH / (kB*T_sub)  
    N_pump_AEP = S_eff_AEP * p_sub_Pa_AEP / (kB*T_sub)  
    
    N_outflux_AEH = B * 0.25 * A_pg_AEH * p_sub_Pa_AEH * vth_outflux / (kB*T_outflux)  
    N_outflux_AEP = B * 0.25 * A_pg_AEP * p_sub_Pa_AEP * vth_outflux / (kB*T_outflux)  
    
    N_leak_AEH = B * 0.25 * A_leakages_AEH * p_sub_Pa_AEH * vth_leak / (kB*T_leak) 
    N_leak_AEP = B * 0.25 * A_leakages_AEP * p_sub_Pa_AEP * vth_leak / (kB*T_leak) 

    #  Convert to mbar
    p_sub_AEH = p_sub_Pa_AEH * 0.01   
    p_sub_AEP = p_sub_Pa_AEP * 0.01   

    # Ratios
    N_out_in_AEH = N_outflux_AEH / N_in_AEH
    N_out_in_AEP = N_outflux_AEP / N_in_AEP
    N_leak_in_AEH = N_leak_AEH / N_in_AEH
    N_leak_in_AEP = N_leak_AEP / N_in_AEP
    N_pump_in_AEH = N_pump_AEH / N_in_AEH
    N_pump_in_AEP = N_pump_AEP / N_in_AEP

    return [p_sub_AEH, p_sub_AEP, N_out_in_AEH, N_out_in_AEP, N_leak_in_AEH, N_leak_in_AEP, N_pump_in_AEH, N_pump_in_AEP]



# AEH
numerator_AEH = (B / 4) * A_pg_AEH * (vth_0)
denominator_AEH = (B / 4) * (A_leak_AEH * vth_leak + A_pg_AEH * vth_0) + S_eff_AEH
print(2-numerator_AEH/denominator_AEH)




## MINE MODEL: MASS BALANCE ##
def ParticleBalance_subdivertor(N_in_AEH, N_in_AEP):

    # We should choose the assumption of in what themperature does the outfluxed gases leave the pumping gap
    T_outflux = T_0   # Choose: T_leak or T_0
    vth_outflux = sqrt(8*kB*T_outflux/(pi*m_H2))  # thermal speed

    # Define the model
    def System_of_Equations(p, N_in_1, N_in_2):
        p1, p2 = p

        ## THROUGHPUTS (Q) ##
        # leakages
        q_leak_1 = B * 0.25 * A_leak_AEH * p1 * vth_leak   # AEH
        q_leak_2 = B * 0.25 * A_leak_AEP * p2 * vth_leak   # AEP
        # pump
        q_pump_1 = S_eff_AEH * p1  # AEH
        q_pump_2 = S_eff_AEP * p2  # AEP
        # outflux
        q_outflux_1 = B * 0.25 * A_pg_AEH * p1 * vth_outflux   # AEH
        q_outflux_2 = B * 0.25 * A_pg_AEP * p2 * vth_outflux   # AEP

        C12 = 0.0
        q_C12 = C12 * (p1 - p2)

        ## PARTICLE BALANCE EQUATIONS 
        # We divide the outgoing throughputs by their respective temperatures to balance PARTICLES.
        eq1 = N_in_1 - (q_outflux_1 / (kB*T_outflux)) - (q_leak_1 / (kB*T_leak)) - (q_C12 / (kB*T_sub)) - (q_pump_1 / (kB*T_sub))
        eq2 = N_in_2 - (q_outflux_2 / (kB*T_outflux)) - (q_leak_2 / (kB*T_leak)) + (q_C12 / (kB*T_sub)) - (q_pump_2 / (kB*T_sub))

        return [eq1, eq2]

    # Initial guess for solving the system
    p_initial = [1e-3, 1e-3]

    # Call the solver to solve the system (Passing N_in directly)
    p_optimal = fsolve(System_of_Equations, p_initial, args=(N_in_AEH, N_in_AEP))

    # Pressures in Pascals
    p_sub_Pa_AEH = p_optimal[0]
    p_sub_Pa_AEP = p_optimal[1]

    # Find the neutral gas flows
    N_pump_AEH = S_eff_AEH * p_sub_Pa_AEH / (kB*T_sub)  
    N_pump_AEP = S_eff_AEP * p_sub_Pa_AEP / (kB*T_sub)  
    
    N_outflux_AEH = B * 0.25 * A_pg_AEH * p_sub_Pa_AEH * vth_outflux / (kB*T_outflux)  
    N_outflux_AEP = B * 0.25 * A_pg_AEP * p_sub_Pa_AEP * vth_outflux / (kB*T_outflux)  
    
    N_leak_AEH = B * 0.25 * A_leakages_AEH * p_sub_Pa_AEH * vth_leak / (kB*T_sub) 
    N_leak_AEP = B * 0.25 * A_leakages_AEP * p_sub_Pa_AEP * vth_leak / (kB*T_sub) 

    # Convert to mbar
    p_sub_AEH = p_sub_Pa_AEH * 0.01   
    p_sub_AEP = p_sub_Pa_AEP * 0.01   

    # Ratios
    N_out_in_AEH = N_outflux_AEH / N_in_AEH
    N_out_in_AEP = N_outflux_AEP / N_in_AEP
    N_leak_in_AEH = N_leak_AEH / N_in_AEH
    N_leak_in_AEP = N_leak_AEP / N_in_AEP
    N_pump_in_AEH = N_pump_AEH / N_in_AEH
    N_pump_in_AEP = N_pump_AEP / N_in_AEP

    return [p_sub_AEH, p_sub_AEP, N_out_in_AEH, N_out_in_AEP, N_leak_in_AEH, N_leak_in_AEP, N_pump_in_AEH, N_pump_in_AEP]











###########################
### PRESSURE COMPARISON ###
###########################

## DEFINE HAAK'S EXPERIMENTAL DATA ##
# AEH Haaks data
N_test_Haak_AEH = np.array([3.23e20 , 1.15e20])  # Discharge [.008 , .031]
p_sub_test_Haak_AEH = np.array([1.4e-4 , 2.4e-5])
# AEP Haaks data
N_test_Haak_AEP = np.array([2.66e19 , 2.78e20])  # Discharge [.008 , .031]
p_sub_test_Haak_AEP = np.array([8.2e-5 , 4.1e-4])



## Compare the Model's pressures ##

# Initialize the vectors
p_sub_test_varout_AEH = np.zeros_like(p_sub_test_Haak_AEH)
p_sub_test_varout_AEP = np.zeros_like(p_sub_test_Haak_AEP)
#
p_sub_test_dirk_AEH = np.zeros_like(p_sub_test_Haak_AEH)
p_sub_test_dirk_AEP = np.zeros_like(p_sub_test_Haak_AEP)
#
p_sub_test_litovoli_AEH = np.zeros_like(p_sub_test_Haak_AEH)
p_sub_test_litovoli_AEP = np.zeros_like(p_sub_test_Haak_AEP)
#
p_sub_test_particle_AEH = np.zeros_like(p_sub_test_Haak_AEH)
p_sub_test_particle_AEP = np.zeros_like(p_sub_test_Haak_AEP)

# Calculate pressures for each discharge
for i in range(len(N_test_Haak_AEH)):
    # Varoutis
    p_sub_test_varout_AEH[i],p_sub_test_varout_AEP[i] = Varoutis_subdivertor(N_test_Haak_AEH[i],N_test_Haak_AEP[i])
    # Dirk
    p_sub_test_dirk_AEH[i],p_sub_test_dirk_AEP[i] = Dirk_subdivertor(N_test_Haak_AEH[i],N_test_Haak_AEP[i])[:2]
    # Conductance Haak's
    p_sub_test_litovoli_AEH[i],p_sub_test_litovoli_AEP[i] = Litovoli_Haak_subdivertor(N_test_Haak_AEH[i],N_test_Haak_AEP[i])[:2]
    # Mine Mass Balance
    p_sub_test_particle_AEH[i],p_sub_test_particle_AEP[i] = ParticleBalance_subdivertor(N_test_Haak_AEH[i],N_test_Haak_AEP[i])[:2]




## Display the table ##

# We combine the AEP and AEH lists into one sequence of 4 values
N_all = np.concatenate((N_test_Haak_AEH, N_test_Haak_AEP))   # Discharge particle flows
p_sub_Haak_all = np.concatenate((p_sub_test_Haak_AEH, p_sub_test_Haak_AEP))   # coresponding experimental pressures
p_sub_Varoutis_all = list(p_sub_test_varout_AEH) + list(p_sub_test_varout_AEP)  # Varoutis pressures
p_sub_Dirk_all = list(p_sub_test_dirk_AEH) + list(p_sub_test_dirk_AEP)    # Dirks pressures
p_sub_Conduct_all = list(p_sub_test_litovoli_AEH) + list(p_sub_test_litovoli_AEP)   # Haak's conductance pressures
p_sub_Mass_all = list(p_sub_test_particle_AEH) + list(p_sub_test_particle_AEP)    # My mass balance pressures

# Create labels for the index
labels = ["AEH (.008)","AEH (.031)", 
          "AEP (.008)", "AEP (.031)"]
# Create DataFrame with Units in the columns for clarity
df = pd.DataFrame({
    "N_in [part/s]": N_all,
    "Experimental [mbar]": p_sub_Haak_all,
    "Varoutis [mbar]": p_sub_Varoutis_all,
    "Dirk [mbar]": p_sub_Dirk_all,
    "Litovoli-Haak [mbar]": p_sub_Conduct_all,
    "Particle Balance [mbar]": p_sub_Mass_all
}, index=labels)

# Formatting options
pd.set_option("display.float_format", "{:.2e}".format)
pd.set_option('display.expand_frame_repr', False) # Prevents wrapping to next line
pd.set_option('display.colheader_justify', 'center') # Centers column headers

print(f"{'MODEL COMPARISON TABLE':^60}")
print("="*60 + "\n")
print(df)




#############################
### KNUDSEN NUMBER CHECK ###
#############################

# Constants for Knudsen calculation
d_H2 = 2.89e-10  # [m] Molecular diameter of H2
sigma_H2 = pi * d_H2**2  # [m^2] Collision cross-section

# Highest experimental pressures in Pascals (worst-case density)
P_max_AEH_Pa = np.max(p_sub_test_Haak_AEH) * 100 
P_max_AEP_Pa = np.max(p_sub_test_Haak_AEP) * 100

# Number density (n = P / kB*T) assuming 300 K
n_AEH = P_max_AEH_Pa / (kB * 300)
n_AEP = P_max_AEP_Pa / (kB * 300)

# Mean free path (lambda)
lambda_AEH = 1.0 / (sqrt(2) * n_AEH * sigma_H2)
lambda_AEP = 1.0 / (sqrt(2) * n_AEP * sigma_H2)

# Characteristic lengths (sqrt of Area)
L_AEH = np.sqrt(A_pg_AEH)
L_AEP = np.sqrt(A_pg_AEP)

# Knudsen Numbers
Kn_AEH = lambda_AEH / L_AEH
Kn_AEP = lambda_AEP / L_AEP

print(f"\n\n{'KNUDSEN NUMBER':^60}")
print("="*60)
print(f"AEH Port: lambda = {lambda_AEH:.2f} m | Kn = {Kn_AEH:.2f}")
print(f"AEP Port: lambda = {lambda_AEP:.2f} m | Kn = {Kn_AEP:.2f}")




## Plot the Models ##

plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14, 
    'legend.fontsize': 14,
    'lines.linewidth': 2.5
})

n_elem = 100
N_in_grid = np.logspace(19, 24, n_elem)

# Initialize the vectors
p_sub_varout_AEH = np.zeros(n_elem)
p_sub_varout_AEP = np.zeros(n_elem)
#
p_sub_dirk_AEH = np.zeros(n_elem)
p_sub_dirk_AEP = np.zeros(n_elem)
#
p_sub_conduct_AEH = np.zeros(n_elem)
p_sub_conduct_AEP = np.zeros(n_elem)
#
p_sub_mass_AEH = np.zeros(n_elem)
p_sub_mass_AEP = np.zeros(n_elem)

# Call the functions for each N_in value
for i in range(n_elem):
    # Varoutis
    p_sub_varout_AEH[i],p_sub_varout_AEP[i] = Varoutis_subdivertor(N_in_grid[i],N_in_grid[i])
    # Dirk
    p_sub_dirk_AEH[i],p_sub_dirk_AEP[i] = Dirk_subdivertor(N_in_grid[i],N_in_grid[i])[:2]
    # Conductance Haak's
    p_sub_conduct_AEH[i],p_sub_conduct_AEP[i] = Litovoli_Haak_subdivertor(N_in_grid[i],N_in_grid[i])[:2]
    # Mine Mass Balance
    p_sub_mass_AEH[i],p_sub_mass_AEP[i] = ParticleBalance_subdivertor(N_in_grid[i],N_in_grid[i])[:2]

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 6), sharey=True)

# Colors
c_varout = '#1565C0'  
c_conduct = '#388E3C' 
c_mass = '#D32F2F'    
c_dirk = '#F57C00'

plot_models = {
    "DIVGAS": True,
    "Dirk": False,
    "Throughput balance": True,
    "Particle balance": True
}

bbox_props = dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8)

## SUBPLOT 1: AEH ##
if plot_models["DIVGAS"]:
    ax1.loglog(N_in_grid, p_sub_varout_AEH, color=c_varout, linestyle='-', label="DIVGAS model")
if plot_models["Dirk"]:
    ax1.loglog(N_in_grid, p_sub_dirk_AEH, color=c_dirk, linestyle='-', label="Dirks model")
if plot_models["Throughput balance"]:
    ax1.loglog(N_in_grid, p_sub_conduct_AEH, color=c_conduct, linestyle='-', label="Throughput balance model")
if plot_models["Particle balance"]:
    ax1.loglog(N_in_grid, p_sub_mass_AEH, color=c_mass, linestyle='-', label="Particle balance model")

ax1.set_title(r'AEH Section (Low iota)')
ax1.set_xlabel(r'$N_{in}$ [s$^{-1}$]')
ax1.set_ylabel(r'$p_{sub}$ [mbar]')

ax1.grid(True, which="major", ls="-", alpha=0.4)
ax1.grid(True, which="minor", ls=":", alpha=0.15)
ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True)

# Plotted exact experimental points as markers and used arrowprops for clean callouts
ax1.plot([3.23e20, 1.15e20], [1.3e-4, 1.8e-5], 'kD', markersize=7, label="Experimental data")
ax1.annotate(".008", (3.23e20, 1.3e-4), textcoords="offset points", xytext=(0,7), ha='right', va='bottom', fontweight='bold', fontsize=12, bbox=bbox_props)
ax1.annotate(".031", (1.15e20, 1.8e-5), textcoords="offset points", xytext=(0,-10), ha='left', va='top', fontweight='bold', fontsize=12, bbox=bbox_props)
#
ax1.legend(loc='upper left')

## SUBPLOT 2: AEP ##
if plot_models["DIVGAS"]:
    ax2.loglog(N_in_grid, p_sub_varout_AEP, color=c_varout, linestyle='-', label="Varouti et al model")
if plot_models["Dirk"]:
    ax2.loglog(N_in_grid, p_sub_dirk_AEP, color=c_dirk, linestyle='-', label="Dirks model")
if plot_models["Throughput balance"]:
    ax2.loglog(N_in_grid, p_sub_conduct_AEP, color=c_conduct, linestyle='-', label="Throughput balance model")
if plot_models["Particle balance"]:
    ax2.loglog(N_in_grid, p_sub_mass_AEP, color=c_mass, linestyle='-', label="Particle balance model")

ax2.set_title(r'AEP Section (High iota)')
ax2.set_xlabel(r'$N_{in}$ [s$^{-1}$]')

ax2.grid(True, which="major", ls="-", alpha=0.4)
ax2.grid(True, which="minor", ls=":", alpha=0.15)
ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True)

# Plotted exact experimental points as markers and used arrowprops for clean callouts
ax2.plot([2.66e19, 2.78e20], [4.1e-5, 1.6e-4], 'kD', markersize=7, label="Experimental data")
ax2.annotate(".008", (2.66e19, 4.1e-5), textcoords="offset points", xytext=(2,7), ha='right', va='bottom', fontweight='bold', fontsize=12, bbox=bbox_props)
ax2.annotate(".031", (2.78e20, 1.6e-4), textcoords="offset points", xytext=(-2,-8), ha='left', va='top', fontweight='bold', fontsize=12, bbox=bbox_props)

plt.tight_layout()


################################
### FLOW FRACTION COMPARISON ###
################################

# We compare the models with Varoutis Fig.11

# Set the parameters 
N_in_AEH = [1e20,1e21,1e22]
N_in_AEP = [1e20,1e21,1e22]

# Initialize the ratio vectors
# Litovoli-Haak Model
N_out_in_litovoli_AEH = np.zeros(np.size(N_in_AEH))
N_out_in_litovoli_AEP = np.zeros(np.size(N_in_AEP))
N_leak_in_litovoli_AEH = []
N_leak_in_litovoli_AEP = []
N_pump_in_litovoli_AEH = np.zeros(np.size(N_in_AEH))
N_pump_in_litovoli_AEP = np.zeros(np.size(N_in_AEP))

# Particle Balance Model
N_out_in_particle_AEH = np.zeros(np.size(N_in_AEH))
N_out_in_particle_AEP = np.zeros(np.size(N_in_AEP))
N_leak_in_particle_AEH = []
N_leak_in_particle_AEP = []
N_pump_in_particle_AEH = np.zeros(np.size(N_in_AEH))
N_pump_in_particle_AEP = np.zeros(np.size(N_in_AEP))

# Run for all the discharges
for i in range(np.size(N_in_AEH)):

    # Litovoli-Haak Model
    results_litovoli = Litovoli_Haak_subdivertor(N_in_AEH[i],N_in_AEP[i])
    N_out_in_litovoli_AEH[i] = results_litovoli[2]
    N_out_in_litovoli_AEP[i] = results_litovoli[3]
    N_leak_in_litovoli_AEH.append(results_litovoli[4])
    N_leak_in_litovoli_AEP.append(results_litovoli[5])
    N_pump_in_litovoli_AEH[i] = results_litovoli[6]
    N_pump_in_litovoli_AEP[i] = results_litovoli[7]
    
    # Particle Balance Model
    results_particle = ParticleBalance_subdivertor(N_in_AEH[i],N_in_AEP[i])
    N_out_in_particle_AEH[i] = results_particle[2]
    N_out_in_particle_AEP[i] = results_particle[3]
    N_leak_in_particle_AEH.append(results_particle[4])
    N_leak_in_particle_AEP.append(results_particle[5])
    N_pump_in_particle_AEH[i] = results_particle[6]
    N_pump_in_particle_AEP[i] = results_particle[7]

# Calculate the fraction sums to see if they add up to 1
k = 0    # Choose for which N_in you plot [0 for 1e20, 1 for 1e21 and 2 for 1e22 ]
frac_sum_litovoli_AEH = N_out_in_litovoli_AEH[k] + np.sum(N_leak_in_litovoli_AEH[k]) + N_pump_in_litovoli_AEH[k]
frac_sum_litovoli_AEP = N_out_in_litovoli_AEP[k] + np.sum(N_leak_in_litovoli_AEP[k]) + N_pump_in_litovoli_AEP[k]
#
frac_sum_particle_AEH = N_out_in_particle_AEH[k] + np.sum(N_leak_in_particle_AEH[k]) + N_pump_in_particle_AEH[k]
frac_sum_particle_AEP = N_out_in_particle_AEP[k] + np.sum(N_leak_in_particle_AEP[k]) + N_pump_in_particle_AEP[k]

## PLOTTING THE BAR CHARTS ##

# Labels for each leakage gap, as Varoutis names them
labels_AEH = ['Outflux','I','D1','D2','D3','H1','H2','H3','H4','H5','C','E','F','G','Pump']
labels_AEP = ['Outflux','J','J1','J2','J3','J4','J5','C','E','F','Pump']

x_AEH = np.arange(len(labels_AEH)) 
x_AEP = np.arange(len(labels_AEP)) 
# Width of the bars
width = 0.25  

# Figures
fig2, ((ax3, ax4), (ax5, ax6)) = plt.subplots(2, 2, figsize=(12, 10))

# Colors
bar_colors = ['#08306B', '#2879B9', '#73B3D8']

for i in range(len(N_in_AEH)):
    # Litovoli-Haak
    values_AEH_litovoli = [N_out_in_litovoli_AEH[i]] + list(N_leak_in_litovoli_AEH[i]) + [N_pump_in_litovoli_AEH[i]]
    values_AEP_litovoli = [N_out_in_litovoli_AEP[i]] + list(N_leak_in_litovoli_AEP[i]) + [N_pump_in_litovoli_AEP[i]]
    # Particle Balance
    values_AEH_particle = [N_out_in_particle_AEH[i]] + list(N_leak_in_particle_AEH[i]) + [N_pump_in_particle_AEH[i]]
    values_AEP_particle = [N_out_in_particle_AEP[i]] + list(N_leak_in_particle_AEP[i]) + [N_pump_in_particle_AEP[i]]

    case_label = f'$N_{{in}}$ = {N_in_AEH[i]:.0e}'
    offset = (i - 1) * width 
    c_bar = bar_colors[i]
    
    # Use distinct x-axes
    ax3.bar(x_AEH + offset, values_AEH_litovoli, width, color=c_bar, linewidth=0.7, label=case_label)
    ax4.bar(x_AEP + offset, values_AEP_litovoli, width, color=c_bar, linewidth=0.7, label=case_label)
    ax5.bar(x_AEH + offset, values_AEH_particle, width, color=c_bar, linewidth=0.7, label=case_label)
    ax6.bar(x_AEP + offset, values_AEP_particle, width, color=c_bar, linewidth=0.7, label=case_label)

# Litovoli-Haak AEH
ax3.set_title('Throughput balance model: AEH')
ax3.set_ylabel(r'Fraction of $N_{in}$')
ax3.set_xticks(x_AEH)
ax3.set_xticklabels(labels_AEH, rotation=45)
ax3.set_ylim(0, 1.1) 
ax3.legend()
ax3.grid(axis='y', alpha=0.3)
ax3.text(0.95, 0.85, f'Sum: {frac_sum_litovoli_AEH:.2f}', transform=ax3.transAxes, ha='right', va='top', fontweight='bold', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Litovoli-Haak AEP
ax4.set_title('Throughput balance model: AEP')
ax4.set_xticks(x_AEP)
ax4.set_xticklabels(labels_AEP, rotation=45)
ax4.set_ylim(0, 1.1)
ax4.legend()
ax4.grid(axis='y', alpha=0.3)
ax4.text(0.95, 0.85, f'Sum: {frac_sum_litovoli_AEP:.2f}', transform=ax4.transAxes, ha='right', va='top', fontweight='bold', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Particle Balance AEH
ax5.set_title('Particle balance model: AEH')
ax5.set_ylabel(r'Fraction of $N_{in}$')
ax5.set_xticks(x_AEH)
ax5.set_xticklabels(labels_AEH, rotation=45)
ax5.set_ylim(0, 1.1) 
ax5.legend()
ax5.grid(axis='y', alpha=0.3)
ax5.text(0.95, 0.85, f'Sum: {frac_sum_particle_AEH:.2f}', transform=ax5.transAxes, ha='right', va='top', fontweight='bold', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Particle Balance AEP
ax6.set_title('Particle balance model: AEP')
ax6.set_xticks(x_AEP)
ax6.set_xticklabels(labels_AEP, rotation=45)
ax6.set_ylim(0, 1.1)
ax6.legend()
ax6.grid(axis='y', alpha=0.3)
ax6.text(0.95, 0.85, f'Sum: {frac_sum_particle_AEP:.2f}', transform=ax6.transAxes, ha='right', va='top', fontweight='bold', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()


# Save figures to .svg
fig1.savefig('Thesis_Figures/Subdiv_models_comparison.svg', format='svg', bbox_inches='tight')
fig2.savefig('Thesis_Figures/Nin_fraction_bar_plots.svg', format='svg', bbox_inches='tight')

plt.show()
