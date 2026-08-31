import numpy as np
from numpy import pi, sqrt, exp, cos, sin
import math
from scipy.optimize import fsolve


###################################################################################
## This script takes the outputs of the Two-point model (target conditions) and 
## uses 
## balance




# Define the known constants and geometric dimentions
m_H = 1.67262e-27  # [kg] mass hydrogen ion
m_H2 = 3.347e-27  # [kg] mass of H2  (= 2 x m_H)
e = 1.60218e-19    # electron charge
kB = 1.380649*1e-23   # [J/K] Boltzmann's constant
a_r = 0.5   # [m] minor radius
R = 5.5   # [m] major radius
a = 1  # angle of incidence
d = 0.4    # [m]  diameter of vacuum cylinders cross-sections
A_pg_AEH = 0.153   # [m^2] punping gap area for AEH (low div)
A_pg_AEP = 0.045   # [m^2] punping gap area for AEP (high div)
A_leak_AEH = 0.03  # [m^2] leakage area for low divertor
A_leak_AEP = 0.022 # [m^2] leakage area for high divertor

# Define constants by assuming typical values
gamma = 8   # Sheath transmission coefficient
koe = 2000  # [W/m eV^7/2] Spitzer-Härm heat conductivity for electrons
koi = 60    # [W/m eV^7/2] Spitzer-Härm heat conductivity for ions
Tvv = 303   # [K] Temperature of the vacuum vessel
T_leak = 300  # [K] Temperature of the leakages
T_div = 300  # [K] Temperature at the Divertor
T_sub = 300  # [K] Temperature at the Sub-divertor
T_0 = 600    # [K] Temperature of the H that enters
ji_AEH = 0.06  # probability of particle pumped out for AEH
ji_AEP = 0.0264  # probability of particle pumped out for AEP

# Effective pumping speed 
S_eff_AEH = (1.0/4)*ji_AEH*(pi*d**2/4)*sqrt(8*kB*Tvv/(pi*m_H2))   # [m^3/s] V.Haak -> 2.35 AEH
S_eff_AEP = (1.0/4)*ji_AEP*(pi*d**2/4)*sqrt(8*kB*Tvv/(pi*m_H2))   # [m^3/s] V.Haak -> 1.18 AEP  

# Imposed bulk flow velocity
s = 0.5/sqrt(pi)
th = 0.0
B = exp(-s**2*cos(th)**2) + sqrt(pi)*s*cos(th)*(1 + math.erf(s*cos(th)))




def Divertor_Subdivertor_Analytical(n_targ,T_targ , theta):

     # Calculates divertor pressure
    def get_p_div(epsilon, a, b, c):
        return epsilon * Gamma_recycl * A_wet * (c / (a * (c - b)))
    
    # Calculates subdivertor pressure
    def get_p_sub(epsilon, b, c):
        return epsilon * Gamma_recycl * A_wet * (1 / (b - c))

    # Calulate the mean thermal speeds for two different temperatures
    vth_0 = sqrt(8*kB*T_0/(pi*m_H2))   # T0
    vth_leak = sqrt(8*kB*T_sub/(pi*m_H2))   # Tsub = T_leak

    # PSI wetting area
    A_wet = 1.5
    # Particle capture efficiency
    epsilon_AEH = 0.03
    epsilon_AEP = 0.03
    # Particle Flux recycling reaching the targets
    Gamma_targ =  n_targ * sqrt(2*e*T_targ/m_H)   # [part/m^2 s]

    # Recycling coefficient
    R = 1
    # Particle Flux recycling from the targets
    Gamma_recycl = R/2 * sin(math.radians(a)) * Gamma_targ   # [part/m^2 s]


    # Form the Divertor-Subdivertor system

    # AEH
    a_AEH = B/4*A_pg_AEH*vth_0/(kB*T_0)
    #
    term_leak_AEH = A_leak_AEH * vth_leak / (kB * T_sub)
    term_pg_AEH = A_pg_AEH * vth_0 / (kB * T_0)
    b_AEH = -epsilon_AEH * B/4 * (term_leak_AEH + term_pg_AEH)
    #
    c_AEH = -B/4 * (term_leak_AEH + term_pg_AEH) - S_eff_AEH/(kB*T_sub)

    # Divertor pressure
    p_div_AEH = get_p_div(epsilon_AEH, a_AEH, b_AEH, c_AEH)  # [Pa]
    # Subdivertor pressure
    p_sub_AEH = get_p_sub(epsilon_AEH, b_AEH, c_AEH)  # [Pa]
    


    # AEP
    a_AEP = B/4*A_pg_AEP*vth_0/(kB*T_0)
    #
    term_leak_AEP = A_leak_AEP * vth_leak / (kB * T_sub)
    term_pg_AEP = A_pg_AEP * vth_0 / (kB * T_0)
    b_AEP = -epsilon_AEP * B/4 * (term_leak_AEP + term_pg_AEP)
    #
    c_AEP = -B/4 * (term_leak_AEP + term_pg_AEP) - S_eff_AEP/(kB*T_sub)

    # Divertor pressure
    p_div_AEP = get_p_div(epsilon_AEP, a_AEP, b_AEP, c_AEP)  # [Pa]
    # Subdivertor pressure
    p_sub_AEP = get_p_sub(epsilon_AEP, b_AEP, c_AEP)  # [Pa]



    # Particle flow throught the pumping gap area
    # AEH
    N_in_AEH = B/4 * A_pg_AEH*vth_0*p_div_AEH / (kB*T_0)  # [particles/s]
    # AEP
    N_in_AEP = B/4 * A_pg_AEP*vth_0*p_div_AEP / (kB*T_0)  # [particles/s]

    # Conver pressures in [mbar]
    # Divertor pressure
    p_div_AEH = 0.01 * p_div_AEH  # [mbar]
    p_div_AEP = 0.01 * p_div_AEP  # [mbar]
    # Subdivertor pressure
    p_sub_AEH = 0.01 * p_sub_AEH  # [mbar]
    p_sub_AEP = 0.01 * p_sub_AEP  # [mbar]


    #-OUTPUT: [0] Gamma_targ, [1] N_in_AEH, [2] N_in_AEP, [3] p_divertor_AEH, [4] p_divertor_AEP, [5] p_subdivertor_AEH, [6] p_subdivertor_AEP

    return (Gamma_targ, N_in_AEH,N_in_AEP, p_div_AEH, p_div_AEP, p_sub_AEH, p_sub_AEP)






### Code that assumes Particle Collection Efficiency  Γ_in/Γ_t = 5% and uses the Dirk's model ###

def Divertor_Subdivertor_Dirk(n_targ,T_targ):

    ## DIVERTOR PARAMETERS ##

    # Attached conditions
    # n_0 = n_targ*sqrt(T_targ/T_0 * e/kB)*sqrt(m_H2/m_H)*sqrt(2*pi)  # neutral gas density for attached case
    # Gamma_0 = (1.0/4)*n_0*sqrt(8*kB*T_0/(pi*m_H2))    # [particles/m^2/s]


    # Particle Collection Efficiency  =(Γ_in/Γ_t)
    PCE = 0.05
    # Particle flux reaching the targets
    Gamma_targ = n_targ * sqrt(2*e*T_targ/m_H)  # [particles/m^2/s]
    # Particle flux that passes through the gaps
    Gamma_in = Gamma_targ * PCE        # [particles/m^2 s]

    # Particle flow that passes through each gap
    N_in_AEH = A_pg_AEH * Gamma_in   # [particles/s]
    N_in_AEP = A_pg_AEP * Gamma_in   # [particles/s]
    

    #print(Gamma_targ,Gamma_in_AEH)

    # Pressure in the divertor
    h_AEH = sqrt(pi*m_H2*kB*T_div/2.0) / (A_leak_AEH + A_pg_AEH*S_eff_AEH/(A_pg_AEH*sqrt(kB*T_sub/(2.0*pi*m_H2)) + S_eff_AEH))
    h_AEP = sqrt(pi*m_H2*kB*T_div/2.0) / (A_leak_AEP + A_pg_AEP*S_eff_AEP/(A_pg_AEP*sqrt(kB*T_sub/(2.0*pi*m_H2)) + S_eff_AEP))
    p_divertor_AEH = (h_AEH * N_in_AEH) * 0.01  # [mbar]
    p_divertor_AEP = (h_AEP * N_in_AEP) * 0.01  # [mbar]
    
    ## SUBDIVERTOR PARAMETERS ##
    # Reduce constants of pressure as the gas flows to the sub-divertor
    w_AEH = A_pg_AEH*sqrt(kB*T_div) / (sqrt(2*pi*m_H2)*S_eff_AEH + A_leak_AEH*sqrt(kB*T_sub))
    w_AEP = A_pg_AEP*sqrt(kB*T_div) / (sqrt(2*pi*m_H2)*S_eff_AEP + A_leak_AEP*sqrt(kB*T_sub))
    
    # Estimating the sub-divertor's pressure 
    p_subdivertor_AEH = p_divertor_AEH * w_AEH   # [mbar]
    p_subdivertor_AEP = p_divertor_AEP * w_AEP   # [mbar]
    

    #-OUTPUT: [0] Gamma_targ, [1] N_divertor_AEH, [2] N_divertor_AEP, [3] p_divertor_AEH, [4] p_divertor_AEP, [5] p_subdivertor_AEH, [6] p_subdivertor_AEP

    return (Gamma_targ, N_in_AEH,N_in_AEP, p_divertor_AEH,p_divertor_AEP, p_subdivertor_AEH,p_subdivertor_AEP)










### Code that assumes Particle Collection Efficiency  Γ_in/Γ_t = 5% and solves for the subdivertor pressure ###



# def Divertor_Subdivertor_Conductance(n_targ,T_targ):

#     # Particle Collection Efficiency  =(Γ_in/Γ_t)
#     PCE = 0.05
#     # Particle flux through 
#     Gamma_targ = n_targ * sqrt(2*e*T_targ/m_H)  # [particles/m^2/s]
#     # Particle flux that passes through the gaps
#     Gamma_in = Gamma_targ * PCE        # [particles/m^2 s]

#     # Particle flow that passes through each gap
#     N_in_AEH = A_pg_AEH * Gamma_in   # [particles/s]
#     N_in_AEP = A_pg_AEP * Gamma_in   # [particles/s]


#     # Calculate the throughputs
#     q_in_1 = N_in_AEH * T_0*kB  # AEH
#     q_in_2 = N_in_AEP * T_0*kB  # AEP

#     # Define the model
#     def System_of_Equations(p , q_in_1,q_in_2):

#         p1,p2 = p

#         ## LEAKAGE THROUGHPUT ##
#         # Leakage Throughputs
#         q_leak_1 = B * 0.25 * A_leak_AEH * p1 * sqrt(8*kB*T_leak/(pi*m_H2))   # AEH
#         q_leak_2 = B * 0.25 * A_leak_AEP * p2 * sqrt(8*kB*T_leak/(pi*m_H2))   # AEP

#         ## PUMPING THROUGHPUT ##
#         # Effective pumping speed 
#         S_eff_AEH = (1.0/4)*ji_AEH*(pi*d**2/4)*sqrt(8*kB*Tvv/(pi*m_H2))   # [m^3/s] V.Haak -> 2.35 AEH
#         S_eff_AEP = (1.0/4)*ji_AEP*(pi*d**2/4)*sqrt(8*kB*Tvv/(pi*m_H2))   # [m^3/s] V.Haak -> 1.18 AEP  
        
#         # Turbomolecular pumping throughput
#         q_pump_1 = S_eff_AEH * p1  # AEH
#         q_pump_2 = S_eff_AEP * p2  # AEP

#         ## OUTFLUX THROUGHPUT ##
#         q_outflux_1 = B * 0.25 * A_pg_AEH * p1 * sqrt(8*kB*T_0/(pi*m_H2))   # AEH
#         q_outflux_2 = B * 0.25 * A_pg_AEP * p2 * sqrt(8*kB*T_0/(pi*m_H2))   # AEP

#         ## CONDUCTANCE for each CHANELL ##
#         C12 = 0.0

#         ## MODELs EQUATIONS: ##
#         eq1 = (q_in_1 - q_outflux_1) - q_leak_1 - C12*(p1-p2) - q_pump_1
#         eq2 = (q_in_2 - q_outflux_2) - q_leak_2 + C12*(p1-p2) - q_pump_2

#         return[eq1,eq2]


#     # Initial guess for solving the system
#     p_initial = np.random.uniform(1e-5 ,1e0, 2)

#     # Call the solver to solve the system
#     p_optimal = fsolve(System_of_Equations,p_initial,args=(q_in_1, q_in_2))

#     # Extract the values and convert to mbar (Assuming solver outputs Pascals: 1 Pa = 0.01 mbar)
#     p_AEH = p_optimal[0] * 0.01   # [mbar]
#     p_AEP = p_optimal[1] * 0.01   # [mbar]

#     # Format as scientific notation directly in the print statement
#     print(f"Subdivertor pressure at AEH: {p_AEH:.2e} [mbar]")
#     print(f"Subdivertor pressure at AEP: {p_AEP:.2e} [mbar]")

#     return[p_AEP,p_AEH]
