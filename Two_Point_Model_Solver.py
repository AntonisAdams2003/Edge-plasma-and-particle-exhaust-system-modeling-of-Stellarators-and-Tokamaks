import numpy as np
from numpy import sqrt , pi
from scipy.optimize import fsolve

#######################################################################################################################
## This script solves the Two-Point Model for:                                                                       ##
## - Tokamak -> Stangeby's Model  call: TTPM(n_upst , Lc,x,q_SOL)  or    TTPM_numerical(n_upst, Lc, x, q_SOL)        ##
## - Stellarator -> Feng's Model  call: Extended_STPM(n_upst,loss_alpha,loss_fcool,loss_fconv , Lc,x,theta,q_SOL)    ##
#######################################################################################################################


## EXTENDED STELLARATOR TWO POINT MODEL ##
def Extended_STPM(n_upst,loss_alpha,loss_fcool,loss_fconv , Lc,x,theta,q_SOL):
    
    # Iterational process parameters
    tolerance = 1e-5   # define the tolerance
    max_iter = 100   # define the maximum number of iterations allowed

    # Define the known constants and geometric dimentions
    m_H = 1.67262e-27  # [kg] mass hydrogen ion
    m_H2 = 3.347e-27  # [kg] mass of H2  (=2 x m_H)
    e = 1.60218e-19    # electron charge
    kB = 1.380649*1e-23   # [J/K] Boltzmann's constant

    # Define constants by assuming typical values
    gamma = 8   # Sheath transmission coefficient
    koe = 2000  # [W/m eV^(7/2)] Spitzer-Härm heat conductivity for electrons

    # Parallel heat flux
    q_par = q_SOL*1e6  # [W/m^2]       #/(4.0*pi**2*a_r*R*theta)
    # Define some helpfull constants
    A = q_par**2*m_H/(gamma**2*e**3)
    B = 7*q_par*Lc/(2.0*koe)
    C = 7*e*x/(4*koe*theta**2)
 
    # Define the equation that we need to solve with respect to n_target 
    def func_val(n_targ,n_upst,A,B,C,loss_alpha,loss_fcool,loss_fconv):
        func = (
            ((4*A*(n_targ/n_upst**3)*(1-loss_fcool)**2)**(1.0/3)* 
            (1 + loss_alpha*(0.5*A*(1/n_targ**2)*(1-loss_fcool)**2)**(-1.0/6)))**(7.0/2) 
            -(0.5*A*(1/n_targ**2)*(1-loss_fcool)**2)**(7.0/6) - B*(1-loss_fconv) + C*(n_upst+n_targ)* 
            ((4*A*(n_targ/n_upst**3)*(1-loss_fcool)**2)**(1.0/3) * 
            (1 + loss_alpha*(0.5*A*(1/n_targ**2)*(1-loss_fcool)**2)**(-1.0/6)) - (0.5*A*(1/n_targ**2)*(1-loss_fcool)**2)**(1.0/3))
            )
        return func
    

    ####### SOLUTION #######
    # Calculating the target particle density
    nt0 = n_upst  # Initiall guess for the algorithm
    n_targ = fsolve(func_val,nt0,(n_upst,A,B,C,loss_alpha,loss_fcool,loss_fconv))
    n_targ = n_targ[0]

    # Upstream Temperature
    T_upst = (4*A*(n_targ/n_upst**3)*(1-loss_fcool)**2)**(1.0/3) * (1 + loss_alpha*(0.5*A*(1/n_targ**2)
        *(1-loss_fcool)**2)**(-1.0/6))   # [eV]

    # Target Temperature
    T_targ = (0.5*A*(1/n_targ**2)*(1-loss_fcool)**2)**(1.0/3)   # [eV]

    #-OUTPUT: [0] n_targ, [1] T_upst, [2] T_targ

    return (n_targ,T_upst,T_targ)




## STANGEBY'S TOKAMAK TPM (ANALYTICAL SOLUTION) ##
def TTPM(n_upst , Lc,x,q_SOL):

    # Define the known constants and geometric dimentions
    m_H = 1.67262e-27  # [kg] mass hydrogen ion
    m_H2 = 2*1.67262e-27  # [kg] mass of H2
    e = 1.60218e-19    # electron charge
    kB = 1.380649*1e-23   # [J/K] Boltzmann's constant

    # Define constants by assuming typical values
    gamma = 8   # Sheath transmission coefficient
    koe = 2000  # [W/m eV^7/2] Spitzer-Härm heat conductivity for electrons

    # The parallel heat flux
    q_par = q_SOL*1e6  # [W/m^2]
    

    ####### SOLUTION #######
    
    ## UPSTREAM ##
    # Upstream Temperature
    T_upst = (7/2*q_par*Lc/koe)**(2/7)   # [eV]

    ## TARGET PARAMETERS ##
    # Target density
    c_n = (7/2*q_par*Lc/koe)**(6/7) * gamma**2*e**3/(4*m_H*q_par**2)
    n_targ = c_n * n_upst**3 
    # Target Temperature
    T_targ = m_H/(2*e) * 4*q_par**2/(gamma*e*n_upst*T_upst)**2   # [eV]

    #-OUTPUT: [0] n_targ, [1] T_upst, [2] T_targ

    return (n_targ,T_upst,T_targ)




## STANGEBY'S TOKAMAK TPM (NUMERICAL SOLUTION) ##
def TTPM_numerical(n_upst, Lc, x, q_SOL):

    # Define the known constants and geometric dimentions
    m_H = 1.67262e-27  # [kg] mass hydrogen ion
    e = 1.60218e-19    # electron charge

    # Define constants by assuming typical values
    gamma = 8   # Sheath transmission coefficient
    koe = 2000  # [W/m eV^7/2] Spitzer-Härm heat conductivity for electrons

    # The parallel heat flux
    q_par = q_SOL*1e6  # [W/m^2]
    
    # Define helpfull constants 
    A = q_par**2 * m_H / (gamma**2 * e**3)
    B = 7 * q_par * Lc / (2.0 * koe)

    ####### SOLUTION #######

    def func_val(n_targ, n_upst, A, B):
        # Algebraically derived from the 3 TTPM equations
        func = (
            ((4 * A * (n_targ / n_upst**3))**(1.0/3))**(7.0/2) -(0.5 * A * (1 / n_targ**2))**(7.0/6) - B
        )
        return func

    # Initiall guess for the algorithm
    nt0 = n_upst  
    n_targ_arr = fsolve(func_val, nt0, args=(n_upst, A, B))
    n_targ = n_targ_arr[0]

    ## UPSTREAM PARAMETERS ##
    T_upst = (4 * A * (n_targ / n_upst**3))**(1.0/3)

    ## TARGET PARAMETERS ##
    T_targ = (0.5 * A * (1 / n_targ**2))**(1.0/3)

    #-OUTPUT: [0] n_targ, [1] T_upst, [2] T_targ
    return (n_targ, T_upst, T_targ)


