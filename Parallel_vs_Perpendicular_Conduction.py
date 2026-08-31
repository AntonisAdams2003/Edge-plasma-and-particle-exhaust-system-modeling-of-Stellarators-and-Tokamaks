import numpy as np
import matplotlib.pyplot as plt



#####################################################################################################################
## In this script we are going to plot the curves where q_parallel_conduction = q_perpaendicular_conduction,       ##
## in a magnetic confined device.                                                                                  ##
##                                                                                                                 ##
## Specifically, we are going to plot hte curves: T^(5/2) = x/(koe Θ^2) * n                                        ##
#####################################################################################################################

# Values for the constans
koe = 2000   # electron spritzer conductivity
koi = 60     # ion spritzer conductivity
kB = 1.38e-23
e = 1.6e-19
x = 3
theta_tok = 0.1    # theta pitch for tokamaks
theta_stel = 1e-3  # theta pitch for stellarators

# Define the particle density vector
N = 200  # points
n = np.linspace(1e19,1e20,N)
n_raised = n**(2/5)




# Tokamak ions
c_tok_ions = (x*e/(koi*theta_tok**2))**(2/5)
T_tok_ions = c_tok_ions*n_raised

# Tokamak electrons
c_tok_elec = (x*e/(koe*theta_tok**2))**(2/5)
T_tok_elec = c_tok_elec*n_raised

# Stellarator ions
c_stel_ions = (x*e/(koi*theta_stel**2))**(2/5)
T_stel_ions = c_stel_ions*n_raised

# Stellarator electrons
c_stel_elec = (x*e/(koe*theta_stel**2))**(2/5)
T_stel_elec = c_stel_elec*n_raised

fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(n, T_tok_ions, color="#830a12", linestyle='dashed')
ax.plot(n, T_tok_elec, color="#830a12", linestyle='solid')
ax.plot(n, T_stel_ions, color="#3b0f70", linestyle='dashed')
ax.plot(n, T_stel_elec, color="#3b0f70", linestyle='solid')

ax.set_xlabel(r'Electron density n$_e$ [m$^{-3}$]',fontsize=16)
ax.set_ylabel(r'Temperature T$_{e,i}$ [eV]',fontsize=16)
ax.set_yscale('log')
ax.set_xlim(n[0],n[-1])
ax.set_ylim(1e0,3.5e2)

x_mid = np.mean(n)
x_lbl = n[int(N / 1.7)]

ax.text(4.0e19, 1.8, 'Tokamak', fontsize=16, bbox=dict(facecolor='white', linewidth=1.5, pad=4), ha='center')
ax.text(4.0e19, 80.0, 'Stellarator', fontsize=16, bbox=dict(facecolor='white', linewidth=1.5, pad=4), ha='center')

ax.text(x_lbl, 1.4, 'electron', fontsize=16)
ax.text(x_lbl, 3.7, 'ion', fontsize=16)
ax.text(x_lbl, 37.0, 'electron', fontsize=16)
ax.text(x_lbl, 156.0, 'ion', fontsize=16)

ax.grid(False)
ax.tick_params(axis='both', which='major', direction='in', length=6, width=1, top=True, right=True,labelsize=14)
ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1, top=True, right=True,labelsize=14)

plt.tight_layout()

# Save the plot
fig.savefig('Thesis_Figures/Conv_Cond_Plot.svg', format='svg', bbox_inches='tight')

plt.show()