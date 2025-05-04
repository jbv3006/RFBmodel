import numpy as np
from sympy import symbols, Eq, solve

# === H2-V RFB Parameters ===
params = {
    #Cell Configuration and Operation parameters
    # Geometry (converted to meters)
    "w_ca": 50e-3,     # Cathode width [m]
    "l_ca": 4.60e-3,   # Cathode length [m]
    "h_ca": 50e-3,     # Cathode height [m]

    "w_an": 50e-3,     # Anode width [m]
    "l_an": 0.19e-3,   # Anode length [m]
    "h_an": 50e-3,     # Anode height [m]

    "l_m": 127e-6,     # Membrane thickness [m]
    "l_ACL": 19e-6,    # Anode catalyst layer thickness [m]
    "l_AGDL": 171e-6,  # Anode gas diffusion layer thickness [m]

    # Flow field dimensions
    "W_ch": 0.8e-3,    # Channel width [m]
    "H_ch": 1e-3,      # Channel depth [m]
    "l_ch": 0.56,      # Channel length [m]

    # Operating conditions
    "j_appl_range": (-5250, 5250),  # Current density range [A/m²]
    "E_min": 0.40,      # Lower voltage limit [V]
    "E_max": 1.40,      # Upper voltage limit [V]
    "p_in_an": 103e3,   # Anode inlet pressure [Pa]
    "RH_an": 1.0,       # Relative humidity (100%)

    # Electrolyte composition
    "c_V_0": 1.03,      # Vanadium conc. [mol/L]
    "c_SO4_0": 6.03,    # Sulfate conc. [mol/L]
 
    # Porosities
    "epsilon_ca": 0.94,        # Cathode porosity
    "epsilon_agdl": 0.6,       # Anode gas diffusion layer porosity
    "epsilon_acl": 0.4,        # Anode catalyst layer porosity
    "epsilon_ionoacl": 0.2,    # Ionomer in ACL porosity

    # Fiber diameter (micrometers converted to meters)
    "d_f": 17.6e-6,            # Fiber diameter [m]

    # Kozeny-Carman constant
    "K_KC": 180,               # Kozeny-Carman constant (unitless)

    #Properties

    # Conductivities (S m^-1)
    "sigma_el_bp": 66000,      # Bipolar plate conductivity
    "sigma_el_cc": 6e7,        # Current Collector conductivity
    "sigma_el_cathode": 4400,  # Cathode electronic conductivity
    "sigma_el_agdl": 8700,     # Anode gas diffusion layer conductivity
    "sigma_el_acl": 240,       # Anode catalyst layer conductivity
    "sigma_VOplus2": 44.7,     # VO^{2+} conductivity
    "sigma_VO2plus": 50.5,     # VO_2^{+} conductivity

    # Density (kg m^-3)
    "rho_ca": 1350,            # Cathode density

    # Ideal gas constant
    "R": 8.314,                # J/(mol·K)

    # Ionomer fixed charge concentration
    "c_f": 1.39,               # Ionomer concentration in M

    # Membrane permeability (m^2)
    "K_m": 2e-18,              # Membrane permeability

    # Dry membrane molar volume (m^3 mol^-1)
    "V_m": 5.56e-4,            # Molar volume of dry membrane

    # Diffusion coefficients (m^2 s^-1)
    "D_wv_ref": 1.055e-4,      # Water vapor diffusion coefficient
    "D_H2_ref": 1.055e-4,      # Hydrogen diffusion coefficient

    # Permeabilities (m^2)
    "K_ACL": 1e-13,            # Anode catalyst layer permeability
    "K_AGDL": 1e-12,           # Anode gas diffusion layer permeability

    # Contact angles (degrees)
    "theta_c_acl": 105,        # Contact angle for ACL
    "theta_c_agdl": 110,       # Contact angle for AGDL

    # Bend resistance (unitless)
    "xi_bend": 1.1,            # Bend resistance

    # Entropy change (J mol^-1 K^-1)
    "delta_S_r": -86.9,        # Entropy change for the reaction

    #Electro-osmotic drag coefficient (Dimensionless)
    "n_dH": 1,

    #Faraday Consant (C/mol)
    "F":96485, 

    #Liquid water density (kg/m3)
    "rho_lw":1000

    #Water molar weight (g/mol)
    "M_w":18,

    #Water viscosity (Pa*s)
    "mu_lw":8.9/1000

})

# Surface tension (N m^-1)
def sigma_s(T):
    return -0.00016767 * T + 0.1218

# Catholyte viscosity as a function of temperature T (in K)
def mu_ca(T):
    R = 8.314  # Ideal gas constant in J/(mol·K)
    return 5.4e-6 * np.exp(16131 / (R * T))

#Half cell potential [V]
E_hc_vo2=0.99 
E_hc_h=0
#Dissociation constant sulfuric acid 
#first dissociation h2so4->h+hso4-
K1=1000
#first dissociation hso4-->h+so4-2
K2=0.0106
def reversible_OCP(T,c1,c2,p_h2,c3,a_wca,F_y):
    R = 8.314  # Ideal gas constant in J/(mol·K)
    n=2
    F=96485 # Faraday constant C/mol
    E_cell= (E_hc_vo2-E_hc_h)+delta_S_r/(n*F)*(T-298.15)
    Eocp= E_cell +R*T/(n*F)*np.log((c1*c2*p_h2^0.5)/(c3*a_wca)*F_y)
    return Eocp



#Mass balance equations. (They need to be solved simultaneously)
c1, c2, c3, c4, c5, c6, c7, c8, SOC= symbols('c1 c2 c3 c4 c5 c6 c7 c8 SOC')

eq1 = Eq(c1+c2+2*c3-c7-2*c8, 0) #Electroneutralidad
eq2 = Eq(SOC-c1/(c1+c3), 0)
eq3 = Eq(z**2 + w, 7)
eq4 = Eq(w**2 + v, 3)
eq5 = Eq(x + y + z + w + v, 20)

solution = solve((eq1, eq2, eq3, eq4, eq5), (c1 c2 c3 c4 c5 c6 c7 c8 SOC))



#Resistance
#R= lenght/conductivity
##FALTAN LAS DIMENSIONES: WIDTH OF BIPOLAR PLATE Y DE CURRENT COLLECTOR 

R_el_bp= l_BP / params["sigma_el_bp"]  # Bipolar Plate Resistance FALTA
R_el_cc=l_el_cc/params["sigma_el_cc"]      # Current Collector Resistance FALTA
R_el_cathode=params["w_ca"]/params["sigma_el_cathode*(1-epsilon_ca)"]*(1-params["epsilon_ca"])^1.5 # Cathode electronic Resistance
R_el_agdl= params["l_AGDL"] / params["sigma_el_agdl"]**(1-params["epsilon_agdl"])^1.5  # Anode gas diffusion layer Resistance
R_el_acl= params["l_ACL"] /  params["sigma_el_acl"]*(1-params["epsilon_acl"])^1.5   # Anode catalyst layer Resistance
R_io_catholyte= params["w_ca"]/sigma_io_catholyte(SOC)*params["epsilon_ca"]^1.5
R_extra=0 #Should be a fitting parameter

#FALTA R_ionomer y sigma_ionomer
R_HFR=(1/R_el_cathode+1/R_io_catholyte)^-1+R_el_acl+R_el_agdl+R_el_bp+R_el_cc+R_extra



#Water transport
J_iw=D_iwm*(c_ca_iw-c_an_iw)/params["l_m"]-params["n_dH"]*params["j_appl_range"]/params["F"]+params["rho_lw"]*params["K_m"]/params["M_w"]/params["mu_lw"]*(params[])