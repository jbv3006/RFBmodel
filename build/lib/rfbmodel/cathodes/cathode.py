from scipy.optimize import fsolve
import numpy as np


class Cathode:
    def __init__(self):

        # Initial guesses for concentrations mol/L
                         # [VO^+2, VO_2^+, H+, HSO4-, SO4-2, H2SO4]
        self.c = np.array([1,      1,      1,  1,     2,     0.01])
        
        # Default parameters
        self.SOC = 0.5          # State of charge
        self.c_sulfate = 6.03   # Total sulfate concentration. 6 M or mol/L
        self.c_vanadium = 1.03  # Total vanadium concentration. M or mol/L
        self.Q_1 = 199.5        # First dissociation quotient of H2SO4 [Adimensional]. Dependence on temperature is negligible
        self.Q_2 = 2.4          # Second dissociation quotient of H2SO4 [Adimensional]. Dependence on temperature is NOT negligible
        self.gamma_1 = 1
        self.gamma_2 = 1
        self.thickness=50/1000 #Cathode thickness [m]
        self.conductivity= 0 #Conductivity of the catholyte / sigma [S m^-1]
        self.resistance=0 #Resistance of the catholyte / m^2/S
        self.T=298 #Temperature of the cell [K]
        print("holaa")
        pass

    def conc_sys(self, x):
        """Returns the residuals eq1…eq6 for the vector x = [c1…c6].
        Name of components
        [c1,    c2,     c3, c4,    c5,    c6]
        [VO^+2, VO_2^+, H+, HSO4-, SO4-2, H2SO4]
        """

        c1, c2, c3, c4, c5, c6 = x

        # Define system of equations
        eq1 = 2*c1+c2+c3-c4-2*c5               # Electroneutrality
        eq2 = self.SOC-c2/(c1+c2)                   # SOC definition
        eq3 = self.c_sulfate - c4+c5+c6             # Total sulfate concentration
        eq4 = self.c_vanadium - c1 - c2             # Total vanadium concentration
        eq5 = self.Q_1 * self.gamma_1 - c3* (c4/c6)  # First dissociation. Falta término gamma de paper de Minnan
        eq6 =  self.Q_2 * self.gamma_2 - c3* (c5/c4) # Second dissociation. Falta término de los gamma
        print("hola hola lele")
        
        # Construct the system as an array
        return np.array([eq1, eq2, eq3, eq4, eq5, eq6])
           
    def conc(self):
        ''' Solves cathode concentration'''
        c = fsolve(self.conc_sys, x0 = self.c)
        # Set the equilibrium concentration in the cathode
        print("lol")
        self.c = c
        

    def cath_resistance(self):
        ''' Solves cathode resistance'''
        # Set the conducttivity in the catholyte
        print("loool")
        self.conductivity= self.SOC*conductivity_funcV5(self.T)+(1-self.SOC)*conductivity_funcV4(self.T)
        self.resistance=self.thickness/self.resistance