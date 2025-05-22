from scipy.optimize import fsolve
import numpy as np
from rfbmodel.data import conductivity_funcV5, conductivity_funcV4

class Cathode:
    def __init__(self, thickness, cath_conductivity, porosity):
        # Initial guesses for concentrations mol/L
        # [VO^+2, VO_2^+, H+, HSO4-, SO4-2, H2SO4]
        self.c = np.array([1,      1,      1,  1,     2,     0.01])

        # Input parameters
        self.thickness = thickness #Cathode thickness [m]
        self.cath_conductivity = cath_conductivity #Conductivity of the cathode / sigma [S m^-1]
        self.porosity = porosity #Cathode porosity (void fraction)

        # Cell parameters
        self.SOC = 0          # State of charge
        self.T = 298 #Temperature of the cell [K]
        
        # Cathode default parameters
        self.c_sulfate = 6.03   # Total sulfate concentration. 6 M or mol/L
        self.c_vanadium = 1.03  # Total vanadium concentration. M or mol/L
        self.Q_1 = 199.5        # First dissociation quotient of H2SO4 [Adimensional]. Dependence on temperature is negligible
        self.Q_2 = 2.4          # Second dissociation quotient of H2SO4 [Adimensional]. Dependence on temperature is NOT negligible
        self.gamma_1 = 1        # Activity coefficient of tank in the first dissociation (?) NOT SURE
        self.gamma_2 = 1
        self.io_conductivity = 0 #Conductivity of the catholyte / sigma [S m^-1]
        self.resistance = 0 #Resistance of the catholyte / 

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
        
        # Construct the system as an array
        return np.array([eq1, eq2, eq3, eq4, eq5, eq6])
           
    def conc(self):
        ''' Solves cathode concentration'''
        c = fsolve(self.conc_sys, x0 = self.c)
        # Set the equilibrium concentration in the cathode
        self.c = c
        
    def ionic_resistance(self):
        ''' Solves catholythe resistance considering void fraction'''
        # Set the conducttivity in the catholyte
        self.io_conductivity= self.SOC*conductivity_funcV5(self.T)+(1-self.SOC)*conductivity_funcV4(self.T)
        # Set the resistance in the catholyte
        self.io_resistance=self.thickness/(self.io_conductivity*self.porosity**1.5)

    def cathode_resistance(self):
        ''' Solves cathode electronic conductivity considering void fraction'''
        # Set the resistance in the cathode
        self.cath_resistance=self.thickness/(self.cath_conductivity*(1-self.porosity)**1.5)

    def total_resistance(self):
        ''' Solves electronic and ionic resistance of the cathode and catholyte'''
        # Set the total resistance in the cathode and catholyte
        self.resistance=(1/self.cath_resistance+1/self.io_resistance)**(-1)
