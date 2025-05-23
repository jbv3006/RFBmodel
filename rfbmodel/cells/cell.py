import numpy as np
from ..acl.acl import ACL
from ..agdl.agdl import AGDL
from ..cathodes.cathode import Cathode
from ..membrane.membrane import Membrane

class Cell:
    #class which contains the 4 classes which make a cell

    F=96485 # Faraday constant (Coulomb/mol)
    R=8.314 # Ideal Gas Constant (Coulomb/mol)

    def __init__(self, 
                 cathode: Cathode, 
                 membrane: Membrane, 
                 acl: ACL, 
                 agdl: AGDL, 
                 T = 0,
                 SOC = 0,
                 Q_v = 0, 
                 j_appl = 0, 
                 reaction_props = dict):
        
        #Cathode initialization
        self.cathode = cathode     # cathode is a class
        #Cell input actualization on cathode
        self.cathode.SOC = SOC
        self.cathode.T=T
        #Cathode functions
        self.cathode.conc()
        self.cathode.ionic_resistance()
        self.cathode.cathode_resistance()
        self.cathode.total_resistance()

        #Membrane initialization
        self.membrane = membrane    # membrane is a class

        #ACL initialization
        self.acl = acl            # acl is a class
        #Cell input actualization on ACL
        self.acl.T = T

        #AGDL initialization
        self.agdl = agdl         # agdl is a class
        #Cell input actualization on AGDL
        self.agdl.T=T



        #Cell operation parameters
        self.T = T                  # Temperature of the cell in Kelvin
        self.SOC = SOC              # State of charge [0 - 1]
        self.Q_v = Q_v              # Catholyte flow rate
        self.j_appl = j_appl        # Applied current density
        self.reaction_props = reaction_props      #reaction_props is a class
        

    
    def E0_cell(self):
        # method for the formal cell potential on temperature
        # SCP: Standard Cell Potential of the two half cell reaction
        n=2
        E0_cell=self.reaction_props['E_SCP'] + self.reaction_props['delta_S_r']/(n*Cell.F)*(self.T-298.15)
        return E0_cell

    def E_OCP(self):
        #method for the open circuit potential of the cell, CON
        # E_OCP=self.E0_cell()+Cell.R*self.T/(self.reaction_props.n*Cell.F)*np.log(##Hay que poner las concentraciones !! en el catodo)
        E_OCP = 1
        return E_OCP

    def E_cell(self):
        """ Calculates total cell voltage based on: Open circuit potential, Membrane ohmic loss, Cathode overpotential, 
        Anode overpotential, Protonic transport in ACL """
        E_OCP=self.E_OCP()
        ohm_resistance=self.membrane.ohm_resistance()
        ca_overpotential=self.cathode.overpotential()
        an_overpotential=self.acl.overpotential()
        H_plus_overpotential=self.acl.Hplus()
        E_cell=E_OCP-ohm_resistance+ca_overpotential-an_overpotential-H_plus_overpotential
        return E_cell
    
    def R_HFR(self):
        """Calculates high frequency resistance, which accounts for the electronic and ionic resistance of the cathode, catholyte, ACL, AGDL, 
         membrane, current collector and bipolar plates"""
        r_hfr=self.cathode.resistance + self.acl.resistance + self.agdl.resistance
        return r_hfr
    
    def calculate_jvw_AGDL(self):
        D_eff = 1.055 * 10**(-4) #Difussion coefficient of water vapor and hydrogen [m2/s]
        self.agdl.calculate_cvw()
        self.acl.calculate_cvw()
        self.agdl.jvw = D_eff * (self.acl.cvw - self.agdl.cvw)/(self.agdl.thickness + self.acl.thickness) + self.agdl.Rw

    def calculate_jvw_ACL(self):
        D_eff = 1.055 * 10**(-4) #Difussion coefficient of water vapor and hydrogen [m2/s]
        self.acl.jvw = D_eff * (self.acl.cvw - self.agdl.cvw)/(self.agdl.thickness + self.acl.thickness) + self.acl.Rw



