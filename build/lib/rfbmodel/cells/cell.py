import numpy as np
from sympy import symbols, Eq, solve
from dataclasses import dataclass
from ..acl.acl import Acl
from ..agdl.agdl import Agdl
from ..cathodes.cathode import Cathode
from ..membrane.membrane import Membrane

@dataclass
class ReactionProperties:
    E_SCP: float  # Standard Cell Potential of the cell (V)
    delta_S_r: float  # Entropy (J/mol·K)
    n: float  # Electron number transfer in reaction ( - )

class Cell:
    #class which contains the 4 classes which make a cell
    F=96485 # Faraday constant (Coulomb/mol)
    R=8.314 # Ideal Gas Constant (Coulomb/mol)

    def __init__(self, 
                 cathode = None, 
                 membrane = None, 
                 acl = None, 
                 agdl = None, 
                 T = 0,
                 SOC = 0,
                 Q_v = 0, 
                 j_appl = 0, 
                 reaction_props = None):
        
        self.membrane = membrane    # membrane is a class
        self.cathode = cathode      # cathode is a class
        self.acl = acl              # acl is a class
        self.agdl = agdl            # agdl is a class
        self.T = T                  # Temperature of the cell in Kelvin
        self.SOC = SOC              # State of charge [0 - 1]
        self.Q_v = Q_v              # Catholyte flow rate
        self.j_appl = j_appl        # Applied current density
        self.reaction_props = reaction_props        #reaction_props is a class

    
    def E0_cell(self):
        # method for the formal cell potential on temperature
        # SCP: Standard Cell Potential of the two half cell reaction
        E0_cell=self.reaction_props.E_SCP+ self.reaction_props.delta_S_r/(n*Cell.F)*(self.T-298.15)
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

