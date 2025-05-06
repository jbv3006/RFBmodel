from dataclasses import dataclass
from sympy import symbols, Eq, solve
@dataclass
class CatholyteProperties:
    K1: float  # Equilibrium constant first dissociation of acid
    K2: float  # Equilibrium constant first dissociation of acid

class Cathode:
    def __init__(self, catholyte_props: CatholyteProperties):
        self.catholyte_props=catholyte_props

    def conc(self,SOC):
        # Define symbolic variables
        c1, c2, c3, c4, c5, c6, c7, c8 = symbols('c1 c2 c3 c4 c5 c6 c7 c8')
        #[VO_2^+, H+, VO^+2, H20, H2, H2SO4, HSO4-, SO4-2, SOC]

        # equations
        eq1 = Eq(c1 + c2 + 2*c3 - c7 - 2*c8, 0)               # electroneutrality
        eq2 = Eq(SOC - c1 / (c1 + c3), 0)                    # SOC definition
        eq3 = Eq(c2 * c1 - self.catholyte_props.K1 * c4, 0)  # First dissociation
        eq4 = Eq(c3 * c1 - self.catholyte_props.K2 * c2, 0)  # Second dissociation
        eq5 = Eq()               # 
        eq6 = Eq()                    # 
        eq7 = Eq()  # 
        eq8 = Eq()  # 
        eq9 = Eq()  # 
        eqs = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9]
        solution = solve(equations, (c1, c2, c3, c4, c5, c6, c7, c8), dict=True)
        return solution