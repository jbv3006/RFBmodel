from sympy import fsolve

class Cathode:
    def __init__(self):
        pass

    def conc(self,SOC):
        # Define symbolic variables
        c1, c2, c3, c4, c5, c6  = vars
        #[ VO^+2, VO_2^+, H+, HSO4-, SO4-2, H2SO4]
        c_sulfate= 6.03  # Total sulfate concentration. 6 M or mol/L
        c_vanadium= 1.03 #Total vanadium concentration. M or mol/L
        Q_1= 199.5 #First dissociation quotient of H2SO4 [Adimensional]. Dependance on temperature is negligible
        Q_2= 2.4 #Second dissociation quotient of H2SO4 [Adimensional]. Dependence on temperature is NOT negligble

        # equations
        eq1 = 2*c1+c2+c3-c4-2*c5               # Electroneutrality
        eq2 = SOC-c2/(c1+c2)                   # SOC definition
        eq3 = c_sulfate - c4+c5+c6             # Total sulfate concentration
        eq4 = c_vanadium - c1 - c2             # Total vanadium concentration
        eq5 = Q_1*gamma_1- c3*c4/c6  # First dissociation. Falta término gamma de paper de Minnan
        eq6 =  Q_2*gamma_2- c3*c5/c4 # Second dissociation. Falta término de los gamma

        eqs = [eq1, eq2, eq3, eq4, eq5, eq6]

        # Initial guesses
        initial_guess = [0.5, 0.53, 3, 2.95, 0.08]
        # Solve
        solution = fsolve(equations, initial_guess)
        return solution