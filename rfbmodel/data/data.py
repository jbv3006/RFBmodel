import numpy as np
import matplotlib.pyplot as plt

#CONDUCTIVITY DATA. FOR T [293, 308] K
#Data from "A multi-stack simulation of shunt currents in vanadium redox flow batteries"
#https://doi.org/10.1016/j.jpowsour.2014.03.054
T = np.array([293, 298, 303, 308])  # Temperature [K]
conductivity_V4 = np.array([27.5, 30.8, 33.6, 37])  # Conductivity [S m^-1]

T = np.array([293, 298, 303, 308])  # Temperature [K]
conductivity_V5 = np.array([41.3, 45.4, 49.2, 53])  # Conductivity [S m^-1]


coeffsV4 = np.polyfit(T, conductivity_V4, deg=2)  # Degree 2 polynomial
coeffsV5 = np.polyfit(T, conductivity_V5, deg=2)  # Degree 2 polynomial


def conductivity_funcV4(T):
    conductivity_funcV4 = np.poly1d(coeffsV4)  # Create a callable function
    sigma=conductivity_funcV4(T)
    return sigma

def conductivity_funcV5(T):
    conductivity_funcV5 = np.poly1d(coeffsV5)  # Create a callable function
    sigma=conductivity_funcV5(T)
    return sigma



