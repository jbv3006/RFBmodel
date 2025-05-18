class ACL:
    def __init__(self):
        self.thickness = 19/10**6 #Anode catalyst layer thickness [m]
        self.conductivity= 240 #ACL electronic conductivity [S m^-1]
        self.porosity = 0.4 #Anode catalyst layer porosity

    def calculate_resistance(self):
        self.resistance = self.thickness/self.conductivity*(1-self.porosity)**1.5
