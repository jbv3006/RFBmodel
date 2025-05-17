class Agdl:
    def __init__(self):
        self.thickness = 171*10**(-6) #[m] 
        self.conductivity = 8700 #[S m^−1]
        self.porosity = 0.6 #AGDL porosity (void fraction)
        self.resistance = float #AGDL resistance in [Ω·m²]

    def agdl_resistance(self):
        self.resistance = self.thickness/self.resistance


    