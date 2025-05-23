class ACL:
    def __init__(self, thickness, conductivity, porosity):
        self.thickness = thickness #Anode catalyst layer thickness [m]
        self.conductivity = conductivity #ACL electronic conductivity [S m^-1]
        self.porosity = porosity #Anode catalyst layer porosity
        self.s = 0.5 # Liquid water saturation. Dimensionless
        self.s_im = 0.12 #Immobile saturation. Water droplets which dont move or undergo a phase transition
        self.s_r = 0 # Reduced saturation, is the accesible saturation
        self.T = 0  #Temperature of ACL, defined on Cell Class [K]
        self.cvw = 0    #Vapor water concentration at the channel interface [mol m^-3]
        self.RH = 0 #Relative humidity [0-1]
        self.jvw = 0 #ACL water flux [mol / m^2 /s]
        self.Rw = 0 #Source term of water vapor [mol /m^2 /s]

    def calculate_resistance(self):
        self.resistance = self.thickness/self.conductivity*(1-self.porosity)**1.5

    def calculate_s_r(self):
        self.s_r = (self.s-self.s_im)/(1-self.s_im)

    def calculate_cvw(self):
        pass

    def calculate_Rw(self):
        pass

