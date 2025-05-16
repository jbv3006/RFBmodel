from dataclasses import dataclass
@dataclass

class Acl:
    def __init__(self, length, electronic_conductivity):
        self.length = length
        self.conductivity=electronic_conductivity

    def resistance(self):
        return self.length/self.conductivity