'''
A particle has :
    a charge (in multiples of the electron charge)
    a mass in eV/c^2
    an energy in eV
    a 2D position vector (x,z) in mm
'''
class Particle:
    def __init__(self,charge,mass,energy,pos):
        self.charge = charge
        self.mass = mass
        self.energy = energy
        self.pos = pos

    def beta(self):
        return (1-(self.mass/self.energy)**2)**0.5
    
    
    
