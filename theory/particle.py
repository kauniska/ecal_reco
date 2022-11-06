'''
A particle has :
    a charge (in multiples of the electron charge)
    a mass in eV/c^2
    an energy in eV
    a 2D position vector (x,z) in mm
    an angle of propagation :
        theta = 0 implies along the x axis
        t = pi/2 implies along the z axis
'''
class Particle:
    def __init__(self,charge,mass,energy,pos,theta):
        self.charge = charge
        self.mass = mass
        self.energy = energy
        self.pos = pos
        self.theta = theta

    def beta(self):
        return (1-(self.mass/self.energy)**2)**0.5
    
    def print(self):
        print("Charge = ", self.charge, "e")
        print("Mass = ", self.mass, "eV")
        print("Energy = ", self.energy, "eV")
        print("Position = ", self.pos, "mm")
        print("Theta = ", self.theta, "rad")
    

    
