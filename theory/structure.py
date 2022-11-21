import numpy as np

'''
A layer has :
    A thickness measured in mm
    An atomic number Z, a mass number A and a density (kg m^-3) and a mean excitation energy I (J) depending on its material
    and a color
'''

class Layer:
    
    def __init__(self,thickness,Z,A,rho,I,color):
        self.thickness = thickness
        self.Z = Z
        self.A = A
        self.rho = rho
        self.I = I
        self.color = color

    #returns the electron number density which is needed in the bethe-bloch function
    def n(self):

        M_u = 1.660538921*10**(-27)
        return self.Z*self.rho/(self.A*M_u)

    def bethe_bloch(self,beta,z):
        m_e = 9.1093837*10**(-31)
        c = 299792458
        e = 1.60217663*10**(-19)
        eps = 8.85418782*10**(-12)
        dE_dx = 4*np.pi/(m_e*c**2) * self.n()*z**2/(beta**2) * (e**2/(4*np.pi*eps))**2 * (np.log((2*m_e*c**2*beta**2)/(self.I*(1-beta**2)))-beta**2)
        return dE_dx/e


class Aluminum(Layer):
    def __init__(self,thickness):
        Z = 13
        A = 26.981539
        rho = 2698.9
        I = 166*1.60217663*10**(-19)
        color = (163/255,163/255,163/255)
        super().__init__(thickness,Z,A,rho,I,color)
    
    def print(self):
        print("Aluminum ", self.thickness, "mm")

class Lead(Layer):
    def __init__(self,thickness):
        Z = 82
        A = 207.2
        rho = 11350
        I = 823*1.60217663*10**(-19)
        color = (76/255,76/255,76/255)
        super().__init__(thickness,Z,A,rho,I,color)
    
    def print(self):
        print("Lead ", self.thickness, "mm")

class Scintillator(Layer):
    def __init__(self,thickness):
        frac_H = 0.085
        frac_C = 0.915
        Z = frac_H + 6*frac_C
        A = 1.00784*frac_H + 12.011*frac_C
        rho = 1032
        I = 64.7*1.60217663*10**(-19)
        color = (225/255,236/255,242/255)
        super().__init__(thickness,Z,A,rho,I,color)
    def print(self):
        print("Scintillator ", self.thickness, "mm")

class Structure:
    '''
    The order of the layers are given from bottom to top
    '''
    def __init__(self,layers):
        self.layers = layers
        self.bound = np.array([[0,layers[0].thickness]])
        for i in range(1,len(layers)):
            self.bound = np.append(self.bound, [[self.bound[-1,1],self.bound[-1,1]+layers[i].thickness]], axis = 0)
    
    #returns the position of the top of the structure
    def top(self):
        return self.bound[-1,1]

    #takes a structure and stacks it on self
    def add_struct(self,new_structure):
        self.layers = self.layers + new_structure.layers
        self.bound = np.concatenate((self.bound,new_structure.bound+self.top()), axis=0)

    #takes a list of layers and adds them to the structure
    def add_layers(self,new_layers):
        new_structure = Structure(new_layers)
        self.add_struct(new_structure)

    #creates and returns a stack of n initial structures
    def mult(self,n):
        new_structure = Structure(self.layers)
        for i in range(1,n):
            new_structure.add_struct(self)
        return new_structure

    #Retruns the layer in which the particle as position pos_z is
    def material(self,pos_z):
        for i in range(len(self.layers)):
            if pos_z >= self.bound[i,0] and pos_z < self.bound[i,1]:
                return self.layers[i]

    def print(self):
        for i in range(len(self.layers)):
            self.layers[i].print()
            print("--- from ", self.bound[i,0], " to ", self.bound[i,1])
    
