import numpy as np

class Layer:
    def __init__(self,thickness,Z,A,rho):
        self.thickness = thickness
        self.top_pos = top_pos
        self.Z = Z
        self.A = A
        self.rho = rho


class Aluminum(Layer):
    def __init__(self,thickness):
        super().__init__(thickness,top_pos,13,26.981539,2698.9)
    
    def print(self):
        print("Aluminum ", thickness, "mm")

class Lead(Layer):
    def __init__(self,thickness):
        super().__init__(thickness,82,207.2,11350)
    
    def print(self):
        print("Lead ", thickness, "mm")

class Scintillator(Layer):
    def __init__(self,thickness):
        super().__init__(thickness,0,0,0)
    def print(self):
        print("Scintillator ", thickness, "mm")

class Structure:
    '''
    The order of the layers are given from bottom to top
    '''
    def __init__(self,layers):
        self.layers = layers
        self.bound = np.array([0,layers[0].thickness])
        for i in range(1,len(layers)):
            np.append(bound, [[bound[-1,1],bound[-1,1]+layers[i].thickness]])
    
    #returns the position of the top of the structure
    def top(self):
        return bound[-1,1]

    #takes a structure and stacks it on self
    def add_struct(self,new_structure):
        self.layers.append(new_structure.layers)
        self.bound = np.concatenate((self.bound,new_structure.bound+self.bound[-1,1]), axis=0)

    #takes a list of layers and adds them to the structure
    def add_layers(self,new_layers):
        new_structure = Structure(new_layers)
        self.add_struct(new_structure)

    #creates and returns a stack of n initial structures
    def mult(self,n):
        new_structure = self
        for i in range(1,n):
            new_structure.add_struct(self)
        return new_structure

    #Retruns the layer in which the particle as position pos_z is
    def material(self,pos_z):
        return None
