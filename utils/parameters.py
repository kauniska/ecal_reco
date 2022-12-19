'''
Geometric properties of the detector
'''

n_strips = 24
n_layers = 8

#Dimensions of the cells expressed in cm
"""
# without the correct geometry implemented
width = 1
thickness = 1
thickness_screen = 0
"""
width = 1.6
thickness = 0.4
thickness_screen = 0.6

#Height of the detector
total_height = 2*n_layers*(thickness+thickness_screen)

#Distance between 2 active layer of the same orientation
Delta_z = 2*(thickness+thickness_screen)




