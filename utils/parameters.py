'''
Geometric properties of the detector
'''
import numpy as np
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

#Distance on the electronic board between each channel's detector and the central tofpet 
distance_Channels_X = np.nan_to_num(np.ndarray(shape=(8,64), dtype=float), nan=0, posinf=0, neginf=0)*0
distance_Channels_Y = np.nan_to_num(np.ndarray(shape=(8,64), dtype=float), nan=0, posinf=0, neginf=0)*0

#General offset on time computed with time alignement algorithm, for each channel
muX = np.nan_to_num(np.ndarray(shape=(8,64), dtype=float), nan=0, posinf=0, neginf=0)*0
muY = np.nan_to_num(np.ndarray(shape=(8,64), dtype=float), nan=0, posinf=0, neginf=0)*0
