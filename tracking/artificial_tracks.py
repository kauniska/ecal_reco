import matplotlib.pyplot as plt
import numpy as np
from track import Track
from hit import Hit
from track_reconstruction import coord_to_pos,coord_to_pos_x,coord_to_pos_z
import sys
sys.path.insert(1, r"C:\Users\nelg\Desktop\Cours\Labo\TP4\Git\utils")
from parameters import *
from physics import overlap_length

def create_artificial_tracks(t,x0,sampling=1):
    """ Takes a tangent and offset as parameters
        returns a list of hits along the track defined by t and x0
        sampling is between 0 and 1 and corresponds to the probability with which a hit 
        trough which the track goes is present in the final list 
    """
    hits = []
    traj = Track()
    traj.t = t
    traj.x0 = x0
    indices = traj.get_tracks()
    tol = 10**(-10)
    for coord in indices:
        xmin = round(coord[0]-0.5*abs(t))
        xmax = round(coord[0]+0.5*abs(t))
        if abs((coord[0]-0.5*abs(t))%1 - 0.5) < tol and xmin == round(coord[0]-0.5*abs(t)-0.1):
            xmin += 1
        if abs((coord[0]+0.5*abs(t))%1 - 0.5) < tol and xmax == round(coord[0]+0.5*abs(t)+0.1):
            xmax -= 1
        for x in np.arange(max(xmin,1),min(xmax,24)+1,1):
            rand = np.random.uniform()
            if sampling>=rand:
                hits.append(Hit([x,coord[1]],True,0,0,0))

    return hits

def fit_artificial_track(t,x0,sampling=1):
    hits = create_artificial_tracks(t,x0,sampling)
    fit_track = Track(hits)
    return fit_track.t, fit_track.x0


def create_artificial_tracks_geom(t,x0,sidex,sampling=1):
    ''' 
    Takes a tangent and offset as parameters
    returns a list of hits on the plane x-z if sidex, y-z otherwise, along the track defined by t and x0.
    sampling is between 0 and 1 and corresponds to the probability with which a hit 
    trough which the track goes is present in the final list 
    '''

    hits = []

    tol = 10**(-10)

    xmin = min(x0,2*t*n_layers*(thickness+thickness_screen)+x0)
    xmin = max(xmin,0)
    xmax = max(x0,2*t*n_layers*(thickness+thickness_screen)+x0)
    xmax = min(xmax,n_strips*width)

    coord_x_min = round(xmin//width + 1)
    coord_x_max = round(xmax//width + 1)

    z_i = []
    for i in range(n_layers):
        z_center = coord_to_pos_z(i+1,sidex)
        z_i.append([z_center-thickness/2,z_center+thickness/2])

    for j in range(coord_x_min,coord_x_max+1):
        if t == 0:
            z_min = 0
            z_max = 2*n_layers*(thickness+thickness_screen)
        else:
            z_max = max((coord_to_pos_x(j)-width/2-x0)/t,(coord_to_pos_x(j)+width/2-x0)/t)
            z_min = min((coord_to_pos_x(j)-width/2-x0)/t,(coord_to_pos_x(j)+width/2-x0)/t)
        for i,z in enumerate(z_i):
            if overlap_length(z,[z_min,z_max]) > tol:
                rand = np.random.uniform()
                if sampling>=rand:
                    hits.append(Hit([j,i+1],sidex,0,0,0))

    return hits





