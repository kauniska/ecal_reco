import matplotlib.pyplot as plt
import numpy as np
from hit import Hit,width,thickness,thickness_screen
from track import Track

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

def create_artificial_tracks_geom(t,x0,sampling=1):
    """ Takes a tangent and offset as parameters
        returns a list of hits along the track defined by t and x0
        sampling is between 0 and 1 and corresponds to the probability with which a hit 
        trough which the track goes is present in the final list 
    """
    hits = []

    z = 
    
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





