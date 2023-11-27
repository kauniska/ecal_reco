import sys
import os 
import fnmatch

sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\utils')
sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\tracking')
from track import Track
from track_reconstruction import *
from track3D import Track3D

def time_correction_fiber(args):
    Speed_In_Fiber = 15 # cm/ns

    #If one argument whcih is Track3D, change the timestamp of each hits of the track and return the track
    if len(args)== 1 : 
        Tx = args[0].x
        Ty = args[0].y
        
        newx = []
        newy = []
        for h in Tx.hits:
            newx.append(h)
            newx[-1].timestamp = h.timestamp - Ty.x(h.get_pos()[1])/Speed_In_Fiber
        for h in Ty.hits:
            newy.append(h)
            newy[-1].timestamp = h.timestamp - Ty.x(h.get_pos()[1])/Speed_In_Fiber
        Txprime = Track(newx)
        Typrime = Track(newy)
        return Track3D(Txprime,Typrime)
    
    # If 2 arguments : a timestamp, a x/y coordinate(in physical units [m]))
    if len(args) == 2 :
        timestamp = args[0] - args[1]/Speed_In_Fiber

        return timestamp
    
    # If 3 arguments : A timestamp, and coordinates (in tofpet id and chanels), change the timestamp and return it
    if len(args) == 3 :
        timestamp = args[0] - mapping_2D(args[1],args[2])[1]/Speed_In_Fiber
        return timestamp