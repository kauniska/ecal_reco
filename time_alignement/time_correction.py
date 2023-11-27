import sys
import os 
import fnmatch
import matplotlib as plt

sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\utils')
sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\tracking')
from track_reconstruction import *
from parameters import *
# from track import Track
# from track3D import Track3D

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
    

    ## Apply a time correction coming from time resultion and light propagation in fibers
def time_correction_offset(args) :

    muX = np.nan_to_num(np.ndarray(shape=(8,64), dtype=float), nan=0, posinf=0, neginf=0)*0
    muY = np.nan_to_num(np.ndarray(shape=(8,64), dtype=float), nan=0, posinf=0, neginf=0)*0

    # If 3 arguments which are a timestamp and coordinate (tofpet id and channel), change the timestamp and returns it
    if len(args) == 3 :
        timestamp = args[0] 
        tofpet_id= args[1] 
        tofpet_channel= args[2] 
    
        
        if is_sidex(tofpet_id) :
            return timestamp - muX[tofpet_id,tofpet_channel]
        else :
            return timestamp - muY[tofpet_id,tofpet_channel]
    
    # If 1 argument which is a track3D, change the timestamp of every hit and return the track3D
    if len(args) == 1 :
        T = args[0]
        for h in T.x.hits:
            h.timestamp = h.timestamp - muX[mapping_inv_2D(1,h.get_pos())]
        for h in T.y.hits:
            h.timestamp = h.timestamp - muY[mapping_inv_2D(0,h.get_pos())]
      
        return T



def time_correction_electronics(args) :

    # If 3 arguments which are a timestamp and coordinate (tofpet id and channel), change the timestamp and returns it
    if len(args) == 3 :
        timestamp = args[0] 
        tofpet_id= args[1] 
        tofpet_channel= args[2] 

        if is_sidex(tofpet_id) :
            return timestamp - mapping_SiPM_delay(tofpet_id,tofpet_channel)
        else :
            return timestamp - mapping_SiPM_delay(tofpet_id,tofpet_channel)
    
    # If 1 argument which is a Track3D, change timestamp of each hit and return the track3D
    if len(args) == 1 :
        T = args[0]
        for h in T.x.hits:
            h.timestamp = h.timestamp - mapping_SiPM_delay[mapping_inv_2D(1,h.get_pos())]
        for h in T.y.hits:
            h.timestamp = h.timestamp - mapping_SiPM_delay[mapping_inv_2D(0,h.get_pos())]
        return T