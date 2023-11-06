import sys
import os 
import fnmatch

sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\utils')
sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\tracking')
from track import Track
from track3D import Track3D

def time_correction(T3D):
    Speed_In_Fiber = 15 # cm/ns
    fac = 6.25 # ns/clock cycle
    Speed_In_Fiber = Speed_In_Fiber*fac # cm/clock cycle
    Tx = T3D.x
    Ty = T3D.y
    
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