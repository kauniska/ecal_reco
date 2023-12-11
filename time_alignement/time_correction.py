import sys
import os 
import fnmatch
import matplotlib as plt
import math
import copy
import pickle

sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\utils')
sys.path.insert(1, r'C:\Users\eliot\EPFL\TP4_ECAL\Code\ecal_reco\tracking')
from track_reconstruction import *
from parameters import *
from track import Track
from track3D import Track3D
# from track import Track
# from track3D import Track3D

def time_correction_fiber(*args):
    Speed_In_Fiber = 15 # cm/ns
    Speed_Of_Light = 30 # cm/ns
    ## clock cycle = 6.25 nanosecond
    Speed_In_Fiber = 1/convert_ns_to_clockcycle(1/Speed_In_Fiber) # cm/clock cycle
    Speed_Of_Light = 1/convert_ns_to_clockcycle(1/Speed_Of_Light) # cm/clock cycle
    zmax =  thickness+thickness_screen + (8-0.5)*thickness + (8-1)*(2*thickness_screen+thickness)

    #If one argument whcih is Track3D, change the timestamp of each hits of the track and return the track
    
    if len(args)== 1 : 
        # print("fiber")
        Tx = copy.deepcopy(args[0].x)
        Ty = copy.deepcopy(args[0].y)
        
        newx = []
        newy = []
        zmax =  thickness+thickness_screen + (8-0.5)*thickness + (8-1)*(2*thickness_screen+thickness)
        for h in Tx.hits:
            newx.append(h)
            heightcorr = math.sqrt(((h.get_pos()[1]-zmax)**2)*(Tx.t**2 + Ty.t**2 + 1))/Speed_Of_Light
            newx[-1].timestamp = h.timestamp - Ty.x(h.get_pos()[1])/Speed_In_Fiber-heightcorr
        for h in Ty.hits:
            newy.append(h)
            heightcorr = math.sqrt(((h.get_pos()[1]-zmax)**2)*(Tx.t**2 + Ty.t**2 + 1))/Speed_Of_Light
            newy[-1].timestamp = h.timestamp - Tx.x(h.get_pos()[1])/Speed_In_Fiber-heightcorr
        Txprime = Track(newx)
        Typrime = Track(newy)
        return Track3D(Txprime,Typrime)
    
    # if len(args)== 5 :
    #     timestamp = args[0]
    #     is_sidex = args[1]
    #     z = args[4]

    #     if is_sidex :
    #         x = args[2]
    #         y = args[3]
    #         # heightcorr = math.sqrt(((z-zmax)**2)*(Tx.t**2 + Ty.t**2 + 1))/Speed_Of_Light
    #         newx[-1].timestamp = h.timestamp - y(z)/Speed_In_Fiber
    #     else : 
    #         y = args[2]
    #         x =args[3]

        

    
def time_correction_electronics(*args) :

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
        # print("PCB")
        T = copy.deepcopy(args[0])
        for h in T.x.hits:
            [x,z] = h.coord
            tofpet = mapping_inv_2D(1,x,z)
            h.timestamp = copy.deepcopy(h.timestamp - mapping_SiPM_delay(tofpet[0], tofpet[1]))
        for h in T.y.hits:
            [y,z] = h.coord
            tofpet = mapping_inv_2D(0,y,z)
            h.timestamp = copy.deepcopy(h.timestamp - mapping_SiPM_delay(tofpet[0], tofpet[1]))
        return T
    

    ## Apply a time correction from general offset and coming from the time alignement procedure
def time_correction_offset(*args) :

    # If 3 arguments which are a timestamp and coordinate (tofpet id and channel), change the timestamp and returns it
    if len(args) == 3 :
        timestamp = args[0] 
        tofpet_id= args[1] 
        tofpet_channel= args[2] 
        [x,z]=mapping_2D(tofpet_id, tofpet_channel)
    
        ##Use the value of muX and muY computed form time alignement procedure and located in parameters.py
        if is_sidex(tofpet_id) :
            return timestamp - muX[x-1,z-1]
        else :
            return timestamp - muY[x-1,z-1]
    
    # If 1 argument which is a track3D, change the timestamp of every hit and return the track3D
    if len(args) == 1 :
        # print("offset")
        T = copy.deepcopy(args[0])
        for h in T.x.hits:
                # tofpet = mapping_inv_2D(1,h.coord[0],h.coord[1])
                h.timestamp = h.timestamp - muX[h.coord[0]-1,h.coord[1]-1]
        for h in T.y.hits:
                # tofpet = mapping_inv_2D(0,h.coord[0],h.coord[1])
                h.timestamp = h.timestamp - muY[h.coord[0]-1,h.coord[1]-1]

        return T

def time_correction_global(*args) :
    # if 1 argument which is a track3D. change the timestamp of every hit 
    # with respect of all the time correction functions.
    if len(args) == 1 :
        # print(args[0].x[0].timestamps)
        T = copy.deepcopy(args[0])
        T2 = copy.deepcopy(time_correction_fiber(T))
        T3 = copy.deepcopy(time_correction_electronics(T2))
        T4= copy.deepcopy(time_correction_offset(T3))
        # print(args[0].x[0].timestamps)
        # print("A")
        return T4
    

    elif len(args)==2 :
        T = args[0] #muon track before decay
        S = args[1] #list of hits of the shower after decay
        # print(S)
        for h in S :
            # if h.is_sidex :
            #     # Find the object with the smallest coord[1]
            #     min_hit_track = min(T.y, key=lambda x: x.coord[1])
            #     h.timestamp = copy.deepcopy(time_correction_fiber(h.timestamp,h.is_sidex,h.coord[0], min_hit_track.coord[0], h.coord[1]))
            #  else :
            #     min_hit_track = min(T.x, key=lambda x: x.coord[1])
            #     h.timestamp = copy.deepcopy(time_correction_fiber(h.timestamp,h.is_sidex,h.coord[0], min_hit_track.coord[0], h.coord[1]))
        
            is_sidex = h.is_sidex
            x = h.coord[0]
            z = h.coord[1]
            # print(x,z)
            tofpet= mapping_inv_2D(is_sidex,x,z)
            # print(tofpet)
            tofpet_id = tofpet[0]
            tofpet_channel = tofpet[1]
            h.timestamp = time_correction_electronics(h.timestamp,tofpet_id, tofpet_channel)
            h.timestamp = time_correction_offset(h.timestamp,tofpet_id, tofpet_channel)
            # print(S)
        return S
    else :
        raise ValueError("Expect one argument (Track3D) or two (Track3D and list of hits)")

