#!/usr/bin/env python3
"""
# Contains the functions used to determine the tracks of a signal. Doesn't contains any filter, chi2 function
# Can be used to determine if the computed track fits with the data 
# Use Hough transfomations on each side of the calorimeter to determine the parameters of the projection of the 
# track. Having the projections of the track make it easy to reconstruct the 3D track 
# Autors: v_0 Irwan Ledorze, Georgios Demetriou
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys
sys.path.insert(1, r"C:\Users\nelg\Desktop\Cours\Labo\TP4\Git\utils")
from parameters import *


# Check that a is tofpet on the x side of the calorimeter
def is_sidex(a):
    if(a==0 or a==1 or a==4 or a==5):
        return 1
    else:
        return 0

# Determine the (X,Z) or (Y,Z) coordinate of a hits, depending on what tofpet_id is entered. Triplet=[channel,t_id,layer]
def mapping_2D(t_id,channel):
    mapping=[[  [10 , 3], [22 , 3], [ 3 , 3], [15 , 3], [ 9 , 3], [21 , 3], [ 4 , 3], [16 , 3], [ 8 , 3], [20 , 3], [ 5 , 3], [17 , 3], [ 7 , 3], [19 , 3], [ 6 , 3], [18 , 3], [ 7 , 2], [24 , 2], [ 1 , 2], [13 , 2], [ 8 , 2], [23 , 2], [ 2 , 2], [14 , 2], [ 9 , 2], [22 , 2], [ 3 , 2], [15 , 2], [10 , 2], [21 , 2], [ 4 , 2], [16 , 2], [11 , 2], [20 , 2], [ 5 , 2], [17 , 2], [12 , 2], [19 , 2], [ 6 , 2], [18  ,2], [ 7 , 1], [24 , 1], [ 1 , 1], [13 , 1], [ 8 , 1], [23 , 1], [ 2 , 1], [14 , 1], [ 9 , 1], [22 , 1], [ 3 , 1], [15 , 1], [10 , 1], [21 , 1], [ 4 , 1], [16 , 1], [11  ,1], [20 , 1], [ 5 , 1], [17 , 1], [12 , 1], [19 , 1], [ 6 , 1], [18 , 1],[ 7,  4], [24,  4], [ 1,  4], [13,  4], [ 8,  4], [23,  4], [ 2,  4], [14,  4], [ 9,  4], [22,  4], [ 3,  4], [15,  4], [10,  4], [21,  4], [ 4,  4], [16,  4], [11,  4], [20,  4], [ 5,  4], [17,  4], [12,  4], [19,  4], [ 6,  4], [18,  4], [12,  3], [24,  3], [ 1,  3], [13,  3], [11,  3], [23,  3], [ 2,  3], [14,  3]],[[10,  7], [22,  7], [ 3,  7], [15,  7], [ 9,  7], [21,  7], [ 4,  7], [16,  7], [ 8,  7], [20,  7], [ 5,  7], [17,  7], [ 7,  7], [19,  7], [ 6,  7], [18,  7], [ 7,  6], [24,  6], [ 1,  6], [13,  6], [ 8,  6], [23,  6], [ 2,  6], [14,  6], [ 9,  6], [22,  6], [ 3,  6], [15,  6], [10 , 6], [21,  6], [ 4,  6], [16,  6], [11,  6], [20,  6], [ 5,  6], [17,  6], [12,  6], [19,  6], [ 6,  6], [18,  6], [ 7,  5], [24,  5], [ 1,  5], [13,  5], [ 8,  5], [23,  5], [ 2,  5], [14,  5], [ 9,  5], [22,  5], [ 3,  5], [15,  5], [10,  5], [21,  5], [ 4,  5], [16 , 5], [11,  5], [20,  5], [ 5,  5], [17,  5], [12,  5], [19,  5], [ 6,  5], [18,  5],[ 7 , 8], [24 , 8], [ 1 , 8], [13 , 8], [ 8 , 8], [23 , 8], [ 2 , 8], [14 , 8], [ 9 , 8], [22 , 8], [ 3 , 8], [15 , 8], [10,  8], [21 , 8], [ 4 , 8], [16 , 8], [11 , 8], [20 , 8], [ 5 , 8], [17 , 8], [12 , 8], [19 , 8], [ 6 , 8], [18 , 8], [12 , 7], [24 , 7], [ 1 , 7], [13 , 7], [11 , 7], [23 , 7], [ 2 , 7], [14 , 7]]] ##
    if is_sidex(t_id):
        return mapping[int(t_id/4)][channel+32*np.mod(t_id,2)]
    else:
        t_id=t_id-2
        return mapping[int(t_id/4)][channel+32*np.mod(t_id,2)]

## Looks how many hits overlap at a certain angle t. Return the the hits index that overlap, the number of overlaping
# and the boundaries, boundaries are the extremal x that belongs to the overlap region
def max_overlap(x0,t):
    nb_hits=len(x0['xu'])
    x0u=[x0['xu'].iloc[i][t] for i in range(nb_hits)]
    x0d=[x0['xd'].iloc[i][t] for i in range(nb_hits)]
    index_sort=np.argsort(np.array(x0u))
    index_sort=index_sort[::-1]
    x0usorted=x0u
    x0usorted.sort(reverse=True)
    index_overlap=max([[j for j in index_sort[i:] if(x0usorted[j]>=x0d[index_sort[i]])] for i in range(nb_hits)],key=len)
    length=len(index_overlap)
    boundaries=[0,0]
    if len(index_overlap)>0:
        boundaries[0]=min([x0u[i] for i in index_overlap])
        boundaries[1]=max([x0d[i] for i in index_overlap])
    return np.sort(index_overlap),length,boundaries


## This function provide the parameters x0 (or y0 depending on which hits we provide) and tx the angle of the trac
def tracks(hits):
    n_points=100
    max=5 # max=5 => angle scanning between [-78.7°,78,7°] 
    tneg=np.linspace(-max,0,n_points) # region over which we want to look for the angle
    tpos=np.linspace(0,max,n_points)
    n_hits=len(hits)
    # x0=pd.DataFrame(columns=['xu','xd'])
    # x0['xu']=[np.zeros(2*n_points) for i in range(n_hits)]
    # x0['xd']=[np.zeros(2*n_points) for i in range(n_hits)]
    x0=pd.DataFrame(data = {
        'xu':[np.zeros(2*n_points) for i in range(n_hits)],
        'xd': [np.zeros(2*n_points) for i in range(n_hits)]
        })
   
    for hit in range(0,n_hits):
        x0['xu'][hit][0:n_points]=(hits[hit][0]+1)*1.6-(hits[hit][1]-8)*2*tneg # Up boundary for x0 for hit number hit in region t<0
        x0['xu'][hit][n_points:2*n_points]=(hits[hit][0]+1)*1.6-(hits[hit][1]-9)*2*tpos # Up boundary for x0 for hit number hit in region t>0
        x0['xd'][hit][0:n_points]=hits[hit][0]*1.6-(hits[hit][1]-9)*2*tneg
        x0['xd'][hit][n_points:2*n_points]=hits[hit][0]*1.6-(hits[hit][1]-8)*2*tpos    
    ##### Now have to find the region of overlap

    T=np.append(tneg,tpos)
    t_overlap=[[] for i in range(2*n_points)]
    boundaries=[[] for i in range(2*n_points)]
    overlap=0
    for t in range(2*n_points):
        a=0
        t_overlap[t],a, boundaries[t]=max_overlap(x0,t)
        if a>overlap:
            overlap=a
    t_max_overlap=[T[t] for t in range(2*n_points) if len(t_overlap[t])==overlap]
    min_max_overlap=[boundaries[t] for t in range(2*n_points) if len(t_overlap[t])==overlap]
    index_=[t_overlap[t] for t in range(2*n_points) if len(t_overlap[t])==overlap][0]
    t=np.mean(t_max_overlap)
    out=np.mean([np.mean(ov) for ov in min_max_overlap])
    return out, t, index_

def chi_2(Hits,track,index_):
    dof=1
    critX=3.841
    expected=np.array([track[Hits[i][1]-1][0] for i in index_])
    X2=sum([(Hits[index_[i]][0]*1.6-expected[i])**2/expected[i] for i in range(len(index_))])
    if X2> critX:
        return False
    else:
        return True

# returns the physical position (in cm) of the center of the cell, on the concerned side of the detector as a function of the (integer) coordinates
def coord_to_pos_x(coord):
    return (coord-0.5)*width

def coord_to_pos_z(coord,x_plane):
    z = (coord-0.5)*thickness + (coord-1)*(2*thickness_screen+thickness)
    if not x_plane:
        z += thickness+thickness_screen
    return z

def coord_to_pos(coord,x_plane):
    return np.array([coord_to_pos_x(coord[0]),coord_to_pos_z(coord[1],x_plane)])

# returns the coordinates and the side of the cell as a function of physical position (in cm)
def pos_to_coord(pos):
    if pos[1] < 0 or pos[1] > 2*n_layers*(thickness+thickness_screen)\
        or pos[0] < 0 or pos[0] > width*n_strips:
        raise ValueError("Position out of bound")
    elif pos[1] % (thickness+thickness_screen) > thickness:
        raise ValueError("z-coordinate cooresponds to passive layer")

    x = pos[0] // width + 1
    z = pos[1] // (thickness+thickness_screen)

    if pos[1] % 2*(thickness+thickness_screen) > thickness+thickness_screen:
        x_plane = False
        z = (z+1)/2
    else :
        x_plane = True
        z = z/2 + 1

    coord = np.array([round(x),round(z)])

    return coord,x_plane


# Plot the fired hits on the chosen side with the right geometry
def plot_hits(hits, x_plane = None, plot_perpendicular = False, scaling = 1, hits_next = None):
    '''
    Arguments :
        -x_plane : if True the plot shows the fired hits on the xz-plane, else on the yz-plane
        -plot_perpendicular : if True the hits on the other plane will be visible (in a darker color)
        -scaling : relative size of the plot
        -hits_next : if given, they will be considered as hits from product of muon decay in a next event. The hits will have another color
    '''
    passive_color = (80/255,80/255,80/255)
    active_color = (22/255,100/255,90/255)
    if x_plane == None:
        x_plane = hits[0].is_sidex
    hits_x = []
    hits_y = []
    for hit in hits:
        if hit.is_sidex:
            hits_x.append(hit)
        else:
            hits_y.append(hit)
    hit_color = (1,232/255,0)
    perpendicular_color = (1,232/255,0,0.25)

    if hits_next != None:
        hits_x_next = []
        hits_y_next = []
        for hit in hits_next:
            if hit.is_sidex:
                hits_x_next.append(hit)
            else:
                hits_y_next.append(hit)   
    hit_color_next = (250/255,100/255,0) 
    perpendicular_color_next = (250/255,100/255,0.25)    

    fig,ax = plt.subplots(figsize = (n_strips*width/2*scaling,2*n_layers*(thickness+thickness_screen)/2*scaling))
    ax.set_facecolor(active_color)
    ax.set_xlim([0,n_strips*width])
    ax.set_ylim([0,2*n_layers*(thickness+thickness_screen)])

    for i in range(2*n_layers):
        ax.axhspan(i*(thickness+thickness_screen)+thickness, (i+1)*(thickness+thickness_screen), facecolor=passive_color, alpha=1)

    if x_plane:
        for hit in hits_x:
            pos = hit.get_pos()
            rectangle = plt.Rectangle((pos[0]-width/2,pos[1]-thickness/2), width, thickness, fc=hit_color, ec='k', lw=1)
            ax.add_patch(rectangle)
            if hits_next != None:
                for hit in hits_x_next:
                    pos = hit.get_pos()
                    rectangle = plt.Rectangle((pos[0]-width/2,pos[1]-thickness/2), width, thickness, fc=hit_color_next, ec='k', lw=1)
                    ax.add_patch(rectangle)
        if plot_perpendicular:
            z_fired = np.zeros((n_layers,))
            for hit in hits_y:
                z_fired[hit.coord[1]-1] = 1
            for i,z in enumerate(z_fired):
                if z:
                    ax.axhspan(coord_to_pos_z(i+1,False)-thickness/2, coord_to_pos_z(i+1,False)+thickness/2, facecolor=perpendicular_color, alpha=0.25)
            if hits_next != None:
                z_next_fired = np.zeros((n_layers,))
                for hit in hits_y_next:
                    z_next_fired[hit.coord[1]-1] = 1
                for i,z in enumerate(z_next_fired):
                    if z:
                        ax.axhspan(coord_to_pos_z(i+1,False)-thickness/2, coord_to_pos_z(i+1,False)+thickness/2, facecolor=perpendicular_color_next, alpha=0.25)
    else:
        for hit in hits_y:
            pos = hit.get_pos()
            rectangle = plt.Rectangle((pos[0]-width/2,pos[1]-thickness/2), width, thickness, fc=hit_color, ec='k', lw=1)
            ax.add_patch(rectangle)
            if hits_next != None:
                for hit in hits_y_next:
                    pos = hit.get_pos()
                    rectangle = plt.Rectangle((pos[0]-width/2,pos[1]-thickness/2), width, thickness, fc=hit_color_next, ec='k', lw=1)
                    ax.add_patch(rectangle)
        if plot_perpendicular:
            z_fired = np.zeros((n_layers,))
            for hit in hits_x:
                z_fired[hit.coord[1]] = 1
            for i,z in enumerate(z_fired):
                if z:
                    ax.axhspan(coord_to_pos_z(i+1,True)-thickness/2, coord_to_pos_z(i+1,True)+thickness/2, facecolor=perpendicular_color, alpha=0.25)
            if hits_next != None:
                z_next_fired = np.zeros((n_layers,))
                for hit in hits_x_next:
                    z_next_fired[hit.coord[1]-1] = 1
                for i,z in enumerate(z_next_fired):
                    if z:
                        ax.axhspan(coord_to_pos_z(i+1,True)-thickness/2, coord_to_pos_z(i+1,True)+thickness/2, facecolor=perpendicular_color_next, alpha=0.25)

    return fig,ax

        

