import numpy as np
import matplotlib.pyplot as plt
from track_reconstruction import max_overlap
import pandas as pd
import sys
sys.path.insert(1, r"C:\Users\nelg\Desktop\Cours\Labo\TP4\Git\utils")
from parameters import *

# event = df_hits.iloc[interesting_events[600]]
def new_method_tracks(hits):
    # hits = []
    # for i in range(event['n_hits']):
    #     h = Hit(event,i)
    #     if h.is_sidex:
    #         hits.append(h)
    n_points = 100
    tmax = 5

    map = np.zeros((2*n_points,2*n_points))
    tneg = np.linspace(-tmax, 0, n_points)
    tpos = np.linspace(0, tmax, n_points)
    T = np.append(tneg,tpos)
    x0_max = n_strips*width + tmax*2*n_layers*(thickness+thickness_screen)
    x0_min = - tmax*2*n_layers*(thickness+thickness_screen)
    x = np.linspace(x0_min,x0_max,2*n_points)
    for h in hits:
        pos = h.get_pos()
        x0u_neg = pos[0] + width/2 - (pos[1] + thickness/2)*tneg
        x0d_neg = pos[0] - width/2 - (pos[1] - thickness/2)*tneg
        x0u_pos = pos[0] + width/2 - (pos[1] - thickness/2)*tpos
        x0d_pos = pos[0] - width/2 - (pos[1] + thickness/2)*tpos
        x0d = np.concatenate((x0d_neg,x0d_pos))
        x0u = np.concatenate((x0u_neg,x0u_pos))

        for j in range(len(T)):
            for i in range(len(x)):
                if x[i] >= x0d[j] and x[i] < x0u[j]:
                    map[i][j] += 1
    # plt.figure()
    # plt.matshow(map)
    index_max_tx = np.where(map == map.max())
    index_max_x = round(np.sum(index_max_tx[0])/len(index_max_tx[0]))
    index_max_t = round(np.sum(index_max_tx[1])/len(index_max_tx[1]))
    # plt.plot(index_max_t,index_max_x,'rx')

    x0 = x[index_max_x]
    t = T[index_max_t]

    # detector_map = np.zeros((8,24))
    # for h in hits:
    #     detector_map[8-h.coord[1]][h.coord[0]-1] += 1

    # z_grid = np.linspace(0,2*8,10)
    # x_track = t*(z_grid-16)+x0-1.6

    # plt.figure()
    # plt.imshow(detector_map, interpolation='nearest', extent=[0,24*1.6,0,8*2])
    # plt.plot(x_track,z_grid,'r-')
    # plt.show
    return t,x0

def old_method_tracks(hits):
    n_points=100
    max=5 # => angle scanning between [-78.7°,78,7°] 
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
   
    for i in range(len(hits)):
        pos = hits[i].get_pos()
        x0['xu'][i][0:n_points] = pos[0] + width/2 - (pos[1] + thickness/2)*tneg # Up boundary for x0 for hit number hit in region t<0
        x0['xu'][i][n_points:2*n_points] = pos[0] + width/2 - (pos[1] - thickness/2)*tpos # Up boundary for x0 for hit number hit in region t>0
        x0['xd'][i][0:n_points] = pos[0] - width/2 - (pos[1] - thickness/2)*tneg
        x0['xd'][i][n_points:2*n_points] = pos[0] - width/2 - (pos[1] + thickness/2)*tpos
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