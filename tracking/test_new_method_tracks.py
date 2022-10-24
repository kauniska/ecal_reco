import numpy as np
import matplotlib.pyplot as plt
# event = df_hits.iloc[interesting_events[600]]
def new_method_tracks(hits):
    # hits = []
    # for i in range(event['n_hits']):
    #     h = Hit(event,i)
    #     if h.is_sidex:
    #         hits.append(h)
    n_points = 100
    tmax = 6.25
    n_strips = 24
    n_layers = 8
    width = 1
    thickness = 1
    map = np.zeros((2*n_points,2*n_points))
    tneg = np.linspace(-tmax, 0, n_points)
    tpos = np.linspace(0, tmax, n_points)
    T = np.append(tneg,tpos)
    x0_max = (n_strips+0.5)*width + tmax*(n_layers+0.5)*thickness
    x0_min = 0.5*width - tmax*(n_layers+0.5)*thickness
    x = np.linspace(x0_min,x0_max,2*n_points)
    for h in hits:
        x0u_neg = (h.coord[0]+0.5)*width + (h.coord[1]+0.5)*tneg*thickness
        x0d_neg = (h.coord[0]-0.5)*width + (h.coord[1]-0.5)*tneg*thickness
        x0u_pos = (h.coord[0]+0.5)*width + (h.coord[1]-0.5)*tpos*thickness
        x0d_pos = (h.coord[0]-0.5)*width + (h.coord[1]+0.5)*tpos*thickness
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
