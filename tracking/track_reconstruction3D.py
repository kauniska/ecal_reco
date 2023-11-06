import sys
import time
tic = time.time() 
sys.path.insert(1, 'C:\\Users\\Pascal\\Desktop\\TP4a\\git\\TP4_ECAL-\\utils')
sys.path.insert(1, 'C:\\Users\\Pascal\\Desktop\\TP4a\\git_final_final\\ecal_reco\\tracking')
from data_loading import *
from tqdm import tqdm
from hit import Hit
from track3D import Track3D
from matplotlib import pyplot as plt
from track_reconstruction import *

import pandas as pd
import uproot
import numpy as np

def dist_from_curve(theta, gamma, x0, x1, isx):
    if(isx):
        return np.abs(-np.tan(theta)*np.sin(gamma)*(x1[1]-x0[2])+x0[0]-x1[0])
    else:
        return np.abs(-np.tan(theta)*np.cos(gamma)*(x1[1]-x0[2])+x0[1]-x1[0])
    
def error_fit(theta,gamma,x0,hits):
    error = 0
    for i in range(len(hits)):
        error += dist_from_curve(theta,gamma,x0,hits[i].get_pos(),hits[i].is_sidex)
    return error

def error_matrix(Hits,theta_vec,gamma_vec,x0,y0,z0):
    Matrix = np.ndarray(shape=(len(theta_vec),len(gamma_vec),len(x0),len(y0)), dtype=float)
    pos = []
    t = time.time()
    for i in range(len(theta_vec)):
        t = time.time()
        for j in range(len(gamma_vec)):
            for k in range(len(x0)):
                for l in range(len(y0)):
                    Matrix[i][j][k][l] = error_fit(theta_vec[i],gamma_vec[j],[x0[k],y0[l],z0],Hits)
        t2 = time.time()-t
    return Matrix

def find_3Dcurve(Hits, Nbr_It = 2, bool_plot = False):
     # Finds the highest cell of index idmax
    import time
    Max = 0
     # number of time the iteration is ran
    for i in range(len(Hits)):
        if(Hits[i].get_pos()[1]>Max):
            idmax = i
            Max = Hits[i].get_pos()[1]

    X0 = Hits[idmax]
    z0 = X0.get_pos()[1]
    # generate the grid
    
    theta = np.linspace(0,0.9*np.pi*0.5,5)
    gamma = np.linspace(0,np.pi*2,5)
    if(X0.is_sidex):
        x0 = np.linspace(X0.get_pos()[0]-0.8,X0.get_pos()[0]+0.8,5)
        y0 = np.linspace(0,38.4,5*24)
    else:
        y0 = np.linspace(X0.get_pos()[0]-0.8,X0.get_pos()[0]+0.8,5)
        x0 = np.linspace(0,38.4,5*24)
    
    # first approx in the whole bar
    
    M = error_matrix(Hits,theta,gamma,x0,y0,z0)
    
    idx=np.unravel_index(np.argmin(M, axis=None), M.shape)

    # iteration

    dt = theta[1]-theta[0]
    dg = gamma[1]-gamma[0]
    dx = x0[1]-x0[0]

    for i in range(Nbr_It):
        dt = dt/5
        dg = dg/5
        dx = dx/5

        #new grid
        theta = np.linspace(theta[idx[0]]-dt,theta[idx[0]]+dt,5)
        gamma = np.linspace(gamma[idx[1]]-dg,gamma[idx[1]]+dg,5)
        x0    = np.linspace(x0[idx[2]]-dx,x0[idx[2]]+dx,5)
        y0    = np.linspace(y0[idx[3]]-dx,y0[idx[3]]+dx,5)

        M = error_matrix(Hits,theta,gamma,x0,y0,z0)
        idx=np.unravel_index(np.argmin(M, axis=None), M.shape)
        ret = [theta[idx[0]], gamma[idx[1]],x0[idx[2]],y0[idx[3]],z0,M[idx]]
    if(bool_plot):
        from mpl_toolkits import mplot3d
        from mpl_toolkits.mplot3d import Axes3D
        theta = ret[0]
        gamma = ret[1]
        X0    = [ret[2], ret[3], ret[4]]

        a = [np.sin(gamma)*np.sin(theta),np.cos(gamma)*np.sin(theta),-np.cos(theta)]
        tt = np.linspace((16-X0[2])/a[2],-X0[2]/a[2],1000)

        xx = X0[0]+a[0]*tt
        yy = X0[1]+a[1]*tt
        zz = X0[2]+a[2]*tt
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot the line in 3D
        ax.plot(xx, yy, zz,'k-', label='3D track')

        # Set labels for the axes
        ax.set_xlabel('X [cm]')
        ax.set_ylabel('Y [cm]')
        ax.set_zlabel('Z [cm]')
        # Set custom axis limits
        ax.axis('equal')
        ax.set_xlim(0, 38.4)  # Set limits for the X axis
        ax.set_ylim(0, 38.4)  # Set limits for the Y axis
        ax.set_zlim(0, 16)  # Set limits for the Z axis
        # Add a legend
        ax.legend()

        # Show the plot
        plt.show()
    return ret