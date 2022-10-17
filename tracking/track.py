"""
Class implementing tracks, contains the list of hits, the physical information, the tangent, the x0.
The methods are the chi^2. 
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import OptimizeWarning
import scipy.optimize as opt


def gaussian2D(xy, x0, y0, sigma_x, sigma_y, A=1, theta=0):
    """ Gaussian in 2D with maximum at (x0, y0) of amplitude A and std deviation (sigma_x, sigma_y) rotated around an angle theta
    : (x,y): position the function is evaluated at
      (x0, y0): center of gaussian
      (sigma_x, sigma_y): std deviation along both axes
      A: amplitude, if -1 then normalized to 1 (-1 by default)
      theta: angle of rotation (radian) (0 by default)
    return: scalar
    """
    (x, y) = np.asarray(xy).reshape(2, int(np.shape(xy)[0]/2))
    a = np.cos(theta)**2/(2*sigma_x**2) + np.sin(theta)**2/(2*sigma_y**2)
    b = -np.sin(2*theta)/(4*sigma_x**2) + np.sin(2*theta)**2/(4*sigma_y**2)
    c = np.sin(theta)**2/(2*sigma_x**2) + np.cos(theta)**2/(2*sigma_y**2)
    r = np.exp(-(a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2))
    if np.sum(r) != 0:
        return A*r/np.sum(r)
    else:
        return np.inf

class Track:
    def __init__(self, *args):
        """Creates a Track from arguments.
            If zero argument, then it's an empty instance
            If one argument (list of Hits), then it computes the best track using Hough transform
            If 3 arguments (list of Hits, t, x0), then it simply writes all the parameters
            
        Args:
            hits (Hit): recorded hits associated to this track
            t (float): tangent angle (0° is vertical)
            x0 (flaot): extrapolated coordinate of the crossing of the top of the box
            hits_index: index of the hits in the event when the data was recorded
        """
        self.steps = 8
        if len(args) == 0:
            self.hits = []
            self.t = None
            self.x0 = None
            self.n_freedom = None
            self.reduced_chi2 = None
            self.mean_time = None
            self.hits_index = None
        elif len(args) == 1:
            self.hits = args[0]
            self.t, self.x0, self.hits_index = self.find_track()
            self.mean_time = None # TODO: implement
            self.n_freedom = len(self.hits) - 1 # two parameters: f(x) = a*x + b, number of data points = len + 1
            if self.n_freedom > 0:
                self.reduced_chi2 = self.chi2()/self.n_freedom
            else:
                self.reduced_chi2 = np.inf
        elif len(args) == 4:
            self.hits = args[0]
            self.t = args[1]
            self.x0 = args[2]
            self.mean_time = args[3]
            self.n_freedom = len(self.hits) - 1  # two parameters: f(x) = a*x + b, number of data points = len + 1
            if self.n_freedom > 0:
                self.reduced_chi2 = self.chi2()/self.n_freedom
            else:
                self.reduced_chi2 = np.inf
            self.hits_index = [i for i in range(len(self.hits))]
        else:
            raise ValueError("not the correct number of arguments given")
        
    def keep_hits_only(self):
        """Keeps only the hits given by the hits_index
        """
        self.hits = [self.hits[i] for i in self.hits_index]
        self.n_freedom = len(self.hits) - 1
        
    def x(self, z):
        if self.t != 0:
            return 8 + (z - self.x0) / self.t
        else:
            return self.x0
            
    def get_timestamps(self):
        """Gets the timestamps for all the hits

        Returns:
            list of float: timestamps of hits
        """
        return [hit.timestamp + hit.timestamp_event for hit in self.hits]

    def print(self):
        """Prints the relevant information (reduced chi^2, tangent angle and x0)
        """
        print("Reduced chi^2 = {:.2f}".format(self.reduced_chi2))
        print("t = {:.2f},\t x0 = {:.2f}".format(self.t, self.x0))

    def chi2(self):
        """Computes the chi^2 between the hits and the track
        Returns:
            float: chi^2
        """
        hits = self.get_hits_coords()
        tracks = self.get_tracks()
        hits.sort(key=lambda x: x[1])
        tracks.sort(key=lambda x: x[1])
        return np.sum([(hit[0] - track[0])**2 / 0.5 for hit, track in zip(hits, tracks)])
    
    def is_good_fit(self):
        return (self.chi2()/self.n_freedom < 5 * 3.841)
    
    def find_track(self, plot = False):
        """Finds the best parameters of a track passing through the hits, can plot the recorded hits and track

        Args:
            hits (list of Hit): recorded hits
            plot (bool, optional): If a figure of the track is to be made or not. Defaults to False.

        Returns:
            x0: coordinate at the top of the box
            t: tan of the angle (0° is vertical)
            indices: indices of the hits considered
        """
        maxi = 6  # maxi=5 => angle scanning between [-78.7°,78,7°]
        sampling = 5
        angle_sampling = 240
        T = np.linspace(-maxi, maxi, angle_sampling)
        x0s = np.empty(len(self.hits) * sampling * sampling * angle_sampling)
        txs = np.empty(len(self.hits) * sampling * sampling * angle_sampling)
        for n, hit in enumerate(self.hits):
            zs = np.linspace(hit.coord[1] - 0.5, hit.coord[1] + 0.5, sampling)
            xs = np.linspace(hit.coord[0] - 0.5, hit.coord[0] + 0.5, sampling)
            for i, z in enumerate(zs):
                for j, x in enumerate(xs):
                    index_prefix = n * sampling * sampling * angle_sampling + i * sampling * angle_sampling + j * angle_sampling
                    x0s[index_prefix:index_prefix + angle_sampling] = z - T * (x - 8)
                    txs[index_prefix:index_prefix + angle_sampling] = T

        H, ts, xs = np.histogram2d(txs, x0s, bins=[100, 100])
        id_t, id_x0 = np.unravel_index(np.argmax(H, axis=None), H.shape)
        self.t = ts[id_t]
        self.x0 = xs[id_x0]
        fit = [[self.x(z), z] for z in np.linspace(1, self.steps, self.steps)]

        if plot:
            fig, axs = plt.subplots(1, 2)
            fig.set_size_inches(12, 3, forward=True)
            fig.set_dpi(140)
            axs[0].hist2d(txs, x0s, bins=angle_sampling, cmap='inferno')
            # axs[0].hist2d(txs, x0s, bins = angle_sampling, cmap = 'inferno')
            axs[0].plot([self.t], [self.x0], 'bx', label='max')
            axs[0].set(xlabel='$t$', ylabel='$x$')
            axs[0].legend()

            hitsX = [hit.coord[0] for hit in self.hits]
            hitsZ = [hit.coord[1] for hit in self.hits]
            axs[1].hist2d(hitsX, hitsZ, bins=[24, 8], range=[[1, 24], [1, 8]], cmap='inferno')
            axs[1].plot([f[0] for f in fit], [f[1] for f in fit], 'b+')
            axs[1].set_xticks(np.linspace(1, 24, 25))
            axs[1].set_yticks(np.linspace(1, 8, 9))
            axs[1].grid(True, which='major')
            axs[1].grid(False, which='minor')
            coords_x = []
            coords_z = []
            for i in self.hits_index:
                if self.hits[i].coord in fit:
                    coords_x.append(self.hits[i].coord[0])
                    coords_z.append(self.hits[i].coord[1])
            axs[1].plot(coords_x, coords_z, 'k*')
            axs[1].set(xlabel='$x$', ylabel='$z$')

        self.hits_index = []
        return self.x0, self.t, self.hits_index

    def precise_track(self, plot = False):
        """Makes a very precise fit of the track, updates the angle of the object
        
        Returns:
            x0: coordinate at the top of the box
            t: tan of the angle (0° is vertical)
            indices: indices of the hits considered
        """
        old_x0 = self.x0
        old_t = self.t
        maxi = 6  # max=5 => angle scanning between [-78.7°,78,7°]
        angle_sampling = 480
        sampling = 20
        T = np.linspace(-maxi, maxi, angle_sampling)
        x0s = np.empty(len(self.hits) * sampling * sampling * angle_sampling)
        txs = np.empty(len(self.hits) * sampling * sampling * angle_sampling)
        for n, hit in enumerate(self.hits):
            zs = np.linspace(hit.coord[1] - 0.5, hit.coord[1] + 0.5, sampling)
            xs = np.linspace(hit.coord[0] - 0.5, hit.coord[0] + 0.5, sampling)
            for i, z in enumerate(zs):
                for j, x in enumerate(xs):
                    index_prefix = n * sampling * sampling * angle_sampling + i * sampling * angle_sampling + j * angle_sampling
                    x0s[index_prefix:index_prefix +
                        angle_sampling] = z - T * (x - 8)
                    txs[index_prefix:index_prefix + angle_sampling] = T

        H, ts, xs = np.histogram2d(txs, x0s, bins = [angle_sampling, angle_sampling], normed = True)
        id_t, id_x0 = np.unravel_index(np.argmax(H, axis=None), H.shape)
        self.t = ts[id_t]
        self.x0 = xs[id_x0]
        
        fit = self.get_tracks()
        self.hits_index = self.get_indices(False)
        
        self.reduced_chi2 = self.chi2() / self.n_freedom
        
        if plot:
            fig, axs = plt.subplots(1, 2)
            fig.set_size_inches(12, 3, forward=True)
            fig.set_dpi(140)
            axs[0].hist2d(txs, x0s, bins = int(2 * angle_sampling / maxi), cmap='inferno', range=[[self.t - 2, self.t+2], [self.x0-10, self.x0+10]])
            # axs[0].hist2d(txs, x0s, bins = angle_sampling, cmap = 'inferno')
            axs[0].plot([self.t], [self.x0], 'bx', label = 'precise')
            axs[0].plot([old_t], [old_x0], 'rx', label = 'rough')
            axs[0].set(xlabel = '$t$', ylabel = '$x$')
            axs[0].legend()

            hitsX = [hit.coord[0] for hit in self.hits]
            hitsZ = [hit.coord[1] for hit in self.hits]
            axs[1].hist2d(hitsX, hitsZ, bins=[24, 8], range=[[1, 24], [1, 8]], cmap='inferno')
            axs[1].plot([f[0] for f in fit], [f[1] for f in fit], 'b+')
            axs[1].set_xticks(np.linspace(1, 24, 25))
            axs[1].set_yticks(np.linspace(1, 8, 9))
            axs[1].grid(True, which = 'major')
            axs[1].grid(False, which = 'minor')
            coords_x = []
            coords_z = []
            for i in self.hits_index:
                if self.hits[i].coord in fit:
                    coords_x.append(self.hits[i].coord[0])
                    coords_z.append(self.hits[i].coord[1])
            axs[1].plot(coords_x, coords_z, 'k*')
            axs[1].set(xlabel = '$x$', ylabel = '$z$')

        return self.t, self.x0, self.hits_index
        
    def get_indices(self, redo_track = True):
        c = self.chi2()
        if redo_track:
            self.precise_track()
        self.hits_index = []
        for t in self.get_tracks():
            for i, hit in enumerate(self.hits):
                if (hit.coord[0] == int(np.round(t[0]))) and (hit.coord[1] == int(np.round(t[1]))):
                    self.hits_index.append(i)
        return self.hits_index
    
    def get_tracks(self):
        """Returns a list of the coordinates of the tracks
        
        Returns:
            list of coordinates (x, z)
        """
        return [[self.x(z), z] for z in np.linspace(1, self.steps, self.steps)]
    
    def get_hits_coords(self):
        """Gets the list of coordinates of the hits

        Returns:
            list of [x, z]: coordinates of all hits
        """
        return [hit.coord for hit in self.hits]

    def get_time_interval(self):
        last_hit = self.hits[np.argmax(self.hits_index, axis=0)]
        if last_hit.coord[1] < 8: # because if last hit is on the last level, we don't know if it left the detector
            ## and if there are hits close to the end point
            distances = self._dr(last_hit, self.hits)
            if np.any(np.array(distances) < 2):
                return np.mean(self.get_timestamps()) - self.hits[0].timestamp_event

    def _dr(self, hit_ref, hits):
        """Computes distance between a reference hit and a list of hits

        Args:
            hit_ref (Hit): reference Hit
            hits (list of Hit): list of Hit

        Returns:
            float : distance
        """
        ## eval array of vectors connecting all hits with ref hit
        dr_u = [hit.coord[0] - hit_ref.coord[0] for hit in hits]
        dr_d = [hit.coord[1] - hit_ref.coord[1] for hit in hits]

        ## eval module
        dr_mod = [np.sqrt(dr_u[i]*dr_u[i] + dr_d[i]*dr_d[i]) for i in range(len(dr_u))]
        return dr_mod
