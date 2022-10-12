"""
Class implementing tracks, contains the list of hits, the physical information, the tangent, the x0.
The methods are the chi^2. 
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
            self.t, self.x0, self.hits_index = self.find_track(self)
            self.t = None
            self.x0 = None
            self.mean_time = None # TODO: implement
            # two parameters: f(x) = a*x + b, number of data points = len + 1
            self.n_freedom = len(self.hits) - 1
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
        return sum([(hit.coord[0] - track[0])**2 / track[0] for hit, track in zip(self.hits, self.get_tracks()) if track[0] != 0])
    
    def is_good_fit(self):
        return (self.chi2() > 3.841)
    
    def tracks(hits, plot=False):
        """Finds the best parameters of a track passing through the hits, can plot the recorded hits and track

        Args:
            hits (list of Hit): recorded hits
            plot (bool, optional): If a figure of the track is to be made or not. Defaults to False.

        Returns:
            x0: coordinate at the top of the box
            t: tan of the angle (0° is vertical)
            indices: indices of the hits considered
        """
        max = 3  # max=5 => angle scanning between [-78.7°,78,7°]
        sampling = 5
        angle_sampling = 150
        T = np.linspace(-max, max, angle_sampling)
        x0s = np.empty(len(hits) * sampling * sampling * angle_sampling)
        txs = np.empty(len(hits) * sampling * sampling * angle_sampling)
        for n, hit in enumerate(hits):
            zs = np.linspace(hit.coord[1] - 0.5, hit[1] + 0.5, sampling)
            xs = np.linspace(hit.coord[0] - 0.5, hit[0] + 0.5, sampling)
            for i, z in enumerate(zs):
                for j, x in enumerate(xs):
                    index_prefix = n * sampling * sampling * angle_sampling + \
                        i * sampling * angle_sampling + j * angle_sampling
                    x0s[index_prefix:index_prefix +
                        angle_sampling] = -(z - 9) * T + x
                    txs[index_prefix:index_prefix + angle_sampling] = T

        H, ts, xs = np.histogram2d(txs, x0s, bins=[100, 100])
        id_t, id_x0 = np.unravel_index(np.argmax(H, axis=None), H.shape)
        t = ts[id_t]
        x0 = xs[id_x0]
        indices = range(len(hits))
        fit = [[np.round(t * (z-9) + x0), z] for z in range(9)]

        if plot:
            plt.figure()
            plt.hist2d(txs, x0s, bins=[len(T), len(T)], cmap='inferno')
            plt.plot([t], [x0], 'rx')
            plt.xlabel('$t$')
            plt.ylabel('$x$')

            plt.figure()
            hitsX = [hit.coord[0] for hit in hits]
            hitsZ = [hit.coord[1] for hit in hits]
            plt.hist2d(hitsX, hitsZ, bins=[24, 8], range=[
                       [1, 24], [1, 8]], cmap='inferno')
            plt.plot([t * (z-9) + x0 for z in range(9)], range(9), 'r-')
            coords_x = []
            coords_z = []
            for hit in hits:
                if hit.coord in fit:
                    coords_x.append(hit.coord[0])
                    coords_z.append(hit.coord[1])
            plt.plot(coords_x, coords_z, 'b+')
            plt.xlabel('$x$')
            plt.ylabel('$z$')

        return x0, t, indices
    
    def precise_track(self):
        """Makes a very precise fit of the track, updates the angle of the object
        
        Returns:
            x0: coordinate at the top of the box
            t: tan of the angle (0° is vertical)
            indices: indices of the hits considered
        """
        max = 6  # max=5 => angle scanning between [-78.7°,78,7°]
        sampling = 20
        angle_sampling = 500
        T = np.linspace(-max, max, angle_sampling)
        x0s = np.empty(len(self.hits) * sampling * sampling * angle_sampling)
        txs = np.empty(len(self.hits) * sampling * sampling * angle_sampling)
        for n, hit in enumerate(self.hits):
            zs = np.linspace(hit.coord[1] - 0.5, hit[1] + 0.5, sampling)
            xs = np.linspace(hit.coord[0] - 0.5, hit[0] + 0.5, sampling)
            for i, z in enumerate(zs):
                for j, x in enumerate(xs):
                    index_prefix = n * sampling * sampling * angle_sampling + \
                        i * sampling * angle_sampling + j * angle_sampling
                    x0s[index_prefix:index_prefix +
                        angle_sampling] = -(z - 9) * T + x
                    txs[index_prefix:index_prefix + angle_sampling] = T

        H, ts, xs = np.histogram2d(txs, x0s, bins=[100, 100])
        id_t, id_x0 = np.unravel_index(np.argmax(H, axis=None), H.shape)
        t = ts[id_t]
        x0 = xs[id_x0]
        fit = [ [-(z - 9) * t + x0, z] for z in range(9)]
        indices = [i for i, hit in enumerate(self.hits) if hit.coord in fit]
        
        self.reduced_chi2 = self.chi2() / self.n_freedom
        
        return t, x0, indices
        
    
    def get_tracks(self):
        """Returns a list of the coordinates of the tracks
        
        Returns:
            list of coordinates
        """
        steps = 9
        zs = np.linspace(1,steps,steps)
        xs = self.t * (zs-9) * 2 + self.x0
        return [[x, (z-1) * 2] for x, z in zip(xs, zs)] # track projection on X
    
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
