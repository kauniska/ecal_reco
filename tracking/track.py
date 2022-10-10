"""
Class implementing tracks, contains the list of hits, the physical information, the tangent, the x0.
The methods are the chi^2. 
"""
import numpy as np
import pandas as pd

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
    
    def find_tracks(self):
        """Finds the best linear track using Hough transform, gives the tangent angle and the coordinate of the crossing of the top of the box

        Returns:
            t: tangent angle (0° is vertical)
            x0: coordinate of the track at the top of the box
        """
        n_points = 100
        max = 5  # max=5 => angle scanning between [-78.7°,78,7°]
        
        # region over which we want to look for the angle
        tneg = np.linspace(-max, 0, n_points)
        tpos = np.linspace(0, max, n_points)
        n_hits = len(self.hits)
        
        x0 = pd.DataFrame(data = {
            'xu' : [np.zeros(2*n_points) for i in range(n_hits)],
            'xd' : [np.zeros(2*n_points) for i in range(n_hits)]
        })
        
        for hit in range(n_hits):
            x0['xu'][hit][:n_points] = (self.hits.coord[hit] + 1) * 1.6 - (self.hits[hit][1] - 8) * 2 * tneg
            x0['xu'][hit][n_points:2*n_points] = (self.hits.coord[hit] + 1) * 1.6 - (self.hits[hit][1] - 9) * 2 * tpos
            x0['xd'][hit][:n_points] = (self.hits.coord[hit] + 1) * 1.6 - (self.hits[hit][1] - 8) * 2 * tneg
            x0['xd'][hit][n_points:2*n_points] = (self.hits.coord[hit] + 1) * 1.6 - (self.hits[hit][1] - 8) * 2 * tpos

        T = np.append(tneg, tpos)
        t_overlap = [[] for i in range(2*n_points)]
        boundaries = [[] for i in range(2*n_points)]
        overlap = 0
        for t in range(2*n_points):
            a = 0
            t_overlap[t], a, boundaries[t] = self.max_overlap(x0, t)
            if a > overlap:
                overlap = a
        t_max_overlap = [T[t] for t in range(2*n_points) if len(t_overlap[t]) == overlap]
        min_max_overlap = [boundaries[t] for t in range(2*n_points) if len(t_overlap[t]) == overlap][0]
        t = np.mean(t_max_overlap)
        x0 = np.mean([np.mean(ov) for ov in min_max_overlap])
        indices = [t_overlap[t] for t in range(2*n_points) if len(t_overlap[t])==overlap][0]
        return t, x0, indices
    
    def max_overlap(self, x0, t):
        nb_hits = len(x0)
        x0u = [x0['xu'].iloc[i][t] for i in range(nb_hits)]
        x0d = [x0['xd'].iloc[i][t] for i in range(nb_hits)]
        index_sort = np.argsort(np.array(x0u))
        index_sort = index_sort[::-1]
        x0u_sorted = x0u
        x0u_sorted.sort(reverse = True)
        index_overlap = max([[j for j in index_sort[i:] if (x0u_sorted[j] >= x0d[index_sort[i]])] for i in range(nb_hits)], key = len)
        boundaries = [0,0]
        if len(index_overlap) > 0:
            boundaries = [min([x0u[i] for i in index_overlap]), max([x0u[i] for i in index_overlap])]
        return np.sort(index_overlap), len(index_overlap), boundaries
    
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
