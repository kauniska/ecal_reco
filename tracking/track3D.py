import numpy as np
import matplotlib.pyplot as plt
from track import Track

class Track3D:
    def __init__(self, *args):
        """Creates a 3D Track class, containing two 2D tracks
        """
        if len(args) == 0:
            self.x = None
            self.y = None
            self.time = None
        elif len(args) == 1:
            hitsX = [hit for hit in args[0] if hit.is_sidex]
            hitsY = [hit for hit in args[0] if not hit.is_sidex]
            self.x = Track(hitsX)
            self.y = Track(hitsY)
            self.time = self.x.get_timestamps() + self.y.get_timestamps()
        elif len(args) == 2:
            self.x = args[0]
            self.y = args[1]
            self.time = self.get_time()

            
    def get_time_interval(self):
        """Get the timings of all hits, returns None if zero hits are associated to the track
        
        Return:
            self.time (List of float): timestamps of hits
        """
        self.time = np.concatenate([self.x.get_timestamps(), self.x.get_timestamps()])
        hits = self.x.hits + self.y.hits
        hits_indices = self.x.hits_index + self.y.hits_index
        if len(hits_indices) > 0:
            last_hit = hits[np.argmax(hits_indices, axis = 0)]
            if last_hit.coord[1] < 8:
                distances = self.x._dr(last_hit, hits)
                if np.any(np.array(distances) < 2):
                    return np.mean(self.time) - hits[0].timestamp_event
        else:
            return None
        
    def precise_track(self):
        self.x.precise_track()
        self.y.precise_track()
        
    def is_good_2D_fit(self):
        return (self.x.is_good_fit() and self.y.is_good_fit())
    
    def reduced_chi2(self):
        return self.x.reduced_chi2 + self.y.reduced_chi2
    
    def print(self, plot = False):
        fig, axs = plt.subplots(1, 2)
        fig.set_size_inches(12, 6)
        self.x.print(plot, axs[0])
        self.y.print(plot, axs[1])
        # fig.show()
        # fig.savefig('test.png')
        
    def kalman_filter(self, sigma = 0.5):
        fig, axs = plt.subplots(1, 2)
        fig.set_size_inches(12, 6)
        self.x.kalman_filter(sigma, axs[0])
        self.y.kalman_filter(sigma, axs[1])