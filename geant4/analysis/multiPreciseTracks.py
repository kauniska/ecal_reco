from multiprocessing import Pool
import pickle

def doPreciseTrack(track):
    track.precise_track()
    return track

def multiPreciseTracks(tracks):
    pool = Pool()
    precise_tracks = pool.map(doPreciseTrack, tracks)
    with open('angle_precise.pkl', 'wb') as oup:
        pickle.dump(precise_tracks, oup, protocol=pickle.HIGHEST_PROTOCOL)