import tqdm
import numpy as np
import pandas as pd
import sys
sys.path.insert(1, r".\utils")
sys.path.insert(1, r".\tracking")
from hit import Hit
from track3D import Track3D
from parameters import *
from physics import dist_line_rect


## This function finds the indices of event which are good candidate for muon decay : good tracks that don't end on a side of the detector
## have enough hits on both planes, and close to the reconstructed track, and for which the next event is not too long after and has hits
## close to the possible decay point.
def find_muon_decay(df, df_total, time_cutoff = 1500, spacial_cutoff = 4, \
                    save_indices = True, save_hits = False, save_stats = False, save_time_intervals = True,\
                    run_name = None, storage_dir = None, \
                    return_stats = True):
    """
    Arguments :
        -df : data frame filtered containing only the events with a certain range of n_hits
        -df_total : data frame containing all events
        -time_cutoff : maximal time interval in clock cycles over which a decay is searched
        -spacial_cutoff : maximal distance between the end of the muon track and the potential electron in the next event 
        -save_indices : if True, the indices of the events candidate for muon decay are stored in files with path {storage_dir"events_indices"run_name.txt}                     
        -save_hits : if True, the lists of hits are stored with pickle in {storage_dir"pickle_events"run_name} for each muon decay event
        -save_stats : if True, the function saves the filtering stats in a dictionary with pickle in {storage_dir"filtering_data"run_name}   
        -save_time_intervals : if True, the time intervals are stored in {storage_dir"time_intervals"run_name.txt}     
        -return_stats : if True, the function returns the number of event filtered out at each step of the algorithm

    Returns :
        -candidate_index : indices of the events considered by the algorithm as muon decay
        -time_intervals : time interval in clock cycle between the muon track and the decay for each decay
       The next return numbers are the stats returned if return_stats = True :       
        -low_number : number of events which contain less that 3 hits in one of the 2 planes
        -bottom_touch : number of events for which the last layer in x or y direction contains a hit
        -side_touch : number of events for which a hit with the lowest z coordinate is a the side of the detector
        -bad_fit : number of events for which the chi square value of the reconstructed track wasn't satisfactory
        -last_event_of_df : number of events for which the muon track is the last event of a run : the decay can't be accessed
        -hits_far_from_track : number of hits for which all hits at a distance higher than spacial_cutoff of the reconstructed track (preventive part of the code, shouldn't happen)
        -too_large_time_interval : number of events for which the next event happend too long after to be considered the product of a decay
        -no_spacial_correlation : number of events for which the hits in the next event are far from the end of the track and can thus not be considered
                                 as caused by an electron coming from a decay
    """
    candidate_index = []
    time_intervals = []
    low_number = 0
    bottom_touch = 0
    side_touch = 0
    bad_fit = 0
    last_event_of_df = 0
    too_large_time_interval = 0
    no_spacial_correlation = 0
    hits_far_from_track = 0

    if save_hits:
        decay_data = {'event_index': [], 'track_x0' : [], 'track_tx' : [], 'track_y0' : [], 'track_ty' : [], 'hits_muon': [], 'hits_electron': [], 'time_interval' : []}

    for index, row in tqdm(df.iterrows(), total = df.shape[0]):
        hits = [Hit(row,i) for i in range(row['n_hits'])]
        hitsX = [h for h in hits if h.is_sidex]
        hitsY = [h for h in hits if not h.is_sidex]
        
        ## Some events don't have three hits on one of the two sides and are thus not considered
        if len(hitsX) > 3 and len(hitsY) > 3:
            #One only considers the events for which the potential track ended inside the detector
            hitsX.sort(key = lambda hit: -hit.coord[1])
            hitsY.sort(key = lambda hit: -hit.coord[1])
            if hitsX[-1].coord[1] > 1 and hitsY[-1].coord[1] > 1:
                hitsX_last = [hit for hit in hitsX if hit.coord[1]==hitsX[-1].coord[1]]
                hitsY_last = [hit for hit in hitsY if hit.coord[1]==hitsY[-1].coord[1]]
                X_ok = True
                Y_ok = True
                for hit in hitsX_last:
                    if hit.coord[0] == 1 or hit.coord[0] == 24:
                        X_ok = False
                        side_touch += 1
                for hit in hitsY_last:
                    if hit.coord[0] == 1 or hit.coord[0] == 24:
                        Y_ok = False
                        if X_ok:
                            side_touch += 1
                if X_ok and Y_ok:
                    # get track parameters
                    track = Track3D(hits)

                    ## check if track has a "good" chi2 value
                    if track.is_good_2D_fit():
                        if index+1 >= len(df_total):
                            last_event_of_df += 1
                        else:
                            next_event = df_total.loc[index+1]
                        
    
                            hits_next_event = [Hit(next_event,i) for i in range(next_event['n_hits'])]
                            hitsX_next_event = [hit for hit in hits_next_event if hit.is_sidex]
                            hitsY_next_event = [hit for hit in hits_next_event if not hit.is_sidex]

                            hitsX = [hit for hit in hitsX if dist_line_rect(track.x.t, track.x.x0, hit.get_pos(), thickness, width) < spacial_cutoff] #Keep only the hits close to the track
                            hitsY = [hit for hit in hitsY if dist_line_rect(track.y.t, track.y.x0, hit.get_pos(), thickness, width) < spacial_cutoff]

                            ## check if there's still hits in the list after removing the ones far from the reconstructed track
                            if len(hitsX) != 0 or len(hitsY) != 0:

                                ## check if the next event happend close enough from the muon track for it to be the product of a decay
                                time_interval = hits_next_event[0].timestamp_event-hits[0].timestamp_event
                                if  time_interval < time_cutoff: 
                                    hitsX.sort(key = lambda hit: -hit.get_pos()[1])
                                    hitsY.sort(key = lambda hit: -hit.get_pos()[1])

                                    X_near = False
                                    Y_near = False

                                    ## Check if the hits in the next event are close to the end of the muon track
                                    if len(hitsX) != 0:
                                        for h in hitsX_next_event:
                                            if np.linalg.norm(h.get_pos()-hitsX[-1].get_pos()) < spacial_cutoff:
                                                X_near = True
                                    if len(hitsY) != 0:
                                        for h in hitsY_next_event:
                                            if np.linalg.norm(h.get_pos()-hitsY[-1].get_pos()) < spacial_cutoff:
                                                Y_near = True
                                    if X_near and Y_near:
                                        candidate_index.append(index)
                                        time_intervals.append(time_interval)
                                        if save_hits:
                                            decay_data['event_index'].append(index)
                                            decay_data['track_x0'].append(track.x.x0)
                                            decay_data['track_tx'].append(track.x.t)
                                            decay_data['track_y0'].append(track.y.x0)
                                            decay_data['track_ty'].append(track.y.t)
                                            decay_data['hits_muon'].append(hits)
                                            decay_data['hits_electron'].append(hits_next_event)
                                            decay_data['time_interval'].append(time_interval)

                                    else:
                                        no_spacial_correlation += 1
                                else:
                                    too_large_time_interval += 1
                            else:
                                hits_far_from_track += 1
                    else:
                        bad_fit += 1
            else:
                bottom_touch += 1
        else:
            low_number += 1
    if save_indices:
        np.savetxt(storage_dir+"events_indices"+run_name+".txt", candidate_index)
    if save_time_intervals:
        np.savetxt(storage_dir+"time_intervals"+run_name+".txt", time_intervals)
    if save_hits:
        decay_data = pd.DataFrame.from_dict(decay_data) # translate the dictionary into a pandas dataframe
        decay_data.to_pickle(storage_dir+"pickle_decay_data"+run_name)
    og_len = len(df_total)
    new_len = len(df)
    filtering = pd.DataFrame({'og_len' : [og_len],
                    'new_len' : [new_len],
                    'low_number' : [low_number],
                    'bottom_touch' : [bottom_touch],
                    'side_touch' : [side_touch],
                    'bad_fit': [bad_fit],
                    'last_event_of_df' : [last_event_of_df],
                    'too_large_time_interval' : [too_large_time_interval],
                    'hits_far_from_track' : [hits_far_from_track],
                    'no_spacial_correlation' : [no_spacial_correlation]})
    if save_stats:
        filtering.to_pickle(storage_dir+"filtering_data"+run_name)
    if return_stats:  
        return candidate_index, filtering
    else:
        return candidate_index, time_intervals