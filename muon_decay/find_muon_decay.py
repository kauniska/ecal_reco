from tqdm import tqdm
import numpy as np
import pandas as pd
import copy
import sys
sys.path.insert(1, ".\\utils")
sys.path.insert(1, ".\\tracking")
sys.path.insert(1, ".\\time_alignement")
from hit import Hit
from track3D import Track3D
from parameters import *
from physics import dist_line_rect
from time_correction import *
from track_reconstruction import mean_timestamp


# This function finds the indices of event which are good candidate for muon decay : good tracks that don't end on a side of the detector
# have enough hits on both planes, and close to the reconstructed track, and for which the next event is not too long after and has hits
# close to the possible decay point.
def find_muon_decay(df, df_total, time_max = 1500, time_min = 0, spacial_cutoff = 4, spacial_min = 0, \
                    spacial_max = 4, n_hits_min = 6, n_hits_max = 50,save_indices = True, \
                    save_hits = False, save_stats = False, save_time_intervals = True,\
                    save_distances = False, run_name = None, storage_dir = None, \
                    return_stats = True, time_corr = False):
    """
    Arguments :
        -df : data frame filtered containing only the events with a certain range of n_hits
        -df_total : data frame containing all events
        -time_max : maximal time interval in clock cycles over which a decay is searched
        -time_min : mininaml time interval in clock cycles over which a decay is searched
        -spacial_cutoff = maximal distance of an hit to the rest of track to be kept 
        -spacial_max : maximal distance between the end of the muon track and the potential electron in the next event 
        -spacial_min : minimal distance between the end of the muon track and the potential electron in the next event 
        -n_hits_min : minimal number of hits in a muon event
        -n_hits_max : maximal number of hits in a muon event
        -save_indices : if True, the indices of the events candidate for muon decay are stored in files with path {storage_dir"events_indices"run_name.txt}                     
        -save_hits : if True, the lists of hits are stored with pickle in {storage_dir"pickle_events"run_name} for each muon decay event
        -save_stats : if True, the function saves the filtering stats in a dictionary with pickle in {storage_dir"filtering_data"run_name}   
        -save_time_intervals : if True, the time intervals are stored in {storage_dir"time_intervals"run_name.txt}     
        -save_distances : if True, the spacial distances between the end of muon track and the electron are stored in {storage_dir"distances"run_name.txt}
        -return_stats : if True, the function returns the number of event filtered out at each step of the algorithm
        -time_corr : if True, apply time correction at the creation of Track3D

    Returns :
        -candidate_index : indices of the events considered by the algorithm as muon decay
        -time_intervals : time interval in clock cycle between the muon track and the decay for each decay
        -distancesX/Y : spatial distance between the end of the muon track and the potential electron in the next event (for X or Y)
       The next return numbers are the stats returned if return_stats = True : 
        -high_number : number of events which contain more than n_hits_max
        -low_number : number of events which contain less than n_hits_min
        -low_number_on_XY : number of events which contain less that 3 hits in one of the 2 planes
        -bottom_touch : number of events for which the last layer in x or y direction contains a hit
        -side_touch : number of events for which a hit with the lowest z coordinate is a the side of the detector
        -bad_fit : number of events for which the chi square value of the reconstructed track wasn't satisfactory
        -last_event_of_df : number of events for which the muon track is the last event of a run : the decay can't be accessed
        -hits_far_from_track : number of hits for which all hits at a distance higher than spacial_cutoff of the reconstructed track (preventive part of the code, shouldn't happen)
        -shower_upward : number events for which the next event have hits "upward" the end of the muon track
        -too_large_time_interval : number of events for which the next event happend too long after to be considered the product of a decay
        -too_small_time_interval : number of events for which the next event happend too short after to be considered the product of a decay
        -no_spacial_correlation : number of events for which the hits in the next event are far from the end of the track and can thus not be considered
                                 as caused by an electron coming from a decay
    """
    candidate_index = []
    time_intervals = []
    distancesX = []
    distancesY = []
    high_number =0
    low_number = 0
    low_number_on_XY = 0
    bottom_touch = 0
    side_touch = 0
    bad_fit = 0
    last_event_of_df = 0
    shower_upward = 0
    too_large_time_interval = 0
    too_small_time_interval = 0
    no_spacial_correlation = 0
    hits_far_from_track = 0
    row_count = 0
    index_count =0


   
    if save_hits:
        decay_data = {'event_index': [], 'track_x0' : [], 'track_tx' : [], 'track_y0' : [], 'track_ty' : [], 'hits_muon': [], 'hits_electron': [], 'time_interval' : [],'distanceX' : [],'distanceY' : []}

    for index, row in tqdm(df.iterrows(), total = df.shape[0]):
        hits = [Hit(row,i) for i in range(row['n_hits'])]
        hitsX = [h for h in hits if h.is_sidex]
        hitsY = [h for h in hits if not h.is_sidex]

        ##Criterion 1 :
        ## Some event have too much hits to be a clean muon track
        if len(hits) < n_hits_max :
            ## Some event don't have enough hits
            if len(hits) > n_hits_min :           
                ## Some events don't have three hits on one of the two sides and are thus not considered
                if len(hitsX) > 3 and len(hitsY) > 3:

                    ##Criterion 2 :
                    #One only considers the events for which the potential track ended inside the detector
                    hitsX.sort(key = lambda hit: -hit.coord[1])
                    hitsY.sort(key = lambda hit: -hit.coord[1])
                    if hitsX[-1].coord[1] > 1 and hitsY[-1].coord[1] > 1:

                        ##Criterion 3 :
                        ## Only keep tracks that doesn't end on the sides of the detector
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

                            ## Time Correction :
                            ## get track parameters and apply time correction
                            if time_corr : 
                                track = time_correction_global(Track3D(hits))
                                hitsX = track.x.hits
                                hitsY = track.y.hits
                                # print("After : "+str(mean_timestamp(hits)))
                            else :
                                track = Track3D(hits)


                            ## Criterion 4 :
                            ## check if track has a "good" chi2 value
                            if track.is_good_2D_fit():

                            ## check if there's still hits in the list after removing the ones far from the reconstructed track
                                hitsX = [hit for hit in hitsX if dist_line_rect(track.x.t, track.x.x0, hit.get_pos(), thickness, width) < spacial_cutoff] #Keep only the hits close to the track
                                hitsY = [hit for hit in hitsY if dist_line_rect(track.y.t, track.y.x0, hit.get_pos(), thickness, width) < spacial_cutoff]
                                if len(hitsX) != 0 and len(hitsY) != 0:

                                    if index+1 <= len(df_total):
                                
                                        next_event = df_total.loc[index+1]

                                        hits_next_event = [Hit(next_event,i) for i in range(next_event['n_hits'])]                       
                                    
                                        # Apply time correction on the hits of the next event
                                        if time_corr :  
                                            hits_next_event = copy.deepcopy(time_correction_global(track,hits_next_event))
                                            hitsX_next_event = [hit for hit in hits_next_event if hit.is_sidex]
                                            hitsY_next_event = [hit for hit in hits_next_event if not hit.is_sidex]
                                        else : 
                                            hitsX_next_event = [hit for hit in hits_next_event if hit.is_sidex]
                                            hitsY_next_event = [hit for hit in hits_next_event if not hit.is_sidex]

                                        ## Criterion 5 :
                                        ## check if the next event happend close enough from the muon track for it to be the product of a decay
                                        if time_corr :
                                        ## The mean value of all singles timestamps of hits in event are computed (and add to the timestamp_event)
                                            # time_interval = hits_next_event[0].timestamp_event + mean_timestamp(hitsX_next_event,hitsY_next_event)- (hits[0].timestamp_event + mean_timestamp(hitsX, hitsY))
                                            # time_interval = hits_next_event[0].timestamp_event + first_timestamp(hits_next_event)- (hits[0].timestamp_event + mean_timestamp(hitsX, hitsY))
                                            time_interval = hits_next_event[0].timestamp_event + first_timestamp(hits_next_event)  - (hits[0].timestamp_event + mean_timestamp(hitsX,hitsY) )
                                        else :
                                            time_interval = hits_next_event[0].timestamp_event  - (hits[0].timestamp_event )
                        
                                        if  time_interval < time_max:
                                            if time_min < time_interval :

                                    
                                                    
                                                hitsX.sort(key = lambda hit: -hit.get_pos()[1])
                                                hitsY.sort(key = lambda hit: -hit.get_pos()[1])

                                                ## Criterion 6 :
                                                ## Check if the electron does not go "upward" in the direction of muon track
                                                X_upward = False
                                                Y_upward = False
                                                for h in hitsX_next_event :
                                                    if h.get_pos()[1] > hitsX[-1].get_pos()[1] :
                                                        X_upward = True
                                                for h in hitsY_next_event :
                                                    if h.get_pos()[1] > hitsY[-1].get_pos()[1] :
                                                        Y_upward = True
                                                if not X_upward and not Y_upward :

                                                    ## Criterion 7 :
                                                    ## Check if the hits in the next event are close to the end of the muon track
                                                    rX = 1e100000
                                                    rY = 1e100000
                                                    if len(hitsX) != 0 :                     
                                                        for h in hitsX_next_event:
                                                            rX_int = np.linalg.norm(h.get_pos()-hitsX[-1].get_pos())
                                                            if rX_int < rX:
                                                                rX = rX_int
                                                    
                                                    if len(hitsY) != 0:
                                                        for h in hitsY_next_event:
                                                            rY_int = np.linalg.norm(h.get_pos()-hitsY[-1].get_pos())
                                                            if rY_int < rY:
                                                                rY = rY_int

                                                    if  rX < spacial_max and rY < spacial_max: 
                                                        if spacial_min <= rX and spacial_min <= rY :
                                                        
                                                    

                                                            ## Saving results :
                                                            candidate_index.append(index)
                                                            time_intervals.append(time_interval)
                                                            if save_distances:
                                                                distancesX.append(rX)
                                                                distancesY.append(rY)
                                                            if save_hits:
                                                                decay_data['event_index'].append(index)
                                                                decay_data['track_x0'].append(track.x.x0)
                                                                decay_data['track_tx'].append(track.x.t)
                                                                decay_data['track_y0'].append(track.y.x0)
                                                                decay_data['track_ty'].append(track.y.t)
                                                                decay_data['hits_muon'].append(hits)
                                                                decay_data['hits_electron'].append(hits_next_event)
                                                                decay_data['time_interval'].append(time_interval)
                                                                if save_distances: 
                                                                    decay_data['distanceX'].append(rX)
                                                                    decay_data['distanceY'].append(rY)
                                                                else : 
                                                                    decay_data['distanceX'].append(0)
                                                                    decay_data['distanceY'].append(0)

                                                        else :
                                                            no_spacial_correlation += 1
                                                    else:
                                                        no_spacial_correlation += 1
                                                else :
                                                    shower_upward +=1
                                            else :
                                                too_small_interval +=1
                                        else :
                                            too_large_time_interval += 1
                                    else :
                                        last_event_of_df += 1
                                else:
                                    hits_far_from_track += 1
                                        
                            else:
                                bad_fit += 1
                            
                    else:
                        bottom_touch += 1
                else:
                    low_number_on_XY += 1
            else : 
                low_number +=1
        else : 
            high_number +=1
 
    if save_indices:
        np.savetxt(storage_dir+"events_indices"+run_name+".txt", candidate_index)
    if save_time_intervals:
        np.savetxt(storage_dir+"time_intervals"+run_name+".txt", time_intervals)
    if save_distances:
        np.savetxt(storage_dir+"distancesX"+run_name+".txt", distancesX)
        np.savetxt(storage_dir+"distancesY"+run_name+".txt", distancesY)
    if save_hits:
        decay_data = pd.DataFrame.from_dict(decay_data) # translate the dictionary into a pandas dataframe
        decay_data.to_pickle(storage_dir+"pickle_decay_data"+run_name)
    og_len = len(df_total)
    new_len = len(df)
    filtering = pd.DataFrame({'og_len' : [og_len],
                    'new_len' : [new_len],
                    'high_number' : [high_number],
                    'low_number' : [low_number],
                    'low_number_on_XY' : [low_number_on_XY],
                    'bottom_touch' : [bottom_touch],
                    'side_touch' : [side_touch],
                    'bad_fit': [bad_fit],
                    'last_event_of_df' : [last_event_of_df],
                    'too_large_time_interval' : [too_large_time_interval],
                    'too_small_time_interval' : [too_small_time_interval],
                    'hits_far_from_track' : [hits_far_from_track],
                    'shower_upward' : [shower_upward],
                    'no_spacial_correlation' : [no_spacial_correlation]})
    if save_stats:
        filtering.to_pickle(storage_dir+"filtering_data"+run_name)
    if return_stats:  
        print(filtering)
        return candidate_index, filtering
    else:
        return candidate_index, time_intervals




# ## Old version of code without modification

# from tqdm import tqdm
# import numpy as np
# import pandas as pd
# import sys
# sys.path.insert(1, ".\\utils")
# sys.path.insert(1, ".\\tracking")
# from hit import Hit
# from track3D import Track3D
# from parameters import *
# from physics import dist_line_rect


# ## This function finds the indices of event which are good candidate for muon decay : good tracks that don't end on a side of the detector
# ## have enough hits on both planes, and close to the reconstructed track, and for which the next event is not too long after and has hits
# ## close to the possible decay point.
# def find_muon_decay(df, df_total, time_cutoff = 1500, spacial_cutoff = 4, \
#                     save_indices = True, save_hits = False, save_stats = False, save_time_intervals = True,\
#                     run_name = None, storage_dir = None, \
#                     return_stats = True, time_corr=False):
#     """
#     Arguments :
#         -df : data frame filtered containing only the events with a certain range of n_hits
#         -df_total : data frame containing all events
#         -time_cutoff : maximal time interval in clock cycles over which a decay is searched
#         -spacial_cutoff : maximal distance between the end of the muon track and the potential electron in the next event 
#         -save_indices : if True, the indices of the events candidate for muon decay are stored in files with path {storage_dir"events_indices"run_name.txt}                     
#         -save_hits : if True, the lists of hits are stored with pickle in {storage_dir"pickle_events"run_name} for each muon decay event
#         -save_stats : if True, the function saves the filtering stats in a dictionary with pickle in {storage_dir"filtering_data"run_name}   
#         -save_time_intervals : if True, the time intervals are stored in {storage_dir"time_intervals"run_name.txt}     
#         -return_stats : if True, the function returns the number of event filtered out at each step of the algorithm

#     Returns :
#         -candidate_index : indices of the events considered by the algorithm as muon decay
#         -time_intervals : time interval in clock cycle between the muon track and the decay for each decay
#        The next return numbers are the stats returned if return_stats = True :       
#         -low_number : number of events which contain less that 3 hits in one of the 2 planes
#         -bottom_touch : number of events for which the last layer in x or y direction contains a hit
#         -side_touch : number of events for which a hit with the lowest z coordinate is a the side of the detector
#         -bad_fit : number of events for which the chi square value of the reconstructed track wasn't satisfactory
#         -last_event_of_df : number of events for which the muon track is the last event of a run : the decay can't be accessed
#         -hits_far_from_track : number of hits for which all hits at a distance higher than spacial_cutoff of the reconstructed track (preventive part of the code, shouldn't happen)
#         -too_large_time_interval : number of events for which the next event happend too long after to be considered the product of a decay
#         -no_spacial_correlation : number of events for which the hits in the next event are far from the end of the track and can thus not be considered
#                                  as caused by an electron coming from a decay
#     """
#     candidate_index = []
#     time_intervals = []
#     low_number = 0
#     bottom_touch = 0
#     side_touch = 0
#     bad_fit = 0
#     last_event_of_df = 0
#     too_large_time_interval = 0
#     no_spacial_correlation = 0
#     hits_far_from_track = 0

#     if save_hits:
#         decay_data = {'event_index': [], 'track_x0' : [], 'track_tx' : [], 'track_y0' : [], 'track_ty' : [], 'hits_muon': [], 'hits_electron': [], 'time_interval' : []}

#     for index, row in tqdm(df.iterrows(), total = df.shape[0]):
#         hits = [Hit(row,i) for i in range(row['n_hits'])]
#         hitsX = [h for h in hits if h.is_sidex]
#         hitsY = [h for h in hits if not h.is_sidex]
        
#         ## Some events don't have three hits on one of the two sides and are thus not considered
#         if len(hitsX) > 3 and len(hitsY) > 3:
#             #One only considers the events for which the potential track ended inside the detector
#             hitsX.sort(key = lambda hit: -hit.coord[1])
#             hitsY.sort(key = lambda hit: -hit.coord[1])
#             if hitsX[-1].coord[1] > 1 and hitsY[-1].coord[1] > 1:
#                 hitsX_last = [hit for hit in hitsX if hit.coord[1]==hitsX[-1].coord[1]]
#                 hitsY_last = [hit for hit in hitsY if hit.coord[1]==hitsY[-1].coord[1]]
#                 X_ok = True
#                 Y_ok = True
#                 for hit in hitsX_last:
#                     if hit.coord[0] == 1 or hit.coord[0] == 24:
#                         X_ok = False
#                         side_touch += 1
#                 for hit in hitsY_last:
#                     if hit.coord[0] == 1 or hit.coord[0] == 24:
#                         Y_ok = False
#                         if X_ok:
#                             side_touch += 1
#                 if X_ok and Y_ok:
#                     # get track parameters
#                     track = Track3D(hits)

#                     ## check if track has a "good" chi2 value
#                     if track.is_good_2D_fit():
#                         if index+1 >= len(df_total):
#                             last_event_of_df += 1
#                         else:
#                             next_event = df_total.loc[index+1]
                        
    
#                             hits_next_event = [Hit(next_event,i) for i in range(next_event['n_hits'])]
#                             hitsX_next_event = [hit for hit in hits_next_event if hit.is_sidex]
#                             hitsY_next_event = [hit for hit in hits_next_event if not hit.is_sidex]

#                             hitsX = [hit for hit in hitsX if dist_line_rect(track.x.t, track.x.x0, hit.get_pos(), thickness, width) < spacial_cutoff] #Keep only the hits close to the track
#                             hitsY = [hit for hit in hitsY if dist_line_rect(track.y.t, track.y.x0, hit.get_pos(), thickness, width) < spacial_cutoff]

#                             ## check if there's still hits in the list after removing the ones far from the reconstructed track
#                             if len(hitsX) != 0 or len(hitsY) != 0:

#                                 ## check if the next event happend close enough from the muon track for it to be the product of a decay
#                                 time_interval = hits_next_event[0].timestamp_event-hits[0].timestamp_event
#                                 if  time_interval < time_cutoff: 
#                                     hitsX.sort(key = lambda hit: -hit.get_pos()[1])
#                                     hitsY.sort(key = lambda hit: -hit.get_pos()[1])

#                                     X_near = False
#                                     Y_near = False

#                                     ## Check if the hits in the next event are close to the end of the muon track
#                                     if len(hitsX) != 0:
#                                         for h in hitsX_next_event:
#                                             if np.linalg.norm(h.get_pos()-hitsX[-1].get_pos()) < spacial_cutoff:
#                                                 X_near = True
#                                     if len(hitsY) != 0:
#                                         for h in hitsY_next_event:
#                                             if np.linalg.norm(h.get_pos()-hitsY[-1].get_pos()) < spacial_cutoff:
#                                                 Y_near = True
#                                     if X_near and Y_near:
#                                         candidate_index.append(index)
#                                         time_intervals.append(time_interval)
#                                         if save_hits:
#                                             decay_data['event_index'].append(index)
#                                             decay_data['track_x0'].append(track.x.x0)
#                                             decay_data['track_tx'].append(track.x.t)
#                                             decay_data['track_y0'].append(track.y.x0)
#                                             decay_data['track_ty'].append(track.y.t)
#                                             decay_data['hits_muon'].append(hits)
#                                             decay_data['hits_electron'].append(hits_next_event)
#                                             decay_data['time_interval'].append(time_interval)

#                                     else:
#                                         no_spacial_correlation += 1
#                                 else:
#                                     too_large_time_interval += 1
#                             else:
#                                 hits_far_from_track += 1
#                     else:
#                         bad_fit += 1
#             else:
#                 bottom_touch += 1
#         else:
#             low_number += 1
#     if save_indices:
#         np.savetxt(storage_dir+"events_indices"+run_name+".txt", candidate_index)
#     if save_time_intervals:
#         np.savetxt(storage_dir+"time_intervals"+run_name+".txt", time_intervals)
#     if save_hits:
#         decay_data = pd.DataFrame.from_dict(decay_data) # translate the dictionary into a pandas dataframe
#         decay_data.to_pickle(storage_dir+"pickle_decay_data"+run_name)
#     og_len = len(df_total)
#     new_len = len(df)
#     filtering = pd.DataFrame({'og_len' : [og_len],
#                     'new_len' : [new_len],
#                     'low_number' : [low_number],
#                     'bottom_touch' : [bottom_touch],
#                     'side_touch' : [side_touch],
#                     'bad_fit': [bad_fit],
#                     'last_event_of_df' : [last_event_of_df],
#                     'too_large_time_interval' : [too_large_time_interval],
#                     'hits_far_from_track' : [hits_far_from_track],
#                     'no_spacial_correlation' : [no_spacial_correlation]})
#     if save_stats:
#         filtering.to_pickle(storage_dir+"filtering_data"+run_name)
#     if return_stats:  
#         return candidate_index, filtering
#     else:
#         return candidate_index, time_intervals